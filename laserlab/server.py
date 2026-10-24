"""HTTP 服务与 JSON API.

提供无第三方依赖的纯标准库 HTTP 服务器，为 Web 前端提供 JSON 接口。

端点
----
GET /                返回静态前端页面
GET /api/health      健康检查
GET /api/cavity      法布里-珀罗腔分析（参数：length, r1, r2, n, loss）
GET /api/beam        稳定腔高斯光束（参数：length, radius1, radius2, wavelength）
GET /api/relax       弛豫振荡仿真（参数：pump_ratio, t_end）
GET /api/modes       纵模与增益（参数：medium, length, g0, nmodes）
"""

import json
import os
import urllib.parse
from http.server import BaseHTTPRequestHandler, HTTPServer

from . import cavity as _cavity
from . import beam as _beam
from . import gain as _gain
from . import dynamics as _dyn
from . import rate as _rate

_WEB_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "web")


def _json_handler(handler, data):
    body = json.dumps(data, ensure_ascii=False).encode("utf-8")
    handler.send_response(200)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Access-Control-Allow-Origin", "*")
    handler.end_headers()
    handler.wfile.write(body)


def _float(params, key, default):
    v = params.get(key, [default])[0]
    if v in ("", "inf", "Infinity"):
        return float("inf")
    return float(v)


def api_cavity(params):
    length = _float(params, "length", 0.2)
    n = _float(params, "n", 1.0)
    r1 = _float(params, "r1", 0.99)
    r2 = _float(params, "r2", 0.98)
    loss = _float(params, "loss", 0.0)
    fp = _cavity.FabryPerot(n, length, r1, r2, loss)
    # 采样透射谱
    center = _cavity.C.wl_to_freq(1064.0)
    freqs = [center + k * fp.fsr_thz / 200.0 for k in range(-200, 201)]
    trans = [fp.transmittance(f) for f in freqs]
    return {
        "fsr_hz": fp.fsr_hz,
        "fsr_thz": fp.fsr_thz,
        "finesse": fp.finesse,
        "linewidth_hz": fp.linewidth_hz(),
        "peak_transmittance": fp.peak_transmittance(),
        "roundtrip_gain_required": fp.roundtrip_gain_required(),
        "spectrum": {"freq_thz": freqs, "transmittance": trans},
    }


def api_beam(params):
    length = _float(params, "length", 0.2)
    r1 = _float(params, "radius1", float("inf"))
    r2 = _float(params, "radius2", 0.5)
    wl = _float(params, "wavelength", 1064.0)
    try:
        gc = _beam.GaussianCavity(length, r1, r2, wl)
    except ValueError as e:
        return {"error": str(e)}
    zs = [length * i / 200.0 for i in range(201)]
    radii = [gc.beam_radius_at(z) * 1e6 for z in zs]
    return {
        "waist_radius_um": gc.waist_radius * 1e6,
        "waist_position_m": gc.waist_position,
        "rayleigh_m": gc.rayleigh,
        "divergence_mrad": _beam.divergence_half_angle(gc.waist_radius, wl) * 1e3,
        "profile": {"z_m": zs, "radius_um": radii},
    }


def api_relax(params):
    pump_ratio = _float(params, "pump_ratio", 2.0)
    t_end = _float(params, "t_end", 1e-3)
    tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1.0e-8
    rp_th = _rate.threshold_pump_rate(tau, sigma, tau_c)
    rp = pump_ratio * rp_th
    sim = _dyn.simulate_relaxation(rp, tau_upper=tau, sigma=sigma, photon_lifetime=tau_c,
                                   t_end=t_end)
    n_ss, phi_ss = _dyn.steady_state(rp, tau, sigma, tau_c)
    f_ro = _dyn.relaxation_frequency(pump_ratio, tau, tau_c)
    # 下采样到 <= 2000 点
    stride = max(1, len(sim.t) // 2000)
    return {
        "threshold_pump_rate": rp_th,
        "steady_phi": phi_ss,
        "relaxation_freq_hz": f_ro,
        "t_us": [t * 1e6 for t in sim.t[::stride]],
        "phi": sim.phi[::stride],
        "N": sim.N[::stride],
        "peak_phi": max(sim.phi),
    }


def api_modes(params):
    medium = params.get("medium", ["nd_yag_1064"])[0]
    length = _float(params, "length", 0.2)
    g0 = _float(params, "g0", 10.0)
    nmodes = int(_float(params, "nmodes", 21))
    fp = _cavity.FabryPerot(1.0, length, 0.99, 0.95)
    med = _gain.make_medium(medium, g0)
    # 纵模对称分布于增益中心两侧
    fsr = fp.fsr_thz
    half = nmodes // 2
    freqs = [med.center_thz + (k - half) * fsr for k in range(nmodes)]
    # 连续增益谱（用于前端画线型）
    span = nmodes * fp.fsr_thz
    cont = [med.center_thz + span * (i / 400.0 - 0.5) for i in range(401)]
    gains = [med.small_signal_gain(f) for f in cont]
    mode_gains = [med.small_signal_gain(f) for f in freqs]
    return {
        "medium": medium,
        "center_thz": med.center_thz,
        "fsr_thz": fp.fsr_thz,
        "spectrum": {"freq_thz": cont, "gain": gains},
        "modes": [{"freq_thz": f, "gain": g} for f, g in zip(freqs, mode_gains)],
    }


_ROUTES = {
    "/api/health": lambda p: {"status": "ok", "version": "1.0.0"},
    "/api/cavity": api_cavity,
    "/api/beam": api_beam,
    "/api/relax": api_relax,
    "/api/modes": api_modes,
}


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        params = urllib.parse.parse_qs(parsed.query)
        if path in _ROUTES:
            try:
                _json_handler(self, _ROUTES[path](params))
            except Exception as e:  # noqa: BLE001
                _json_handler(self, {"error": str(e)})
            return
        if path in ("/", "/index.html"):
            self._serve_static("index.html")
            return
        # 静态资源
        fname = os.path.basename(path)
        if fname.endswith((".js", ".css")) and os.path.exists(os.path.join(_WEB_DIR, fname)):
            self._serve_static(fname)
            return
        self.send_response(404)
        self.end_headers()
        self.wfile.write(b"not found")

    def _serve_static(self, fname):
        fp = os.path.join(_WEB_DIR, fname)
        ctype = "text/html" if fname.endswith(".html") else (
            "application/javascript" if fname.endswith(".js") else "text/css")
        with open(fp, "rb") as f:
            body = f.read()
        self.send_response(200)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        pass  # 静默日志


def run(host="127.0.0.1", port=8779):
    srv = HTTPServer((host, port), Handler)
    print(f"LaserLab server running at http://{host}:{port}")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


if __name__ == "__main__":
    run()
