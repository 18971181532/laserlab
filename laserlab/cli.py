"""命令行接口.

用法示例
--------
.. code-block:: bash

    python run.py cavity --length 0.2 --r1 0.99 --r2 0.98
    python run.py beam --length 0.2 --r1 inf --r2 0.5 --wavelength 1064
    python run.py relax --pump-ratio 2.0 --t-end 0.001
    python run.py modes --medium nd_yag_1064 --length 0.2 --g0 10
    python run.py serve --port 8779
"""

import argparse
import json
import sys


def _cmd_cavity(args):
    from .cavity import FabryPerot, StableCavity

    fp = FabryPerot(args.n, args.length, args.r1, args.r2, args.loss)
    sc = StableCavity(args.length, args.radius1, args.radius2) if args.radius1 else None
    out = {
        "fsr_hz": fp.fsr_hz,
        "fsr_thz": fp.fsr_thz,
        "finesse": fp.finesse,
        "linewidth_hz": fp.linewidth_hz(),
        "peak_transmittance": fp.peak_transmittance(),
        "roundtrip_gain_required": fp.roundtrip_gain_required(),
    }
    if sc is not None:
        out["g1"] = sc.g1
        out["g2"] = sc.g2
        out["g_product"] = sc.g_product
        out["stability"] = sc.classify()
    print(json.dumps(out, indent=2, ensure_ascii=False))


def _cmd_beam(args):
    from .beam import GaussianCavity

    r1 = args.radius1 if args.radius1 != float("inf") else float("inf")
    r2 = args.radius2 if args.radius2 != float("inf") else float("inf")
    gc = GaussianCavity(args.length, r1, r2, args.wavelength)
    out = {
        "waist_radius_um": gc.waist_radius * 1e6,
        "waist_position_m": gc.waist_position,
        "rayleigh_range_m": gc.rayleigh,
        "w_at_m1_um": gc.beam_radius_at(0.0) * 1e6,
        "w_at_m2_um": gc.beam_radius_at(args.length) * 1e6,
        "divergence_mrad": __import__("laserlab.beam", fromlist=["divergence_half_angle"])
        .divergence_half_angle(gc.waist_radius, args.wavelength) * 1e3,
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))


def _cmd_relax(args):
    from .rate import threshold_pump_rate
    from .dynamics import simulate_relaxation, relaxation_frequency

    # 使用 Nd:YAG 典型参数
    tau = 2.3e-4
    sigma = 2.8e-23
    tau_c = 1.0e-8
    rp_th = threshold_pump_rate(tau, sigma, tau_c)
    rp = args.pump_ratio * rp_th
    sim = simulate_relaxation(rp, tau_upper=tau, sigma=sigma, photon_lifetime=tau_c,
                              t_end=args.t_end)
    # 求稳态
    from .dynamics import steady_state
    n_ss, phi_ss = steady_state(rp, tau, sigma, tau_c)
    f_ro = relaxation_frequency(args.pump_ratio, tau, tau_c)
    out = {
        "threshold_pump_rate": rp_th,
        "pump_rate": rp,
        "steady_N": n_ss,
        "steady_phi": phi_ss,
        "relaxation_freq_hz": f_ro,
        "final_phi": sim.phi[-1],
        "peak_phi": max(sim.phi),
        "n_samples": len(sim.t),
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))


def _cmd_modes(args):
    from .cavity import FabryPerot
    from .gain import make_medium

    fp = FabryPerot(1.0, args.length, 0.99, 0.95)
    med = make_medium(args.medium, args.g0)
    freqs = fp.mode_frequencies(med.center_thz, 0, args.nmodes)
    gains = [med.small_signal_gain(f) for f in freqs]
    out = {
        "medium": args.medium,
        "center_thz": med.center_thz,
        "fsr_thz": fp.fsr_thz,
        "modes": [{"freq_thz": f, "gain_1pm": g} for f, g in zip(freqs, gains)],
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))


def _cmd_serve(args):
    from .server import run
    run(args.host, args.port)


def main(argv=None):
    p = argparse.ArgumentParser(prog="laserlab", description="激光器物理仿真工作台")
    sub = p.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("cavity", help="法布里-珀罗腔分析")
    pc.add_argument("--length", type=float, default=0.2)
    pc.add_argument("--n", type=float, default=1.0)
    pc.add_argument("--r1", type=float, default=0.99)
    pc.add_argument("--r2", type=float, default=0.98)
    pc.add_argument("--loss", type=float, default=0.0)
    pc.add_argument("--radius1", type=float, default=None)
    pc.add_argument("--radius2", type=float, default=None)
    pc.set_defaults(func=_cmd_cavity)

    pb = sub.add_parser("beam", help="稳定腔高斯光束")
    pb.add_argument("--length", type=float, default=0.2)
    pb.add_argument("--radius1", type=float, default=float("inf"))
    pb.add_argument("--radius2", type=float, default=0.5)
    pb.add_argument("--wavelength", type=float, default=1064.0)
    pb.set_defaults(func=_cmd_beam)

    pr = sub.add_parser("relax", help="弛豫振荡仿真")
    pr.add_argument("--pump-ratio", type=float, default=2.0)
    pr.add_argument("--t-end", type=float, default=1e-3)
    pr.set_defaults(func=_cmd_relax)

    pm = sub.add_parser("modes", help="纵模频率与增益")
    pm.add_argument("--medium", default="nd_yag_1064")
    pm.add_argument("--length", type=float, default=0.2)
    pm.add_argument("--g0", type=float, default=10.0)
    pm.add_argument("--nmodes", type=int, default=11)
    pm.set_defaults(func=_cmd_modes)

    ps = sub.add_parser("serve", help="启动 Web 界面")
    ps.add_argument("--host", default="127.0.0.1")
    ps.add_argument("--port", type=int, default=8779)
    ps.set_defaults(func=_cmd_serve)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
