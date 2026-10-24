// LaserLab 前端逻辑（零依赖，手写 Canvas 绘图）
(function () {
  "use strict";

  const API = "/api";

  async function get(path, params) {
    const qs = new URLSearchParams(params).toString();
    const res = await fetch(API + path + (qs ? "?" + qs : ""));
    return res.json();
  }

  // ---------- 通用绘图工具 ----------
  function clearCanvas(cv) {
    const ctx = cv.getContext("2d");
    ctx.clearRect(0, 0, cv.width, cv.height);
    ctx.fillStyle = "#0f172a";
    ctx.fillRect(0, 0, cv.width, cv.height);
    return ctx;
  }

  function drawAxes(ctx, w, h, pad, yLabel) {
    ctx.strokeStyle = "#334155";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(pad.l, h - pad.b);
    ctx.lineTo(w - pad.r, h - pad.b);
    ctx.moveTo(pad.l, pad.t);
    ctx.lineTo(pad.l, h - pad.b);
    ctx.stroke();
    ctx.fillStyle = "#94a3b8";
    ctx.font = "11px monospace";
    ctx.textAlign = "center";
    ctx.fillText(yLabel, pad.l - 8, pad.t - 8);
  }

  function niceTicks(min, max, n) {
    const span = max - min;
    const step = span / n;
    const mag = Math.pow(10, Math.floor(Math.log10(step)));
    const norm = step / mag;
    const sn = norm >= 5 ? 10 : norm >= 2 ? 5 : norm >= 1 ? 2 : 1;
    const s = sn * mag;
    const start = Math.ceil(min / s) * s;
    const ticks = [];
    for (let v = start; v <= max; v += s) ticks.push(v);
    return ticks;
  }

  function fmt(v, digits) {
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1e9 || a < 1e-3) return v.toExponential(digits == null ? 2 : digits);
    return v.toPrecision(digits == null ? 4 : digits);
  }

  // ---------- 谐振腔 ----------
  async function renderCavity() {
    const params = {
      length: +document.getElementById("cav-length").value,
      r1: +document.getElementById("cav-r1").value,
      r2: +document.getElementById("cav-r2").value,
      n: +document.getElementById("cav-n").value,
      loss: +document.getElementById("cav-loss").value,
    };
    const d = await get("/cavity", params);
    const m = document.getElementById("cav-metrics");
    m.innerHTML = `
      <div class="metric"><div class="k">自由光谱范围 FSR</div><div class="v">${fmt(d.fsr_thz)} <span class="unit">THz</span></div></div>
      <div class="metric"><div class="k">精细度</div><div class="v">${fmt(d.finesse)}</div></div>
      <div class="metric"><div class="k">线宽 FWHM</div><div class="v">${fmt(d.linewidth_hz)} <span class="unit">Hz</span></div></div>
      <div class="metric"><div class="k">峰值透射率</div><div class="v">${fmt(d.peak_transmittance)}</div></div>
      <div class="metric"><div class="k">阈值往返增益</div><div class="v">${fmt(d.roundtrip_gain_required)}</div></div>`;
    drawCavityCanvas(d);
  }

  function drawCavityCanvas(d) {
    const cv = document.getElementById("cav-canvas");
    const ctx = clearCanvas(cv);
    const w = cv.width, h = cv.height, pad = { l: 60, r: 20, t: 24, b: 40 };
    const freqs = d.spectrum.freq_thz, trans = d.spectrum.transmittance;
    const f0 = freqs[0], f1 = freqs[freqs.length - 1];
    const x = f => pad.l + (f - f0) / (f1 - f0) * (w - pad.l - pad.r);
    const y = t => h - pad.b - t * (h - pad.t - pad.b);
    drawAxes(ctx, w, h, pad, "透射率");
    // 网格
    ctx.strokeStyle = "#1e293b";
    for (const tv of [0, 0.25, 0.5, 0.75, 1]) {
      ctx.beginPath(); ctx.moveTo(pad.l, y(tv)); ctx.lineTo(w - pad.r, y(tv)); ctx.stroke();
      ctx.fillStyle = "#64748b"; ctx.textAlign = "right";
      ctx.fillText(tv.toFixed(2), pad.l - 6, y(tv) + 4);
    }
    // 透射谱
    ctx.strokeStyle = "#e63946";
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < freqs.length; i++) {
      const px = x(freqs[i]), py = y(trans[i]);
      i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
    }
    ctx.stroke();
    // 标注 FSR
    const fsr = d.fsr_thz;
    const centerIdx = Math.floor(freqs.length / 2);
    const cf = freqs[centerIdx];
    ctx.strokeStyle = "#2a9d8f";
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(x(cf), pad.t); ctx.lineTo(x(cf), h - pad.b); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x(cf + fsr), pad.t); ctx.lineTo(x(cf + fsr), h - pad.b); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#2a9d8f"; ctx.textAlign = "center";
    ctx.fillText("FSR = " + fmt(fsr) + " THz", (x(cf) + x(cf + fsr)) / 2, pad.t + 12);
  }

  // ---------- 高斯光束 ----------
  async function renderBeam() {
    const r1v = document.getElementById("beam-r1").value.trim();
    const params = {
      length: +document.getElementById("beam-length").value,
      radius1: r1v === "inf" || r1v === "" ? "inf" : +r1v,
      radius2: +document.getElementById("beam-r2").value,
      wavelength: +document.getElementById("beam-wl").value,
    };
    const d = await get("/beam", params);
    const m = document.getElementById("beam-metrics");
    if (d.error) {
      m.innerHTML = `<div class="metric"><div class="k">提示</div><div class="v">${d.error}</div></div>`;
      clearCanvas(document.getElementById("beam-canvas"));
      return;
    }
    m.innerHTML = `
      <div class="metric"><div class="k">束腰半径 w₀</div><div class="v">${fmt(d.waist_radius_um)} <span class="unit">μm</span></div></div>
      <div class="metric"><div class="k">束腰位置</div><div class="v">${fmt(d.waist_position_m)} <span class="unit">m</span></div></div>
      <div class="metric"><div class="k">瑞利长度 z_R</div><div class="v">${fmt(d.rayleigh_m)} <span class="unit">m</span></div></div>
      <div class="metric"><div class="k">远场发散半角</div><div class="v">${fmt(d.divergence_mrad)} <span class="unit">mrad</span></div></div>`;
    drawBeamCanvas(d);
  }

  function drawBeamCanvas(d) {
    const cv = document.getElementById("beam-canvas");
    const ctx = clearCanvas(cv);
    const w = cv.width, h = cv.height, pad = { l: 60, r: 20, t: 24, b: 40 };
    const zs = d.profile.z_m, radii = d.profile.radius_um;
    const z0 = zs[0], z1 = zs[zs.length - 1];
    const rmax = Math.max(...radii) * 1.15 || 1;
    const x = z => pad.l + (z - z0) / (z1 - z0) * (w - pad.l - pad.r);
    const y = r => h / 2 - (r / rmax) * ((h - pad.t - pad.b) / 2);
    drawAxes(ctx, w, h, pad, "半径");
    // 光束包络（上下对称）
    ctx.fillStyle = "rgba(230, 57, 70, 0.18)";
    ctx.beginPath();
    for (let i = 0; i < zs.length; i++) ctx.lineTo(x(zs[i]), y(radii[i]));
    for (let i = zs.length - 1; i >= 0; i--) ctx.lineTo(x(zs[i]), y(-radii[i]));
    ctx.closePath(); ctx.fill();
    ctx.strokeStyle = "#e63946"; ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < zs.length; i++) { const px = x(zs[i]), py = y(radii[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
    ctx.stroke();
    ctx.beginPath();
    for (let i = 0; i < zs.length; i++) { const px = x(zs[i]), py = y(-radii[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
    ctx.stroke();
    // 标注束腰
    const wz = d.waist_position_m;
    ctx.strokeStyle = "#2a9d8f"; ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(x(wz), y(d.waist_radius_um)); ctx.lineTo(x(wz), y(-d.waist_radius_um)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#2a9d8f"; ctx.textAlign = "center";
    ctx.fillText("w₀ = " + fmt(d.waist_radius_um) + " μm", x(wz), h - pad.b + 16);
    // 镜面位置标注
    ctx.fillStyle = "#94a3b8";
    ctx.fillText("M1", x(z0), h - pad.b + 16);
    ctx.fillText("M2", x(z1), h - pad.b + 16);
  }

  // ---------- 弛豫振荡 ----------
  async function renderRelax() {
    const params = {
      pump_ratio: +document.getElementById("relax-ratio").value,
      t_end: +document.getElementById("relax-tend").value / 1000,
    };
    const d = await get("/relax", params);
    const m = document.getElementById("relax-metrics");
    m.innerHTML = `
      <div class="metric"><div class="k">阈值泵浦速率</div><div class="v">${fmt(d.threshold_pump_rate)}</div></div>
      <div class="metric"><div class="k">稳态光子密度 φ_ss</div><div class="v">${fmt(d.steady_phi)}</div></div>
      <div class="metric"><div class="k">弛豫振荡频率</div><div class="v">${fmt(d.relaxation_freq_hz)} <span class="unit">Hz</span></div></div>
      <div class="metric"><div class="k">峰值光子密度</div><div class="v">${fmt(d.peak_phi)}</div></div>`;
    drawRelaxCanvas(d);
  }

  function drawRelaxCanvas(d) {
    const cv = document.getElementById("relax-canvas");
    const ctx = clearCanvas(cv);
    const w = cv.width, h = cv.height, pad = { l: 70, r: 20, t: 24, b: 40 };
    const ts = d.t_us, phis = d.phi;
    const t0 = ts[0], t1 = ts[ts.length - 1];
    const pmax = Math.max(...phis) * 1.1 || 1;
    const x = t => pad.l + (t - t0) / (t1 - t0) * (w - pad.l - pad.r);
    const y = v => h - pad.b - (v / pmax) * (h - pad.t - pad.b);
    drawAxes(ctx, w, h, pad, "光子密度 φ");
    ctx.fillStyle = "#64748b"; ctx.textAlign = "right";
    ctx.fillText(fmt(pmax), pad.l - 6, y(pmax) + 4);
    ctx.fillText("0", pad.l - 6, y(0) + 4);
    ctx.textAlign = "center";
    for (const tv of niceTicks(t0, t1, 5)) ctx.fillText(tv.toFixed(0) + "μs", x(tv), h - pad.b + 16);
    // 光子密度曲线
    ctx.strokeStyle = "#e63946"; ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < ts.length; i++) { const px = x(ts[i]), py = y(phis[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
    ctx.stroke();
    // 稳态线
    ctx.strokeStyle = "#2a9d8f"; ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(pad.l, y(d.steady_phi)); ctx.lineTo(w - pad.r, y(d.steady_phi)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#2a9d8f"; ctx.textAlign = "left";
    ctx.fillText("稳态 φ_ss", pad.l + 8, y(d.steady_phi) - 6);
  }

  // ---------- 纵模谱 ----------
  async function renderModes() {
    const params = {
      medium: document.getElementById("modes-medium").value,
      length: +document.getElementById("modes-length").value,
      g0: +document.getElementById("modes-g0").value,
      nmodes: 31,
    };
    const d = await get("/modes", params);
    const m = document.getElementById("modes-metrics");
    const nm = d.modes.length;
    m.innerHTML = `
      <div class="metric"><div class="k">介质</div><div class="v">${d.medium}</div></div>
      <div class="metric"><div class="k">增益中心频率</div><div class="v">${fmt(d.center_thz)} <span class="unit">THz</span></div></div>
      <div class="metric"><div class="k">纵模间距 FSR</div><div class="v">${fmt(d.fsr_thz)} <span class="unit">THz</span></div></div>
      <div class="metric"><div class="k">纵模数</div><div class="v">${nm}</div></div>`;
    drawModesCanvas(d);
  }

  function drawModesCanvas(d) {
    const cv = document.getElementById("modes-canvas");
    const ctx = clearCanvas(cv);
    const w = cv.width, h = cv.height, pad = { l: 60, r: 20, t: 24, b: 40 };
    const freqs = d.spectrum.freq_thz, gains = d.spectrum.gain;
    const f0 = freqs[0], f1 = freqs[freqs.length - 1];
    const gmax = Math.max(...gains) * 1.15 || 1;
    const x = f => pad.l + (f - f0) / (f1 - f0) * (w - pad.l - pad.r);
    const y = g => h - pad.b - (g / gmax) * (h - pad.t - pad.b);
    drawAxes(ctx, w, h, pad, "增益 g(ν)");
    // 连续增益谱（洛伦兹线型）
    ctx.strokeStyle = "#94a3b8"; ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < freqs.length; i++) { const px = x(freqs[i]), py = y(gains[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
    ctx.stroke();
    // 纵模竖线
    ctx.strokeStyle = "#e63946"; ctx.lineWidth = 2;
    for (const mode of d.modes) {
      const px = x(mode.freq_thz);
      ctx.beginPath(); ctx.moveTo(px, y(0)); ctx.lineTo(px, y(mode.gain)); ctx.stroke();
    }
    // FSR 标注
    ctx.fillStyle = "#2a9d8f"; ctx.textAlign = "center";
    const c = d.center_thz;
    ctx.fillText("FSR = " + fmt(d.fsr_thz) + " THz", x(c), pad.t + 12);
    ctx.fillStyle = "#94a3b8";
    ctx.fillText("增益线型（灰）与纵模（红）", (pad.l + w - pad.r) / 2, pad.t + 28);
  }

  // ---------- 绑定事件 ----------
  document.querySelectorAll(".tab").forEach(btn => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".panel").forEach(p => p.classList.remove("active"));
      btn.classList.add("active");
      document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
    });
  });

  document.getElementById("cav-run").addEventListener("click", renderCavity);
  document.getElementById("beam-run").addEventListener("click", renderBeam);
  document.getElementById("relax-run").addEventListener("click", renderRelax);
  document.getElementById("modes-run").addEventListener("click", renderModes);

  // 初始渲染
  renderCavity();
  renderBeam();
  renderRelax();
  renderModes();
})();
