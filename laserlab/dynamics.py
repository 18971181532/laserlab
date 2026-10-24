"""激光动力学：弛豫振荡、稳态输出与多模竞争.

弛豫振荡（relaxation oscillation）是激光器在泵浦或腔内损耗突变后
光子数与反转粒子数的阻尼振荡，是 Q 开关、调制响应的基础。其
（角）频率在小信号近似下约为

    Ω_R ≈ sqrt( (r - 1) / (τ τc) )

其中 r = Rp/Rp_th 为泵浦超阈值比，τ 为上能级寿命，τc 为光子腔寿命。
"""

from math import pi, sqrt

from . import rate as _rate


class RelaxationOscillation:
    """四能级单模激光的时域仿真结果封装."""

    def __init__(self, t, N, phi):
        self.t = t
        self.N = N
        self.phi = phi

    def steady_state_power(self, photon_energy_j):
        """稳态输出功率估计 (W)：φ_ss · V_cavity · hν / τc（相对量）.

        这里返回腔内光子数密度的稳态值（用户可自行乘上几何/能量因子）。
        """
        return self.phi[-1]


def simulate_relaxation(pump_rate, tau_upper=2.3e-4, sigma=2.8e-23,
                        photon_lifetime=1.0e-8, beta=1e-6, t_end=1e-3,
                        c=3e8, adaptive=True, n0=0.0, phi0=0.0):
    """仿真四能级激光的弛豫振荡.

    返回 RelaxationOscillation 对象，含时间轴、反转密度 N、光子密度 φ。
    初始反转 N0=0 表示从冷启动泵浦，观察起振与弛豫振荡全过程。
    """
    f = _rate.four_level_rates(pump_rate, tau_upper, sigma, photon_lifetime, beta, c=c)
    if adaptive:
        t, y = _rate.solve_rk45(f, 0.0, [n0, phi0], t_end,
                                h0=1e-10, rtol=1e-6, atol=1e-9, hmax=t_end / 500.0)
    else:
        t, y = _rate.solve_rk4(f, 0.0, [n0, phi0], t_end, t_end / 20000.0)
    N = [v[0] for v in y]
    phi = [v[1] for v in y]
    return RelaxationOscillation(t, N, phi)


def relaxation_frequency(pump_ratio, tau_upper, photon_lifetime):
    """小信号弛豫振荡频率 (Hz).

    pump_ratio r = Rp/Rp_th（>1）。Ω_R/(2π) = sqrt((r-1)/(τ·τc))/(2π)。
    """
    if pump_ratio <= 1.0:
        return 0.0
    omega = sqrt((pump_ratio - 1.0) / (tau_upper * photon_lifetime))
    return omega / (2.0 * pi)


def steady_state(pump_rate, tau_upper, sigma, photon_lifetime, beta=0.0, c=3e8):
    """四能级激光的稳态解析解（忽略自发项或给定 β）.

    稳态：dN/dt=0, dφ/dt=0。忽略 β 时：

        N_ss = 1/(c σ τc)
        φ_ss = τc (Rp - N_ss/τ)

    返回 (N_ss, φ_ss)。若泵浦低于阈值，φ_ss 为 0。
    """
    n_ss = 1.0 / (c * sigma * photon_lifetime)
    phi_ss = photon_lifetime * (pump_rate - n_ss / tau_upper)
    if phi_ss < 0.0:
        phi_ss = 0.0
    return n_ss, phi_ss


class MultiMode:
    """多纵模竞争简化模型.

    均匀加宽增益下，各纵模共享同一反转粒子数，增益被总强度饱和，
    因此最终往往只剩中心模（模式竞争/烧孔）。本模型离散化增益谱，
    对每个纵模独立积分光子数，反转由所有模的总受激辐射共同消耗。
    """

    def __init__(self, gain_medium, cavity, pump_ratio, tau_upper=2.3e-4,
                 sigma=2.8e-23, beta=1e-8, c=3e8):
        self.gain = gain_medium
        self.cavity = cavity
        self.pump_ratio = pump_ratio
        self.tau_upper = tau_upper
        self.sigma = sigma
        self.beta = beta
        self.c = c
        # 纵模频率：增益中心附近取若干 FSR 间隔
        self.modes = cavity.mode_frequencies(gain_medium.center_thz, 0, 0)  # placeholder

    def mode_frequencies(self, n_each_side):
        """返回中心 ±n_each_side 个 FSR 范围内的纵模频率 (THz)."""
        fsr = self.cavity.fsr_thz
        center = self.gain.center_thz
        return [center + k * fsr for k in range(-n_each_side, n_each_side + 1)]

    def simulate(self, n_each_side=5, t_end=1e-3, steps=40000):
        """积分多模光子数，返回 (t, N, modes_phi) 其中 modes_phi 为 dict.

        反转 N 被所有模的总受激辐射消耗，导致增益竞争（烧孔）。
        泵浦项取使中心模小信号增益恰能补偿腔损耗的 1.5 倍（超阈值），
        使系统起振并观察到模式竞争稳态。
        """
        freqs = self.mode_frequencies(n_each_side)
        tau_c = 5e-8  # 光子腔寿命
        f = self._make_rhs(freqs, tau_c)
        t, y = _rate.solve_rk45(f, 0.0, [0.0] + [0.0] * len(freqs), t_end,
                                h0=1e-10, rtol=1e-5, atol=1e-9, hmax=t_end / 500.0)
        N = [v[0] for v in y]
        modes_phi = {freqs[i]: [v[i + 1] for v in y] for i in range(len(freqs))}
        return t, N, modes_phi

    def _make_rhs(self, freqs, tau_c):
        gain = self.gain
        c, sigma, tau, beta = self.c, self.sigma, self.tau_upper, self.beta
        # 泵浦项：使中心模超阈值 1.5 倍
        g_center = gain.small_signal_gain(gain.center_thz)
        pump = 1.5 / (c * sigma * g_center * tau_c * tau)

        def f(t, y):
            N = y[0]
            phis = y[1:]
            dN = pump - N / tau
            dphis = []
            for i, nu in enumerate(freqs):
                g = gain.small_signal_gain(nu)
                stim = c * sigma * g * N * phis[i]
                dphis.append(stim - phis[i] / tau_c + beta * N / tau)
                dN -= stim
            return [dN] + dphis

        return f
