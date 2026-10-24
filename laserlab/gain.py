"""增益介质与增益饱和模型.

覆盖内容：

* 均匀加宽（homogeneous broadening）：洛伦兹线型
* 非均匀加宽（inhomogeneous broadening）：高斯线型
* 小信号增益谱 g(ν) 与饱和增益 g(ν)/(1 + I/Isat)
* 常用激光介质参数（Nd:YAG、Er:YAG、He-Ne、半导体等）

增益饱和是激光动力学与多模竞争的核心：当腔内强度接近饱和强度
I_sat 时，增益被"烧孔"（均匀加宽）或在谱上打洞（非均匀加宽）。
"""

from math import pi, exp, log, sqrt

from . import constants as C


class GainMedium:
    """均匀加宽增益介质.

    Parameters
    ----------
    g0 : float
        小信号增益系数（每单位长度，1/m）。
    center_thz : float
        增益中心频率 (THz)。
    linewidth_thz : float
        增益谱半高全宽 FWHM (THz)。
    isat : float
        饱和强度 (W/m²)。
    """

    def __init__(self, g0, center_thz, linewidth_thz, isat):
        self.g0 = g0
        self.center_thz = center_thz
        self.linewidth_thz = linewidth_thz
        self.isat = isat

    def small_signal_gain(self, freq_thz):
        """洛伦兹线型小信号增益 g(ν)."""
        dv = freq_thz - self.center_thz
        hw = self.linewidth_thz / 2.0
        return self.g0 * (hw ** 2) / (hw ** 2 + dv ** 2)

    def saturated_gain(self, freq_thz, intensity):
        """饱和增益 g(ν)/(1 + I/I_sat)."""
        return self.small_signal_gain(freq_thz) / (1.0 + intensity / self.isat)


class InhomogeneousMedium:
    """非均匀加宽增益介质（高斯线型，如气体多普勒加宽）."""

    def __init__(self, g0, center_thz, linewidth_thz, isat):
        self.g0 = g0
        self.center_thz = center_thz
        self.linewidth_thz = linewidth_thz
        self.isat = isat

    def small_signal_gain(self, freq_thz):
        dv = freq_thz - self.center_thz
        sigma = self.linewidth_thz / (2.0 * sqrt(2.0 * log(2.0)))
        return self.g0 * exp(-dv ** 2 / (2.0 * sigma ** 2))


# ----------------------------------------------------------------------
# 常用介质参数（中心波长 nm, 线宽, 饱和强度）
# ----------------------------------------------------------------------
MEDIA = {
    # 名称: (中心波长 nm, 线宽 nm(近似), 饱和强度 W/m²)
    "nd_yag_1064": (1064.0, 0.6, 2.9e7),
    "nd_yag_1319": (1319.0, 0.8, 3.0e7),
    "er_glass_1540": (1540.0, 30.0, 4.0e7),
    "hene_633": (632.8, 0.0015, 1.0e3),
    "ti_sapphire_800": (800.0, 125.0, 1.0e10),
    "co2_10600": (10600.0, 1.0, 1.0e6),
}


def make_medium(name, g0):
    """按名称构造 GainMedium，g0 为小信号增益系数 (1/m)."""
    if name not in MEDIA:
        raise KeyError(f"未知介质 {name!r}，可选: {sorted(MEDIA)}")
    wl_nm, linewidth_nm, isat = MEDIA[name]
    center_thz = C.wl_to_freq(wl_nm)
    linewidth_thz = abs(C.wl_to_freq(wl_nm) - C.wl_to_freq(wl_nm + linewidth_nm))
    return GainMedium(g0, center_thz, linewidth_thz, isat)
