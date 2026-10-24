"""光学谐振腔分析.

覆盖内容：

* 法布里-珀罗（Fabry-Pérot）标准具的透射率/反射率/精细度
* 自由光谱范围（FSR）、纵模间距、纵模频率梳
* 谐振腔稳定性判据（g 因子、稳定/临界/非稳定腔）
* 往返增益与损耗、激光出光阈值条件

物理背景
--------
两片反射率 R1、R2 的平行平面镜构成法布里-珀罗腔。腔内单程相位
``δ = 2 π n L / λ``，标准具透射率满足 Airy 公式：

    T = T_max / (1 + F sin²(δ))

其中 T_max = T1 T2 / (1 - √(R1 R2))² 为峰值透射率，精细度系数
``F = 4 √(R1 R2) / (1 - √(R1 R2))²``。精细度（finesse）定义为
FSR 与透射峰半高全宽（FWHM）之比，描述腔的谐振锐度。
"""

from math import pi, sqrt, sin, cos

from . import constants as C


class FabryPerot:
    """法布里-珀罗标准具.

    Parameters
    ----------
    n : float
        腔内介质折射率。
    length : float
        腔长（m）。
    r1 : float
        左端镜面强度反射率（0~1）。
    r2 : float
        右端镜面强度反射率（0~1）。
    loss : float, optional
        单程强度损耗系数（吸收/散射，0~1）。
    """

    def __init__(self, n, length, r1, r2, loss=0.0):
        if not (0.0 <= r1 < 1.0 and 0.0 <= r2 < 1.0):
            raise ValueError("镜面反射率必须落在 [0, 1) 区间")
        self.n = n
        self.length = length
        self.r1 = r1
        self.r2 = r2
        self.loss = loss
        # 有效场反射率幅度（含单程损耗）
        self._sqrt_r = sqrt(r1 * r2) * (1.0 - loss)

    # ------------------------------------------------------------------
    # 频率/波长特征
    # ------------------------------------------------------------------
    @property
    def optical_length(self):
        """光程长度 n·L (m)."""
        return self.n * self.length

    @property
    def fsr_hz(self):
        """自由光谱范围 (Hz)：c / (2 n L)."""
        return C.C / (2.0 * self.optical_length)

    @property
    def fsr_thz(self):
        return self.fsr_hz / 1e12

    def mode_frequencies(self, center_thz, order, nmodes):
        """返回围绕中心频率的纵模频率梳 (THz).

        纵模满足 ``q·c/(2nL) = ν_q``，相邻纵模间距即 FSR。
        """
        f0 = center_thz - order * self.fsr_thz
        return [f0 + k * self.fsr_thz for k in range(nmodes)]

    # ------------------------------------------------------------------
    # 精细度与线宽
    # ------------------------------------------------------------------
    @property
    def finesse_coeff(self):
        """Airy 精细度系数 F = 4 √(R1 R2) / (1 - √(R1 R2))²."""
        sq = self._sqrt_r
        return 4.0 * sq / (1.0 - sq) ** 2

    @property
    def finesse(self):
        """精细度（finesse）."""
        sq = self._sqrt_r
        return pi * sqrt(sq) / (1.0 - sq)

    def linewidth_hz(self):
        """谐振峰半高全宽 (Hz) = FSR / finesse."""
        return self.fsr_hz / self.finesse

    # ------------------------------------------------------------------
    # 透射率 / 反射率
    # ------------------------------------------------------------------
    def transmittance(self, freq_thz):
        """给定真空频率 (THz) 的标准具强度透射率（含损耗）. """
        delta = 2.0 * pi * self.optical_length * freq_thz * 1e12 / C.C
        return 1.0 / (1.0 + self.finesse_coeff * sin(delta) ** 2)

    def reflectivity(self, freq_thz):
        """强度反射率 = 1 - 透射率（单程损耗假设已并入）."""
        return 1.0 - self.transmittance(freq_thz)

    def transmittance_vs_wavelength(self, wl_nm):
        """给定真空波长 (nm) 的透射率."""
        return self.transmittance(C.wl_to_freq(wl_nm))

    def peak_transmittance(self):
        """峰值透射率 T_max."""
        sq = self._sqrt_r
        return (1.0 - self.r1 - self.r2 + self.r1 * self.r2) * (1.0 - self.loss) ** 2 / (1.0 - sq) ** 2

    # ------------------------------------------------------------------
    # 往返损耗与阈值
    # ------------------------------------------------------------------
    @property
    def roundtrip_loss(self):
        """往返强度损耗（不含镜面透射）."""
        return 1.0 - (1.0 - self.loss) ** 2

    def roundtrip_gain_required(self):
        """达到阈值所需的往返增益系数（>1 出光）.

        阈值条件：R1 R2 (1-loss)² · g_rt ≥ 1，即往返净增益补偿所有损耗。
        """
        return 1.0 / (self.r1 * self.r2 * (1.0 - self.loss) ** 2)


class StableCavity:
    """两镜稳定腔的 g 因子稳定性判据.

    对曲率半径 R1、R2（凹面为正）的腔，定义

        g1 = 1 - L/R1,  g2 = 1 - L/R2

    稳定条件为 ``0 ≤ g1·g2 ≤ 1``。g1·g2 == 0 或 == 1 为临界腔。
    """

    def __init__(self, length, radius1, radius2):
        """radius 取 +∞ 表示平面镜（R=∞）。传入 float('inf') 即可."""
        self.length = length
        self.radius1 = radius1
        self.radius2 = radius2
        self.g1 = 1.0 - length / radius1
        self.g2 = 1.0 - length / radius2

    @property
    def g_product(self):
        return self.g1 * self.g2

    @property
    def is_stable(self):
        """是否严格稳定腔."""
        return 0.0 < self.g_product < 1.0

    @property
    def is_critical(self):
        """是否临界腔（边界）."""
        gp = self.g_product
        return abs(gp) < 1e-12 or abs(gp - 1.0) < 1e-12

    def classify(self):
        """返回字符串分类：stable / critical / unstable."""
        if self.is_stable:
            return "stable"
        if self.is_critical:
            return "critical"
        return "unstable"


def finesse_from_reflectance(r1, r2, loss=0.0):
    """由两镜反射率直接求精细度（便捷函数）."""
    fp = FabryPerot(1.0, 1.0, r1, r2, loss)
    return fp.finesse
