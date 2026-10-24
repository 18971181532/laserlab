"""高斯光束传播与 ABCD 矩阵光学.

高斯光束的复 q 参数完全描述其传播行为：

    1/q = 1/R - i·λ/(π w²)

其中 R 为等相面曲率半径，w 为光束半径（1/e² 强度）。经过一个由
2×2 ABCD 矩阵描述的光学元件后，q 参数按双线性变换更新：

    q' = (A q + B) / (C q + D)

本模块实现：

* 常见元件的 ABCD 矩阵（自由空间、薄透镜、平面/球面镜、介质界面）
* q 参数随 z 的传播、束腰与瑞利长度
* 稳定腔内的自再现解（self-consistent mode）：往返 ABCD 矩阵求高斯模
* 模匹配（两个高斯光束之间的耦合效率）
"""

from cmath import sqrt as csqrt
from math import pi, sqrt, atan

from . import constants as C


def q_parameter(w0, z=0.0, wavelength_nm=1064.0):
    """由束腰半径 w0 与位置 z（相对束腰）构造 q 参数 (m).

    返回复 q 参数。z=0 处 R=∞，q 为纯虚数 i·zR。
    """
    lam = wavelength_nm * 1e-9
    zr = pi * w0 ** 2 / lam
    return z + 1j * zr


def beam_radius(q, wavelength_nm=1064.0):
    """从 q 参数求 1/e² 光束半径 w (m).

    w² = -λ / (π · Im(1/q))。
    """
    lam = wavelength_nm * 1e-9
    return sqrt(-lam / (pi * (1.0 / q).imag))


def curvature_radius(q):
    """从 q 参数求等相面曲率半径 R (m).

    等相面为平面时曲率半径为无穷大，返回 ``float('inf')``。
    """
    re = (1.0 / q).real
    if abs(re) < 1e-18:
        return float("inf")
    return 1.0 / re


def rayleigh_range(w0, wavelength_nm=1064.0):
    """束腰 w0 对应的瑞利长度 zR = π w0² / λ."""
    return pi * w0 ** 2 / (wavelength_nm * 1e-9)


def divergence_half_angle(w0, wavelength_nm=1064.0):
    """远场发散半角 (rad) = λ / (π w0)."""
    return wavelength_nm * 1e-9 / (pi * w0)


# ----------------------------------------------------------------------
# ABCD 矩阵
# ----------------------------------------------------------------------
def free_space(d):
    """自由空间传播距离 d (m)."""
    return [[1.0, d], [0.0, 1.0]]


def thin_lens(f):
    """薄透镜，焦距 f (m，正为会聚)."""
    return [[1.0, 0.0], [-1.0 / f, 1.0]]


def curved_mirror(radius):
    """球面镜，曲率半径 radius (m，凹面正). 等效焦距 f = R/2."""
    f = radius / 2.0
    return [[1.0, 0.0], [-1.0 / f, 1.0]]


def flat_mirror():
    """平面镜."""
    return [[1.0, 0.0], [0.0, 1.0]]


def refraction(n1, n2):
    """介质界面（平面），折射率 n1 -> n2."""
    return [[1.0, 0.0], [0.0, n1 / n2]]


def matmul(a, b):
    """2×2 矩阵乘法 a·b."""
    return [
        [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
        [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
    ]


def compose(*matrices):
    """按顺序合成多个矩阵（从左到右依次作用）."""
    result = [[1.0, 0.0], [0.0, 1.0]]
    for m in matrices:
        result = matmul(m, result)
    return result


def apply_matrix(q, matrix):
    """对 q 参数施加 ABCD 矩阵的双线性变换."""
    a, b = matrix[0][0], matrix[0][1]
    c, d = matrix[1][0], matrix[1][1]
    return (a * q + b) / (c * q + d)


# ----------------------------------------------------------------------
# 稳定腔自再现模
# ----------------------------------------------------------------------
class GaussianCavity:
    """两镜稳定腔内基模高斯光束的自再现解.

    设腔镜 M1、M2 曲率半径 R1、R2，腔长 L。在镜 M1 处，自再现要求
    等相面曲率半径与镜面匹配：R(z) = R1（镜面）。由往返 ABCD 矩阵
    的稳定性条件 0 < (A+D)/2 < 1（本约定下为 -1 < (A+D)/2 < 1 的形式
    视符号约定），可解出镜面处的束腰，再回溯腔内束腰位置与尺寸。

    这里采用直接方法：在 M1 处令 R(z1) 等于镜面曲率半径，配合往返
    矩阵求出自洽的 q。
    """

    def __init__(self, length, radius1, radius2, wavelength_nm=1064.0):
        self.length = length
        self.radius1 = radius1
        self.radius2 = radius2
        self.wavelength_nm = wavelength_nm
        self._solve()

    def _solve(self):
        L, R1, R2 = self.length, self.radius1, self.radius2
        lam = self.wavelength_nm * 1e-9
        g1 = 1.0 - L / R1
        g2 = 1.0 - L / R2
        gp = g1 * g2
        if not (0.0 < gp < 1.0):
            raise ValueError("腔为非稳定或临界腔，不存在高斯自再现模")

        # 标准两镜稳定腔闭式解（Kogelnik & Li 结果）：
        # 束腰距 M1 的距离与束腰半径
        denom = g1 + g2 - 2.0 * g1 * g2
        z1 = L * g2 * (1.0 - g1) / denom
        w0 = sqrt(lam * L / pi * sqrt(gp * (1.0 - gp)) / denom)

        self.waist_position = z1
        self.waist_radius = w0
        self.rayleigh = rayleigh_range(w0, self.wavelength_nm)

        # 镜面处的 q 参数：束腰处 q = i·zR，传播 (z - z1) 后平移
        zr = self.rayleigh
        self.q_at_m1 = (0.0 - z1) + 1j * zr
        self.q_at_m2 = (L - z1) + 1j * zr

    def beam_radius_at(self, z):
        """腔内距 M1 为 z 处的光束半径 (m).

        高斯光束半径 w(z) = w0 · sqrt(1 + ((z - z1)/zR)²)。
        """
        zr = self.rayleigh
        dz = z - self.waist_position
        return self.waist_radius * sqrt(1.0 + (dz / zr) ** 2)

    def curvature_at(self, z):
        """腔内距 M1 为 z 处的等相面曲率半径 (m).

        R(z) = (z - z1) · (1 + (zR/(z-z1))²)。
        """
        dz = z - self.waist_position
        if abs(dz) < 1e-18:
            return float("inf")
        zr = self.rayleigh
        return dz * (1.0 + (zr / dz) ** 2)

    @property
    def is_stable(self):
        return True  # 构造成功即稳定


# ----------------------------------------------------------------------
# 模匹配
# ----------------------------------------------------------------------
def mode_overlap(w0a, w0b, wavelength_nm=1064.0):
    """两个同轴基模高斯光束（束腰 w0a、w0b 重叠于同一位置）的耦合效率.

    返回 0~1 的功率重叠因子：
        η = 4 w0a² w0b² / (w0a² + w0b²)²
    """
    return 4.0 * w0a ** 2 * w0b ** 2 / (w0a ** 2 + w0b ** 2) ** 2


def coupling_efficiency(q1, q2, wavelength_nm=1064.0):
    """两个任意 q 参数高斯光束在同一点的功率耦合效率.

    η = 4 / ( (w1/w2 + w2/w1)² + (π w1 w2 / λ)² (1/R1 - 1/R2)² )
    曲率失配与尺寸失配都会降低耦合。
    """
    w1 = beam_radius(q1, wavelength_nm)
    w2 = beam_radius(q2, wavelength_nm)
    r1 = curvature_radius(q1)
    r2 = curvature_radius(q2)
    inv_r1 = 0.0 if r1 == float("inf") else 1.0 / r1
    inv_r2 = 0.0 if r2 == float("inf") else 1.0 / r2
    lam = wavelength_nm * 1e-9
    size_term = (w1 / w2 + w2 / w1) ** 2
    curv_term = (pi * w1 * w2 / lam) ** 2 * (inv_r1 - inv_r2) ** 2
    return 4.0 / (size_term + curv_term)
