"""物理常量.

数值取自 CODATA 2018 推荐值，精度足以满足本工作台的工程仿真需求。
"""

from math import pi

# 光速 (m/s)
C = 299_792_458.0

# 普朗克常数 (J·s)
H = 6.62607015e-34

# 玻尔兹曼常数 (J/K)
KB = 1.380649e-23

# 电子电荷 (C)
E = 1.602176634e-19

# 真空介电常数 (F/m)
EPS0 = 8.8541878128e-12

# 约化普朗克常数 (J·s)
HBAR = H / (2.0 * pi)


def freq_to_wl(freq_thz):
    """真空频率 (THz) -> 真空波长 (nm)."""
    return C / (freq_thz * 1e12) * 1e9


def wl_to_freq(wl_nm):
    """真空波长 (nm) -> 真空频率 (THz)."""
    return C / (wl_nm * 1e-9) / 1e12
