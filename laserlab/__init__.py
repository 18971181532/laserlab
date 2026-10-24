"""LaserLab — 激光器物理与谐振腔仿真工作台.

一个零第三方依赖的纯 Python 库，用于：

* 光学谐振腔分析（自由光谱范围 FSR、纵模间距、g 因子稳定性判据）
* 高斯光束传播（ABCD 矩阵、q 参数、束腰/瑞利长度演化）
* 增益介质与速率方程（三能级/四能级、增益饱和、泵浦阈值）
* 激光动力学（弛豫振荡、单模/多模竞争、小信号增益）

物理约定
--------
* 时间谐波约定为 ``exp(-i ω t)``。
* 波长 λ 使用真空波长（单位 nm），频率用 THz。
* 所有角度单位均为弧度，长度单位均为米（m），除非注释另行说明。
"""

__version__ = "1.0.0"
__author__ = "LaserLab contributors"

from . import constants  # noqa: F401
from . import cavity    # noqa: F401
from . import beam      # noqa: F401
from . import gain      # noqa: F401
from . import rate      # noqa: F401
from . import dynamics  # noqa: F401
