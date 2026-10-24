"""激光速率方程与数值积分器.

速率方程是激光动力学的核心。对四能级系统，反转粒子数密度 N 与
腔内光子数密度 φ 的演化满足（单模近似）：

    dN/dt = Rp - N/τ - c·σ·N·φ
    dφ/dt = c·σ·N·φ - φ/τc + β·N/τ

其中 Rp 为泵浦速率，τ 为上能级寿命，σ 为受激发射截面，τc 为光子
腔寿命，β 为自发辐射进入激光模的比例因子。

本模块实现：

* RK4 定步长积分器（通用，适合刚性问题需配合小步长）
* 自适应步长 RK45（Dormand–Prince，误差控制，适合弛豫振荡这类
  快慢时间尺度分离的刚性问题）
* 常见右端项的工厂函数
"""


def rk4_step(f, t, y, h):
    """经典四阶 Runge–Kutta 单步.

    Parameters
    ----------
    f : callable(t, y) -> array
        右端函数。
    t, y : float, list/array
        当前时间与状态。
    h : float
        步长。

    Returns
    -------
    new_y : list
        下一步状态。
    """
    k1 = f(t, y)
    k2 = f(t + h / 2, [y[i] + h / 2 * k1[i] for i in range(len(y))])
    k3 = f(t + h / 2, [y[i] + h / 2 * k2[i] for i in range(len(y))])
    k4 = f(t + h, [y[i] + h * k3[i] for i in range(len(y))])
    return [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(len(y))]


def solve_rk4(f, t0, y0, t1, dt):
    """定步长 RK4 积分，返回 (t_list, y_list)."""
    n = int(round((t1 - t0) / dt))
    ts = [t0]
    ys = [list(y0)]
    t, y = t0, list(y0)
    for _ in range(n):
        y = rk4_step(f, t, y, dt)
        t += dt
        ts.append(t)
        ys.append(list(y))
    return ts, ys


# Dormand–Prince 5(4) 系数（经典 RK45 表）
_A = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0]
_B = [
    [0.2],
    [3.0 / 40.0, 9.0 / 40.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
]
_CV = [0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0]
_C5 = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0]
_C4 = [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0,
       -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0]


def rk45_step(f, t, y, h):
    """Dormand–Prince 5(4) 单步，返回 (y5, y4, err)."""
    dim = len(y)
    k = [[0.0] * dim for _ in range(7)]
    k[0] = f(t, y)
    for i in range(1, 7):
        yi = [y[j] + h * sum(_B[i - 1][m] * k[m][j] for m in range(i)) for j in range(dim)]
        k[i] = f(t + _CV[i] * h, yi)
    y5 = [y[j] + h * sum(_C5[i] * k[i][j] for i in range(7)) for j in range(dim)]
    y4 = [y[j] + h * sum(_C4[i] * k[i][j] for i in range(7)) for j in range(dim)]
    err = max(abs(y5[j] - y4[j]) for j in range(dim))
    return y5, y4, err


def solve_rk45(f, t0, y0, t1, h0=1e-9, rtol=1e-6, atol=1e-12, hmin=1e-14, hmax=None):
    """自适应步长 RK45 积分，返回 (t_list, y_list).

    误差控制：|err| ≤ atol + rtol·max(|y|)。步长按 5 阶与 4 阶解之差
    自适应缩放，适合弛豫振荡这类快慢时间尺度分离的刚性问题。
    """
    if hmax is None:
        hmax = (t1 - t0) / 10.0
    ts = [t0]
    ys = [list(y0)]
    t, y, h = t0, list(y0), h0
    safety = 0.9
    while t < t1 - 1e-16:
        if h > hmax:
            h = hmax
        if t + h > t1:
            h = t1 - t
        try:
            y5, y4, err = rk45_step(f, t, y, h)
        except OverflowError:
            h *= 0.5
            continue
        scale = atol + rtol * max(max(abs(v) for v in y), 1.0)
        if err <= scale or h <= hmin:
            y = y5
            t += h
            ts.append(t)
            ys.append(list(y))
        if err == 0.0:
            h *= 2.0
        else:
            h *= safety * (scale / err) ** 0.2
        if h < hmin:
            h = hmin
    return ts, ys


# ----------------------------------------------------------------------
# 四能级速率方程工厂
# ----------------------------------------------------------------------
def four_level_rates(pump_rate, tau_upper, sigma, photon_lifetime, beta,
                     c=3e8, roundtrip_factor=1.0):
    """构造四能级单模速率方程右端函数.

    Parameters
    ----------
    pump_rate : float
        泵浦速率 Rp (1/(m³·s))。
    tau_upper : float
        上能级寿命 τ (s)。
    sigma : float
        受激发射截面 σ (m²)。
    photon_lifetime : float
        光子腔寿命 τc (s)。
    beta : float
        自发辐射进入激光模的比例因子。
    c : float
        光速（默认真空）。
    roundtrip_factor : float
        往返净增益缩放因子（默认 1）。

    Returns
    -------
    f(t, [N, phi]) -> [dN/dt, dφ/dt]
    """

    def f(t, y):
        N, phi = y
        stim = c * sigma * N * phi * roundtrip_factor
        dN = pump_rate - N / tau_upper - stim
        dphi = stim - phi / photon_lifetime + beta * N / tau_upper
        return [dN, dphi]

    return f


def threshold_pump_rate(tau_upper, sigma, photon_lifetime, c=3e8):
    """四能级激光的阈值泵浦速率.

    稳态、无自发项时 dφ/dt=0 给出 c·σ·N_th = 1/τc，配合 dN/dt=0
    （忽略受激项）得 Rp_th = N_th/τ = 1/(c·σ·τ·τc)。
    """
    return 1.0 / (c * sigma * tau_upper * photon_lifetime)
