# LaserLab — 激光器物理与谐振腔仿真工作台

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![CI](https://github.com/18971181532/laserlab/actions/workflows/tests.yml/badge.svg)](https://github.com/18971181532/laserlab/actions/workflows/tests.yml)

一个**零第三方依赖**、纯标准库实现的激光器物理仿真库，覆盖从谐振腔到增益介质、从速率方程到激光动力学的完整物理链路。既可作为光电/物理专业的学习工具，也可用于快速估算实际激光系统（如 Nd:YAG、He-Ne、Ti:Sapphire）的关键参数。

## 项目简介

激光器是光电专业最具代表性的器件之一。理解其工作原理需要同时掌握多个物理层面的知识：

1. **光学谐振腔** —— 决定纵模频率、线宽与稳定性；
2. **高斯光束** —— 决定腔内模式的空间分布与束腰；
3. **增益介质** —— 决定可增益的频谱范围与饱和行为；
4. **速率方程** —— 决定激光从起振到稳态的动态过程。

LaserLab 用清晰的模块划分与手写的数值方法，将这些物理过程逐一实现，并提供一个交互式 Web 界面直观展示仿真结果。

### 特色

- **纯标准库**：不依赖 NumPy/SciPy，物理与数值算法全部手写，便于理解底层原理；
- **解析解交叉验证**：法布里-珀罗透射率、精细度、高斯光束束腰、稳态速率方程均有解析基准，测试直接比对；
- **自适应步长积分器**：手写 Dormand–Prince RK45，正确求解弛豫振荡这类快慢时间尺度分离的刚性问题；
- **手写 Canvas 可视化**：Web 界面零 CDN，前端绘图逻辑完全自包含。

## 核心功能

| 模块 | 功能 |
|------|------|
| `cavity` | 法布里-珀罗标准具（透射率/精细度/线宽/FSR）、g 因子稳定性判据、纵模频率梳 |
| `beam` | 高斯光束 q 参数、ABCD 矩阵光学、稳定腔自再现解（Kogelnik & Li）、模匹配耦合效率 |
| `gain` | 均匀加宽（洛伦兹）/非均匀加宽（高斯）线型、增益饱和、常用介质参数库 |
| `rate` | 四能级速率方程、RK4 与自适应 RK45 积分器、阈值泵浦速率 |
| `dynamics` | 弛豫振荡仿真、稳态解析解、多纵模竞争（烧孔效应） |

## 环境与依赖

- **Python 3.9+**（无任何第三方依赖，仅用标准库 `math`/`cmath`/`json`/`http` 等）
- 可选：现代浏览器（用于 Web 界面）

## 安装与运行

```bash
git clone https://github.com/18971181532/laserlab.git
cd laserlab
```

无需 `pip install`，直接运行：

### 1. 命令行

```bash
# 法布里-珀罗腔分析
python run.py cavity --length 0.2 --r1 0.99 --r2 0.98

# 稳定腔高斯光束
python run.py beam --length 0.2 --radius1 inf --radius2 0.5 --wavelength 1064

# 弛豫振荡仿真
python run.py relax --pump-ratio 2.0 --t-end 0.001

# 纵模频率与增益
python run.py modes --medium nd_yag_1064 --length 0.2 --g0 10
```

### 2. Web 界面

```bash
python run.py serve --port 8779
```

浏览器打开 http://127.0.0.1:8779 ，四个标签页对应谐振腔、高斯光束、弛豫振荡、纵模谱。

### 3. 演示脚本

```bash
python run_demos.py
```

### 4. 运行测试

```bash
python run_tests.py
```

## 使用示例

### 法布里-珀罗腔

```python
from laserlab import cavity

fp = cavity.FabryPerot(n=1.0, length=0.2, r1=0.99, r2=0.98)
print(f"FSR = {fp.fsr_thz:.4f} THz")
print(f"精细度 = {fp.finesse:.2f}")
print(f"线宽 = {fp.linewidth_hz():.3e} Hz")
```

### 稳定腔高斯光束

```python
from laserlab import beam

gc = beam.GaussianCavity(length=0.2, radius1=float("inf"), radius2=0.5,
                         wavelength_nm=1064.0)
print(f"束腰半径 = {gc.waist_radius * 1e6:.2f} μm")
print(f"束腰位置 = {gc.waist_position:.4f} m")
print(f"瑞利长度 = {gc.rayleigh:.4f} m")
```

### 弛豫振荡

```python
from laserlab import rate, dynamics

tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1.0e-8   # Nd:YAG 典型参数
rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
sim = dynamics.simulate_relaxation(2.0 * rp_th, tau_upper=tau, sigma=sigma,
                                   photon_lifetime=tau_c, t_end=1e-3)
print(f"弛豫振荡频率 = {dynamics.relaxation_frequency(2.0, tau, tau_c):.3e} Hz")
print(f"稳态光子密度 = {sim.phi[-1]:.3e}")
```

### 多纵模竞争

```python
from laserlab import cavity, gain, dynamics

fp = cavity.FabryPerot(1.0, 0.3, 0.99, 0.95)
med = gain.make_medium("nd_yag_1064", g0=50.0)
mm = dynamics.MultiMode(med, fp, pump_ratio=2.0)
t, N, modes = mm.simulate(n_each_side=3, t_end=1e-3)
```

## 物理背景与实现细节

### 法布里-珀罗标准具

两片反射率 R₁、R₂ 的平行镜构成 FP 腔，Airy 透射率：

```
T(δ) = 1 / (1 + F·sin²δ)，  δ = 2π·n·L·ν/c
```

精细度系数 `F = 4√(R₁R₂)/(1-√(R₁R₂))²`，精细度 `𝓕 = π√(R₁R₂)/(1-√(R₁R₂))`，
谐振峰线宽 `Δν = FSR/𝓕`。FSR = `c/(2nL)`。

### 稳定腔判据

对曲率半径 R₁、R₂ 的两镜腔，定义 `g₁ = 1 - L/R₁`、`g₂ = 1 - L/R₂`，
稳定条件为 `0 < g₁g₂ < 1`。束腰位置与尺寸由 Kogelnik & Li 闭式解给出：

```
z₁ = L·g₂(1-g₁) / (g₁+g₂-2g₁g₂)
w₀² = (λL/π)·√(g₁g₂(1-g₁g₂)) / (g₁+g₂-2g₁g₂)
```

### 速率方程与弛豫振荡

四能级单模速率方程（见 `rate.py`），弛豫振荡角频率：

```
Ω_R = √((r-1)/(τ·τc))，  r = Rp/Rp_th 为超阈值比
```

弛豫振荡是快慢时间尺度分离的刚性问题，故实现 Dormand–Prince 5(4)
自适应步长积分器（RK45），步长按 5 阶与 4 阶解之差自适应缩放。

## 目录结构

```
laserlab/
├── laserlab/               # 核心库
│   ├── __init__.py         # 包入口
│   ├── constants.py        # 物理常量（CODATA 2018）
│   ├── cavity.py           # 法布里-珀罗腔、稳定腔判据
│   ├── beam.py             # 高斯光束、ABCD 矩阵、自再现解
│   ├── gain.py             # 增益介质、线型、饱和、介质库
│   ├── rate.py             # 速率方程、RK4/RK45 积分器
│   ├── dynamics.py         # 弛豫振荡、稳态、多模竞争
│   ├── cli.py              # 命令行接口
│   └── server.py           # HTTP 服务与 JSON API
├── web/                    # 前端（手写 Canvas，零 CDN）
│   ├── index.html
│   ├── style.css
│   └── app.js
├── tests/
│   └── test_laserlab.py    # 40 个单元测试
├── examples/               # 示例配置
├── run.py                  # 引导脚本（CLI 入口）
├── run_demos.py            # 6 个演示
├── run_tests.py            # 测试入口
├── README.md
├── LICENSE
├── pyproject.toml
└── .github/workflows/tests.yml
```

## 局限性

- 速率方程采用单模/简化多模近似，未包含横模空间烧孔、高阶模等精细效应；
- 增益介质参数为典型值，实际器件需按具体掺杂浓度/温度标定；
- 未包含量子噪声与相位噪声模型（Schawlow-Townes 线宽）——可作为后续扩展。

## License

MIT License，详见 [LICENSE](LICENSE)。
