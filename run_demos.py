"""可运行的演示脚本.

每个 demo 独立打印一段结果，展示库的核心能力。运行方式：

    python run_demos.py
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from laserlab import cavity, beam, gain, dynamics, rate  # noqa: E402


def demo_fabry_perot():
    print("=" * 60)
    print("Demo 1: 法布里-珀罗标准具")
    print("=" * 60)
    fp = cavity.FabryPerot(n=1.0, length=0.2, r1=0.99, r2=0.98, loss=0.0)
    print(f"  腔长 L = {fp.length} m, 光程 = {fp.optical_length:.4f} m")
    print(f"  FSR = {fp.fsr_thz:.4f} THz ({fp.fsr_hz:.3e} Hz)")
    print(f"  精细度 = {fp.finesse:.2f}")
    print(f"  谐振峰线宽 = {fp.linewidth_hz():.3e} Hz")
    print(f"  峰值透射率 = {fp.peak_transmittance():.4f}")
    print(f"  阈值往返增益 = {fp.roundtrip_gain_required():.4f}")
    # 透射率采样
    center = 1064.0  # nm -> THz
    f0 = 3e8 / 1064e-9 / 1e12
    print(f"  中心频率附近透射率: T(ν0) = {fp.transmittance(f0):.4f}")


def demo_stable_cavity():
    print()
    print("=" * 60)
    print("Demo 2: 谐振腔稳定性 (g 因子判据)")
    print("=" * 60)
    cases = [
        ("平面平行腔", 0.2, float("inf"), float("inf")),
        ("共焦腔", 0.2, 0.2, 0.2),
        ("半球腔", 0.2, float("inf"), 0.2),
        ("非稳定腔", 0.2, 0.05, 0.05),
    ]
    for name, L, R1, R2 in cases:
        sc = cavity.StableCavity(L, R1, R2)
        print(f"  {name:10s} L={L} R1={R1} R2={R2} -> g1·g2 = {sc.g_product:+.4f} "
              f"({sc.classify()})")


def demo_gaussian_beam():
    print()
    print("=" * 60)
    print("Demo 3: 稳定腔内高斯光束自再现解")
    print("=" * 60)
    gc = beam.GaussianCavity(length=0.2, radius1=float("inf"), radius2=0.5,
                             wavelength_nm=1064.0)
    print(f"  束腰半径 w0 = {gc.waist_radius * 1e6:.2f} μm")
    print(f"  束腰位置（距 M1）= {gc.waist_position:.4f} m")
    print(f"  瑞利长度 zR = {gc.rayleigh:.4f} m")
    print(f"  M1 处光束半径 = {gc.beam_radius_at(0.0) * 1e6:.2f} μm")
    print(f"  M2 处光束半径 = {gc.beam_radius_at(0.2) * 1e6:.2f} μm")
    print(f"  远场发散半角 = {beam.divergence_half_angle(gc.waist_radius, 1064.0) * 1e3:.2f} mrad")


def demo_gain():
    print()
    print("=" * 60)
    print("Demo 4: 增益介质线型")
    print("=" * 60)
    for name in ["nd_yag_1064", "hene_633", "ti_sapphire_800"]:
        med = gain.make_medium(name, g0=10.0)
        g_center = med.small_signal_gain(med.center_thz)
        g_half = med.small_signal_gain(med.center_thz + med.linewidth_thz / 2)
        print(f"  {name:16s} 中心增益 = {g_center:.3f} 1/m, "
              f"半宽处增益 = {g_half:.3f} 1/m (应为中心的一半: {g_center / 2:.3f})")


def demo_relaxation():
    print()
    print("=" * 60)
    print("Demo 5: 弛豫振荡 (Nd:YAG 四能级)")
    print("=" * 60)
    tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1.0e-8
    rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
    print(f"  阈值泵浦速率 Rp_th = {rp_th:.3e} 1/(m³·s)")
    for ratio in [1.5, 2.0, 4.0]:
        f_ro = dynamics.relaxation_frequency(ratio, tau, tau_c)
        sim = dynamics.simulate_relaxation(ratio * rp_th, tau_upper=tau, sigma=sigma,
                                           photon_lifetime=tau_c, t_end=5e-4)
        n_ss, phi_ss = dynamics.steady_state(ratio * rp_th, tau, sigma, tau_c)
        print(f"  r = {ratio:.1f}: 弛豫频率 = {f_ro:.3e} Hz, "
              f"稳态 φ_ss = {phi_ss:.3e}, 仿真末值 = {sim.phi[-1]:.3e}")


def demo_multimode():
    print()
    print("=" * 60)
    print("Demo 6: 多纵模竞争")
    print("=" * 60)
    fp = cavity.FabryPerot(1.0, 0.3, 0.99, 0.95)
    med = gain.make_medium("nd_yag_1064", g0=50.0)
    mm = dynamics.MultiMode(med, fp, pump_ratio=2.0)
    t, N, modes = mm.simulate(n_each_side=3, t_end=1e-3)
    print(f"  纵模数 = {len(modes)}, 采样点数 = {len(t)}")
    # 输出各模最终光子数，观察中心模占优
    finals = sorted(((nu, phi[-1]) for nu, phi in modes.items()), key=lambda x: -x[1])
    print("  稳态光子数分布（从高到低）:")
    for nu, val in finals:
        print(f"    ν = {nu:.4f} THz -> φ = {val:.3e}")


def main():
    demo_fabry_perot()
    demo_stable_cavity()
    demo_gaussian_beam()
    demo_gain()
    demo_relaxation()
    demo_multimode()
    print()
    print("所有演示完成。")


if __name__ == "__main__":
    main()
