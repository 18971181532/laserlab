"""LaserLab 测试套件.

覆盖解析解验证、物理一致性契约、数值积分器正确性。
运行：python run_tests.py
"""

import sys
import os
import unittest
from math import pi, sqrt, isclose

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from laserlab import constants as C  # noqa: E402
from laserlab import cavity  # noqa: E402
from laserlab import beam  # noqa: E402
from laserlab import gain  # noqa: E402
from laserlab import rate  # noqa: E402
from laserlab import dynamics  # noqa: E402


class TestConstants(unittest.TestCase):
    def test_freq_wl_roundtrip(self):
        self.assertAlmostEqual(C.wl_to_freq(1064.0), 281.760, places=3)
        self.assertAlmostEqual(C.freq_to_wl(C.wl_to_freq(1064.0)), 1064.0, places=8)

    def test_wl_to_freq_relation(self):
        # c = λ·ν 精确成立
        lam = 532e-9
        self.assertAlmostEqual(C.C / lam / 1e12, C.wl_to_freq(532.0), places=8)


class TestFabryPerot(unittest.TestCase):
    def setUp(self):
        self.fp = cavity.FabryPerot(n=1.0, length=0.1, r1=0.9, r2=0.9)

    def test_fsr(self):
        # FSR = c/(2nL) = 2.9979e8/0.2 = 1.49896e9 Hz
        expected = C.C / (2 * 0.1)
        self.assertAlmostEqual(self.fp.fsr_hz, expected, places=4)

    def test_fsr_thz(self):
        self.assertAlmostEqual(self.fp.fsr_thz, C.C / 0.2 / 1e12, places=6)

    def test_finesse_known_value(self):
        # R=0.9 -> finesse = π√R/(1-R) = π·0.9487/0.1 = 29.8
        expected = pi * sqrt(0.9) / 0.1
        self.assertAlmostEqual(self.fp.finesse, expected, places=2)

    def test_transmittance_resonance(self):
        # 谐振频率处 δ=0，透射率应接近峰值
        f0 = C.C / (2 * self.fp.optical_length)  # 第一个纵模
        # 单程相位 δ = 2π n L ν / c = π（对第一个谐振需整数倍 π）
        # 用中心频率对应 δ 为 π 的整数倍来构造谐振
        nu = 1.5e9  # 恰好一个 FSR，δ = π·整数
        t = self.fp.transmittance(nu / 1e12)  # 转 THz
        self.assertGreater(t, 0.99)  # 接近峰值透射

    def test_transmittance_antiresonance(self):
        # 半 FSR 处透射率最小
        nu = 0.75e9
        t = self.fp.transmittance(nu / 1e12)
        self.assertLess(t, 0.1)

    def test_reflectivity_complement(self):
        nu = 1.2e9
        t = self.fp.transmittance(nu / 1e12)
        r = self.fp.reflectivity(nu / 1e12)
        self.assertAlmostEqual(t + r, 1.0, places=10)

    def test_linewidth_equals_fsr_over_finesse(self):
        self.assertAlmostEqual(self.fp.linewidth_hz(),
                               self.fp.fsr_hz / self.fp.finesse, places=6)

    def test_invalid_reflectance(self):
        with self.assertRaises(ValueError):
            cavity.FabryPerot(1.0, 0.1, 1.0, 0.9)

    def test_roundtrip_gain(self):
        # 阈值往返增益 = 1/(R1 R2)
        self.assertAlmostEqual(self.fp.roundtrip_gain_required(), 1.0 / (0.9 * 0.9), places=10)


class TestStableCavity(unittest.TestCase):
    def test_plane_parallel_critical(self):
        sc = cavity.StableCavity(0.2, float("inf"), float("inf"))
        self.assertEqual(sc.classify(), "critical")

    def test_confocal_critical(self):
        sc = cavity.StableCavity(0.2, 0.2, 0.2)
        self.assertEqual(sc.classify(), "critical")

    def test_unstable(self):
        sc = cavity.StableCavity(0.2, 0.05, 0.05)
        self.assertEqual(sc.classify(), "unstable")
        self.assertGreater(sc.g_product, 1.0)

    def test_stable(self):
        # R1=R2=1.0, L=0.2 -> g1=g2=0.8, g1g2=0.64
        sc = cavity.StableCavity(0.2, 1.0, 1.0)
        self.assertEqual(sc.classify(), "stable")
        self.assertAlmostEqual(sc.g_product, 0.64, places=10)


class TestBeam(unittest.TestCase):
    def test_rayleigh_range(self):
        # zR = π w0² / λ
        w0 = 1e-3
        lam = 1064e-9
        self.assertAlmostEqual(beam.rayleigh_range(w0, 1064.0), pi * w0 ** 2 / lam, places=6)

    def test_beam_radius_at_rayleigh(self):
        # z = zR 处 w(zR) = √2 · w0
        w0 = 1e-3
        zr = beam.rayleigh_range(w0, 1064.0)
        q = beam.q_parameter(w0, 0.0, 1064.0)
        qz = q + zr
        w = beam.beam_radius(qz, 1064.0)
        self.assertAlmostEqual(w, sqrt(2) * w0, places=6)

    def test_divergence(self):
        w0 = 1e-3
        self.assertAlmostEqual(beam.divergence_half_angle(w0, 1064.0),
                               1064e-9 / (pi * w0), places=8)

    def test_thin_lens_q(self):
        # 平行光束过薄透镜，焦点处曲率半径为 f
        q = 1e10j  # 近似平面波（大束腰）
        f = 0.1
        q2 = beam.apply_matrix(q, beam.thin_lens(f))
        # 出射后 q = 1/(-1/f + 1/q) ≈ -f + i0 附近
        self.assertAlmostEqual(q2.real, -f, places=2)

    def test_matrix_compose_identity(self):
        m = beam.compose(beam.free_space(0.1), beam.free_space(-0.1))
        self.assertAlmostEqual(m[0][0], 1.0)
        self.assertAlmostEqual(m[1][1], 1.0)
        self.assertAlmostEqual(m[0][1], 0.0)

    def test_gaussian_cavity_symmetric(self):
        # 对称共心稳定腔（R1=R2=R, R>L/2）：束腰应在腔中心 z=L/2
        # 取 R=1.0, L=0.2 -> g1=g2=0.8 稳定
        gc = beam.GaussianCavity(0.2, 1.0, 1.0, 1064.0)
        self.assertAlmostEqual(gc.waist_position, 0.1, places=6)

    def test_gaussian_cavity_unstable_raises(self):
        with self.assertRaises(ValueError):
            beam.GaussianCavity(0.2, 0.05, 0.05, 1064.0)

    def test_mode_overlap_identical(self):
        self.assertAlmostEqual(beam.mode_overlap(1e-3, 1e-3), 1.0, places=10)

    def test_mode_overlap_mismatch(self):
        eta = beam.mode_overlap(1e-3, 2e-3)
        self.assertLess(eta, 1.0)
        self.assertGreater(eta, 0.0)

    def test_coupling_efficiency_perfect(self):
        w0 = 1e-3
        q = beam.q_parameter(w0, 0.0, 1064.0)
        self.assertAlmostEqual(beam.coupling_efficiency(q, q), 1.0, places=8)


class TestGain(unittest.TestCase):
    def test_lorentzian_half_max(self):
        med = gain.make_medium("nd_yag_1064", 10.0)
        g0 = med.small_signal_gain(med.center_thz)
        ghalf = med.small_signal_gain(med.center_thz + med.linewidth_thz / 2)
        self.assertAlmostEqual(ghalf, g0 / 2.0, places=6)

    def test_saturation_reduces_gain(self):
        med = gain.make_medium("nd_yag_1064", 10.0)
        g_unsat = med.saturated_gain(med.center_thz, 0.0)
        g_sat = med.saturated_gain(med.center_thz, med.isat)
        self.assertAlmostEqual(g_sat, g_unsat / 2.0, places=6)

    def test_gaussian_lineshape(self):
        med = gain.InhomogeneousMedium(10.0, 281.7, 1.0, 1e7)
        g0 = med.small_signal_gain(281.7)
        ghalf = med.small_signal_gain(281.7 + 1.0 / 2)
        self.assertAlmostEqual(ghalf, g0 / 2.0, places=4)

    def test_unknown_medium(self):
        with self.assertRaises(KeyError):
            gain.make_medium("nonexistent", 10.0)


class TestRate(unittest.TestCase):
    def test_rk4_exponential(self):
        # dy/dt = -y, y(0)=1 -> y(t) = e^-t
        f = lambda t, y: [-y[0]]
        t, y = rate.solve_rk4(f, 0.0, [1.0], 1.0, 0.001)
        self.assertAlmostEqual(y[-1][0], 1.0 / 2.718281828, places=3)

    def test_rk45_exponential_accuracy(self):
        f = lambda t, y: [-y[0]]
        t, y = rate.solve_rk45(f, 0.0, [1.0], 1.0)
        self.assertAlmostEqual(y[-1][0], 1.0 / 2.718281828, places=5)

    def test_rk45_harmonic(self):
        # d²y/dt² = -y -> y = cos(t), dy/dt = -sin(t)
        f = lambda t, y: [y[1], -y[0]]
        t, y = rate.solve_rk45(f, 0.0, [1.0, 0.0], 2 * pi)
        self.assertAlmostEqual(y[-1][0], 1.0, places=4)
        self.assertAlmostEqual(y[-1][1], 0.0, places=4)

    def test_threshold_pump_rate(self):
        # Rp_th = 1/(c σ τ τc)
        tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1e-8
        expected = 1.0 / (3e8 * sigma * tau * tau_c)
        self.assertAlmostEqual(rate.threshold_pump_rate(tau, sigma, tau_c), expected, places=12)


class TestDynamics(unittest.TestCase):
    def test_steady_state_below_threshold(self):
        tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1e-8
        rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
        n_ss, phi_ss = dynamics.steady_state(0.5 * rp_th, tau, sigma, tau_c)
        self.assertEqual(phi_ss, 0.0)

    def test_steady_state_above_threshold(self):
        tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1e-8
        rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
        n_ss, phi_ss = dynamics.steady_state(2.0 * rp_th, tau, sigma, tau_c)
        self.assertGreater(phi_ss, 0.0)
        # N_ss 应等于 1/(c σ τc)（增益钳制）
        self.assertAlmostEqual(n_ss, 1.0 / (3e8 * sigma * tau_c), places=12)

    def test_relaxation_frequency_scaling(self):
        # Ω ∝ sqrt(r-1)，r 增大频率升高
        tau, tau_c = 2.3e-4, 1e-8
        f1 = dynamics.relaxation_frequency(2.0, tau, tau_c)
        f2 = dynamics.relaxation_frequency(4.0, tau, tau_c)
        self.assertGreater(f2, f1)
        self.assertAlmostEqual(f2 / f1, sqrt(3.0 / 1.0), places=4)

    def test_relaxation_frequency_below_threshold(self):
        self.assertEqual(dynamics.relaxation_frequency(0.5, 2.3e-4, 1e-8), 0.0)

    def test_simulate_converges_to_steady(self):
        tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1e-8
        rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
        sim = dynamics.simulate_relaxation(2.0 * rp_th, tau_upper=tau, sigma=sigma,
                                           photon_lifetime=tau_c, t_end=1e-3)
        n_ss, phi_ss = dynamics.steady_state(2.0 * rp_th, tau, sigma, tau_c)
        # 仿真末值应收敛到稳态（相对误差 < 5%）
        self.assertAlmostEqual(sim.phi[-1], phi_ss, delta=phi_ss * 0.05)

    def test_simulate_relaxation_has_peak(self):
        # 弛豫振荡应存在超过稳态的峰值（超调）
        tau, sigma, tau_c = 2.3e-4, 2.8e-23, 1e-8
        rp_th = rate.threshold_pump_rate(tau, sigma, tau_c)
        sim = dynamics.simulate_relaxation(2.0 * rp_th, tau_upper=tau, sigma=sigma,
                                           photon_lifetime=tau_c, t_end=5e-4)
        n_ss, phi_ss = dynamics.steady_state(2.0 * rp_th, tau, sigma, tau_c)
        self.assertGreater(max(sim.phi), phi_ss * 1.1)


class TestMultiMode(unittest.TestCase):
    def test_center_mode_dominates(self):
        fp = cavity.FabryPerot(1.0, 0.3, 0.99, 0.95)
        med = gain.make_medium("nd_yag_1064", 50.0)
        mm = dynamics.MultiMode(med, fp, pump_ratio=2.0)
        t, N, modes = mm.simulate(n_each_side=3, t_end=5e-4)
        finals = {nu: phi[-1] for nu, phi in modes.items()}
        center = max(finals, key=lambda nu: -abs(nu - med.center_thz))
        # 中心模最终光子数应最高（模式竞争结果）
        self.assertEqual(max(finals, key=finals.get), center)


if __name__ == "__main__":
    unittest.main(verbosity=2)
