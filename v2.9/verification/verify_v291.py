#!/usr/bin/env python3
"""Verification manifest for the claims retained in DS V2.9.1.

The suite verifies consequences of the published equations.  It does not test
whether the postulated equations describe nature.  Historical V2.7 checks that
target withdrawn claims are deliberately excluded.
"""

from __future__ import annotations

import importlib.util
import itertools
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.constants import alpha as alpha_codata
from scipy.constants import e, epsilon_0, h, hbar, m_e, mu_0, pi, c

ROOT = Path(__file__).resolve().parents[1]
SUBSTRATE_PATH = ROOT / "scripts" / "substrate" / "substrate_stars_v291.py"
JOHNS_DIR = ROOT / "scripts" / "johns"
SCHRODINGER_DIR = ROOT / "scripts" / "schrodinger"
LEGACY_DIR = ROOT / "scripts" / "legacy"

# The reconstructed Appendix-E modules retain their historical flat import layout.
sys.path.insert(0, str(SCHRODINGER_DIR))
import calc1 as sch_calc1  # noqa: E402
import calc1_common as sch_common  # noqa: E402
import calc1b as sch_calc1b  # noqa: E402
import calc3 as sch_calc3  # noqa: E402

spec = importlib.util.spec_from_file_location("substrate_v291", SUBSTRATE_PATH)
substrate = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(substrate)

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    RESULTS.append((status, name, detail))
    print(f"[{status}] {name}: {detail}")


def close(value: float, target: float, relative: float) -> bool:
    return abs(value / target - 1.0) <= relative


def block_geometry_phase() -> None:
    scales = substrate.physical_scales()
    check("Electron reduced-Compton anchor", close(scales["R3_m"], 3.86159e-13, 2e-5), f"R3={scales['R3_m']:.7e} m")
    check("String length anchor", close(scales["ell1_m"], 8.0877e-13, 2e-5), f"ell1={scales['ell1_m']:.7e} m")
    check("Cell energy density", close(scales["u0_J_m3"], 4.6428e23, 3e-4), f"u0={scales['u0_J_m3']:.7e} J/m^3")

    ratio = 2.0 * np.cosh(pi)
    gamma_pole = (1.0 - 0.5) / (1.0 + 0.5)
    check("Bifilar geometry", close(ratio, 23.1839, 2e-5), f"D/r=2cosh(pi)={ratio:.7f}")
    check("Pole reflection coefficient", abs(gamma_pole - 1.0 / 3.0) < 1e-15, f"Gamma={gamma_pole:.12g}")

    z0 = np.sqrt(mu_0 / epsilon_0)
    alpha_imp = e**2 * z0 / (4.0 * pi * hbar)
    q = e / 3.0
    winding = 1.0 / 3.0
    coefficient = 3.0 * pi * hbar * winding**2 / (q**2 * z0)
    expected = 3.0 / (4.0 * alpha_imp)
    check("Impedance form of alpha", close(alpha_imp, alpha_codata, 5e-10), f"alpha_Z={alpha_imp:.12g}")
    check("Conditional Z3 coefficient", abs(coefficient / expected - 1.0) < 2e-15,
          f"Delta_Z3={coefficient:.9f}=3/(4alpha_Z)")


def block_substrate() -> None:
    # Dense enough to reproduce the printed values while keeping the manifest quick.
    theta = np.geomspace(1.01, 100.0, 360)
    points = np.array([substrate.integrate_model(value) for value in theta])
    x_surface, mass_dimensionless = points[:, 0], points[:, 1]
    scales = substrate.physical_scales()
    mass_solar = mass_dimensionless * scales["b_kg"] / substrate.M_SUN
    compactness = mass_dimensionless / (2.0 * x_surface)
    index = int(np.argmax(mass_solar))
    theta_critical = float(theta[index])
    mass_critical = float(mass_solar[index])
    compact_critical = float(compactness[index])

    check("A8 critical density", abs(theta_critical - 14.5) < 0.5, f"uc/u0={theta_critical:.4f}")
    check("A8 critical mass", close(mass_critical, 1.29e6, 0.006), f"Mcrit={mass_critical:.7g} Msun")
    check("A8 critical compactness", abs(compact_critical - 0.81) < 0.02, f"rs/R={compact_critical:.5f}")

    distribution = substrate.stable_distribution(theta_critical, samples=650)
    log_mass = np.log10(distribution)
    fraction = float(np.mean(distribution >= mass_critical / np.sqrt(10.0)))
    check("Mass-function median", abs(np.median(log_mass) - 6.02) < 0.02, f"median={np.median(log_mass):.5f}")
    check("Mass-function 10-percent quantile", abs(np.quantile(log_mass, 0.10) - 5.49) < 0.03,
          f"q10={np.quantile(log_mass, 0.10):.5f}")
    check("Registered half-decade fraction", abs(fraction - 0.86) < 0.02, f"f1/2={fraction:.5f}")

    local_theta = np.linspace(13.8, 15.1, 27)
    local_mass = np.array([substrate.integrate_model(value)[1] for value in local_theta])
    poly = np.polyfit(local_theta, local_mass, 2)
    fitted = np.polyval(poly, local_theta)
    r2 = 1.0 - np.sum((local_mass - fitted) ** 2) / np.sum((local_mass - local_mass.mean()) ** 2)
    vertex = -poly[1] / (2.0 * poly[0])
    check("Simple quadratic turning point", poly[0] < 0 and r2 > 0.999 and abs(vertex - 14.5) < 0.4,
          f"vertex={vertex:.4f}, R2={r2:.7f}; implies exponent -1/2")

    # From u0 proportional to ell1^-4 and b proportional to u0^-1/2.
    ell_nominal = scales["ell1_m"]
    scale_factors = np.array([0.5, 1.63, 4.0])
    def mass_scale(ell: float) -> float:
        u0 = h * c / ell**4
        return c**4 / (16.0 * substrate.G * np.sqrt(pi * substrate.G * u0))
    normalized = np.array([mass_scale(factor * ell_nominal) / mass_scale(ell_nominal) / factor**2
                           for factor in scale_factors])
    check("Exact ell1-squared mass scaling", np.allclose(normalized, 1.0, atol=2e-15),
          f"normalized={normalized}")


def johns_ports() -> tuple[list[tuple[int, int, int]], dict[tuple[int, int, int], int]]:
    ports = [(d, s, p) for d in range(3) for s in (+1, -1) for p in range(3) if p != d]
    return ports, {port: index for index, port in enumerate(ports)}


def reconstruct_johns_candidates() -> list[np.ndarray]:
    ports, index = johns_ports()

    def third(d: int, p: int) -> int:
        return 3 - d - p

    candidates = []
    for choice in itertools.product((0, 1), repeat=12):
        matrix = np.zeros((12, 12))
        for row, (d, _s, p) in enumerate(ports):
            t = third(d, p)
            matrix[row, index[(t, +1, p)]] = 0.5
            matrix[row, index[(t, -1, p)]] = 0.5
            matrix[row, index[(p, +1, d)]] = -0.5 if choice[row] == 0 else 0.5
            matrix[row, index[(p, -1, d)]] = -0.5 if choice[row] == 1 else 0.5
        if np.allclose(matrix, matrix.T) and np.allclose(matrix @ matrix.T, np.eye(12)):
            candidates.append(matrix)
    return candidates


def bloch_frequencies(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    ports, index = johns_ports()
    propagation = np.zeros((12, 12), dtype=complex)
    for column, (direction, side, polarization) in enumerate(ports):
        propagation[index[(direction, -side, polarization)], column] = np.exp(1j * vector[direction] * side)
    return np.sort(np.angle(np.linalg.eigvals(propagation @ matrix)))


def acoustic_speeds(matrix: np.ndarray, direction: np.ndarray, magnitude: float = 0.05) -> np.ndarray:
    omega = bloch_frequencies(matrix, direction / np.linalg.norm(direction) * magnitude)
    positive = omega[(omega > 1e-8) & (omega < 1.0)]
    return positive[:2] / magnitude if positive.size >= 2 else np.array([])


def block_johns() -> None:
    matrix = np.load(JOHNS_DIR / "S_johns.npy")
    check("Johns matrix shape", matrix.shape == (12, 12), f"shape={matrix.shape}")
    check("Johns symmetry", np.allclose(matrix, matrix.T, atol=1e-14), "S=S^T")
    check("Johns unitarity", np.allclose(matrix.T @ matrix, np.eye(12), atol=1e-14), "S^T S=I")
    check("Johns entry alphabet", set(np.unique(matrix)).issubset({-0.5, 0.0, 0.5}), f"entries={np.unique(matrix)}")
    check("No self-reflection", np.allclose(np.diag(matrix), 0.0), "diag(S)=0")

    candidates = reconstruct_johns_candidates()
    check("Constraint scan", len(candidates) == 8, f"symmetric unitary candidates={len(candidates)}")
    ports, index = johns_ports()
    def rotation_matrix(mapping: dict[int, tuple[int, int]]) -> np.ndarray:
        rotation = np.zeros((12, 12))
        for column, (direction, side, polarization) in enumerate(ports):
            direction_new, direction_sign = mapping[direction]
            polarization_new, polarization_sign = mapping[polarization]
            rotation[index[(direction_new, side * direction_sign, polarization_new)], column] = polarization_sign
        return rotation
    rz = rotation_matrix({0: (1, +1), 1: (0, -1), 2: (2, +1)})
    rx = rotation_matrix({0: (0, +1), 1: (2, +1), 2: (1, -1)})
    cubic = [candidate for candidate in candidates
             if np.allclose(rz @ candidate @ rz.T, candidate)
             and np.allclose(rx @ candidate @ rx.T, candidate)]
    check("Cubic-equivariant candidates", len(cubic) == 2, f"rotation-equivariant candidates={len(cubic)}")

    isotropic = []
    for candidate in cubic:
        speed_x = acoustic_speeds(candidate, np.array([1.0, 0.0, 0.0]))
        speed_diag = acoustic_speeds(candidate, np.array([1.0, 1.0, 1.0]))
        if speed_x.size == 2 and speed_diag.size == 2 and np.allclose(speed_x, speed_diag, atol=5e-4):
            isotropic.append(candidate)
    check("Low-q isotropic selection", len(isotropic) == 1, f"isotropic candidates={len(isotropic)}")
    if isotropic:
        check("Archived Johns node selected", np.allclose(matrix, isotropic[0]), "S_johns equals selected candidate")
        speed = acoustic_speeds(matrix, np.array([1.0, 1.0, 1.0]))
        check("Two transverse acoustic modes", speed.size == 2 and np.allclose(speed, 0.5, atol=5e-4),
              f"v/c_link={speed}")
    generic = bloch_frequencies(matrix, np.array([0.31, 0.17, 0.53]))
    check("Two zero-frequency constraint modes", np.sum(np.abs(generic) < 1e-9) == 2,
          f"count={np.sum(np.abs(generic) < 1e-9)}")


def block_manuscript_hygiene() -> None:
    source = (ROOT / "manuscript" / "DS_model_V2_9_1_Zenodo.tex").read_text(encoding="utf-8")
    check("Low-impedance equation withdrawn", "eq:vsync" not in source, "no synchronization-speed equation label")
    check("A8 source printed", "u - u_0" in source, "corrected source present")
    check("Conditional phase equation printed", "eq:phase_action" in source and "eq:z3_coefficient" in source,
          "phase normalization and Z3 coefficient present")
    check("Instanton-to-G gap disclosed", "The two missing bridges to $G$" in source,
          "missing gravitational map stated")


def block_schrodinger_reconstruction() -> None:
    expected_files = ["calc1.py", "calc1b.py", "calc1_common.py", "calc2.py", "calc3.py"]
    present = [name for name in expected_files if (SCHRODINGER_DIR / name).is_file()]
    check("Appendix-E reconstruction inventory", len(present) == 5, f"files={len(present)}/5")

    band = sch_common.trace_bound_band()
    cache_ok = band.k.shape == (30,) and band.modes.shape == (30, 12 * 11 * 11)
    check("Reconstructed 30-point band cache", cache_ok, f"k={band.k.shape}, modes={band.modes.shape}")

    near = int(np.argmin(np.abs(band.k - 1.5)))
    anchor_ok = abs(band.omega[near] - 0.36558142) < 2e-8 and band.core_weight[near] > 0.49
    check("Braided-defect band anchor", anchor_ok,
          f"omega={band.omega[near]:.8f}, core={band.core_weight[near]:.5f}")

    free = sch_calc1.run()
    free_ok = abs(float(free["variance_T_spectral"]) - 23.455) < 0.003 and abs(float(free["variance_T_quadratic"]) - 22.87) < 0.02
    check("Free-envelope anchors", free_ok,
          f"variance={free['variance_T_spectral']:.6f}, quadratic={free['variance_T_quadratic']:.6f}")

    interference = sch_calc1b.run()
    interference_ok = abs(float(interference["measured_peak"]) - 0.3926990817) < 1e-9 and float(interference["correlation_with_complex_envelope"]) > 0.9995
    check("Two-packet interference", interference_ok,
          f"peak={interference['measured_peak']:.10f}, corr={interference['correlation_with_complex_envelope']:.6f}")

    action = sch_calc3.run()
    action_ok = abs(float(action["E_over_N"]) - 0.36713) < 1e-5 and abs(float(action["P_over_N"]) - 1.5024) < 1e-4 and abs(float(action["mF_times_omega2_mean_single_band"]) - 1.0) < 1e-12
    check("Single-action-quantum closures", action_ok,
          f"E/N={action['E_over_N']:.7f}, P/N={action['P_over_N']:.7f}, m*omega2={action['mF_times_omega2_mean_single_band']:.6f}")


def optional_symbolic_checks() -> None:
    try:
        import sympy  # noqa: F401
    except ModuleNotFoundError:
        print("[SKIP] Optional symbolic checks: install requirements.txt (sympy missing)")
        return
    for script, label in [
        ("theorem3_covariance.py", "Scalar-action covariance identity"),
        ("conformal_identity.py", "Conformal representative identity"),
        ("boosted_residual.py", "Flat linear wave residual only"),
        ("ppn_2pn.py", "Scalar-profile PPN diagnostic only"),
    ]:
        result = subprocess.run([sys.executable, str(LEGACY_DIR / script)], cwd=LEGACY_DIR, capture_output=True, text=True)
        check(label, result.returncode == 0, f"{script}, exit={result.returncode}")


def main() -> int:
    start = time.time()
    print("DS V2.9.1 CURRENT-CLAIM VERIFICATION")
    print("=" * 72)
    block_geometry_phase()
    block_substrate()
    block_johns()
    block_manuscript_hygiene()
    block_schrodinger_reconstruction()
    optional_symbolic_checks()
    failed = [item for item in RESULTS if item[0] == "FAIL"]
    print("=" * 72)
    print(f"RESULT: {len(RESULTS) - len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL; {time.time()-start:.1f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
