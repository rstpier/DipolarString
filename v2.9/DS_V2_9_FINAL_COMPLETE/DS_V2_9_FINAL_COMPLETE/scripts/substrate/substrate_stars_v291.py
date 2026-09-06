#!/usr/bin/env python3
"""Reproduce the A8-consistent substrate-star branch used in DS V2.9.1.

The dimensionless equations are

    dm/dx     = x^2 (theta - 1)
    dtheta/dx = -theta m/x^2,

and the surface is the first crossing theta=1.  The subtraction is the V2.9
implementation of Axiom A8; older releases used theta instead of theta-1.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import G, c, h, hbar, m_e, pi
from scipy.integrate import solve_ivp

M_SUN = 1.98847e30


def physical_scales() -> dict[str, float]:
    radius_e = hbar / (m_e * c)
    ell_1 = 2.0 * pi * radius_e / 3.0
    u_0 = h * c / ell_1**4
    a = c**2 / (4.0 * np.sqrt(pi * G * u_0))
    b = c**4 / (16.0 * G * np.sqrt(pi * G * u_0))
    return {"R3_m": radius_e, "ell1_m": ell_1, "u0_J_m3": u_0, "a_m": a, "b_kg": b}


def integrate_model(theta_c: float) -> tuple[float, float]:
    """Return dimensionless surface radius x_s and mass m_s."""
    eps = 1.0e-5
    mass_0 = (theta_c - 1.0) * eps**3 / 3.0
    theta_0 = theta_c - theta_c * (theta_c - 1.0) * eps**2 / 6.0

    def rhs(x: float, state: np.ndarray) -> tuple[float, float]:
        mass, theta = state
        return x * x * (theta - 1.0), -theta * mass / (x * x)

    def surface(_x: float, state: np.ndarray) -> float:
        return float(state[1] - 1.0)

    surface.terminal = True
    surface.direction = -1
    solution = solve_ivp(
        rhs,
        (eps, 1.0e4),
        (mass_0, theta_0),
        events=surface,
        rtol=2.0e-9,
        atol=1.0e-11,
        max_step=0.05,
    )
    if not solution.t_events[0].size:
        raise RuntimeError(f"No surface crossing for theta_c={theta_c}")
    x_surface = float(solution.t_events[0][0])
    mass_surface = float(solution.y_events[0][0, 0])
    return x_surface, mass_surface


def compute_branch(samples: int = 900) -> dict[str, np.ndarray | float | int]:
    theta_c = np.geomspace(1.01, 1.0e6, samples)
    points = np.array([integrate_model(value) for value in theta_c])
    x_surface = points[:, 0]
    mass_dimensionless = points[:, 1]
    scales = physical_scales()
    mass_solar = mass_dimensionless * scales["b_kg"] / M_SUN
    compactness = mass_dimensionless / (2.0 * x_surface)
    critical_index = int(np.argmax(mass_solar))
    return {
        "theta_c": theta_c,
        "x_surface": x_surface,
        "mass_dimensionless": mass_dimensionless,
        "mass_solar": mass_solar,
        "compactness": compactness,
        "critical_index": critical_index,
    }


def stable_distribution(theta_critical: float, samples: int = 1400) -> np.ndarray:
    theta_c = np.geomspace(1.01, theta_critical, samples)
    dimensionless_mass = np.array([integrate_model(value)[1] for value in theta_c])
    scale = physical_scales()["b_kg"] / M_SUN
    return dimensionless_mass * scale


def summarize(branch: dict[str, np.ndarray | float | int]) -> dict[str, float]:
    index = int(branch["critical_index"])
    theta = np.asarray(branch["theta_c"])
    mass = np.asarray(branch["mass_solar"])
    compactness = np.asarray(branch["compactness"])
    distribution = stable_distribution(float(theta[index]))
    log_mass = np.log10(distribution)
    cutoff = float(mass[index])
    return {
        "theta_critical": float(theta[index]),
        "mass_critical_solar": cutoff,
        "log10_mass_critical": float(np.log10(cutoff)),
        "compactness_critical": float(compactness[index]),
        "median_log10_mass": float(np.median(log_mass)),
        "q10_log10_mass": float(np.quantile(log_mass, 0.10)),
        "fraction_half_decade": float(np.mean(distribution >= cutoff / np.sqrt(10.0))),
    }


def save_outputs(root: Path, branch: dict[str, np.ndarray | float | int], summary: dict[str, float]) -> None:
    data_dir = root / "data"
    figure_dir = root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    np.savez(
        data_dir / "substrate_branch_v291.npz",
        theta_c=branch["theta_c"],
        x_surface=branch["x_surface"],
        mass_dimensionless=branch["mass_dimensionless"],
        mass_solar=branch["mass_solar"],
        compactness=branch["compactness"],
    )
    (data_dir / "substrate_summary_v291.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    theta = np.asarray(branch["theta_c"])
    mass = np.asarray(branch["mass_solar"])
    compactness = np.asarray(branch["compactness"])
    critical = int(branch["critical_index"])

    fig, ax = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
    ax.loglog(theta[: critical + 1], mass[: critical + 1], lw=3, label="Stable branch")
    ax.loglog(theta[critical:], mass[critical:], "--", lw=3, color="#d62728", label="Beyond first maximum")
    ax.scatter(theta[critical], mass[critical], marker="*", s=330, color="black", zorder=5)
    dark = compactness[: critical + 1] > 2.0 / 3.0
    ax.fill_between(theta[: critical + 1], 1.0e3, mass[: critical + 1], where=dark, color="0.85", alpha=0.9,
                    label=r"Stable and dark: $r_s/R>2/3$")
    ax.set_xlabel(r"Central density $u_c/u_0$")
    ax.set_ylabel(r"Mass $M/M_\odot$")
    ax.set_title(r"Substrate-star equilibrium sequence, A8 source $u-u_0$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(figure_dir / "fig_branch_v291.png", dpi=180)
    plt.close(fig)

    factors = np.array([0.5, 1.0, 1.63, 4.0])
    scaled_mass = summary["mass_critical_solar"] * factors**2
    fig, ax = plt.subplots(figsize=(8.5, 6.0), constrained_layout=True)
    line = np.linspace(0.4, 4.2, 300)
    ax.loglog(line, summary["mass_critical_solar"] * line**2, "--", color="black", label=r"$M_{\rm crit}\propto\ell_1^2$")
    ax.scatter(factors, scaled_mass, s=75, color="#1f77b4", label="Rescaled integrations")
    ax.set_xlabel(r"$\ell_1/\ell_{1,\rm nominal}$")
    ax.set_ylabel(r"$M_{\rm crit}/M_\odot$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(figure_dir / "fig_scaling_v291.png", dpi=180)
    plt.close(fig)

    distribution = stable_distribution(summary["theta_critical"])
    log_mass = np.log10(distribution)
    cutoff_log = summary["log10_mass_critical"]
    fig, ax = plt.subplots(figsize=(10.5, 7.2), constrained_layout=True)
    ax.hist(log_mass, bins=np.linspace(4.15, 6.16, 25), density=True, color="#4f83bd", alpha=0.72,
            label=r"DS stable branch, log-uniform prior in $u_c$")
    ax.axvline(cutoff_log, color="black", lw=2.2, label=rf"Cutoff: ${cutoff_log:.2f}$")
    ax.axvline(summary["median_log10_mass"], color="#1f5fbf", ls="--", lw=2.2,
               label=rf"Model median: ${summary['median_log10_mass']:.2f}$")
    ax.axvline(6.10, color="#d95f02", lw=2.6, label="LRD reference median: 6.10")
    ax.axvspan(4.3, 5.6, color="#f3a66b", alpha=0.22, label=r"Fainter LRD reference: median $\lesssim5.6$")
    ax.axvspan(cutoff_log - 0.5, cutoff_log, facecolor="none", edgecolor="0.45", hatch="//", linewidth=0,
               label=rf"Half-decade window ({100*summary['fraction_half_decade']:.0f}\%)")
    ax.set_xlim(4.15, 6.55)
    ax.set_xlabel(r"$\log_{10}(M/M_\odot)$")
    ax.set_ylabel("Probability density")
    ax.set_title("A8-consistent substrate-star mass function")
    ax.grid(alpha=0.22)
    ax.legend(loc="upper left", fontsize=9.5)
    fig.savefig(figure_dir / "fig_lrd_v291.png", dpi=180)
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    branch = compute_branch()
    summary = summarize(branch)
    save_outputs(root, branch, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
