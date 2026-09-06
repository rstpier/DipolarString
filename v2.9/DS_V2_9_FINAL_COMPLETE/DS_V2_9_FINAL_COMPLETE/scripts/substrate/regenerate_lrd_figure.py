"""Regenerate the A8-consistent LRD comparison figure for DS V2.9.1."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def surface_mass(theta_c: float) -> float:
    """Integrate dm/dx=x^2(theta-1), dtheta/dx=-theta*m/x^2 to theta=1."""
    eps = 1.0e-5
    m0 = (theta_c - 1.0) * eps**3 / 3.0
    theta0 = theta_c - theta_c * (theta_c - 1.0) * eps**2 / 6.0

    def rhs(x, state):
        mass, theta = state
        return x * x * (theta - 1.0), -theta * mass / (x * x)

    def surface(_x, state):
        return state[1] - 1.0

    surface.terminal = True
    surface.direction = -1
    solution = solve_ivp(
        rhs,
        (eps, 1.0e4),
        (m0, theta0),
        events=surface,
        rtol=2.0e-9,
        atol=1.0e-11,
        max_step=0.05,
    )
    if not solution.t_events[0].size:
        raise RuntimeError(f"No surface crossing for theta_c={theta_c}")
    return float(solution.y_events[0][0, 0])


theta_c = np.geomspace(1.01, 14.5, 1200)
dimensionless_mass = np.array([surface_mass(value) for value in theta_c])
mass = dimensionless_mass / dimensionless_mass.max() * 1.29e6
log_mass = np.log10(mass)

critical_log_mass = np.log10(1.29e6)
model_median = float(np.median(log_mass))
model_q10 = float(np.quantile(log_mass, 0.10))
registered_lower = critical_log_mass - 0.5

fig, ax = plt.subplots(figsize=(10.5, 7.2), constrained_layout=True)
ax.hist(
    log_mass,
    bins=np.linspace(4.15, 6.16, 25),
    density=True,
    color="#4f83bd",
    alpha=0.72,
    label=r"DS prediction (stable branch, log-uniform prior in $u_c$)",
)
ax.axvline(
    critical_log_mass,
    color="black",
    linewidth=2.2,
    label=rf"Model cutoff: $\log_{{10}} M_{{\rm crit}}={critical_log_mass:.2f}$",
)
ax.axvline(
    model_median,
    color="#1f5fbf",
    linestyle="--",
    linewidth=2.2,
    label=rf"Model median: ${model_median:.2f}$",
)
ax.axvline(
    6.10,
    color="#d95f02",
    linewidth=2.6,
    label=r"LRD median, 7 firm masses: $6.10$ (Rusakov et al. 2026)",
)
ax.axvspan(
    4.3,
    5.6,
    color="#f3a66b",
    alpha=0.22,
    label=r"18 fainter LRDs: reported median $\log M\lesssim5.6$",
)
ax.axvspan(
    registered_lower,
    critical_log_mass,
    facecolor="none",
    edgecolor="#777777",
    hatch="//",
    linewidth=0.0,
    label=r"Registered half-decade window ($86\%$ on corrected branch)",
)
ax.set_xlim(4.15, 6.55)
ax.set_xlabel(r"$\log_{10}(M/M_\odot)$", fontsize=14)
ax.set_ylabel("Probability density", fontsize=14)
ax.set_title("A8-consistent substrate-star mass function and LRD reference values", fontsize=16)
ax.grid(alpha=0.22)
ax.legend(loc="upper left", fontsize=10.5, framealpha=0.96)
ax.text(
    0.985,
    0.035,
    rf"10th percentile ${model_q10:.2f}$; dynamical-mass tension discussed in text",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=10,
    color="#444444",
)
fig.savefig("fig_lrd.png", dpi=180)
