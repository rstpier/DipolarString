#!/usr/bin/env python3
"""Reconstructed Appendix-E calculation 3: one action quantum and inertia."""

from __future__ import annotations

import argparse
import json

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

from calc1_common import (
    K0,
    assemble_full_state,
    circular_moments,
    full_map_step,
    gaussian_amplitudes,
    load_scattering,
    spatial_density,
    trace_bound_band,
)


def run(alpha: float = 4e-4, steps: int = 1200, full_map: bool = False) -> dict[str, float | list[float] | str]:
    band = trace_bound_band()
    amplitude = gaussian_amplitudes(band.k)
    modal_energy = np.abs(amplitude) ** 2
    count = float(np.sum(modal_energy / band.omega))
    energy = float(modal_energy.sum())
    momentum = float(np.sum(modal_energy * band.k / band.omega))

    spline = CubicSpline(band.k, band.omega)
    times = np.linspace(0.0, steps, 13)
    swept_k = K0 + alpha * times
    velocity = spline(swept_k, 1)
    curvature = spline(swept_k, 2)
    acceleration = alpha * curvature
    force_mass = np.divide(alpha, acceleration, out=np.full_like(acceleration, np.nan), where=np.abs(acceleration) > 1e-14)
    closure = force_mass * curvature

    result: dict[str, float | list[float] | str] = {
        "status": "RECONSTRUCTED; historical finite-difference mass convention not recovered",
        "hbar_lattice": 1.0,
        "N": count,
        "E_over_N": energy / count,
        "omega_k0": float(spline(K0)),
        "P_over_N": momentum / count,
        "weighted_k": float(np.sum(modal_energy * band.k / band.omega) / np.sum(modal_energy / band.omega)),
        "sweep_k": [float(swept_k[0]), float(swept_k[-1])],
        "velocity_minmax": [float(np.min(velocity)), float(np.max(velocity))],
        "mF_times_omega2_mean_single_band": float(np.nanmean(closure)),
        "reference_full_map_mean": 0.975,
        "reference_full_map_range": [0.92, 0.995],
        "norm": float(modal_energy.sum()),
        "interpretation": (
            "single-band normalization consistency; not an independent derivation of the "
            "physical Planck constant or the Born rule"
        ),
    }
    if full_map:
        z0 = 96.0
        state = assemble_full_state(band, amplitude, z0)
        scattering = load_scattering()
        stride = 5
        sampled_times = [0.0]
        centers = [z0]
        center_hint = z0
        for time_index in range(steps):
            state = full_map_step(state, scattering=scattering, gauge_phase=alpha * time_index)
            if (time_index + 1) % stride == 0:
                center, _ = circular_moments(spatial_density(state), center_hint)
                center_hint = center
                sampled_times.append(float(time_index + 1))
                centers.append(center)
        sampled_times_array = np.asarray(sampled_times)
        centers_array = np.asarray(centers)
        dt = float(stride)
        window = min(31, len(centers_array) // 2 * 2 - 1)
        velocity_measured = savgol_filter(centers_array, window, 3, deriv=1, delta=dt)
        acceleration_measured = savgol_filter(centers_array, window, 3, deriv=2, delta=dt)
        k_sweep = K0 + alpha * sampled_times_array
        velocity_predicted = spline(k_sweep, 1)
        curvature_predicted = spline(k_sweep, 2)
        interior = slice(window // 2, -window // 2)
        relative_velocity = np.abs(
            (velocity_measured[interior] - velocity_predicted[interior]) / velocity_predicted[interior]
        )
        inertial_closure = alpha * curvature_predicted[interior] / acceleration_measured[interior]
        finite = np.isfinite(inertial_closure) & (np.abs(acceleration_measured[interior]) > 1e-7)
        result.update(
            full_map_dimension=int(np.prod(state.shape)),
            full_map_norm=float(np.vdot(state, state).real),
            velocity_relative_error_max=float(np.max(relative_velocity)),
            mF_times_omega2_mean_full_map=float(np.mean(inertial_closure[finite])),
            mF_times_omega2_range_full_map=[
                float(np.min(inertial_closure[finite])),
                float(np.max(inertial_closure[finite])),
            ],
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=4e-4)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--full-map", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.alpha, args.steps, args.full_map), indent=2))


if __name__ == "__main__":
    main()
