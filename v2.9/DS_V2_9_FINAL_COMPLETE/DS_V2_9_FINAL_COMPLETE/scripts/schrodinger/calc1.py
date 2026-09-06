#!/usr/bin/env python3
"""Reconstructed Appendix-E calculation 1: free single-band envelope."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from calc1_common import (
    K0,
    SIGMA_K,
    assemble_full_state,
    circular_moments,
    full_map_step,
    gaussian_amplitudes,
    local_derivatives,
    load_scattering,
    project_band,
    spatial_density,
    trace_bound_band,
)


REFERENCE = {
    "variance_T": 23.455,
    "variance_exact_predictor": 23.452,
    "variance_quadratic_predictor": 22.87,
    "modal_amplitude_error": 1.9e-8,
    "modal_phase_error_rad": 1.5e-9,
    "quadratic_phase_residual": 0.084,
    "cubic_phase_residual": 0.016,
}


def run(steps: int = 480, full_map: bool = False) -> dict[str, float | str]:
    band = trace_bound_band()
    amplitude = gaussian_amplitudes(band.k)
    velocity, curvature = local_derivatives(band.k, band.omega)
    mean_velocity = float(np.sum(np.abs(amplitude) ** 2 * velocity))
    state0 = assemble_full_state(band, amplitude, z0=64.0)
    state_t = assemble_full_state(band, amplitude * np.exp(1j * band.omega * steps), z0=64.0)
    _, variance0 = circular_moments(spatial_density(state0), 64.0)
    _, variance_t = circular_moments(spatial_density(state_t), 64.0 + mean_velocity * steps)
    weights = np.abs(amplitude) ** 2
    exact_predictor = variance0 + steps**2 * float(weights @ (velocity - weights @ velocity) ** 2)

    fit = np.poly1d(np.polyfit(band.k - K0, band.omega, 4, w=np.sqrt(weights)))
    omega0 = float(fit(0.0))
    v0 = float(np.polyder(fit, 1)(0.0))
    omega2 = float(np.polyder(fit, 2)(0.0))
    quadratic = omega0 + v0 * (band.k - K0) + 0.5 * omega2 * (band.k - K0) ** 2
    v_quadratic = v0 + omega2 * (band.k - K0)
    quadratic_state = assemble_full_state(
        band, amplitude * np.exp(1j * quadratic * steps), z0=64.0
    )
    _, quadratic_predictor = circular_moments(
        spatial_density(quadratic_state), 64.0 + float(weights @ v_quadratic) * steps
    )
    phase_quadratic = float(np.sqrt(weights @ ((band.omega - quadratic) * steps) ** 2))
    omega3 = float(np.polyder(fit, 3)(0.0))
    cubic = quadratic + omega3 * (band.k - K0) ** 3 / 6.0
    phase_cubic = float(np.sqrt(weights @ ((band.omega - cubic) * steps) ** 2))

    result: dict[str, float | str] = {
        "status": "RECONSTRUCTED, not original source recovery",
        "band_points": int(len(band.k)),
        "omega_at_k0_fit": omega0,
        "group_velocity_at_k0": v0,
        "curvature_at_k0": omega2,
        "initial_variance": variance0,
        "variance_T_spectral": variance_t,
        "variance_T_exact_predictor": exact_predictor,
        "variance_T_quadratic": quadratic_predictor,
        "quadratic_phase_rms_rad": phase_quadratic,
        "cubic_fit_phase_rms_rad": phase_cubic,
        "reference_variance_T": REFERENCE["variance_T"],
    }
    if full_map:
        direct = state0.copy()
        scattering = load_scattering()
        initial_coefficients = project_band(direct, band)
        for _ in range(steps):
            direct = full_map_step(direct, scattering=scattering)
        final_coefficients = project_band(direct, band)
        expected = initial_coefficients * np.exp(1j * band.omega * steps)
        supported = np.abs(initial_coefficients) > 1e-7
        amplitude_error = np.max(
            np.abs(np.abs(final_coefficients[supported]) / np.abs(initial_coefficients[supported]) - 1.0)
        )
        phase_error = np.angle(final_coefficients[supported] / expected[supported])
        result.update(
            full_map_dimension=int(np.prod(direct.shape)),
            full_map_norm=float(np.vdot(direct, direct).real),
            modal_amplitude_error=float(amplitude_error),
            modal_phase_error_rad=float(np.max(np.abs(phase_error))),
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=480)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--full-map", action="store_true")
    args = parser.parse_args()
    result = run(args.steps, args.full_map)
    print(json.dumps(result, indent=2))
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
