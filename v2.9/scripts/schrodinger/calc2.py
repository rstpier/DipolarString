#!/usr/bin/env python3
"""Reconstructed Appendix-E calculation 2: a local potential on the defect band.

The quick path evolves the exact 30-mode projected operator.  ``--full-map``
also evolves the 185,856-component Johns state and measures leakage out of the
continued band.  Runtime for the full path is machine-dependent.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from calc1_common import (
    assemble_full_state,
    circular_moments,
    correlation,
    effective_band_map,
    envelope,
    full_map_step,
    gaussian_amplitudes,
    gaussian_profile,
    project_band,
    spatial_density,
    trace_bound_band,
)


def evolve_projected(operator: np.ndarray, amplitude: np.ndarray, steps: int) -> np.ndarray:
    state = amplitude.copy()
    for _ in range(steps):
        state = operator @ state
    return state


def delay_from_centroids(free_density: np.ndarray, perturbed_density: np.ndarray, center_hint: float, velocity: float) -> float:
    free_center, _ = circular_moments(free_density, center_hint)
    perturbed_center, _ = circular_moments(perturbed_density, free_center)
    return float(-(perturbed_center - free_center) / velocity)


def run(height: float, steps: int = 480, full_map: bool = False) -> dict[str, float | str]:
    band = trace_bound_band()
    z0 = 92.0
    profile = gaussian_profile(band.nz, center=64.0, width=6.0, height=height)
    amplitude0 = gaussian_amplitudes(band.k) * np.exp(1j * band.k * z0)
    free_amplitude = amplitude0 * np.exp(1j * band.omega * steps)
    operator = effective_band_map(band, profile)
    projected_amplitude = evolve_projected(operator, amplitude0, steps)
    free_env = envelope(free_amplitude, band.k, np.zeros_like(band.omega), band.nz)
    projected_env = envelope(projected_amplitude, band.k, np.zeros_like(band.omega), band.nz)
    velocity = -0.076
    center_hint = z0 + velocity * steps
    result: dict[str, float | str] = {
        "status": "RECONSTRUCTED, not original source recovery",
        "V0": height,
        "projected_norm": float(np.vdot(projected_amplitude, projected_amplitude).real),
        "projected_delay_steps": delay_from_centroids(
            np.abs(free_env) ** 2, np.abs(projected_env) ** 2, center_hint, velocity
        ),
    }

    if full_map:
        free_state = assemble_full_state(band, amplitude0, z0=0.0)
        potential_state = free_state.copy()
        for _ in range(steps):
            free_state = full_map_step(free_state)
            potential_state = full_map_step(potential_state, potential=profile)
        full_coefficients = project_band(potential_state, band)
        on_band = float(np.vdot(full_coefficients, full_coefficients).real / np.vdot(potential_state, potential_state).real)
        full_envelope = envelope(full_coefficients, band.k, np.zeros_like(band.omega), band.nz)
        free_density = spatial_density(free_state)
        perturbed_density = spatial_density(potential_state)
        result.update(
            on_band_metric_definition=(
                "projection onto the reconstructed 30-mode continued bound-band subspace; "
                "historical full-band projector not recovered"
            ),
            full_map_dimension=int(np.prod(potential_state.shape)),
            full_map_norm=float(np.vdot(potential_state, potential_state).real),
            on_band_fraction=on_band,
            full_map_delay_steps=delay_from_centroids(free_density, perturbed_density, center_hint, velocity),
            density_correlation_with_projected_model=correlation(
                np.abs(full_envelope) ** 2, np.abs(projected_env) ** 2
            ),
        )
    return result


def localized_state(height: float = 0.030) -> dict[str, float]:
    band = trace_bound_band()
    profile = gaussian_profile(band.nz, center=64.0, width=6.0, height=height)
    operator = effective_band_map(band, profile)
    eigenvalues, eigenvectors = np.linalg.eig(operator)
    z = np.arange(band.nz)
    distance = (z - 64.0 + band.nz / 2.0) % band.nz - band.nz / 2.0
    # The historical appendix reports its "weight in the well" without defining
    # the integration window.  We record the reconstructed window explicitly.
    well = np.abs(distance) <= 11.0
    best = None
    for index in range(len(eigenvalues)):
        amplitude = eigenvectors[:, index]
        psi = envelope(amplitude, band.k, np.zeros_like(band.omega), band.nz)
        weight = float(np.sum(np.abs(psi[well]) ** 2))
        if best is None or weight > best[0]:
            best = (weight, float(np.angle(eigenvalues[index])))
    assert best is not None
    return {
        "well_weight_radius_11": best[0],
        "projected_eigenfrequency": best[1],
        "reference_weight": 0.88,
        "reference_projected_frequency": 0.41823,
        "reference_full_map_frequency": 0.41819,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=float, default=0.015)
    parser.add_argument("--steps", type=int, default=480)
    parser.add_argument("--full-map", action="store_true")
    parser.add_argument("--localized", action="store_true")
    args = parser.parse_args()
    result = run(args.height, args.steps, args.full_map)
    if args.localized:
        result["localized_state_V0_0.030"] = localized_state()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
