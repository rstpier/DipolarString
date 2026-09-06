#!/usr/bin/env python3
"""Reconstructed Appendix-E calculation 1b: two-packet interference."""

from __future__ import annotations

import argparse
import json

import numpy as np

from calc1_common import correlation, gaussian_amplitudes, trace_bound_band


def packet_state(band, center: float, sigma: float, t: float, z0: float) -> tuple[np.ndarray, np.ndarray]:
    amplitude = gaussian_amplitudes(band.k, center=center, sigma=sigma)
    amplitude = amplitude * np.exp(1j * band.omega * t)
    z = np.arange(band.nz)
    state = np.zeros((band.nx, band.nx, band.nz, 12), dtype=complex)
    for coefficient, k, mode in zip(amplitude, band.k, band.modes):
        phase = np.exp(-1j * k * (z - z0)) / np.sqrt(band.nz)
        state += coefficient * mode.reshape(band.nx, band.nx, 12)[:, :, None, :] * phase[None, None, :, None]
    envelope = amplitude @ np.exp(-1j * np.outer(band.k, z - z0))
    return state, envelope


def run(steps: int = 480, sigma: float = 0.10) -> dict[str, float | str]:
    band = trace_bound_band()
    state1, env1 = packet_state(band, 1.3, sigma, steps, 64.0)
    state2, env2 = packet_state(band, 1.7, sigma, steps, 64.0)
    p12 = np.sum(np.abs(state1 + state2) ** 2, axis=(0, 1, 3))
    p1 = np.sum(np.abs(state1) ** 2, axis=(0, 1, 3))
    p2 = np.sum(np.abs(state2) ** 2, axis=(0, 1, 3))
    cross = p12 - p1 - p2
    predicted = 2.0 * np.real(env1 * np.conj(env2)) / band.nz
    spectrum = np.abs(np.fft.fft(cross - cross.mean()))
    index = int(np.argmax(spectrum[1 : band.nz // 2]) + 1)
    peak = 2.0 * np.pi * index / band.nz
    return {
        "status": "RECONSTRUCTED, not original source recovery",
        "measured_peak": peak,
        "expected_delta_k": 0.4,
        "ring_bin": 2.0 * np.pi / band.nz,
        "correlation_with_complex_envelope": correlation(cross, predicted),
        "reference_peak": 0.3926990817,
        "reference_correlation": 0.99962,
        "interpretation": (
            "linear-superposition/phase/Fourier-bin consistency check; the existence of the "
            "cross term itself is algebraic and is not a Born-rule test"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=480)
    parser.add_argument("--sigma", type=float, default=0.10)
    args = parser.parse_args()
    print(json.dumps(run(args.steps, args.sigma), indent=2))


if __name__ == "__main__":
    main()
