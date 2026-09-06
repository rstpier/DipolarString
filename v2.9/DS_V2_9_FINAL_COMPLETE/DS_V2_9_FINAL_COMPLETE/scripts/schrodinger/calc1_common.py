#!/usr/bin/env python3
"""Shared Johns-lattice machinery for the reconstructed Appendix-E calculations.

This file was reconstructed on 2026-08-31 from the archived V2.8 session code
(`defect.py`, `track.py`) and from the numerical protocol printed in Appendix E.
It is not represented as a byte-for-byte recovery of the lost original.

Conventions follow the archived code: a Bloch eigenvalue is exp(+i*omega), and
a positive-z hop carries exp(+i*k).  A real-space Bloch state is therefore
proportional to exp(-i*k*z).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
S_PATH = HERE.parent / "johns" / "S_johns.npy"
BAND_CACHE = HERE / "band_nx11_nz128_reconstructed.npz"

NX = 11
NZ = 128
K0 = 1.5
SIGMA_K = 0.15
K_INDICES = np.arange(19, 49)  # exactly 30 ring momenta, 0.933 <= k <= 2.356
PORTS = [(d, s, p) for d in range(3) for s in (+1, -1) for p in range(3) if p != d]
PIDX = {port: index for index, port in enumerate(PORTS)}


@dataclass(frozen=True)
class BandData:
    k: np.ndarray
    omega: np.ndarray
    core_weight: np.ndarray
    modes: np.ndarray
    nx: int
    nz: int


def load_scattering() -> np.ndarray:
    matrix = np.load(S_PATH)
    if matrix.shape != (12, 12) or not np.allclose(matrix.T @ matrix, np.eye(12), atol=1e-13):
        raise ValueError(f"invalid Johns matrix: {S_PATH}")
    return matrix


def core_columns(nx: int = NX, handedness: int = +1) -> list[tuple[int, int]]:
    c = nx // 2
    positive = [(c, c), (c + 1, c), (c, c + 1)]
    return positive if handedness >= 0 else [positive[0], positive[2], positive[1]]


def core_mask(nx: int, core: list[tuple[int, int]]) -> np.ndarray:
    mask = np.zeros(12 * nx * nx, dtype=bool)
    for ix, iy in core:
        start = 12 * ((ix % nx) * nx + (iy % nx))
        mask[start : start + 12] = True
    return mask


def bloch_operator(nx: int, k: float, handedness: int = +1) -> sp.csc_matrix:
    """Return the sparse transverse supercell map U(k) for the D3 braid."""
    scattering = load_scattering()
    core = core_columns(nx, handedness)
    core_set = set(core)
    n_nodes = nx * nx
    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []

    def gi(ix: int, iy: int, port: tuple[int, int, int]) -> int:
        return 12 * ((ix % nx) * nx + (iy % nx)) + PIDX[port]

    def add(row: int, col: int, value: complex) -> None:
        rows.append(row)
        cols.append(col)
        vals.append(value)

    theta = handedness * 2.0 * np.pi / 3.0
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]], dtype=float
    )
    for ix in range(nx):
        for iy in range(nx):
            for polarization in (1, 2):
                add(gi(ix + 1, iy, (0, -1, polarization)), gi(ix, iy, (0, +1, polarization)), 1.0)
                add(gi(ix - 1, iy, (0, +1, polarization)), gi(ix, iy, (0, -1, polarization)), 1.0)
            for polarization in (0, 2):
                add(gi(ix, iy + 1, (1, -1, polarization)), gi(ix, iy, (1, +1, polarization)), 1.0)
                add(gi(ix, iy - 1, (1, +1, polarization)), gi(ix, iy, (1, -1, polarization)), 1.0)
            if (ix, iy) in core_set:
                j = core.index((ix, iy))
                up = core[(j + 1) % 3]
                down = core[(j - 1) % 3]
                for arrival_pol in range(2):
                    for departure_pol in range(2):
                        add(
                            gi(up[0], up[1], (2, -1, arrival_pol)),
                            gi(ix, iy, (2, +1, departure_pol)),
                            rotation[arrival_pol, departure_pol] * np.exp(1j * k),
                        )
                        add(
                            gi(down[0], down[1], (2, +1, arrival_pol)),
                            gi(ix, iy, (2, -1, departure_pol)),
                            rotation.T[arrival_pol, departure_pol] * np.exp(-1j * k),
                        )
            else:
                for polarization in (0, 1):
                    add(
                        gi(ix, iy, (2, -1, polarization)),
                        gi(ix, iy, (2, +1, polarization)),
                        np.exp(1j * k),
                    )
                    add(
                        gi(ix, iy, (2, +1, polarization)),
                        gi(ix, iy, (2, -1, polarization)),
                        np.exp(-1j * k),
                    )
    propagation = sp.csr_matrix((vals, (rows, cols)), shape=(12 * n_nodes, 12 * n_nodes))
    scatter = sp.block_diag([scattering] * n_nodes, format="csr")
    return (propagation @ scatter).tocsc()


def _align_phase(vector: np.ndarray, reference: np.ndarray | None) -> np.ndarray:
    vector = vector / np.linalg.norm(vector)
    if reference is None:
        pivot = int(np.argmax(np.abs(vector)))
        return vector * np.exp(-1j * np.angle(vector[pivot]))
    overlap = np.vdot(reference, vector)
    return vector * np.exp(-1j * np.angle(overlap))


def trace_bound_band(
    nx: int = NX,
    nz: int = NZ,
    cache: Path | None = BAND_CACHE,
    force: bool = False,
) -> BandData:
    """Trace the localized band by eigenvector continuation on 30 ring momenta."""
    if cache is not None and cache.exists() and not force:
        data = np.load(cache)
        return BandData(data["k"], data["omega"], data["core_weight"], data["modes"], int(data["nx"]), int(data["nz"]))

    indices = K_INDICES if nz == NZ else np.arange(np.ceil(0.9 * nz / (2 * np.pi)), np.floor(2.4 * nz / (2 * np.pi)) + 1, dtype=int)
    momenta = 2.0 * np.pi * indices / nz
    core = core_columns(nx)
    mask = core_mask(nx, core)
    records: list[tuple[float, float, float, np.ndarray]] = []
    previous: np.ndarray | None = None
    deterministic_v0 = np.ones(12 * nx * nx, dtype=complex)

    for k in momenta[::-1]:
        guess = float(np.interp(k, [0.9, 1.5, 2.0, 2.4], [0.4037, 0.3673, 0.33246, 0.31424]))
        operator = bloch_operator(nx, float(k))
        eigenvalues, eigenvectors = spla.eigs(
            operator,
            k=8,
            sigma=np.exp(1j * guess),
            which="LM",
            tol=1e-11,
            v0=deterministic_v0,
        )
        probabilities = np.abs(eigenvectors) ** 2
        probabilities /= probabilities.sum(axis=0, keepdims=True)
        weights = probabilities[mask].sum(axis=0)
        phases = np.angle(eigenvalues)
        if previous is None:
            eligible = (phases > 0.02) & (phases < k / 2.0)
            chosen = int(np.argmax(weights * eligible))
        else:
            chosen = int(np.argmax(np.abs(previous.conj() @ eigenvectors)))
        mode = _align_phase(eigenvectors[:, chosen], previous)
        previous = mode
        records.append((float(k), float(phases[chosen]), float(weights[chosen]), mode))

    records.reverse()
    result = BandData(
        k=np.array([item[0] for item in records]),
        omega=np.unwrap(np.array([item[1] for item in records])),
        core_weight=np.array([item[2] for item in records]),
        modes=np.stack([item[3] for item in records]),
        nx=nx,
        nz=nz,
    )
    if cache is not None:
        np.savez_compressed(
            cache,
            k=result.k,
            omega=result.omega,
            core_weight=result.core_weight,
            modes=result.modes,
            nx=result.nx,
            nz=result.nz,
            provenance="reconstructed from archived V2.8 session code on 2026-08-31",
        )
    return result


def local_derivatives(k: np.ndarray, omega: np.ndarray, degree: int = 4, half_window: int = 3) -> tuple[np.ndarray, np.ndarray]:
    velocity = np.empty_like(omega)
    curvature = np.empty_like(omega)
    for index in range(len(k)):
        lo = max(0, index - half_window)
        hi = min(len(k), index + half_window + 1)
        order = min(degree, hi - lo - 1)
        polynomial = np.poly1d(np.polyfit(k[lo:hi] - k[index], omega[lo:hi], order))
        velocity[index] = np.polyder(polynomial, 1)(0.0)
        curvature[index] = np.polyder(polynomial, 2)(0.0)
    return velocity, curvature


def gaussian_amplitudes(k: np.ndarray, center: float = K0, sigma: float = SIGMA_K) -> np.ndarray:
    """Gaussian amplitude; sigma is the standard deviation of amplitude, as in Appendix E."""
    amplitude = np.exp(-0.5 * ((k - center) / sigma) ** 2).astype(complex)
    return amplitude / np.linalg.norm(amplitude)


def envelope(amplitude: np.ndarray, k: np.ndarray, omega: np.ndarray, nz: int, t: float = 0.0, z0: float = 0.0) -> np.ndarray:
    z = np.arange(nz)
    wave = np.exp(-1j * np.outer(k, z - z0) + 1j * omega[:, None] * t)
    psi = amplitude @ wave
    return psi / np.linalg.norm(psi)


def circular_moments(probability: np.ndarray, center_hint: float) -> tuple[float, float]:
    probability = np.asarray(probability, dtype=float)
    probability = probability / probability.sum()
    n = probability.size
    z = np.arange(n, dtype=float)
    displacement = (z - center_hint + n / 2.0) % n - n / 2.0
    mean_displacement = float(probability @ displacement)
    center = center_hint + mean_displacement
    centered = displacement - mean_displacement
    variance = float(probability @ centered**2)
    return center, variance


def assemble_full_state(band: BandData, amplitude: np.ndarray, z0: float) -> np.ndarray:
    state = np.zeros((band.nx, band.nx, band.nz, 12), dtype=complex)
    z = np.arange(band.nz)
    for coefficient, k, mode in zip(amplitude, band.k, band.modes):
        transverse = mode.reshape(band.nx, band.nx, 12)
        phase = np.exp(-1j * k * (z - z0)) / np.sqrt(band.nz)
        state += coefficient * transverse[:, :, None, :] * phase[None, None, :, None]
    return state / np.linalg.norm(state)


def full_map_step(
    state: np.ndarray,
    scattering: np.ndarray | None = None,
    potential: np.ndarray | None = None,
    gauge_phase: float = 0.0,
) -> np.ndarray:
    """One exact scatter-propagate step on the Nx x Nx x Nz braided lattice."""
    nx, ny, nz, channels = state.shape
    if nx != ny or channels != 12:
        raise ValueError("state must have shape (Nx, Nx, Nz, 12)")
    scattering = load_scattering() if scattering is None else scattering
    outgoing = np.einsum("ij,xyzj->xyzi", scattering, state, optimize=True)
    incoming = np.zeros_like(outgoing)

    for direction in (0, 1):
        axis = direction
        polarizations = (1, 2) if direction == 0 else (0, 2)
        for side in (+1, -1):
            for polarization in polarizations:
                source = PIDX[(direction, side, polarization)]
                target = PIDX[(direction, -side, polarization)]
                incoming[..., target] = np.roll(outgoing[..., source], shift=side, axis=axis)

    core = core_columns(nx)
    core_set = set(core)
    normal = [(ix, iy) for ix in range(nx) for iy in range(nx) if (ix, iy) not in core_set]
    for ix, iy in normal:
        for side in (+1, -1):
            phase = np.exp(1j * side * gauge_phase)
            for polarization in (0, 1):
                source = PIDX[(2, side, polarization)]
                target = PIDX[(2, -side, polarization)]
                incoming[ix, iy, :, target] = phase * np.roll(outgoing[ix, iy, :, source], shift=side)

    theta = 2.0 * np.pi / 3.0
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    for j, (ix, iy) in enumerate(core):
        up = core[(j + 1) % 3]
        down = core[(j - 1) % 3]
        plus = np.stack([outgoing[ix, iy, :, PIDX[(2, +1, p)]] for p in (0, 1)], axis=-1)
        minus = np.stack([outgoing[ix, iy, :, PIDX[(2, -1, p)]] for p in (0, 1)], axis=-1)
        plus_arrival = np.roll(plus @ rotation.T, shift=+1, axis=0) * np.exp(1j * gauge_phase)
        minus_arrival = np.roll(minus @ rotation, shift=-1, axis=0) * np.exp(-1j * gauge_phase)
        for p in (0, 1):
            incoming[up[0], up[1], :, PIDX[(2, -1, p)]] = plus_arrival[:, p]
            incoming[down[0], down[1], :, PIDX[(2, +1, p)]] = minus_arrival[:, p]

    if potential is not None:
        if potential.shape != (nz,):
            raise ValueError("potential must have shape (Nz,)")
        incoming *= np.exp(1j * potential)[None, None, :, None]
    return incoming


def project_band(state: np.ndarray, band: BandData) -> np.ndarray:
    """Project a full 3-D state onto the 30 continued Bloch modes."""
    z = np.arange(band.nz)
    coefficients = np.empty(len(band.k), dtype=complex)
    for i, (k, mode) in enumerate(zip(band.k, band.modes)):
        transverse = np.tensordot(state, np.exp(1j * k * z) / np.sqrt(band.nz), axes=([2], [0]))
        coefficients[i] = np.vdot(mode, transverse.reshape(-1))
    return coefficients


def spatial_density(state: np.ndarray) -> np.ndarray:
    density = np.sum(np.abs(state) ** 2, axis=(0, 1, 3))
    return density / density.sum()


def gaussian_profile(nz: int, center: float, width: float, height: float) -> np.ndarray:
    z = np.arange(nz, dtype=float)
    distance = (z - center + nz / 2.0) % nz - nz / 2.0
    return height * np.exp(-0.5 * (distance / width) ** 2)


def effective_band_map(band: BandData, potential: np.ndarray) -> np.ndarray:
    """Projected one-step map P_band exp(iV) U P_band in the continued band basis."""
    z = np.arange(band.nz)
    phase = np.exp(1j * potential)
    overlap = band.modes.conj() @ band.modes.T
    fourier = np.empty((len(band.k), len(band.k)), dtype=complex)
    for i, ki in enumerate(band.k):
        for j, kj in enumerate(band.k):
            fourier[i, j] = np.mean(phase * np.exp(1j * (ki - kj) * z))
    return overlap * fourier * np.exp(1j * band.omega)[None, :]


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    a0 = np.asarray(a, float) - np.mean(a)
    b0 = np.asarray(b, float) - np.mean(b)
    return float(np.dot(a0, b0) / np.sqrt(np.dot(a0, a0) * np.dot(b0, b0)))


def cli_band() -> None:
    parser = argparse.ArgumentParser(description="Rebuild/cache the Appendix-E bound band")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    band = trace_bound_band(force=args.force)
    print(f"30-point band: k=[{band.k[0]:.6f}, {band.k[-1]:.6f}]")
    print(f"omega(k nearest 1.5)={band.omega[np.argmin(abs(band.k-K0))]:.8f}")
    print(f"core weight range=[{band.core_weight.min():.4f}, {band.core_weight.max():.4f}]")
    print(f"cache: {BAND_CACHE.name}")


if __name__ == "__main__":
    cli_band()
