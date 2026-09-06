#!/usr/bin/env python3
"""Constructive reconstruction of the Johns symmetric condensed node (SCN).

This supplementary script closes the artifact named in Theorem 4 of the V2.9.1
manuscript.  It is a reconstruction dated 31 August 2026, not recovery of the
historical source file with the same name.

It performs, in one deterministic run:
  1. exhaustive enumeration of the 2^12 allowed sign motifs of the structural
     SCN ansatz;
  2. selection of symmetric, unitary candidates;
  3. full proper-cubic (24-rotation) equivariance filtering;
  4. a deterministic 50-direction low-q Bloch scan;
  5. comparison of the unique isotropic candidate with S_johns.npy.

The 50-direction criterion used here is the one asserted in the manuscript:
two degenerate transverse acoustic modes with v = omega/|q| = 1/2 to absolute
accuracy better than 1e-6 in cell units.  The scan is performed at |q|=0.005,
where finite-lattice corrections remain below that tolerance.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PORTS = [(d, s, p) for d in range(3) for s in (+1, -1) for p in range(3) if p != d]
INDEX = {port: index for index, port in enumerate(PORTS)}
N = len(PORTS)
Q_MAG = 0.005
N_DIRECTIONS = 50
SEED = 291
ISOTROPY_ATOL = 1.0e-6
DEGENERACY_ATOL = 1.0e-9


def third_axis(d: int, p: int) -> int:
    return 3 - d - p


def build_candidate(choice: tuple[int, ...]) -> np.ndarray:
    """Build one sign branch of the structural Johns-node ansatz."""
    matrix = np.zeros((N, N), dtype=float)
    for row, (d, _side, p) in enumerate(PORTS):
        t = third_axis(d, p)
        matrix[row, INDEX[(t, +1, p)]] = 0.5
        matrix[row, INDEX[(t, -1, p)]] = 0.5
        matrix[row, INDEX[(p, +1, d)]] = -0.5 if choice[row] == 0 else 0.5
        matrix[row, INDEX[(p, -1, d)]] = -0.5 if choice[row] == 1 else 0.5
    return matrix


def enumerate_symmetric_unitary() -> list[np.ndarray]:
    candidates: list[np.ndarray] = []
    for choice in itertools.product((0, 1), repeat=N):
        matrix = build_candidate(choice)
        if np.allclose(matrix, matrix.T, atol=1e-14) and np.allclose(
            matrix.T @ matrix, np.eye(N), atol=1e-14
        ):
            candidates.append(matrix)
    return candidates


def proper_cubic_rotations() -> list[np.ndarray]:
    """Return the 24 proper signed-permutation rotations of the cube."""
    rotations: list[np.ndarray] = []
    for permutation in itertools.permutations(range(3)):
        perm = np.eye(3, dtype=int)[:, permutation]
        for signs in itertools.product((-1, +1), repeat=3):
            rotation = perm @ np.diag(signs)
            if round(np.linalg.det(rotation)) == 1:
                rotations.append(rotation.astype(int))
    assert len(rotations) == 24
    return rotations


def port_representation(rotation: np.ndarray) -> np.ndarray:
    """Representation of a spatial rotation on directional/polarization ports."""
    transform = np.zeros((N, N), dtype=float)
    for column, (direction, side, polarization) in enumerate(PORTS):
        direction_vector = np.zeros(3, dtype=int)
        direction_vector[direction] = side
        polarization_vector = np.zeros(3, dtype=int)
        polarization_vector[polarization] = 1

        direction_rotated = rotation @ direction_vector
        polarization_rotated = rotation @ polarization_vector

        direction_new = int(np.argmax(np.abs(direction_rotated)))
        side_new = int(direction_rotated[direction_new])
        polarization_new = int(np.argmax(np.abs(polarization_rotated)))
        polarization_sign = int(polarization_rotated[polarization_new])

        transform[INDEX[(direction_new, side_new, polarization_new)], column] = polarization_sign
    return transform


def cubic_equivariant(matrix: np.ndarray, representations: list[np.ndarray]) -> bool:
    return all(
        np.allclose(rep @ matrix @ rep.T, matrix, atol=1e-14)
        for rep in representations
    )


def bloch_frequencies(matrix: np.ndarray, q_vector: np.ndarray) -> np.ndarray:
    propagation = np.zeros((N, N), dtype=complex)
    for column, (direction, side, polarization) in enumerate(PORTS):
        propagation[INDEX[(direction, -side, polarization)], column] = np.exp(
            1j * q_vector[direction] * side
        )
    return np.sort(np.angle(np.linalg.eigvals(propagation @ matrix)))


def acoustic_speeds(matrix: np.ndarray, direction: np.ndarray) -> np.ndarray:
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    omega = bloch_frequencies(matrix, Q_MAG * direction)
    positive = omega[(omega > 1e-10) & (omega < 1.0)]
    if positive.size < 2:
        return np.array([], dtype=float)
    return positive[:2] / Q_MAG


def direction_set() -> np.ndarray:
    rng = np.random.default_rng(SEED)
    directions = rng.normal(size=(N_DIRECTIONS, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    return directions


def scan_isotropy(matrix: np.ndarray, directions: np.ndarray) -> dict[str, float | bool]:
    speeds = []
    for direction in directions:
        value = acoustic_speeds(matrix, direction)
        if value.size != 2:
            return {
                "pass": False,
                "max_abs_speed_error": float("inf"),
                "max_polarization_split": float("inf"),
                "v_min": float("nan"),
                "v_max": float("nan"),
            }
        speeds.append(value)

    array = np.asarray(speeds)
    max_error = float(np.max(np.abs(array - 0.5)))
    max_split = float(np.max(np.abs(array[:, 0] - array[:, 1])))
    return {
        "pass": max_error < ISOTROPY_ATOL and max_split < DEGENERACY_ATOL,
        "max_abs_speed_error": max_error,
        "max_polarization_split": max_split,
        "v_min": float(np.min(array)),
        "v_max": float(np.max(array)),
    }


def main() -> int:
    candidates = enumerate_symmetric_unitary()
    rotations = proper_cubic_rotations()
    representations = [port_representation(rotation) for rotation in rotations]
    cubic = [matrix for matrix in candidates if cubic_equivariant(matrix, representations)]

    directions = direction_set()
    scans = [scan_isotropy(matrix, directions) for matrix in cubic]
    isotropic = [matrix for matrix, result in zip(cubic, scans) if bool(result["pass"])]

    archived = np.load(ROOT / "S_johns.npy")
    archived_match = len(isotropic) == 1 and np.allclose(isotropic[0], archived, atol=1e-14)

    print("DS V2.9.1 — reconstructed scn_from_dqd.py")
    print(f"sign motifs scanned               : {2**N}")
    print(f"symmetric + unitary candidates    : {len(candidates)}")
    print(f"proper cubic rotations checked    : {len(rotations)}")
    print(f"cubic-equivariant candidates      : {len(cubic)}")
    print(f"Bloch directions                  : {N_DIRECTIONS} (seed={SEED})")
    print(f"|q|                               : {Q_MAG}")
    for index, result in enumerate(scans, start=1):
        print(
            f"cubic candidate {index}: pass={result['pass']}  "
            f"v_range=[{result['v_min']:.12f}, {result['v_max']:.12f}]  "
            f"max|v-1/2|={result['max_abs_speed_error']:.3e}  "
            f"max polarization split={result['max_polarization_split']:.3e}"
        )
    print(f"unique 50-direction isotropic node: {len(isotropic) == 1}")
    print(f"equals archived S_johns.npy       : {archived_match}")

    success = (
        len(candidates) == 8
        and len(cubic) == 2
        and len(isotropic) == 1
        and archived_match
    )
    print("RESULT: " + ("PASS" if success else "FAIL"))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
