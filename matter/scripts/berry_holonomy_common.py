"""Shared machinery: regularised TEM modes of n indistinguishable conductors on an n-gon.

E = sum_k q_k (rho - p_k) / (|rho - p_k|^2 + eps^2), eps = conductor radius r.  Smooth in the
n-gon angle psi, so Berry connections are well-defined; symmetries hold to grid precision.
D/r = 2 cosh(pi) is the DS matching condition; the Berry structure below does not depend on it.
"""

from __future__ import annotations

import math
import numpy as np

R = 1.0
EPS = R
D = 2 * math.cosh(math.pi) * R
N_GRID = 1401
L = 5 * D
XS = np.linspace(-L, L, N_GRID)
X, Y = np.meshgrid(XS, XS, indexing="ij")
DA = (XS[1] - XS[0]) ** 2

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    RESULTS.append(("PASS" if condition else "FAIL", name, detail))


def report(conclusion: str) -> int:
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print()
    print(f"RESULT: {len(RESULTS) - len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    print(conclusion)
    return 1 if failed else 0


def positions(n: int, psi: float):
    rho0 = D / 2 if n == 2 else D / (2 * math.sin(math.pi / n))   # n-gon of side D
    return [rho0 * np.array([math.cos(psi + 2 * math.pi * k / n), math.sin(psi + 2 * math.pi * k / n)])
            for k in range(n)]


def field(q, pts):
    ex, ey = np.zeros_like(X), np.zeros_like(X)
    for qi, p in zip(q, pts):
        dx, dy = X - p[0], Y - p[1]
        g = 1.0 / (dx * dx + dy * dy + EPS * EPS)
        ex += qi * dx * g
        ey += qi * dy * g
    return ex, ey


def ip(a, b) -> float:
    return float(np.sum(a[0] * b[0] + a[1] * b[1]) * DA)


def normalise(f):
    nrm = math.sqrt(ip(f, f))
    return (f[0] / nrm, f[1] / nrm)


def differential_charges(n: int):
    """Orthonormal charge patterns with zero total charge (the differential modes)."""
    if n == 2:
        return [np.array([1.0, -1.0]) / math.sqrt(2)]
    if n == 3:
        return [np.array([1.0, -1.0, 0.0]) / math.sqrt(2), np.array([1.0, 1.0, -2.0]) / math.sqrt(6)]
    raise ValueError(n)


def modes(n: int, psi: float):
    pts = positions(n, psi)
    return [normalise(field(q, pts)) for q in differential_charges(n)]


def common_mode(n: int, psi: float):
    return normalise(field(np.ones(n) / math.sqrt(n), positions(n, psi)))


def berry_connection(n: int, psi: float, h: float = 1e-4) -> np.ndarray:
    """A_ij(psi) = <E_i | d_psi E_j>, central difference on smooth fields."""
    ep, em, e0 = modes(n, psi + h), modes(n, psi - h), modes(n, psi)
    k = len(e0)
    return np.array([[ip(e0[i], ((ep[j][0] - em[j][0]) / (2 * h), (ep[j][1] - em[j][1]) / (2 * h)))
                      for j in range(k)] for i in range(k)])


def transport(n: int, psi: float) -> np.ndarray:
    """Overlap matrix U_ij = <E_i(0) | E_j(psi)>: exact transport of the differential subspace."""
    e0, e1 = modes(n, 0.0), modes(n, psi)
    k = len(e0)
    return np.array([[ip(e0[i], e1[j]) for j in range(k)] for i in range(k)])
