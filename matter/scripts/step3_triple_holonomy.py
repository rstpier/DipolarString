#!/usr/bin/env python3
"""Step 3a: three indistinguishable conductors (N_DS = 3, the electron's string count).

Unlike the pair, the differential modes of an equilateral triple form a DEGENERATE DOUBLET
(the E irrep of C3v).  A degenerate doublet is exactly what a continuous Berry connection
needs.  Checks: degeneracy, orthogonality, a non-zero antisymmetric (so(2)) connection
constant in psi, the exact holonomy over the closed loop psi = 2pi/3 of the unordered triple,
and the contrast with the pair.

CORRECTION (step 3b, scripts/step3b_triple_curvature.py): the "holonomy over 2pi/3" measured
here is the overlap of the doublet with itself after a third of a turn -- the PERMUTATION part
of the holonomy (the three conductors exchanged).  The full holonomy of a closed string
composes it with parallel transport by the connection A measured separately below; for a
planar ring twisted by 1/3 it is (1 - a) x 120 deg = 97.7 deg, and exactly 120 deg only when the
triple closes by writhe with zero twist.  The labels below are kept as run on 15 September;
read "holonomy" as "permutation part".

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import (berry_connection, check, common_mode, differential_charges,
                                   field, ip, positions, report, transport)


def block_degeneracy() -> None:
    raw = [field(q, positions(3, 0.0)) for q in differential_charges(3)]
    g = np.array([[ip(a, b) for b in raw] for a in raw])
    rel = abs(g[0, 0] - g[1, 1]) / g[0, 0]
    off = abs(g[0, 1]) / g[0, 0]
    check("Differential doublet is degenerate", rel < 1e-6, f"W_aa = W_bb to {rel:.1e} relative (E irrep of C3v)")
    check("Differential doublet is orthogonal", off < 1e-6, f"|W_ab|/W_aa = {off:.1e}")
    wc = ip(field(np.ones(3) / math.sqrt(3), positions(3, 0.0)), field(np.ones(3) / math.sqrt(3), positions(3, 0.0)))
    check("Common mode is non-degenerate with the doublet", abs(wc - g[0, 0]) / g[0, 0] > 0.5,
          f"W_common / W_diff = {wc/g[0,0]:.3f}")


def block_connection() -> None:
    mats = [berry_connection(3, psi) for psi in (0.0, 0.5, 1.1)]
    diag = max(np.max(np.abs(np.diag(m))) for m in mats)
    asym = max(abs(m[0, 1] + m[1, 0]) for m in mats)
    vals = [m[0, 1] for m in mats]
    const = max(abs(v - vals[0]) for v in vals)
    check("Non-abelian Berry connection is antisymmetric (so(2))", diag < 1e-9 and asym < 1e-5,
          f"max |diag| = {diag:.1e}, max |A_ab + A_ba| = {asym:.1e}")
    check("Berry connection is NON-ZERO (unlike the pair)", abs(vals[0]) > 0.1,
          f"A_ab = {vals[0]:+.6f} per unit psi")
    check("Berry connection is constant in psi (C3 symmetry)", const < 1e-5,
          f"spread over three angles = {const:.1e}")


def block_z3_holonomy() -> None:
    u = transport(3, 2 * math.pi / 3)
    ev = np.linalg.eigvals(u)
    target = np.exp(2j * math.pi / 3)
    dev = min(abs(ev[0] - target), abs(ev[1] - target))
    ang = math.degrees(math.atan2(u[1, 0], u[0, 0]))
    check("Overlap (permutation part of the holonomy) over the closed loop 2pi/3 is a rotation by 120 deg", abs(abs(ang) - 120) < 1e-2 and abs(np.linalg.det(u) - 1) < 1e-6,
          f"rotation = {ang:+.3f} deg, det U = {np.linalg.det(u):.6f}")
    check("Its eigenvalues are exp(+/- 2 pi i / 3): fractional winding 1/3 on e+/- (permutation part)", dev < 1e-5,
          f"eigenvalues {np.round(ev, 6)}; deviation from exp(2pi i/3) = {dev:.1e}")
    u1 = transport(3, 2 * math.pi)
    check("Full turn is the identity", np.allclose(u1, np.eye(2), atol=1e-6),
          f"U(2pi) = {np.round(u1, 6).tolist()}")
    cm = ip(common_mode(3, 0.0), common_mode(3, 2 * math.pi / 3))
    check("Common mode carries no phase", abs(cm - 1) < 1e-9, f"<Ec(0)|Ec(2pi/3)> = {cm:+.9f}")


def block_contrast() -> None:
    a2 = abs(berry_connection(2, 0.7)[0, 0])
    a3 = abs(berry_connection(3, 0.7)[0, 1])
    check("Pair: A = 0; triple: A != 0 -- the degenerate doublet makes the difference", a2 < 1e-9 and a3 > 0.1,
          f"|A| pair = {a2:.1e}, |A_ab| triple = {a3:.4f}")
    z2 = transport(2, math.pi)[0, 0]
    check("Z_n holonomy: n = 2 gives -1 (1/2), n = 3 gives exp(2pi i/3) (1/3)", abs(z2 + 1) < 1e-9,
          "n indistinguishable conductors on an n-gon -> fractional winding 1/n; the G conjecture takes 1/3 as input")


def main() -> int:
    block_degeneracy()
    block_connection()
    block_z3_holonomy()
    block_contrast()
    return report("Conclusion (n = 3): the differential doublet carries a NON-ZERO Berry connection, and the "
                  "permutation of the conductors over 2pi/3 acts on it as an exact Z3 rotation. The two are "
                  "composed, and the curvature computed, in step 3b (step3b_triple_curvature.py): full holonomy "
                  "2 pi [(1 - a) Lk + a Wr], curvature a x (area form of the tangent), a = 0.186.")


if __name__ == "__main__":
    raise SystemExit(main())
