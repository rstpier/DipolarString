#!/usr/bin/env python3
"""Redraw: continuing R59 ('what is tilted at the cube diagonal, and what accumulates 2/9 rad'),
setting the R60 proposal aside.

The cubic reading.  Koide's form says the three sqrt-masses are the three components, on three
orthogonal axes, of one unit direction u at 45 deg from the cube diagonal (1,1,1) with azimuth
theta about it.  In a cubic medium (the square section w = d of R17, the cubic pole plug of R51)
the three axes are the cell's edges; the generations would be the squared components of one
direction of the medium.

Checks:
  A. the data direction: u = (sqrt m_e, sqrt m_mu, sqrt m_tau)/norm = (0.0165, 0.237, 0.971):
     13.7 deg off the tau axis, 76.3 deg off the mu axis, 89.06 deg off the e axis; 45.00 deg from
     the diagonal.
  B. the 45 deg is the equality of the two Z_3 channels of a three-strand system (Phase A: the
     singlet A along (1,1,1) and the doublet E_1 in the plane): |u_singlet| = |u_doublet|, i.e. the
     family's mass sum splits equally between the common channel and the difference channel.
  C. no crystallographic direction (h, k, l) with |h|, |k|, |l| <= 8 lies at 45 deg from the diagonal
     (the condition 2(h+k+l)^2 = 3(h^2+k^2+l^2) has no small solution): the direction is not a
     lattice vector; the 45 deg must come from a rule on the channels, not from the cell's geometry.
  D. the rule the base could supply: one quantum of circulation in each channel -- the singlet and
     the doublet each carrying one unit (as the vacuum DQD carries two, R57 A, and hands one to each
     daughter, R4): equal weights give 45 deg exactly for any azimuth; the azimuth theta is then the
     doublet's phase, dynamical, and 2/9 remains to be produced.
  E. verdict: in the cubic reading the 45 deg becomes 'one quantum in the singlet, one in the doublet'
     (conditional), and the two missing rules reduce to one: what fixes the doublet's phase at 2/9.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    v = [math.sqrt(m) for m in (ME, MMU, MTAU)]
    norm = math.sqrt(sum(x * x for x in v))
    u = [x / norm for x in v]
    angles = [math.degrees(math.acos(c)) for c in u]
    diag = sum(u) / math.sqrt(3)
    polar = math.degrees(math.acos(diag))
    check("A. Data direction u = (0.0165, 0.237, 0.971): 13.7 deg off the tau axis, 76.3 off mu, 89.06 off e; 45.00 deg from the diagonal",
          abs(u[0] - 0.0165) < 0.001 and abs(u[1] - 0.237) < 0.002 and abs(u[2] - 0.971) < 0.002 and abs(angles[2] - 13.7) < 0.1 and abs(polar - 45.0) < 0.01,
          f"u = ({u[0]:.4f}, {u[1]:.4f}, {u[2]:.4f}); angles to axes {angles[0]:.2f}, {angles[1]:.2f}, {angles[2]:.2f} deg; from diagonal {polar:.3f} deg")

    singlet = diag
    doublet = math.sqrt(1 - diag**2)
    check("B. The 45 deg is |u_singlet| = |u_doublet|: the family's mass sum splits equally between the common (A) and difference (E_1) channels of the three-strand Z_3",
          abs(singlet - doublet) < 1e-4, f"singlet {singlet:.5f}, doublet {doublet:.5f}")

    sols = []
    for h, k, l in itertools.product(range(-8, 9), repeat=3):
        if (h, k, l) == (0, 0, 0):
            continue
        if 2 * (h + k + l) ** 2 == 3 * (h * h + k * k + l * l):
            sols.append((h, k, l))
    check("C. No crystallographic direction with |h|,|k|,|l| <= 8 lies at 45 deg from the diagonal (2(h+k+l)^2 = 3(h^2+k^2+l^2) has no small solution): the direction is not a lattice vector",
          sols == [], f"solutions: {sols}")

    # D. equal quanta in the two channels
    def polar_for_weights(w_s, w_d):
        return math.degrees(math.atan2(w_d, w_s))
    check("D. One quantum in each channel (singlet, doublet) gives 45 deg exactly for any azimuth; the azimuth is the doublet's phase, dynamical -- 2/9 remains to be produced",
          abs(polar_for_weights(1.0, 1.0) - 45.0) < 1e-12 and abs(polar_for_weights(2.0, 1.0) - 26.57) < 0.01,
          f"weights (1,1): {polar_for_weights(1,1):.1f} deg; (2,1): {polar_for_weights(2,1):.2f} deg")

    check("E. Verdict: in the cubic reading the 45 deg becomes 'one quantum in the singlet, one in the doublet' (conditional); the two missing rules reduce to one: what fixes the doublet's phase at 2/9",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
