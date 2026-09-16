#!/usr/bin/env python3
"""Redraw: 'it is the junction mass: the more connections a junction has, the more inertia space
grants it' -- the author's law for why the same strings weigh 140 MeV as a ring and 938 MeV as a
star.  Tested in the base (BASE.md R2, R15, R19, R20, R22).

Checks:
  A. the base's own Coulomb law does give inertia growing with connections, linearly: a balanced
     k-junction (unit poles at a common distance) binds like floor(k/2) pairs; the odd pole adds
     nothing (k = 3 binds like k = 2, k = 5 like k = 4).  The law is DERIVED, but linear, and
     capped by frustration.
  B. single-centre reading (proton = one 5-junction, neutron = one 4-junction): the floor law makes
     the two centres equal, so p = n + one string + its Coulomb part: the proton would be heavier.
     Observed: the neutron is heavier by 1.29 MeV.  Excluded by the sign.
  C. the content that gives the same graph: u = pair (2+, one k = 2 junction), d = Y (1+, 2-, one
     k = 3 centre).  p = uud = 7 strings (5+, 2-), n = udd = 8 (4+, 4-): charges and parities right,
     and the pion rings of R15 recovered: pi+ = u dbar = (4+, 1-), pi- = d ubar = (1+, 4-),
     pi0 = u ubar = (2+, 2-).
  D. sign of n - p: n has one more string and a k = 3 centre in place of a k = 2 junction (equal
     under the floor law): heavier, as observed; the extra string weighs 1.29 + 1.0 (Coulomb) =
     2.3 MeV at the nucleon's scale.
  E. where the GeV sits: the pion ring has only k = 2 junctions and weighs its mode (140 MeV); the
     baryon's ~930 MeV must be its centre.  The d's own Y-centre has the same valence and weighs
     < 2.3 MeV: same k, 400 x different mass.  By the linear law at 2 MeV per pair the centre would
     need ~930 connections; with unit poles at fm scale it is impossible, the Coulomb junction reaches
     930 MeV only at s = 1.5 am (matched poles: 21 am).  Valence is not the variable, scale is.
  F. the pattern in the data: light mesons reach the proton's mass without any centre (eta' 958,
     a0 980) and rho/pi = 5.55 at equal content and equal ring topology, close to p/pi+ = 6.72:
     a factor 5-7 arises inside one topology, so the ring-star ratio does not measure a junction.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction as F

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
MP, MN, ME = 938.27209, 939.56542, 0.51099895
MPIP, MPI0, MRHO, META_P, MA0 = 139.57039, 134.9768, 775.26, 957.78, 980.0
DELTA_MATCHED = 1 / (math.pi * math.sqrt(ALPHA))
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def best_junction(k: int) -> float:
    """Most negative Coulomb energy of k unit poles at mutual distance 1 fm, over all sign choices."""
    best = math.inf
    for signs in itertools.product((+1, -1), repeat=k):
        e = sum(K * a * b for a, b in itertools.combinations(signs, 2))
        best = min(best, e)
    return best


def main() -> int:
    # A. frustration
    ej = {k: best_junction(k) for k in range(2, 8)}
    check("A. Base Coulomb law: a balanced k-junction binds like floor(k/2) pairs (k = 2..7): inertia does grow with connections, linearly, and the odd pole adds nothing (k = 3 = k = 2, k = 5 = k = 4). Derived, linear, capped by frustration",
          all(abs(ej[k] - (k // 2) * ej[2]) < 1e-9 for k in ej),
          "; ".join(f"k={k}: {e:.3f} MeV = {e/ej[2]:.0f} pair(s)" for k, e in ej.items()) + " (unit poles, 1 fm)")

    # B. single centre
    check("B. Single-centre reading (p = one 5-junction, n = one 4-junction): the floor law makes the two centres equal, so p = n + one string + its Coulomb part, and the proton would be heavier; observed n - p = +1.29 MeV: excluded by the sign",
          abs(ej[5] - ej[4]) < 1e-9 and MN > MP, f"m_j(5) = m_j(4) = {ej[4]:.3f} MeV; m_n - m_p = {MN-MP:+.3f} MeV")

    # C. content
    u, d = (F(2), F(0)), (F(1), F(2))                 # (n+, n-)
    ubar, dbar = (u[1], u[0]), (d[1], d[0])
    def add(*qs): return (sum(q[0] for q in qs), sum(q[1] for q in qs))
    def q(c): return (c[0] - c[1]) / 3
    p, n = add(u, u, d), add(u, d, d)
    pip, pim, pi0 = add(u, dbar), add(d, ubar), add(u, ubar)
    ok = (q(p) == 1 and q(n) == 0 and sum(p) % 2 == 1 and sum(n) % 2 == 0
          and pip == (4, 1) and pim == (1, 4) and pi0 == (2, 2) and q(u) == F(2, 3) and q(d) == F(-1, 3))
    check("C. u = pair (2+), d = Y (1+, 2-): p = uud = 7 strings (5+, 2-), n = udd = 8 (4+, 4-); pi+ = (4+, 1-), pi- = (1+, 4-), pi0 = (2+, 2-) -- the rings of R15 recovered from the quark content",
          ok, f"p = ({int(p[0])}+, {int(p[1])}-) q = {q(p)}, n = ({int(n[0])}+, {int(n[1])}-) q = {q(n)}; pi+ = ({int(pip[0])}+, {int(pip[1])}-), pi- = ({int(pim[0])}+, {int(pim[1])}-), pi0 = ({int(pi0[0])}+, {int(pi0[1])}-)")

    # D. sign and budget of n - p
    budget = (MN - MP) + 1.00
    check("D. n has one string more than p and a k = 3 centre in place of a k = 2 junction (equal under the floor law): heavier, as observed (+1.29 MeV); with the proton's Coulomb part (~1.0 MeV) the extra string weighs 2.3 MeV at the nucleon's scale",
          MN > MP and abs(ej[3] - ej[2]) < 1e-9 and abs(budget - 2.29) < 0.01, f"m_n - m_p = {MN-MP:.3f} MeV; + 1.00 MeV Coulomb = {budget:.2f} MeV per string")

    # E. where the GeV sits
    centre = MP - 2.0 * 2 - 2.3                        # rough: all but small junctions and one string
    k_needed = 2 * centre / 2.0
    s_unit = K / centre
    s_matched = DELTA_MATCHED**2 * K / centre
    check("E. Pion ring (k = 2 only) = its mode, 140 MeV; the baryon's ~930 MeV must sit in its centre, while the d's own Y-centre (same valence) weighs < 2.3 MeV: same k, 400 x the mass. The linear law at 2 MeV per pair needs ~930 connections; a Coulomb junction reaches 930 MeV only at s = 1.5 am (unit poles) or 21 am (matched): valence is not the variable, scale is",
          centre / 2.3 > 400 and abs(k_needed - 932) < 2 and abs(s_unit * 1e3 - 1.55) < 0.05 and abs(s_matched * 1e3 - 21.5) < 0.5,
          f"centre ~ {centre:.0f} MeV = {centre/2.3:.0f} x 2.3 MeV; connections needed at 2 MeV/pair: {k_needed:.0f}; s(unit) = {s_unit*1e3:.2f} am, s(matched) = {s_matched*1e3:.1f} am")

    # F. pattern in the data
    check("F. Light mesons reach the proton's mass with no centre (eta' 958, a0 980); rho/pi = 5.55 at equal content and equal ring topology vs p/pi+ = 6.72: a factor 5-7 arises inside one topology, so the ring-star ratio does not measure a junction",
          META_P > MP and MA0 > MP and abs(MRHO / MPIP - 5.55) < 0.01 and abs(MP / MPIP - 6.72) < 0.01,
          f"eta' = {META_P:.0f}, a0 = {MA0:.0f} > m_p; rho/pi+ = {MRHO/MPIP:.2f}, p/pi+ = {MP/MPIP:.2f}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
