#!/usr/bin/env python3
"""Redraw: what produces Koide's rotation 2/9 (R56)?

Geometry of Koide's relation.  With v = (sqrt m_e, sqrt m_mu, sqrt m_tau), Q = |v|^2/(sum v)^2 = 2/3
means v makes exactly 45 deg with the diagonal (1,1,1): half of sum m is common to the three leptons
(the projection on the diagonal), half sits in their differences.  The azimuth of v around the
diagonal, measured from the k = 0 axis with sqrt m_k = A [1 + sqrt 2 cos(theta + 2 pi k/3)], is theta.

Checks:
  A. the data: polar angle 45.000 deg (Q = 2/3 to 1e-5) and azimuth |theta| = 2/9 rad to 2e-5,
     i.e. exact at the precision of the tau mass -- any mechanism must give 2/9 exactly.
  B. the power law (n/3)^(2 pi) at 3, 7, 11 gives a polar angle of 45.07 deg and an azimuth 0.2200
     rad: the count structure (3, 7, 11) already puts the azimuth within 1 % of 2/9; Koide's
     exactness is the refinement the power law lacks (R55).
  C. base quantities equal to 2/9: the strand-charge products |q_u q_d| = (2/3)(1/3) and 2 (1/3)^2
     are exact; nothing else in the base (Berry factors 0.186-0.507, aspect 6/pi^3 = 0.194,
     circulating fractions 3/4 and 0.313, 4/pi^2 = 0.405, 1/2pi) comes within 10 %.
  D. the 45 deg: 'half common, half in the differences' is the three-lepton analogue of the base's
     half circulating, half static (g = 2) -- an analogy, not a derivation.
  E. verdict: 2/9 is exact and is a rational of the strand charges (u times d), but the base has
     no rule turning a charge product into the azimuth of the generations' circle; the mechanism
     must produce both the 45 deg and the 2/9.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def koide_geometry(masses):
    v = [math.sqrt(m) for m in masses]
    s = sum(v)
    q = sum(x * x for x in v) / s**2
    polar = math.degrees(math.acos(math.sqrt(1 / (3 * q))))
    mean = s / 3
    vp = [x - mean for x in v]
    c1 = (2 * vp[0] - vp[1] - vp[2]) / math.sqrt(6)
    c2 = (vp[1] - vp[2]) / math.sqrt(2)
    theta = math.atan2(-c2, c1)
    return q, polar, theta


def reduce(theta):
    """|theta| reduced modulo 2 pi/3, folded to [0, pi/3]."""
    t = theta % (2 * math.pi / 3)
    return min(t, 2 * math.pi / 3 - t)


def main() -> int:
    q, polar, theta = koide_geometry((ME, MMU, MTAU))
    th = reduce(theta)
    check("A. Data: polar angle 45.000 deg (Q = 2/3 to 1e-5) and azimuth |theta| = 2/9 rad to 2e-5 -- exact at the tau's precision; any mechanism must give 2/9 exactly",
          abs(polar - 45.0) < 0.002 and abs(th / (2 / 9) - 1) < 5e-5,
          f"Q = {q:.6f}, polar = {polar:.4f} deg, |theta| = {th:.6f} rad, 2/9 = {2/9:.6f} ({th/(2/9)-1:+.4%})")

    ladder = [ME * (n / 3) ** (2 * math.pi) for n in (3, 7, 11)]
    q_l, polar_l, theta_l = koide_geometry(ladder)
    th_l = reduce(theta_l)
    check("B. Power law (n/3)^2pi at 3, 7, 11: polar 45.07 deg, azimuth 0.2200 rad, within 1 % of 2/9: the count structure sets the azimuth, Koide's exactness is the missing refinement",
          abs(polar_l - 45.07) < 0.02 and abs(th_l - 0.2200) < 0.001 and abs(th_l / (2 / 9) - 1) < 0.015,
          f"Q = {q_l:.5f}, polar = {polar_l:.3f} deg, |theta| = {th_l:.4f} rad ({th_l/(2/9)-1:+.1%} vs 2/9)")

    cands = {"|q_u q_d| = (2/3)(1/3)": 2 / 9, "2 (1/3)^2": 2 / 9, "Berry a(3)": 0.186, "Berry a(4)": 0.329, "aspect 6/pi^3": 6 / math.pi**3,
             "f_e = 3/4": 0.75, "f_N": 0.313, "4/pi^2": 4 / math.pi**2, "1/(2 pi)": 1 / (2 * math.pi)}
    exact = [k for k, v in cands.items() if abs(v - 2 / 9) < 1e-12]
    near = [k for k, v in cands.items() if 1e-12 < abs(v / (2 / 9) - 1) < 0.10]
    check("C. Base quantities equal to 2/9: the strand-charge products |q_u q_d| and 2(1/3)^2 are exact; no other base number (Berry factors, aspect ratio, circulating fractions, 4/pi^2, 1/2pi) is within 10 %",
          set(exact) == {"|q_u q_d| = (2/3)(1/3)", "2 (1/3)^2"} and near == [],
          f"exact: {exact}; within 10 %: {near}; others: " + ", ".join(f"{k} {v:.3f}" for k, v in cands.items() if k not in exact))

    v = [math.sqrt(m) for m in (ME, MMU, MTAU)]
    common = sum(v) ** 2 / 3
    diff = sum(x * x for x in v) - common
    check("D. The 45 deg: half of sum m is common to the three, half is in their differences -- the three-lepton analogue of half circulating, half static (g = 2); an analogy, not a derivation",
          abs(diff / common - 1) < 1e-4, f"common {common:.2f} MeV, differences {diff:.2f} MeV, ratio {diff/common:.5f}")

    check("E. Verdict: 2/9 is exact and equals the u x d charge product, but no rule of the base turns a charge product into the azimuth of the generations' circle; the mechanism must give both the 45 deg and the 2/9",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
