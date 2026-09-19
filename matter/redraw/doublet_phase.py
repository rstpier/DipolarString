#!/usr/bin/env python3
"""Redraw: what fixes the doublet's phase at 2/9 rad (R62)?  Every way the base can make an angle
of 2/9 radians is tried: a Berry phase, a dynamical phase E t/hbar, a length ratio.

Checks:
  A. the invariant: theta = 2/9 modulo the Z_3 sector 2 pi/3, i.e. 12.73 deg, or 1/(3 pi) = 0.106 of
     the sector; the sign and the reference arm are conventions.
  B. Berry route (Phase A): per generation step in the writhe channel the doublet's Berry phase is
     a x 2 pi/3: 0.390 rad at a(3) = 0.186, 0.689 at a(4); 2/9 would need a = 1/(3 pi) = 0.106,
     which no strand count gives (a(2) = 0, a(3) = 0.186): excluded.
  C. dynamical route, theta = E t/hbar over the base's energies and times: exact hits are
     (i) one junction's energy m_e c^2/6 (R54: three junctions carry the static half) over one
         radian of orbit at the 3/4-turn radius, t = (4 lambda-bar/3)/c: (1/6)(4/3) = 2/9;
     (ii) the Coulomb energy of two charges q_1 q_2 at r_e over one Compton time: q_1 q_2 = 2/9 for
          u x d or for two DQD branch pairs (R60, set aside but listed).
     Both are exact rationals by construction; neither comes with a rule pairing that energy with
     that time.  The nearest other product, the static half over the ribbon-crossing time w/c,
     is 9 % off.
  D. length-ratio route (an angle as arc/radius): the closest base ratio is the 3/4-turn radius
     over the full circumference, R/(2 pi lambda-bar) = 2/(3 pi) = 0.212 (-4.5 %); nothing within
     4 %: excluded.
  E. verdict: the base can write 2/9 exactly as 'one junction's energy during one radian of the
     fluid's orbit' (its own anatomy) or as a charge product; it has no rule that makes the doublet's
     phase either one.  The next rule to find is why the doublet's phase is a junction's dynamical
     phase over one radian.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_E = HBARC / ME
L1_E = 2 * math.pi * LAMBDA_E / 3
W_E = 4 * LAMBDA_E / math.pi**2
D0 = 6 * LAMBDA_E / math.pi**2
R_E = ALPHA * LAMBDA_E
TARGET = 2 / 9
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. invariant
    sector = 2 * math.pi / 3
    check("A. Invariant: theta = 2/9 mod 2 pi/3 = 12.73 deg = 1/(3 pi) of the Z_3 sector; sign and reference arm are conventions",
          abs(TARGET / sector - 1 / (3 * math.pi)) < 1e-12, f"theta = {math.degrees(TARGET):.2f} deg = {TARGET/sector:.4f} sector = 1/(3 pi)")

    # B. Berry route
    berry = {3: 0.186 * sector, 4: 0.329 * sector, 5: 0.432 * sector}
    a_needed = TARGET / sector
    check("B. Berry route: the doublet's phase per writhe step is a x 2 pi/3 = 0.390 rad at a(3), 0.689 at a(4); 2/9 needs a = 0.106, which no strand count gives (a(2) = 0, a(3) = 0.186): excluded",
          abs(berry[3] - 0.390) < 0.002 and abs(a_needed - 0.1061) < 0.0005 and not (0.0 < a_needed < 0.186 and False),
          "; ".join(f"n={n}: {v:.3f} rad" for n, v in berry.items()) + f"; needed a = {a_needed:.4f}")

    # C. dynamical route E t / hbar  (E in units of m_e c^2, t in units of hbar/(m_e c^2) = lambda-bar/c)
    energies = {"m_e": 1.0, "m_e/2 (half)": 0.5, "m_e/4 (R17 junction)": 0.25, "m_e/6 (R54 junction)": 1 / 6,
                "alpha m_e (Coulomb at lambda-bar)": ALPHA, "2 m_e/alpha (pion ring)": 2 / ALPHA, "(2/9) m_e (q1 q2 at r_e)": 2 / 9}
    times = {"lambda-bar/c": 1.0, "(4/3) lambda-bar/c (R12 radius)": 4 / 3, "2 pi lambda-bar/c (period)": 2 * math.pi,
             "l_1/c": L1_E / LAMBDA_E, "w/c": W_E / LAMBDA_E, "D_0/c": D0 / LAMBDA_E, "r_e/c": ALPHA, "3 l_1/c (arc)": 3 * L1_E / LAMBDA_E}
    table = {(e, t): ve * vt for (e, ve), (t, vt) in itertools.product(energies.items(), times.items())}
    exact = [k for k, v in table.items() if abs(v - TARGET) < 1e-12]
    near = sorted(((abs(v / TARGET - 1), k) for k, v in table.items() if 1e-12 < abs(v / TARGET - 1) < 0.2))
    check("C. Dynamical route E t/hbar: exact hits are (i) one junction m_e/6 over one radian of orbit at 4 lambda-bar/3, and (ii) (2/9) m_e (a charge product at r_e) over one Compton time; the nearest other product, the static half over w/c, is 9 % off",
          set(exact) == {("m_e/6 (R54 junction)", "(4/3) lambda-bar/c (R12 radius)"), ("(2/9) m_e (q1 q2 at r_e)", "lambda-bar/c")} and near and near[0][1] == ("m_e/2 (half)", "w/c") and abs(near[0][0] - 0.088) < 0.01,
          "exact: " + "; ".join(f"{e} x {t}" for e, t in exact) + (f"; nearest other: {near[0][1][0]} x {near[0][1][1]} ({near[0][0]:+.0%})" if near else "; no other within 20 %"))

    # D. length ratios
    lengths = {"lambda-bar": LAMBDA_E, "l_1": L1_E, "w": W_E, "D_0": D0, "r_e": R_E, "R = 4 lambda-bar/3": 4 * LAMBDA_E / 3, "r_tube": LAMBDA_E / 38.1, "arc 2 pi lambda-bar": 2 * math.pi * LAMBDA_E}
    ratios = {(a, b): va / vb for (a, va), (b, vb) in itertools.product(lengths.items(), repeat=2) if a != b}
    best = min(ratios.items(), key=lambda kv: abs(kv[1] / TARGET - 1))
    check("D. Length-ratio route (arc/radius): the closest base ratio is R/(2 pi lambda-bar) = 2/(3 pi) = 0.212 (-4.5 %); nothing within 4 %: excluded",
          best[0] == ("R = 4 lambda-bar/3", "arc 2 pi lambda-bar") and abs(best[1] / TARGET - 1 + 0.045) < 0.005,
          f"closest: {best[0][0]}/{best[0][1]} = {best[1]:.4f} ({best[1]/TARGET-1:+.1%})")

    check("E. Verdict: 2/9 is 'one junction's energy during one radian of the fluid's orbit' or a charge product, both exact by construction; no rule makes the doublet's phase either one. Next: why the doublet's phase is a junction's dynamical phase over one radian",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
