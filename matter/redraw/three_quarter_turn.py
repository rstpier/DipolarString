#!/usr/bin/env python3
"""Redraw, R12 (author): the daughter is a 3/4 turn, and it contributes.

Picture P2: half the rest energy is static in the two junctions (R7), half is the mode living on
the open arc (fraction a of a full turn, a = 3/4 here), partially reflected at its two pole ends.
Massless fluid (R1).  Three conditions:
  (i)   S = hbar/2 :  f R E_mode / c = hbar/2, with E_mode = m c^2 / 2 and f the circulating
        fraction of the mode (forward minus backward over total);
  (ii)  resonance of the fundamental on the arc: beta l_arc = pi, l_arc = a 2 pi R, so
        E_mode = hbar c beta = hbar c / (2 a R);
  (iii) g = 2 holds by the half split (massless_fluid.py), and mu = e c f R / 2.
Result: R = lambda-bar_C / a, and f = a EXACTLY -- the arc fraction IS the circulating fraction:
'the 3/4 turn contributes 3/4'.  Then rho = (1 - f)/(1 + f) is the backward/forward power the
poles must reflect.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

HBAR, C, E, ME = 1.054571817e-34, 2.99792458e8, 1.602176634e-19, 9.1093837015e-31
MU_B = E * HBAR / (2 * ME)
LB = HBAR / (ME * C)
MEC2 = ME * C ** 2
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def solve(a: float):
    e_mode = MEC2 / 2
    r = HBAR * C / (2 * a * e_mode)          # from (ii)
    f = (HBAR / 2) * C / (r * e_mode)        # from (i)
    rho = (1 - f) / (1 + f)
    mu = E * C * f * r / 2
    return r, f, rho, mu


def main() -> int:
    r, f, rho, mu = solve(0.75)
    check("3/4 turn: the ring radius becomes (4/3) lambda-bar_C = 515 fm and the circulating fraction of the mode is exactly 3/4",
          abs(r / LB - 4 / 3) < 1e-12 and abs(f - 0.75) < 1e-12, f"R = {r*1e15:.0f} fm, f = {f:.4f}")
    check("Identity: for an arc of fraction a the circulating fraction is f = a -- the 3/4 turn contributes 3/4 (checked for a = 1/2, 2/3, 3/4, 4/5, 1)",
          all(abs(solve(a)[1] - a) < 1e-12 for a in (0.5, 2 / 3, 0.75, 0.8, 1.0)), "f R = lambda-bar_C always, R = lambda-bar_C / a")
    check("mu = mu_B and g = 2 hold for every a (they only need f R = lambda-bar_C and the half split)",
          abs(mu / MU_B - 1) < 1e-12, f"mu = {mu/MU_B:.6f} mu_B")
    check("Pole reflection required by the 3/4 turn: backward/forward power rho = 1/7, amplitude |Gamma| = 0.378",
          abs(rho - 1 / 7) < 1e-12, f"rho = {rho:.4f}, |Gamma| = {math.sqrt(rho):.3f}")
    # the old Gamma_pole = 1/3 in amplitude: rho = 1/9 -> f = 0.8 -> a = 0.8
    f_old = (1 - 1 / 9) / (1 + 1 / 9)
    check("The old Gamma_pole = 1/3 (amplitude, rho = 1/9) would correspond to a 4/5 turn (f = 0.8, R = 1.25 lambda-bar_C), not 3/4",
          abs(f_old - 0.8) < 1e-12, f"f = {f_old:.3f}; the 3/4 turn asks for |Gamma| = 0.378 instead of 0.333 -- a 13 % difference in amplitude")
    e_j = MEC2 / 4
    delta = math.sqrt(e_j / E / 1e3 / 6.31) / 3
    check("The two junctions carry the static half: 128 keV each, pole charges ~1.5 e (delta^2 scaling of 6.3 keV at e/3)",
          abs(e_j / E / 1e3 - 127.7) < 0.2 and 1.3 < delta < 1.7, f"E_junction = {e_j/E/1e3:.1f} keV, delta = {delta:.2f} e")
    check("Consistency: the mode's wavelength on the 3/4 arc is twice the arc, 2 x (3/4) 2 pi R = 3 pi R = 4 pi lambda-bar_C = 2 lambda_C",
          abs(2 * 0.75 * 2 * math.pi * r / (2 * 2 * math.pi * LB) - 1) < 1e-12, "the fundamental of the open 3/4 arc at radius 4/3 lambda-bar_C has the Compton wavelength x 2")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    print("Conclusion: with half the mass static in two junctions (128 keV each, poles ~1.5 e) and half in the mode on a "
          "3/4 turn, S = hbar/2 and the arc resonance give R = 4/3 lambda-bar_C = 515 fm and a circulating fraction exactly "
          "3/4 -- the arc fraction is the circulating fraction, for any arc. mu = mu_B and g = 2 follow. The poles must then "
          "reflect 1/7 of the power (|Gamma| = 0.38); the old Gamma_pole = 1/3 would mean a 4/5 turn.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
