#!/usr/bin/env python3
"""Redraw, after R24 ('what fixes the size of the centre'): apply the base's own size-fixing
mechanism -- the fluid at c on a ring, S = R E_circ / c = hbar/2, mu = q c R / 2, g = m c^2 / E_circ
(rod_z0, massless_fluid, g_bookkeeping) -- to the nucleons, with the R24 content p = (5+, 2-),
n = (4+, 4-).

Checks:
  A. proton with its whole charge e on one ring: mu_p = 2.793 mu_N gives R = 2.793 lambda-bar_p =
     0.587 fm (0.70 r_p); g_p = 5.586 = m c^2 / E_circ means 18 % of the proton's energy circulates,
     82 % is static in the centre (electron: 50 / 50).
  B. R24 content on two rings (positive strings at R+, negative at R-), same graph for p and n so
     the same radii: mu_p and mu_n fix R+ = 3.75 lambda-bar_p = 0.789 fm (-6 % vs r_p) and
     R- = 5.19 lambda-bar_p = 1.09 fm: the negative strings circulate outside, as the neutron's
     negative skin requires.  Two data, two unknowns: R+ near r_p is the only content.
  C. the same rings as static charge fail the charge radii: r_p would be 0.50 fm (-41 %) and
     <r^2>_n = -0.75 fm^2 (6.5 x too negative): in the base the charge sits at the poles and the
     current on the arc (R5), so the moments locate the fluid, not the charge.
  D. spin: the + ring alone carrying S = hbar/2 circulates E = hbar c / (2 R+) = 125 MeV, 13 % of m_p.
  E. what is NOT fixed: nothing in the base sets the 82-87 % static energy of the centre (~800 MeV);
     the electron's static half was fixed by g = 2, the nucleon's g only tells how much circulates.
Exit status is zero only if every check passes.  PDG 2024 values.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
MP, MN = 938.27209, 939.56542
RP, R2N = 0.8409, -0.1155
MU_P, MU_N = 2.79284734, -1.91304273
LAMBDA_P = HBARC / MP
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. whole charge on one ring
    R1 = MU_P * LAMBDA_P
    g_p = 2 * MU_P
    frac = 1 / g_p
    check("A. Whole charge e on one ring: mu_p fixes R = 2.793 lambda-bar_p = 0.587 fm (0.70 r_p); g_p = 5.586 = m c^2/E_circ: 18 % of the energy circulates, 82 % is static in the centre (electron: 50/50)",
          abs(R1 - 0.587) < 0.002 and abs(R1 / RP - 0.698) < 0.003 and abs(frac - 0.179) < 0.001,
          f"R = {R1:.3f} fm = {R1/RP:.3f} r_p; E_circ/m_p = {frac:.3f}")

    # B. two rings, R24 content
    # mu_p/mu_N = (5/3) R+ - (2/3) R-  (in lambda_p);  mu_n/mu_N = (4/3)(R+ - R-)
    diff = MU_N / (4 / 3)                      # R+ - R-
    Rplus = (MU_P + (2 / 3) * (-diff)) / (5 / 3 - 2 / 3)
    Rminus = Rplus - diff
    check("B. R24 content on two rings with the same radii for p and n: mu_p and mu_n give R+ = 3.75 lambda-bar_p = 0.789 fm (-6 % vs r_p) and R- = 5.19 lambda-bar_p = 1.09 fm: the negative strings circulate outside (the neutron's negative skin). Two data, two unknowns; R+ near r_p is the only content",
          abs(Rplus - 3.750) < 0.005 and abs(Rminus - 5.185) < 0.005 and abs(Rplus * LAMBDA_P / RP - 0.938) < 0.003,
          f"R+ = {Rplus:.3f} lambda_p = {Rplus*LAMBDA_P:.3f} fm ({Rplus*LAMBDA_P/RP-1:+.1%} vs r_p); R- = {Rminus:.3f} lambda_p = {Rminus*LAMBDA_P:.3f} fm")

    # C. the same rings as static charge
    rp2 = (5 / 3) * (Rplus * LAMBDA_P)**2 - (2 / 3) * (Rminus * LAMBDA_P)**2
    rn2 = (4 / 3) * ((Rplus * LAMBDA_P)**2 - (Rminus * LAMBDA_P)**2)
    check("C. The same rings as static charge fail: r_p would be 0.50 fm (-41 %), <r^2>_n = -0.75 fm^2 (6.5 x too negative): the moments locate the fluid, not the charge (R5: charge at the poles, current on the arc)",
          abs(math.sqrt(rp2) - 0.495) < 0.005 and abs(rn2 / R2N - 6.5) < 0.2,
          f"sqrt<r^2>_p = {math.sqrt(rp2):.3f} fm ({math.sqrt(rp2)/RP-1:+.0%}); <r^2>_n = {rn2:.3f} fm^2 ({rn2/R2N:.1f} x measured)")

    # D. spin on the + ring
    E_circ = HBARC / (2 * Rplus * LAMBDA_P)
    check("D. S = hbar/2 on the + ring: E_circ = hbar c/(2 R+) = 125 MeV, 13 % of the proton",
          abs(E_circ - 125.1) < 0.5 and abs(E_circ / MP - 0.133) < 0.002, f"E_circ = {E_circ:.1f} MeV = {E_circ/MP:.3f} m_p")

    # E. what is not fixed
    static = MP - E_circ
    check("E. Not fixed: the static 87 % (~810 MeV) of the centre; the electron's static half came from g = 2, the nucleon's g only says how much circulates",
          static > 0.8 * MP, f"static = {static:.0f} MeV = {static/MP:.0%} of m_p")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
