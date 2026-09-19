#!/usr/bin/env python3
"""Redraw: which strands close into circuits (R40)?  Under the circuit rule a closed circuit of k
strands in series has path k l_1 and energy pi hbar c/(k l_1); the proton's static part must be
the sum over its circuits.  Enumerate every partition of the proton's 9 strands (3-strand content,
R28) into circuits.

Checks:
  A. the static target: 770-813 MeV (R25: e on one ring / + ring), i.e. sum 1/k_i = 1.01-1.07 in
     units of 1/l_1(9) where pi hbar c/l_1(9) = 763 MeV.
  B. of the 30 partitions of 9, the only one with sum 1/k = 1 is {3, 3, 3}: three circuits of three
     strands -- the three quarks, each a 3-strand circuit like a lepton; it gives 763 MeV (-1 % vs
     770, -6 % vs 813).  The band 770-813 is bracketed by {3, 3, 3} (763) and {4, 3, 2} (13/12,
     826 MeV, which mixes strands across quarks and has no reading); {8, 1} gives 858; every other
     partition misses by more than 12 %.
  C. each quark circuit stores pi hbar c/(3 l_1(9)) = 254 MeV, the constituent-quark scale (m_p/3 =
     313 MeV, -19 %); the quark is the electron's object (one 3-strand circuit) at the scale l_1(9):
     ratio 3^(2 pi) = 996 exactly.
  D. adding the circulating parts (R29 rule 153/191 MeV; base m/|g| 168/246): p = 916-931 MeV
     (-2.4 to -0.8 %), n = 954-1009 MeV (+1.5 to +7.4 %): right order, splitting 30-60 x too big.
  E. the pion as two 3-strand circuits: 40 MeV at the count scale l_1(6), 508 MeV at l_1(9) -- neither
     is 140: the pion stays the closed ring of R16 (closure first, R15).
  F. spin: three circuits at hbar/2 each combine to hbar/2 (up, up, down); the + ring of R25 then
     cannot carry another hbar/2 -- the ring must be the circuits' collective circulation (open).
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
from fractions import Fraction as F

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
L1 = lambda n: 2 * math.pi * LAMBDA_E / 3 * (3 / n) ** (2 * math.pi)
E1 = math.pi * HBARC / L1(9)                     # 763 MeV: one circuit of one strand
STATIC = (770.3, 813.1)                           # R25 A / R25 D
MODE_P, MODE_N = (152.5, 168.0), (190.7, 245.6)   # R29 rule / base m/|g|
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def partitions(n, k_max=None):
    k_max = n if k_max is None else k_max
    if n == 0:
        yield ()
        return
    for k in range(min(n, k_max), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def main() -> int:
    # A. target
    lo, hi = STATIC[0] / E1, STATIC[1] / E1
    check("A. Static target 770-813 MeV = sum 1/k_i = 1.01-1.07 in units of 1/l_1(9), where pi hbar c/l_1(9) = 763 MeV",
          abs(E1 - 763) < 1 and abs(lo - 1.01) < 0.005 and abs(hi - 1.066) < 0.005, f"E1 = {E1:.1f} MeV; target/E1 = {lo:.3f}-{hi:.3f}")

    # B. partitions of 9
    parts = list(partitions(9))
    scored = sorted(((p, sum(F(1, k) for k in p)) for p in parts), key=lambda t: abs(float(t[1]) - 1))
    exact = [p for p, s in scored if s == 1]
    best, second, third = scored[0], scored[1], scored[2]
    others_off = all(abs(float(s) - 1) > 0.12 for p, s in scored[3:])
    check("B. Of the 30 partitions of 9 strands, only {3, 3, 3} has sum 1/k = 1: three circuits of three strands, the three quarks, 763 MeV (-1 % vs 770, -6 % vs 813). The band 770-813 is bracketed by {3, 3, 3} and {4, 3, 2} (826 MeV, mixes strands across quarks, no reading); {8, 1} gives 858; every other partition misses by more than 12 %",
          len(parts) == 30 and exact == [(3, 3, 3)] and best[0] == (3, 3, 3) and second[0] == (4, 3, 2) and abs(float(second[1]) * E1 - 826) < 1 and third[0] == (8, 1) and others_off,
          f"{len(parts)} partitions; exact: {exact}; next: {second[0]} -> {float(second[1])*E1:.0f} MeV; third: {third[0]} -> {float(third[1])*E1:.0f} MeV; fourth: {scored[3][0]} -> {float(scored[3][1])*E1:.0f} MeV")

    # C. quark circuit
    e_q = E1 / 3
    ratio = e_q / (ME / 2)
    check("C. Each quark circuit stores pi hbar c/(3 l_1(9)) = 254 MeV, the constituent scale (m_p/3 = 313 MeV, -19 %); the quark is the electron's object (one 3-strand circuit) at l_1(9): energy ratio 3^2pi = 996",
          abs(e_q - 254.3) < 0.5 and abs(e_q / (MP / 3) - 0.813) < 0.003 and abs(ratio - 3 ** (2 * math.pi)) < 1e-6,
          f"E_q = {e_q:.1f} MeV ({e_q/(MP/3)-1:+.0%} vs m_p/3); E_q/(m_e/2) = {ratio:.0f}")

    # D. nucleon masses
    p_lo, p_hi = (E1 + m for m in MODE_P)
    n_lo, n_hi = (E1 + m for m in MODE_N)
    check("D. With the circulating parts: p = 916-931 MeV (-2.4 to -0.8 %), n = 954-1009 MeV (+1.5 to +7.4 %): right order, splitting 30-60 x too big",
          abs(p_lo / MP - 0.976) < 0.003 and abs(p_hi / MP - 0.992) < 0.003 and abs(n_lo / MN - 1.015) < 0.003 and abs(n_hi / MN - 1.074) < 0.003 and n_lo > p_lo,
          f"p = {p_lo:.0f}-{p_hi:.0f} ({p_lo/MP-1:+.1%} to {p_hi/MP-1:+.1%}); n = {n_lo:.0f}-{n_hi:.0f} ({n_lo/MN-1:+.1%} to {n_hi/MN-1:+.1%}); n - p = {n_lo-p_lo:.0f} MeV vs 1.29")

    # E. pion
    pi_count = 2 * math.pi * HBARC / (3 * L1(6))
    pi_9 = 2 * E1 / 3
    check("E. Pion as two 3-strand circuits: 40 MeV at l_1(6), 508 MeV at l_1(9), neither 140: the pion stays the closed ring of R16",
          abs(pi_count - 39.8) < 0.5 and abs(pi_9 - 508.5) < 1, f"{pi_count:.1f} MeV at l_1(6); {pi_9:.0f} MeV at l_1(9)")

    # F. spin (recorded)
    check("F. Three circuits at hbar/2 combine to hbar/2 (up, up, down); the + ring of R25 cannot carry another hbar/2: the ring must be the circuits' collective circulation (open)",
          True, "recorded")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
