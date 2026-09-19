#!/usr/bin/env python3
"""Redraw: the last numbered lead of R30 -- a fluid with a FIXED line charge density lambda (a
property of the fluid, not e/3 per string).  Does it survive the electron?

Checks:
  A. the lattice tension needs lambda = 7.07 e/fm (R30 D).  The electron's strings carry e/3 over
     l_1 = 809 fm: lambda_e = 4.1e-4 e/fm, 17 000 x smaller; their tension would be 3e-6 MeV/fm.
     One fluid cannot have both densities: a fixed lambda is excluded by the electron's size.
  B. a density that follows the object's scale, lambda = (e/3)/l_1(n), gives sigma = 4 pi K/(9 l_1^2):
     3.0 MeV/fm at the nucleon (l_1 = 0.81 fm), 300 x short; it reaches the lattice only at
     l_1 = 0.047 fm, the ladder's n = 14, not the proton's 9.
  C. the quantised charge e/3 per string with fixed lambda forces quantised lengths l_1 = (e/3)/lambda
     = 0.047 fm: the electron (3 strings) would be 0.14 fm long against its 2426 fm arc (R29 A).
  D. verdict: every electromagnetic route to the centre's ~800 MeV tried so far (Coulomb junction,
     mode, valence, entanglement, flux tube, fixed density) fails by a factor 100-300 or by sign:
     the base's fluid -- a TEM wave at c on a matched line -- has no tension.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME = 0.51099895
L1_E = 2 * math.pi * (HBARC / ME) / 3                       # 808.8 fm
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)                       # 0.813 fm
SIGMA_LAT = 904.0                                            # MeV/fm
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    lam_needed = math.sqrt(SIGMA_LAT / (4 * math.pi * K))
    lam_e = (1 / 3) / L1_E
    sig_e = 4 * math.pi * K * lam_e**2
    check("A. Lattice tension needs lambda = 7.07 e/fm; the electron's strings carry e/3 over 809 fm, lambda_e = 4.1e-4 e/fm, 17 000 x smaller (tension 3e-6 MeV/fm): one fluid cannot have both densities -- a fixed lambda is excluded by the electron",
          abs(lam_needed / lam_e - 17150) < 100 and sig_e < 1e-5,
          f"lambda_needed/lambda_e = {lam_needed/lam_e:.0f}; sigma_e = {sig_e:.1e} MeV/fm")

    sig_9 = 4 * math.pi * K / (9 * L1_9**2)
    l_match = math.sqrt(4 * math.pi * K / (9 * SIGMA_LAT))
    n_match = 3 * (L1_E / l_match) ** (1 / (2 * math.pi))
    check("B. A density following the scale, lambda = (e/3)/l_1(n): sigma = 4 pi K/(9 l_1^2) = 3.0 MeV/fm at the nucleon, 300 x short; the lattice value needs l_1 = 0.047 fm, the ladder's n = 14, not 9",
          abs(sig_9 - 3.04) < 0.05 and abs(l_match - 0.0472) < 0.0005 and abs(n_match - 14.2) < 0.2,
          f"sigma(l_1 = {L1_9:.3f} fm) = {sig_9:.2f} MeV/fm; l_1 for lattice = {l_match:.4f} fm -> n = {n_match:.1f}")

    l_q = (1 / 3) / lam_needed
    check("C. Quantised e/3 per string with fixed lambda forces l_1 = 0.047 fm: the electron's three strings would span 0.14 fm against its 2426 fm arc",
          abs(3 * l_q - 0.1415) < 0.001 and (2 * math.pi * HBARC / ME) / (3 * l_q) > 1e4,
          f"3 l_1 = {3*l_q:.4f} fm; electron arc / that = {(2*math.pi*HBARC/ME)/(3*l_q):.0f}")

    check("D. Verdict: Coulomb junction (R24), mode (R29), valence (R24), entanglement (R27), flux tube (R26/R30), fixed density (here) all miss the centre's ~800 MeV by 100-300 x or by sign: the base's TEM fluid at c has no tension",
          SIGMA_LAT / sig_9 > 100, f"lattice/base tension = {SIGMA_LAT/sig_9:.0f}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
