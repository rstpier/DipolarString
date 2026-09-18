#!/usr/bin/env python3
"""Redraw: 'stowing = intertwining' -- the ends of two strings are wound around each other, which
prevents their separation.  Tested in the base (BASE.md R1, R14, R17, R22, R26).

Checks:
  A. topology: the linking number of two intertwined ends is conserved only while all four ends are
     held (closed network); an isolated pair unwinds by free rotation.  So the lock is a property of
     the network, not of the pair: an isolated quark holds no twist (recorded).
  B. base mode energy of a twisted pair at fixed axial length l = 0.90 fm (R26 tube) and helix radius
     lambda-bar_p: the arc length grows with the number of turns and the mode hbar c pi/s FALLS:
     N = 1 gains 302 MeV, N = 2 gains 467 MeV, the limit is hbar c pi/l = 690 MeV.  Winding is
     spontaneous and unwinding (needed to separate) costs up to 690 MeV: the right order, but a
     BINDING -- the twist holds the ends, it does not store the 813 MeV.
  C. magnetic dual: the wound fluid is a solenoid of a perfect conductor, its flux is conserved, and
     stretching at fixed flux per turn costs sigma_B = Phi_1^2/(2 mu0 A).  With Phi_1 = h/e and
     A = pi lambda-bar_p^2: sigma_B = (pi^2/4 alpha) x sigma_E = 338 x 904 MeV/fm = 306 GeV/fm.
     Matching the lattice needs a tube of radius 18.4 lambda-bar_p = 3.9 fm (bigger than the proton)
     or a flux of 0.054 h/e.  Right sign, wrong scale: R26's electric tube is the one that fits.
  D. duality: sigma_E = sigma_B when Phi = Z_0 delta; the base's matched pole delta = e/(pi sqrt alpha)
     corresponds to Phi = (2 sqrt alpha/pi) h/e = 0.054 h/e -- not a flux quantum.
  E. intertwined arms are 3-strand braids (Bilson-Thompson).  BT's strands carry {-1/3, 0, +1/3} and
     make u = (+,+,0), d = (-,0,0), e = (-,-,-) with three strands; the base's rule (every string
     +-1/3) cannot make +2/3 with three strings (R22).  One neutral string -- R1's open question,
     a fluid with poles but no net charge -- would let the base build every quark as a 3-strand
     braid, as BT does.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction as F

HBARC = 197.3269804
ALPHA = 1 / 137.035999
MP = 938.27209
LAMBDA_P = HBARC / MP
L_TUBE = 0.899                       # fm, R26
SIGMA_E = 2 * MP**2 / (math.pi**2 * HBARC)   # 904 MeV/fm, R26
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. topology (recorded)
    check("A. Linking number is conserved only for held ends (closed network); an isolated intertwined pair unwinds by free rotation: the lock belongs to the network, an isolated quark holds no twist (recorded)",
          True, "Calugareanu-White: Lk = Tw + Wr is invariant for closed curves; open ends let twist flow out")

    # B. mode energy of a twisted pair
    def arc(N, l=L_TUBE, R=LAMBDA_P):
        return l * math.sqrt(1 + (2 * math.pi * R * N / l) ** 2)
    e_straight = HBARC * math.pi / L_TUBE
    gains = {N: e_straight - HBARC * math.pi / arc(N) for N in (1, 2, 3, 10)}
    check("B. Twisted pair at fixed axial length 0.90 fm, helix radius lambda-bar_p: the mode hbar c pi/s falls with the turns (N = 1: -302 MeV, N = 2: -467, limit -690 = hbar c pi/l): winding is spontaneous, unwinding costs up to 690 MeV -- a binding of the right order, not the stored 813 MeV",
          abs(gains[1] - 302) < 3 and abs(gains[2] - 467) < 3 and abs(e_straight - 690) < 2 and gains[10] < e_straight,
          "; ".join(f"N={N}: s = {arc(N):.2f} fm, gain {g:.0f} MeV" for N, g in gains.items()) + f"; hbar c pi/l = {e_straight:.0f} MeV")

    # C. magnetic dual
    ratio = math.pi**2 / (4 * ALPHA)                     # sigma_B(h/e) / sigma_E(delta matched)
    sigma_B = ratio * SIGMA_E
    r_match = math.sqrt(ratio) * LAMBDA_P
    phi_match = 1 / math.sqrt(ratio)                     # in h/e
    check("C. Wound fluid as a flux-conserving solenoid: with one flux quantum h/e per turn in a tube of radius lambda-bar_p, sigma_B = (pi^2/4 alpha) sigma_E = 338 x 904 MeV/fm = 306 GeV/fm; the lattice needs a 3.9 fm tube or 0.054 h/e: right sign, wrong scale",
          abs(ratio - 338.3) < 0.5 and abs(sigma_B / 1e3 - 306) < 1 and abs(r_match - 3.87) < 0.02 and abs(phi_match - 0.0544) < 0.0005,
          f"ratio = {ratio:.1f}; sigma_B = {sigma_B/1e3:.0f} GeV/fm; radius for lattice sigma = {r_match/LAMBDA_P:.1f} lambda_p = {r_match:.2f} fm; flux for lattice sigma = {phi_match:.4f} h/e")

    # D. duality
    phi_of_delta = 2 * math.sqrt(ALPHA) / math.pi        # Phi = Z0 delta in units of h/e
    check("D. sigma_E = sigma_B when Phi = Z_0 delta: the matched pole e/(pi sqrt alpha) corresponds to Phi = (2 sqrt alpha/pi) h/e = 0.054 h/e, not a flux quantum",
          abs(phi_of_delta - phi_match) < 1e-12, f"Phi/(h/e) = {phi_of_delta:.4f} = 0.109 x h/2e")

    # E. braids and neutral strands
    base_3 = {sum(c) for c in itertools.product((F(-1, 3), F(1, 3)), repeat=3)}
    bt_3 = {sum(c) for c in itertools.product((F(-1, 3), F(0), F(1, 3)), repeat=3)}
    check("E. Three intertwined strands (Bilson-Thompson braid): with strands in {-1/3, 0, +1/3} the charges include +2/3 (u = (+,+,0)); with the base's +-1/3 only they do not. One neutral string (R1's open question) lets the base build every quark as a 3-strand braid",
          F(2, 3) not in base_3 and F(2, 3) in bt_3 and F(-1, 3) in bt_3 and F(-1) in bt_3,
          "base: " + ", ".join(str(q) for q in sorted(base_3)) + "; BT: " + ", ".join(str(q) for q in sorted(bt_3)))

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
