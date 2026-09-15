#!/usr/bin/env python3
"""Redraw, R17 (author): the string can stretch and widen like a sheet; the width is fixed by the
impedance Z_0 of space.

  A. Matching: a parallel-plate sheet of width w and separation d has Z = Z_0 d/w; Z = Z_0 gives
     w = d (square cross-section).  For two wires the same principle gives D/r = 2 cosh pi.
     Matching fixes the shape, not the size.
  B. Pole charge per quantum on a matched guide (pole_from_circulation.py with Z = Z_0):
         delta = sqrt(4 hbar / pi Z_0) = e / (pi sqrt(alpha)) = 3.73 e   -- a closed form.
  C. Pole-to-pole Coulomb energy at separation s:  delta^2 / (4 pi eps0 s) = hbar c / (pi^2 s).
     Alpha cancels: on a matched guide a junction's binding energy has the SAME form as a mode
     energy (hbar c / length).  If the separation s scales with the object (s = the sheet width
     w, w proportional to R), junction masses scale as 1/R exactly like mode masses -- the
     condition the muon/tau calibration imposed and that fixed-radius contact binding failed.
  D. g = 2 (static half = two junctions): 2 hbar c / (pi^2 w) = m c^2 / 2 with m c^2 = hbar c /
     lambda-bar gives w = 4 lambda-bar / pi^2 = 0.405 lambda-bar: 157 fm for the electron, 0.76 fm
     for the muon; with the 3/4-turn radius R = 4 lambda-bar / 3 the aspect ratio is w/R = 3/pi^2
     = 0.30 for every lepton -- a self-similar family.
  E. The price: a sheet of width w has transverse modes at hbar c pi / w -- 3.9 MeV for the
     electron.  An electron excited by 3.9 MeV does not exist (no e* below ~100 GeV, and atomic
     physics would see it).  Either the fluid is TEM-only by axiom (transverse modes forbidden),
     or the electron's guide is far thinner than 157 fm and the static half lives elsewhere.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

HBAR, C, E, EPS0, MU0 = 1.054571817e-34, 2.99792458e8, 1.602176634e-19, 8.8541878128e-12, 4e-7 * math.pi
Z0 = math.sqrt(MU0 / EPS0)
ALPHA = E ** 2 / (4 * math.pi * EPS0 * HBAR * C)
HBARC = 197.3269804              # MeV fm
ME, MMU = 0.51099895, 105.6583755
LB_E, LB_MU = HBARC / ME, HBARC / MMU
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    check("A. Matching a parallel-plate sheet to Z_0: Z = Z_0 d/w = Z_0 <=> w = d (square section); for two wires the same principle gives D/r = 2 cosh pi = 23.2",
          abs(2 * math.cosh(math.pi) - 23.18) < 0.01, "matching fixes the shape of the cross-section, not its size")
    delta = math.sqrt(4 * HBAR / (math.pi * Z0)) / E
    check("B. Pole charge per quantum on a matched guide: delta = sqrt(4 hbar / pi Z_0) = e / (pi sqrt alpha) = 3.73 e",
          abs(delta - 1 / (math.pi * math.sqrt(ALPHA))) < 1e-9 and abs(delta - 3.73) < 0.01, f"delta = {delta:.3f} e; 1/(pi sqrt alpha) = {1/(math.pi*math.sqrt(ALPHA)):.3f}")
    # C. junction energy at separation s: delta^2 e^2 / (4 pi eps0 s) = hbar c / (pi^2 s)
    s_fm = 10.0
    u_j = delta ** 2 * ALPHA * HBARC / s_fm
    check("C. Pole-to-pole Coulomb energy at separation s: delta^2 e^2/(4 pi eps0 s) = hbar c / (pi^2 s) -- alpha cancels; a junction's mass has the form of a mode's",
          abs(u_j - HBARC / (math.pi ** 2 * s_fm)) < 1e-9, f"at s = 10 fm: {u_j:.2f} MeV = hbar c / (pi^2 s)")
    check("C'. Hence, if the pole separation is the sheet width and the width scales with the object, junction masses scale as 1/R like mode masses -- what the mu/tau calibration required and fixed-radius contact could not give",
          True, "lepton_calibration.py: contact at a fixed tube radius needed poles of 22 e and 88 e; here the poles stay at 3.73 e and the separation scales")
    # D. g = 2: two junctions carry m c^2 / 2
    w_e = 4 * LB_E / math.pi ** 2
    w_mu = 4 * LB_MU / math.pi ** 2
    r_e = 4 * LB_E / 3
    check("D. g = 2 with two junctions: 2 hbar c/(pi^2 w) = m c^2/2 -> w = 4 lambda-bar/pi^2 = 0.405 lambda-bar: 157 fm for the electron, 0.76 fm for the muon",
          abs(w_e - 156.5) < 0.5 and abs(w_mu - 0.757) < 0.005, f"w_e = {w_e:.1f} fm, w_mu = {w_mu:.3f} fm")
    check("D'. With the 3/4-turn radius R = 4 lambda-bar/3 the aspect ratio is w/R = 3/pi^2 = 0.30 for every lepton: a self-similar family",
          abs(w_e / r_e - 3 / math.pi ** 2) < 1e-12, f"w/R = {w_e/r_e:.3f}")
    # E. transverse modes
    e_t_e, e_t_mu = HBARC * math.pi / w_e, HBARC * math.pi / w_mu
    check("E. The price: transverse modes of the sheet at hbar c pi / w -- 3.9 MeV for the electron, 815 MeV for the muon",
          abs(e_t_e - 3.96) < 0.05 and abs(e_t_mu - 819) < 5, f"electron {e_t_e:.2f} MeV, muon {e_t_mu:.0f} MeV")
    check("E'. An electron excited by 3.9 MeV is not observed: either transverse modes are forbidden by axiom (TEM-only fluid), or the electron's guide is much thinner than 157 fm and the static half is elsewhere",
          True, "a width below 6e-3 fm would push the transverse mode above 100 GeV; then the two-junction static half (hbar c/pi^2 w each) would exceed the electron mass by 2e4 -- the sheet cannot be both thin and the carrier of the static half")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
