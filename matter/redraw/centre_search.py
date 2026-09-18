#!/usr/bin/env python3
"""Redraw: what fixes the ~810 MeV static energy of the nucleon's centre?  The author's lead: 'the
stowing of the ends of two strings, which prevents their separation'.  Every energy the base offers
is tried against E_s = m_p - E_circ = 813 MeV (R25: 125 MeV circulate on the + ring).

Candidates:
  A. the base's own energies all fall as 1/L (mode hbar c pi/L, Coulomb K/s, cavity between two
     stowed ends).  A cavity gives 813 MeV at l = 0.763 fm = R+ (-3 %): the sheet's transverse mode
     (R17) at the nucleon's width.  Right scale, wrong sign: a 1/L energy DROPS when the ends
     separate, so it cannot be what prevents separation.
  B. the cost of cutting (two new matched poles, 2 hbar c/(pi^2 w)): 813 MeV at w = 0.049 fm; the
     lepton rule w = 4 lambda-bar/pi^2 gives 469 MeV = m_p/2 by construction of g = 2.  Nothing fixes w.
  C. the only energy that GROWS with separation: the pole's flux confined to the sheet (flux tube).
     With the matched pole delta = e/(pi sqrt alpha) the tension is sigma = delta^2/(2 eps0 A) =
     2 hbar c/(pi A).  For a tube of radius lambda-bar_p: sigma = 2 (m_p c^2)^2/(pi^2 hbar c) =
     904 MeV/fm, i.e. sqrt(sigma) = (sqrt 2/pi) m_p c^2 = 422 MeV -- inside the lattice range
     (420-440 MeV; 0.89-0.98 GeV/fm).  813 MeV is then 0.90 fm of tube (r_p: 760 MeV, -6.5 %;
     a Y of three arms of 0.30 fm).  CONDITIONAL: the base's conducting line spreads a pole's
     charge along itself (energy ~ 1/L); keeping the flux in the sheet is a new rule (Wen / dual
     superconductor), and the tube radius lambda-bar_p is the proton's own length, so sigma from m_p
     is a consistency relation, not a derivation of m_p.
  D. sign test: pulling the stowed ends 1 fm apart releases 461 MeV under A
     and costs +904 MeV under C.
  E. caveat under C: a pion ring (140 MeV) pays for 0.15 fm of tube, so a 0.9 fm tube stores six
     pions' worth: the base needs a rule that forbids breaking (string breaking in QCD happens
     near 1.2 fm because a pair costs ~1 GeV, not 140 MeV).
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
MP, MPI = 938.27209, 139.57039
RP = 0.8409
LAMBDA_P = HBARC / MP
R_PLUS = 3.749 * LAMBDA_P                 # R25
E_CIRC = HBARC / (2 * R_PLUS)             # 125.1 MeV
E_S = MP - E_CIRC                         # 813 MeV
SQRT_SIGMA_LATTICE = (420.0, 440.0)       # MeV
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    print(f"target: E_s = m_p - E_circ = {E_S:.1f} MeV (E_circ = {E_CIRC:.1f} MeV on R+ = {R_PLUS:.3f} fm)\n")

    # A. 1/L energies
    l_cav = HBARC * math.pi / E_S
    check("A. Cavity between two stowed ends (mode hbar c pi/l): 813 MeV at l = 0.763 fm = R+ (-3 %), the sheet's transverse mode at the nucleon's width; but a 1/L energy drops when the ends separate: right scale, wrong sign",
          abs(l_cav - 0.763) < 0.002 and abs(l_cav / R_PLUS - 0.967) < 0.005,
          f"l = {l_cav:.3f} fm = {l_cav/R_PLUS:.3f} R+; dE/dl = {-HBARC*math.pi/l_cav**2:.0f} MeV/fm (negative: separation releases energy)")

    # B. cutting cost
    w_cut = 2 * HBARC / (math.pi**2 * E_S)
    w_lep = 4 * LAMBDA_P / math.pi**2
    e_lep = 2 * HBARC / (math.pi**2 * w_lep)
    check("B. Cost of a cut (two new matched poles, 2 hbar c/(pi^2 w)): 813 MeV at w = 0.049 fm; the lepton rule w = 4 lambda-bar/pi^2 gives 469 MeV = m_p/2 by construction: nothing fixes w",
          abs(w_cut - 0.0492) < 0.0005 and abs(e_lep / MP - 0.5) < 1e-9,
          f"w(813 MeV) = {w_cut:.4f} fm; w_lep = {w_lep:.4f} fm -> {e_lep:.0f} MeV = {e_lep/MP:.2f} m_p")

    # C. confined flux
    delta2 = 1 / (math.pi**2 * ALPHA)                       # (delta/e)^2
    def sigma(area_fm2):                                     # MeV/fm: delta^2 e^2/(2 eps0 A) = 2 K delta^2 ... = 2 hbar c/(pi A)
        return delta2 * 4 * math.pi * K / (2 * area_fm2)
    sig_p = sigma(math.pi * LAMBDA_P**2)
    sig_formula = 2 * MP**2 / (math.pi**2 * HBARC)
    sqrt_sig = math.sqrt(sig_p * HBARC)
    lat_lo, lat_hi = (s**2 / HBARC for s in SQRT_SIGMA_LATTICE)
    L_tube = E_S / sig_p
    w_square = math.sqrt(2 * HBARC / (math.pi * sig_p))
    check("C. Flux of the matched pole confined to a tube of radius lambda-bar_p: sigma = 2 (m_p c^2)^2/(pi^2 hbar c) = 904 MeV/fm, sqrt(sigma) = (sqrt 2/pi) m_p = 422 MeV, inside the lattice range 420-440 MeV; 813 MeV = 0.90 fm of tube (r_p: -6.5 %). CONDITIONAL: the base's line spreads a pole's charge (1/L); keeping the flux in the sheet is a new rule, and lambda-bar_p is the proton's own length",
          abs(sig_p - sig_formula) < 1e-6 and abs(sqrt_sig - 422.4) < 0.5 and lat_lo <= sig_p <= lat_hi and abs(L_tube - 0.899) < 0.003 and abs(sig_p * RP / E_S - 0.935) < 0.005,
          f"sigma = {sig_p:.0f} MeV/fm (lattice {lat_lo:.0f}-{lat_hi:.0f}); sqrt(sigma) = {sqrt_sig:.0f} MeV; L = {L_tube:.3f} fm; sigma r_p = {sig_p*RP:.0f} MeV ({sig_p*RP/E_S-1:+.1%}); square sheet w = {w_square:.3f} fm; three arms of {L_tube/3:.2f} fm")

    # D. sign test
    dE_A = HBARC * math.pi * (1 / (l_cav + 1) - 1 / l_cav)
    dE_C = sig_p * 1.0
    check("D. Pulling the stowed ends 1 fm apart: -461 MeV under A (released), +904 MeV under C (paid): only the confined flux prevents separation",
          dE_A < 0 and abs(dE_A + 461) < 5 and abs(dE_C - 904) < 5, f"dE(A) = {dE_A:+.0f} MeV, dE(C) = {dE_C:+.0f} MeV")

    # E. breaking caveat
    L_break = MPI / sig_p
    check("E. Caveat under C: one pion ring (140 MeV) pays for 0.15 fm of tube, so a 0.9 fm tube holds six pions' worth -- the base needs a rule forbidding the break (QCD: pair cost ~1 GeV, breaking near 1.2 fm)",
          abs(L_break - 0.154) < 0.002 and abs(L_tube / L_break - 5.8) < 0.1,
          f"L_break = {L_break:.3f} fm; L_tube/L_break = {L_tube/L_break:.1f}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
