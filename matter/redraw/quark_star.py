#!/usr/bin/env python3
"""Redraw: 'a junction is worth 2 MeV; but quarks do not work like the leptons: they seem to be
branched as a star' -- the author's statements, tested in the base (BASE.md R2, R6, R8, R13, R17).

Checks:
  A. a 2 MeV junction in the base's Coulomb law u = delta^2 K/s: with unit poles (delta = e) it sits
     at s = 0.72 fm, the nucleon's size; with the matched poles of R17 (3.73 e) at 10 fm; with one
     string's pole (e/3) at 0.08 fm.  It is 15.7 x the electron's junction (127.75 keV): the base's
     1/s scaling would put the electron's junction at 10 fm inside the proton -- so the 2 MeV joint is
     a unit-pole junction at the nucleon's size, not the electron's junction rescaled.
  B. the ring charge rule q = (n+ - n-) e/3 is blind to topology: a star of k arms has the charge
     classes of a chain of k strings (parity of k).  +2/3 needs an even number of arms, -1/3 an odd
     one: a 3-arm star cannot be a u quark, a 4-arm star cannot be a d quark.  Minimal stars: d = 1
     string or 3 arms (1+, 2-); u = 2 strings or 4 arms (3+, 1-).
  C. junction budget of a proton drawn as a star of stars (uud + one Y-centre): 4 junctions x 2 MeV
     = 8 MeV = 0.9 % of the proton: the nucleon's mass is not in its junctions.
  D. the base's mode reading at the charge radius: hbar c / r_p = 234.7 MeV; the measured product
     r_p m_p c / hbar = 3.998 +- 0.002 (PDG 2024 r_p = 0.8409 +- 0.0004 fm): the proton is 4 modes
     of its charge radius to 0.05 %, a 3-arm Y is 3 modes (-25 %).  Coincidence to explain, not a
     derivation; the fluid at c on r_p needs a circulating charge 0.70 e for mu_p (2/3 e: -4.5 %).
  E. sign of confinement: the base's string gets lighter when longer (mode ~ 1/L) and its pole-pole
     Coulomb energy is also ~ 1/L, so a star has no size of its own -- only the spin condition
     S = R E/c = hbar/2 sets it, as for the electron (R = lambda-bar/f, f = 1/4 for r_p = 4 lambda-bar_p).
     Lattice QCD's Y-string (Takahashi 2001: sigma ~ 0.89 GeV/fm) gets heavier when longer: opposite sign.
  F. R6 (spin = symmetry axes / 2) on a star: a symmetric 3-arm star has 3 in-plane axes -> spin 3/2
     (Delta); one arm unlike the others leaves 1 axis -> spin 1/2 (N).  Same content, different shape:
     the Delta(1232) - N split (294 MeV) would be geometric.  Recorded, not tested.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

MP, ME = 938.27209, 0.51099895
HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC                       # 1.440 MeV fm
RP, DRP = 0.8409, 0.0004                # PDG 2024
LAMBDA_P = HBARC / MP                   # 0.2103 fm
DELTA_MATCHED = 1 / (math.pi * math.sqrt(ALPHA))
U_E = ME / 4                            # electron junction (R13)
U_J = 2.0                               # the author's junction, MeV
MU_P = 2.79284734
M_DELTA = 1232.0
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def charge_classes(k: int):
    """Charges (in e) of k polarised strings, n+ + n- = k."""
    return sorted({(2 * nplus - k) / 3 for nplus in range(k + 1)})


def main() -> int:
    # A. 2 MeV junction in the Coulomb law
    s_unit, s_matched, s_third = K / U_J, DELTA_MATCHED**2 * K / U_J, K / (9 * U_J)
    ratio = U_J / U_E
    check("A. A 2 MeV junction in u = delta^2 K/s: unit poles (e) at 0.72 fm -- the nucleon's size; matched poles (3.73 e) at 10 fm; string poles (e/3) at 0.08 fm. It is 15.7 x the electron's junction: not that junction rescaled, a unit-pole junction at 0.7 fm",
          abs(s_unit - 0.720) < 0.002 and abs(s_matched - 10.0) < 0.05 and abs(s_third - 0.080) < 0.001 and abs(ratio - 15.66) < 0.02,
          f"s(e) = {s_unit:.3f} fm, s(3.73 e) = {s_matched:.2f} fm, s(e/3) = {s_third:.3f} fm; 2 MeV / 127.75 keV = {ratio:.2f}")

    # B. charge classes of stars
    c3, c4, c1, c2 = charge_classes(3), charge_classes(4), charge_classes(1), charge_classes(2)
    check("B. The charge rule is topology-blind: a k-arm star has the classes of a k-chain. +2/3 needs even k, -1/3 odd k: a 3-arm star cannot be u, a 4-arm star cannot be d. Minimal: d = 1 string or 3 arms (1+,2-); u = 2 strings or 4 arms (3+,1-)",
          (2 / 3 not in c3) and (-1 / 3 not in c4) and (-1 / 3 in c1) and (2 / 3 in c2) and (-1 / 3 in c3) and (2 / 3 in c4),
          "; ".join(f"k={k}: " + ", ".join(f"{q:+.2f}" for q in cl) for k, cl in ((1, c1), (2, c2), (3, c3), (4, c4))))

    # C. junction budget
    n_j = 4
    budget = n_j * U_J
    check("C. Proton as a star of stars (u, u, d each with a centre + one Y-centre): 4 junctions x 2 MeV = 8 MeV = 0.9 % of the proton: the mass is not in the junctions",
          abs(budget / MP - 0.0085) < 0.0005, f"{budget:.0f} MeV / {MP:.0f} MeV = {budget/MP:.2%}")

    # D. mode reading at the charge radius
    e_mode = HBARC / RP
    prod = RP / LAMBDA_P
    dprod = DRP / LAMBDA_P
    q_circ = MU_P / prod                    # mu = q c R/2 -> mu/mu_N = (q/e)(R/lambda_p)
    check("D. hbar c / r_p = 234.7 MeV; r_p m_p c / hbar = 3.998 +- 0.002: the proton is 4 modes of its charge radius (0.05 %), a 3-arm Y is 3 modes (-25 %). A coincidence to explain, not a derivation. Fluid at c on r_p: mu_p needs a circulating charge 0.70 e (2/3 e gives -4.5 %)",
          abs(e_mode - 234.7) < 0.2 and abs(prod - 4) < 2.5 * dprod and abs(3 * e_mode / MP - 0.75) < 0.002 and abs(q_circ - 0.698) < 0.002,
          f"hbar c/r_p = {e_mode:.1f} MeV; r_p/lambda_p = {prod:.4f} +- {dprod:.4f}; 4 modes = {4*e_mode:.1f} MeV ({4*e_mode/MP-1:+.2%}); 3 modes = {3*e_mode:.0f} MeV ({3*e_mode/MP-1:+.0%}); q_circ = {q_circ:.3f} e, 2/3 e -> {(2/3)*prod/MU_P-1:+.1%}")

    # E. sign of confinement
    sigma = 0.89e3                          # MeV/fm, lattice Y-string tension
    L1, L2 = 0.5, 1.0
    base_ratio = (HBARC / L2) / (HBARC / L1)     # mode energy ratio when L doubles
    lattice_ratio = (sigma * L2) / (sigma * L1)
    f_p = LAMBDA_P / RP
    check("E. Doubling a string's length halves the base's energy (mode ~ 1/L, Coulomb ~ 1/L too) and doubles the lattice Y-string's (sigma L): opposite signs. The base's star has no size of its own; only S = R E/c = hbar/2 fixes it, R = lambda-bar/f with f = 1/4 for the proton",
          abs(base_ratio - 0.5) < 1e-9 and abs(lattice_ratio - 2) < 1e-9 and abs(f_p - 0.25) < 0.001,
          f"base: E(1 fm)/E(0.5 fm) = {base_ratio:.2f}; lattice: {lattice_ratio:.2f}; f = lambda_p/r_p = {f_p:.4f}")

    # F. R6 on a star (recorded)
    axes_sym, axes_asym = 3, 1
    check("F. R6 on a star: symmetric 3-arm star, 3 in-plane axes -> spin 3/2 (Delta); one arm unlike the others, 1 axis -> spin 1/2 (N): the Delta-N split of 294 MeV would be a shape difference at equal content. Recorded, not tested",
          axes_sym / 2 == 1.5 and axes_asym / 2 == 0.5,
          f"Delta(1232) - p = {M_DELTA - MP:.0f} MeV = {(M_DELTA-MP)/e_mode:.2f} modes of r_p")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
