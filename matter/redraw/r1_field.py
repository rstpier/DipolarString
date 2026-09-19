#!/usr/bin/env python3
"""Redraw: dig into R1 -- what does the field of the polarised fluid, and of a pole, do in the
base's line?  Is there a flux tube (R26) in the base as written?

The base's string is a bifilar line: two anti-parallel branches +-lambda at spacing D, wire radius
r, D/r = 2 cosh pi (matched to Z_0).  Checks:
  A. the transverse field of the bifilar line is confined by geometry: the energy density falls as
     1/rho^4 beyond D; numerically, 92 % of the field energy per unit length lies within rho = D and
     99 % within 3 D.  The base's line is its own sheet: no new rule is needed to confine the
     POLARISED fluid's field.
  B. but that confined field gives a tension only if the charge per unit length is fixed.  With the
     ring rule (e/3 per string) and free string length, U = (e/3)^2/(eps0 l) falls as 1/l: no
     tension.  With the string length fixed at the ladder's l_1(9) = 0.81 fm, each string stores
     4 pi K/(9 l_1) = 2.47 MeV (the 2.3 MeV 'extra string' of R24/R29, +8 %), a tension of only
     3.0 MeV/fm: 300 x below the lattice.
  C. a pole's monopole field inside a conducting tube is screened, not channelled: it decays as
     exp(-2.405 z/a) along the tube, within 0.09 fm for a = lambda-bar_p, with a finite local energy
     ~ K delta^2/a ~ 95 MeV.  R26's tube (flux carried the whole length) is not what a conducting
     fluid does.
  D. what the lattice tension would need: a fixed line charge of 7.1 e/fm on the branches, i.e.
     strings of e/3 only 0.047 fm long (19 per 0.9 fm), which the ladder places at n = 14 strands,
     not the 9 of the three-strand proton; or a fixed current of 4.8e5 A, 13 x the base's e c/(2 pi lambda-bar_p).
  E. verdict: the base cannot hold both 'e/3 per string of free length' and a flux-tube tension;
     R26 is not a mechanism of the base as written.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC                       # MeV fm
MP = 938.27209
LAMBDA_P = HBARC / MP
D_OVER_R = 2 * math.cosh(math.pi)
L1_9 = 2 * math.pi * (HBARC / 0.51099895) / 3 * (3 / 9) ** (2 * math.pi)   # 0.813 fm
SIGMA_LAT = 904.0                       # MeV/fm (R26)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def field_energy_fraction(D=1.0, rho_max=1.0, n=1400):
    """Fraction of the bifilar line's transverse field energy (per unit length) inside rho < rho_max,
    for two line charges +-1 at x = +-D/2, wire radius r = D / D_OVER_R (interiors excluded)."""
    r = D / D_OVER_R
    R_big = 400 * D
    # log-polar grid
    rho = np.exp(np.linspace(math.log(r * 0.2), math.log(R_big), n))
    phi = np.linspace(0, 2 * math.pi, 720, endpoint=False)
    RR, PP = np.meshgrid(rho, phi, indexing="ij")
    X, Y = RR * np.cos(PP), RR * np.sin(PP)
    def E(x0):
        dx, dy = X - x0, Y
        d2 = dx**2 + dy**2
        return dx / d2, dy / d2
    Ex1, Ey1 = E(+D / 2)
    Ex2, Ey2 = E(-D / 2)
    u = (Ex1 - Ex2) ** 2 + (Ey1 - Ey2) ** 2
    inside_wire = ((X - D / 2) ** 2 + Y**2 < r**2) | ((X + D / 2) ** 2 + Y**2 < r**2)
    u = np.where(inside_wire, 0.0, u)
    dA = RR * RR * math.log(rho[1] / rho[0]) * (phi[1] - phi[0])
    total = np.sum(u * dA)
    inside = np.sum(np.where(RR < rho_max, u * dA, 0.0))
    return inside / total


def main() -> int:
    # A. confinement by geometry
    f1, f3, f10 = (field_energy_fraction(rho_max=m) for m in (1.0, 3.0, 10.0))
    check("A. Bifilar line (D/r = 2 cosh pi): the transverse field energy is confined by geometry -- 92 % within rho = D, 99 % within 3 D, energy density ~ 1/rho^4 beyond: the polarised fluid's field needs no new rule to stay on the line",
          abs(f1 - 0.92) < 0.01 and f3 > 0.98 and f10 > 0.995,
          f"fraction inside D: {f1:.3f}, 3D: {f3:.4f}, 10D: {f10:.5f}")

    # B. tension only with fixed line charge
    u_string_free = lambda l: 4 * math.pi * K / (9 * l)          # (e/3)^2/(eps0 l) = 4 pi K/(9 l)
    U_09 = u_string_free(L1_9)
    tension_09 = U_09 / L1_9
    check("B. With e/3 per string and free length, U = (e/3)^2/(eps0 l) = 4 pi K/(9 l) falls as 1/l: no tension. At the fixed l_1(9) = 0.81 fm each string stores 2.47 MeV (R24's 2.3 MeV extra string, +8 %), a tension of 3.0 MeV/fm, 300 x below the lattice",
          abs(U_09 - 2.474) < 0.01 and abs(tension_09 - 3.04) < 0.05 and SIGMA_LAT / tension_09 > 250,
          f"U(l_1 = {L1_9:.3f} fm) = {U_09:.3f} MeV; sigma = {tension_09:.2f} MeV/fm = lattice/{SIGMA_LAT/tension_09:.0f}; U(2 l) / U(l) = {u_string_free(2*L1_9)/U_09:.2f}")

    # C. screening of a monopole in a conducting tube
    a = LAMBDA_P
    decay = a / 2.405
    E_local = K * DELTA**2 / a
    check("C. A pole inside a conducting tube of radius lambda-bar_p is screened, not channelled: its field decays as exp(-2.405 z/a), within 0.09 fm, with a local energy ~ K delta^2/a ~ 95 MeV; the flux never travels the 0.9 fm of R26's tube",
          abs(decay - 0.0874) < 0.001 and abs(E_local - 95.2) < 1 and decay < 0.1 * 0.9,
          f"decay length = {decay:.4f} fm; K delta^2/a = {E_local:.1f} MeV")

    # D. what the lattice tension needs
    lam_needed = math.sqrt(SIGMA_LAT / (4 * math.pi * K))          # e/fm
    l_string = (1 / 3) / lam_needed
    n_per_tube = 0.9 / l_string
    n_ladder = 3 * (2 * math.pi * (HBARC / 0.51099895) / 3 / l_string) ** (1 / (2 * math.pi))
    I_base = 1.602176634e-19 * 2.99792458e8 / (2 * math.pi * LAMBDA_P * 1e-15)          # A
    sigma_SI = SIGMA_LAT * 1.602176634e-13 / 1e-15                                       # J/m
    I_needed = math.sqrt(2 * sigma_SI / (4 * math.pi * 1e-7))
    check("D. The lattice tension needs a fixed line charge of 7.1 e/fm: strings of e/3 only 0.047 fm long (19 per 0.9 fm), which the ladder puts at n = 14 strands, not 9; or a fixed current of 4.8e5 A, 13 x the base's e c/(2 pi lambda-bar_p)",
          abs(lam_needed - 7.07) < 0.05 and abs(l_string - 0.0471) < 0.0005 and abs(n_ladder - 14.2) < 0.2 and abs(I_needed / I_base - 13.2) < 0.2,
          f"lambda = {lam_needed:.2f} e/fm, string length {l_string:.4f} fm, {n_per_tube:.0f} per 0.9 fm, ladder n = {n_ladder:.1f}; I = {I_needed:.2e} A = {I_needed/I_base:.1f} x {I_base:.2e} A")

    # E. verdict
    check("E. Verdict: the base cannot hold both 'e/3 per string of free length' and a flux-tube tension; a conducting fluid screens a pole instead of channelling it: R26's tube is not a mechanism of the base as written",
          SIGMA_LAT / tension_09 > 100 and decay < 0.1, "see B, C, D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
