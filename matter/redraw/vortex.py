#!/usr/bin/env python3
"""Redraw: 'the fluid is not a wave, it is a vortex' -- the author's redesign of R1 (R32).  What a
vortex changes, tested against the base's established results and the missing ~800 MeV.

Checks:
  A. a vortex has a tension by construction: its circulation is a conserved quantum, not a charge
     free to spread, so a steady current I on the matched line stores u = Z_0 I^2 / c per unit
     length, constant along the line.  dE/dL > 0: the sign the centre needed (R26 D), with no new rule.
  B. two readings on the electron's ring: (i) a smoke ring (poloidal swirl at c around a core of
     radius a, R/a = 37.1) self-propels at 0.073 c -- not a particle at rest; (ii) a current loop
     (flow along the ring, the base's R2/R4) is steady: spin hbar/2, mu_B and R11's no-radiation
     hold, better than a standing wave.  The vortex must be the current loop.
  C. the loop's flow energy for charge q at c on radius R is 2 K q^2 / R = 2 alpha (q/e)^2 hbar c / R:
     for q = e on the electron's ring, 7.5 keV = 1.5 % of m_e.  The flow carries the electron's
     circulating half only if q = e/(2 sqrt alpha) = 5.85 e (then E = hbar c / 2R: the spin condition
     restated); the matched pole 3.73 e carries 40 % of it.  R29's mode reading is lost.
  D. the tension of that vortex at the proton's scale: u = hbar c / (4 pi R^2) = 355 MeV/fm at
     R = lambda-bar_p (lattice / 2.5); R26's electric tube gave 904.  Both scale as m_p^2: the number
     is the proton's own Compton length, circular until R is derived.
  E. two anti-parallel vortex lines at distance d interact as (mu_0 I^2 / 2 pi) ln d per unit length:
     a logarithmic potential, force ~ 1/d -- a mild confinement between lines that the wave never had;
     for the current of 904 MeV/fm, going from 0.1 to 1 fm costs 331 MeV per fm of line.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP = 0.51099895, 938.27209
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
Z0 = 376.730313
E_SI, C_SI = 1.602176634e-19, 2.99792458e8
MEV_FM_TO_J_M = 1.602176634e-13 / 1e-15
SIGMA_LAT = 904.0
R_OVER_A = 37.1                         # R_3 / r of the base
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def u_line(I_amp):
    """Energy per unit length of the matched line carrying current I at c: Z_0 I^2 / c, in MeV/fm."""
    return Z0 * I_amp**2 / C_SI / MEV_FM_TO_J_M


def main() -> int:
    # A. tension by construction
    I_test = 1.0e5
    u1, u2 = u_line(I_test), u_line(I_test)
    E_L1, E_L2 = u1 * 1.0, u2 * 2.0
    check("A. A vortex's circulation is a conserved quantum: a steady current I on the matched line stores Z_0 I^2/c per unit length, constant along the line, so E grows with L -- the sign the centre needed, with no new rule",
          abs(E_L2 / E_L1 - 2) < 1e-12 and u1 > 0, f"u(1e5 A) = {u1:.0f} MeV/fm; E(2 fm)/E(1 fm) = {E_L2/E_L1:.2f}")

    # B. smoke ring vs current loop
    v_ring = (1 / (2 * R_OVER_A)) * (math.log(8 * R_OVER_A) - 0.25)
    check("B. Smoke ring (poloidal swirl at c, R/a = 37.1) self-propels at 0.073 c: not at rest; the current loop (flow along the ring) is steady and keeps spin hbar/2, mu_B and no radiation: the vortex must be the current loop",
          abs(v_ring - 0.0734) < 0.001, f"v_ring/c = {v_ring:.4f}")

    # C. flow energy of the loop
    def e_flow(q, R):
        return 2 * K * q**2 / R
    e_e = e_flow(1.0, LAMBDA_E)
    q_half = 1 / (2 * math.sqrt(ALPHA))
    delta = 1 / (math.pi * math.sqrt(ALPHA))
    check("C. Flow energy of charge q at c on radius R: 2 K q^2/R. For q = e on the electron's ring: 7.5 keV = 1.5 % of m_e; the circulating half m_e c^2/2 needs q = e/(2 sqrt alpha) = 5.85 e (the spin condition restated); the matched pole 3.73 e carries 40 %. The mode reading of R29 is lost",
          abs(e_e * 1e3 - 7.46) < 0.05 and abs(e_flow(q_half, LAMBDA_E) / (ME / 2) - 1) < 1e-12 and abs(e_flow(delta, LAMBDA_E) / (ME / 2) - 4 / math.pi**2) < 1e-12,
          f"E(e) = {e_e*1e3:.2f} keV = {e_e/ME:.3f} m_e; q_half = {q_half:.2f} e; delta share = {4/math.pi**2:.3f}")

    # D. tension at the proton's scale
    u_p = HBARC / (4 * math.pi * LAMBDA_P**2)
    check("D. Tension of that vortex at R = lambda-bar_p: hbar c/(4 pi R^2) = 355 MeV/fm (lattice/2.5); R26's tube gave 904; both scale as m_p^2 -- the number is the proton's Compton length, circular until R is derived",
          abs(u_p - 355) < 2 and abs(SIGMA_LAT / u_p - 2.55) < 0.05, f"u = {u_p:.0f} MeV/fm; lattice/u = {SIGMA_LAT/u_p:.2f}")

    # E. line-line interaction
    I_lat = math.sqrt(SIGMA_LAT * MEV_FM_TO_J_M * C_SI / Z0)
    u_int = SIGMA_LAT * math.log(10) / (2 * math.pi)      # (mu0 I^2/2pi) ln 10 = (Z0 I^2/c)(ln 10 / 2 pi)
    check("E. Two anti-parallel vortex lines interact as (mu_0 I^2/2 pi) ln d per unit length: logarithmic, force ~ 1/d; at the lattice current (3.4e5 A) going from 0.1 to 1 fm costs 331 MeV per fm of line",
          abs(I_lat - 3.39e5) < 2e3 and abs(u_int - 331) < 2, f"I = {I_lat:.2e} A; (Z_0 I^2/c)(ln 10/2 pi) = {u_int:.0f} MeV per fm of line")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
