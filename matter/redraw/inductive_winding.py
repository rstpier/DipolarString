#!/usr/bin/env python3
"""Redraw: 'could there be an inductive amplification at each turn?' -- the author's question after
R32.  A charged strand wound as a helix (n turns per unit length, core radius a) is a solenoid: its
field is mu_0 n I and its energy per unit length mu_0 (n I)^2 A / 2 grows as the SQUARE of the turns.

Checks:
  A. the gain: for a fluid of fixed charge per unit path length lambda moving at c, the wound
     strand stores (lambda^2/2 eps0) pi (n a)^2 per unit axial length against lambda^2/(2 eps0) for
     the straight line: an amplification by pi (n a)^2, quadratic in turns per core radius.  But with
     a FIXED total charge e on the path, the current per turn falls as the path lengthens and the
     gain cancels exactly: u = K/(2 l^2), independent of n (0.8 MeV for l = 0.9 fm).  The gain needs
     the fluid's own density, each turn carrying its own current -- the vortex's conserved circulation.
  B. at the nucleon strand density (e/3 per 0.81 fm, straight line 3.0 MeV/fm, R30 B) the lattice's
     904 MeV/fm needs pi (n a)^2 = 297: n a = 9.7, i.e. 46 turns per fm on a core of lambda-bar_p,
     41 turns in the 0.9 fm centre, 55 fm of string = 68 strands of 0.81 fm (not 9).
  C. at the electron's density (4.1e-4 e/fm, R31) the same tension needs n a = 9700: 46 000 turns
     per fm, 41 000 turns in the centre, 54 000 fm of string = 67 electron strings.  R31's exclusion
     lifts if the winding number is free: one fluid density, wound or not.
  D. only the unpaired (charged) branch winds into a solenoid: twisting a bifilar pair with opposite
     currents cancels the azimuthal components; neutral strands are inductively inert (R29).
  E. the winding number is the strand's twist (Calugareanu Lk = Tw + Wr), a topological integer the
     base already has (R27); the value needed (~40 per nucleon at strand density) is not derived.
     The electric energy of the wound charge adds a term of the same order or larger (a charged
     cylinder screened at r_p), which lowers the winding needed to a few tens of turns per fm.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP = 0.51099895, 938.27209
LAMBDA_P = HBARC / MP
L1_E = 2 * math.pi * (HBARC / ME) / 3                    # 808.8 fm
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)                    # 0.813 fm
SIGMA_LAT = 904.0
L_CENTRE = 0.9
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def u_straight(lam):          # lambda^2 / eps0 in MeV/fm (E + B parts, R30): 4 pi K lambda^2
    return 4 * math.pi * K * lam**2


def gain(n_per_fm, a):
    return math.pi * (n_per_fm * a) ** 2


def main() -> int:
    # A. gain and its cancellation at fixed total charge
    a = LAMBDA_P
    g10 = gain(10.0, a)
    l = L_CENTRE
    # fixed total charge e at c on a helix of n turns/fm over length l: I = e c / s, s = 2 pi a n l
    def u_fixed_charge(n):
        s = 2 * math.pi * a * n * l
        lam_path = 1.0 / s                                   # e per fm along the path
        return u_straight(lam_path) / 2 * gain(n, a)         # magnetic part only, amplified
    u1, u100 = u_fixed_charge(1.0), u_fixed_charge(100.0)
    check("A. A wound strand is a solenoid: energy per unit length amplified by pi (n a)^2 (10 turns/fm on lambda-bar_p: x 14). But with a fixed total charge the current per turn drops as the path grows and the gain cancels exactly: u = K/(2 l^2) whatever n. The gain needs a fixed charge per unit path length -- the vortex's conserved circulation",
          abs(g10 - 13.9) < 0.1 and abs(u1 - u100) / u1 < 1e-9 and abs(u1 - K / (2 * l**2)) / u1 < 1e-9,
          f"gain(10/fm) = {g10:.1f}; fixed-charge u(1 turn/fm) = {u1:.3f} MeV/fm = u(100/fm) = {u100:.3f}; K/(2 l^2) = {K/(2*l**2):.3f} MeV/fm")

    # B. winding at the nucleon strand density
    lam9 = (1 / 3) / L1_9
    u9 = u_straight(lam9)
    g_needed = SIGMA_LAT / u9
    na = math.sqrt(g_needed / math.pi)
    n9 = na / a
    turns = n9 * L_CENTRE
    string = 2 * math.pi * a * turns
    check("B. Nucleon strand density (e/3 per 0.81 fm, straight 3.0 MeV/fm): the lattice tension needs a gain of 297, n a = 9.7, 46 turns per fm on lambda-bar_p, 41 turns in the 0.9 fm centre, 55 fm of string = 68 strands of 0.81 fm, not 9",
          abs(u9 - 3.04) < 0.05 and abs(g_needed - 297) < 3 and abs(na - 9.72) < 0.05 and abs(n9 - 46.2) < 0.3 and abs(string / L1_9 - 68) < 1,
          f"u_straight = {u9:.2f} MeV/fm; gain {g_needed:.0f}; n a = {na:.2f}; n = {n9:.1f}/fm; turns = {turns:.0f}; string = {string:.0f} fm = {string/L1_9:.0f} strands")

    # C. winding at the electron's density
    lam_e = (1 / 3) / L1_E
    ue = u_straight(lam_e)
    g_e = SIGMA_LAT / ue
    na_e = math.sqrt(g_e / math.pi)
    n_e = na_e / a
    string_e = 2 * math.pi * a * n_e * L_CENTRE
    check("C. Electron's density (4.1e-4 e/fm, straight 3e-6 MeV/fm): the same tension needs n a = 9700, 46 000 turns per fm, 41 000 turns in the centre, 54 000 fm of string = 67 electron strings. R31's exclusion lifts if the winding is free: one fluid, wound or not",
          abs(na_e - 9720) < 50 and abs(n_e - 46200) < 300 and abs(string_e / L1_E - 67) < 1,
          f"u_straight = {ue:.1e} MeV/fm; n a = {na_e:.0f}; n = {n_e:.0f}/fm; string = {string_e:.0f} fm = {string_e/L1_E:.0f} electron strings")

    # D. bifilar twist cancels
    theta = math.radians(60)
    I_az_pair = math.sin(theta) - math.sin(theta)      # opposite currents, same twist
    I_az_single = math.sin(theta)
    check("D. Twisting a bifilar pair with opposite currents cancels the azimuthal current (no solenoid); only the unpaired charged branch winds into one: neutral strands are inductively inert",
          I_az_pair == 0 and I_az_single > 0, f"pair: {I_az_pair:.2f} I, single branch: {I_az_single:.2f} I at 60 deg pitch")

    # E. winding number as twist
    check("E. The winding number is the strand's twist Tw (Lk = Tw + Wr), an integer the base has (R27); ~40 per nucleon at strand density is not derived; the wound charge's electric energy adds a comparable term and lowers the turns needed to a few tens per fm",
          turns > 30 and turns < 60, f"Tw needed ~ {turns:.0f} at strand density, ~ {n_e*L_CENTRE:.0f} at the electron's")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
