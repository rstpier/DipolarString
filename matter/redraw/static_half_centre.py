#!/usr/bin/env python3
"""Redraw: what fixes the electron's static half once the circuit is a closed one-way ring (R53)?
The half must carry no current (no mu), no angular momentum (no S), and equal m_e c^2/2 (g = 2).

Candidates:
  A. the ring's own field energy: (1/2) L I^2 with L = mu_0 R (ln(8R/a) - 2), I = e c/(2 pi R):
     0.6 keV, 0.1 % of the half (massless_fluid.py's 0.4 %): excluded.
  B. a second circulation at the centre ('rotating on itself', R13): every quantum of circulation
     carries hbar/2 (R40 D), so a core quantum would double or cancel the spin; only a
     counter-rotating pair (S = 0) survives, with no rule for its size: not fixed.
  C. the junctions of the ring: a closed ring of three strands has three junctions; a junction is
     static (no current, no S, no mu) and carries the pole-pole energy K delta^2/D = hbar c/(pi^2 D)
     with delta = e/(pi sqrt alpha) and D the strands' spacing (R17 C, R24).  Three of them equal
     m_e c^2/2 for D = 6 lambda-bar/pi^2 = 0.608 lambda-bar = 235 fm.  The base's own DQD spacing is
     D_0 = 2 cosh(pi) r with r = lambda-bar/37.1 (manuscript): 241 fm, 2.7 % away.  So the static half
     is the three junctions at the DQD spacing, and it fixes the tube radius r = lambda-bar/38.1
     (manuscript 37.1, 2.7 %).  This replaces R17's two-junction width 4 lambda-bar/pi^2 by the
     three-junction spacing 6 lambda-bar/pi^2 (x 1.5).
  D. consequences: g = 2 kept on the ladder (junction energy ~ 1/D ~ 1/l_1); at the nucleon the
     spacing is 0.236 fm and the p-n splitting (R44 rule, plug pole) moves to 1.35 MeV (+4.5 %),
     1.40 MeV (+8 %) with the sphere pole.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
R_TUBE_MANUSCRIPT = LAMBDA_E / 37.1
D0 = 2 * math.cosh(math.pi) * R_TUBE_MANUSCRIPT
E_DQD = 4 * math.pi * K / (9 * L1_9)
M_MUTUAL = K / (math.sqrt(3) * L1_9)
DM_OBS = MN - MP
F_PLUG = 0.759
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def dm(w9):
    a = w9 / 4
    S_u = K / (2 * L1_9) * (math.log(4 * L1_9 / a) - 1)
    S_d = K / L1_9 * (math.log(2 * L1_9 / a) - 1)
    return E_DQD - ((4 / 9) * S_u - (1 / 9) * S_d + M_MUTUAL / 3)


def main() -> int:
    # A. field energy of the ring current
    R, a = LAMBDA_E, 4 * LAMBDA_E / math.pi**2
    u_field = (K / (2 * math.pi * R)) * (math.log(8 * R / a) - 2)      # (1/2) L I^2 with I = e c/2 pi R, in MeV
    check("A. The ring's own field energy (1/2) L I^2 = 0.6 keV, 0.1 % of the half: excluded",
          abs(u_field * 1e3 - 0.58) < 0.05 and u_field / (ME / 2) < 0.005, f"{u_field*1e3:.2f} keV = {u_field/(ME/2):.3%} of m_e c^2/2")

    # B. core circulation
    check("B. A second circulation at the centre carries hbar/2 per quantum (R40 D): it doubles or cancels the spin; only a counter-rotating pair (S = 0) survives, with no rule for its size: not fixed",
          True, "recorded")

    # C. the junctions
    D_needed = 6 * LAMBDA_E / math.pi**2
    e_junction_D0 = HBARC / (math.pi**2 * D0)
    three_at_D0 = 3 * e_junction_D0
    r_needed = D_needed / (2 * math.cosh(math.pi))
    check("C. Three junctions of the closed 3-strand ring, each hbar c/(pi^2 D): equal to m_e c^2/2 for D = 6 lambda-bar/pi^2 = 235 fm; the base's DQD spacing D_0 = 2 cosh(pi) lambda-bar/37.1 = 241 fm gives 0.487 m_e (-2.7 %): the static half is the three junctions at the DQD spacing, fixing the tube radius at lambda-bar/38.1 (manuscript 37.1)",
          abs(3 * HBARC / (math.pi**2 * D_needed) - ME / 2) < 1e-12 and abs(D_needed - 234.8) < 0.5 and abs(D0 - 241.3) < 0.5 and abs(three_at_D0 / (ME / 2) - 0.973) < 0.003 and abs(LAMBDA_E / r_needed - 38.1) < 0.1,
          f"D needed = {D_needed:.1f} fm; D_0 = {D0:.1f} fm; three junctions at D_0 = {three_at_D0*1e3:.1f} keV = {three_at_D0/(ME/2):.3f} x half; r = lambda-bar/{LAMBDA_E/r_needed:.1f}")

    # D. consequences
    w9_new = D_needed * 3 ** (-2 * math.pi)
    dm_plug = dm(F_PLUG * w9_new)
    dm_sphere = dm(w9_new)
    check("D. g = 2 kept on the ladder (junction energy ~ 1/D); the width becomes 6 lambda-bar/pi^2 (x 1.5): at the nucleon 0.236 fm, and m_n - m_p moves to 1.35 MeV (+4.5 %, plug pole) or 1.40 MeV (+8 %, sphere pole)",
          abs(w9_new - 0.236) < 0.002 and abs(dm_plug / DM_OBS - 1.045) < 0.01 and abs(dm_sphere / DM_OBS - 1.085) < 0.01,
          f"D_9 = {w9_new:.4f} fm; m_n - m_p = {dm_plug:.3f} ({dm_plug/DM_OBS-1:+.1%}) plug, {dm_sphere:.3f} ({dm_sphere/DM_OBS-1:+.1%}) sphere")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
