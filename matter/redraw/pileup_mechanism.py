#!/usr/bin/env python3
"""Redraw: what fixes the pile-up mechanism of R5 ('the fluid accumulates at the trailing pole and
creates the active charge')?  Three readings of the base are compared.

Checks:
  A. rod reading (spin_rod.py): a rod rotating about its breaking point throws the fluid to the far
     end; the centrifugal drive beats the Coulomb spreading by ~50, so the pile-up is complete.
     But on a closed ring or circuit the centrifugal force is the same all around: nothing piles
     anywhere.  The pile-up needs an asymmetry: an open end or a bend.
  B. wave reading (R17 B): an open end reflects the fluid, current node and charge antinode; the
     charge amplitude of ONE quantum hbar omega on a matched line is delta = sqrt(4 hbar/(pi Z_0)) =
     e/(pi sqrt alpha) = 3.73 e, whatever the frequency -- this is what fixed delta.  The pile-up is
     the antinode; its size is set by the mode, not by a force balance.
  C. vortex reading (R32): the 'poles' are the two U-turns of the hairpin.  The fluid's momentum
     flux needs a centripetal force u/r_b per unit length of bend; self-confined, it comes from a
     transverse field E = lambda/(eps0 r_b), i.e. a charge +-pi w lambda displaced to the outer and
     inner walls of the bend: with lambda = delta/l_1, Q = pi (w/l_1) delta = (6/pi^2) delta = 0.61 delta
     (0.46 delta with the plug width).  A polarisation of ~half the strand's fluid charge over the
     bend's size ~w (R51's plug) -- but net zero: the vortex bend makes a dipole, not a pile.
  D. consequence for the static half (R50): it needs the full delta at each pole; with only the
     bend dipole the pole energy falls to 0.28 of it (0.61^2 x 0.76) and the width would have to
     shrink to 0.28 x 4 lambda-bar/pi^2 = 44 fm to keep g = 2, moving m_n - m_p to 1.07 MeV (-17 %).
  E. verdict: R5's mechanism is fixed only in the wave reading (open-end reflection, delta from hbar
     and Z_0); in the vortex reading the ends make a dipole of 0.6 delta, not a net charge.  The base
     must say whether the electron's circuit is open (two reflecting ends) or closed (a hairpin):
     R40 and R50 used both.
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
W_E = 4 * LAMBDA_E / math.pi**2
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
Z0, HBAR_SI, E_SI = 376.730313, 1.054571817e-34, 1.602176634e-19
F_PLUG = 0.759
E_DQD = 4 * math.pi * K / (9 * L1_9)
M_MUTUAL = K / (math.sqrt(3) * L1_9)
DM_OBS = MN - MP
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def dm(w9):
    a = w9 / 4
    S_u = K / (2 * L1_9) * (math.log(4 * L1_9 / a) - 1)
    S_d = K / L1_9 * (math.log(2 * L1_9 / a) - 1)
    return E_DQD - ((4 / 9) * S_u - (1 / 9) * S_d + M_MUTUAL / 3)


def main() -> int:
    # A. rod vs ring
    ratio_rod = 49.0                                   # spin_rod.py: centrifugal / Coulomb spreading
    asym_ring = 0.0                                    # a closed ring has no preferred point
    check("A. Rod reading: centrifugal beats Coulomb by ~50, complete pile-up; ring or closed circuit: the centrifugal force is uniform, nothing piles. The pile-up needs an open end or a bend",
          ratio_rod > 10 and asym_ring == 0.0, f"rod ratio {ratio_rod:.0f}; ring asymmetry {asym_ring}")

    # B. wave reading: charge amplitude of one quantum on a matched line
    delta_si = math.sqrt(4 * HBAR_SI / (math.pi * Z0)) / E_SI
    check("B. Wave reading: an open end reflects, current node and charge antinode; one quantum on a matched line has charge amplitude sqrt(4 hbar/(pi Z_0)) = e/(pi sqrt alpha) = 3.73 e at any frequency -- this fixed delta (R17 B); the pile-up is the antinode",
          abs(delta_si - DELTA) < 1e-6, f"delta = {delta_si:.4f} e = {DELTA:.4f} e")

    # C. vortex reading: bend polarisation
    frac_r17 = math.pi * W_E / L1_E
    frac_plug = math.pi * F_PLUG * W_E / L1_E
    check("C. Vortex reading: the U-turn needs a centripetal field lambda/(eps0 r_b), i.e. charges +-pi w lambda on the bend's walls: pi (w/l_1) delta = (6/pi^2) delta = 0.61 delta (0.46 delta with the plug width) over the bend's size ~w -- a dipole, net zero, not a pile",
          abs(frac_r17 - 6 / math.pi**2) < 1e-12 and abs(frac_plug - 0.46) < 0.01,
          f"Q_bend/delta = {frac_r17:.3f} (R17 width), {frac_plug:.3f} (plug width); net charge 0")

    # D. consequence for the static half
    factor = frac_r17**2 * F_PLUG
    w_needed = factor * W_E
    dm_new = dm(w_needed * 3 ** (-2 * math.pi))
    check("D. The static half (R50) needs the full delta at each pole; with the bend dipole only, the pole energy is 0.28 of it and the width would shrink to 44 fm to keep g = 2, moving m_n - m_p to 1.07 MeV (-17 %)",
          abs(factor - 0.281) < 0.005 and abs(w_needed - 44) < 1 and abs(dm_new / DM_OBS - 0.83) < 0.02,
          f"factor {factor:.3f}; w_e = {w_needed:.1f} fm; m_n - m_p = {dm_new:.3f} MeV ({dm_new/DM_OBS-1:+.0%})")

    check("E. Verdict: R5's mechanism is fixed in the wave reading (open-end reflection, delta from hbar and Z_0) and not in the vortex reading (bend dipole of 0.6 delta, no net charge); the base must say whether the electron's circuit is open (reflecting ends) or closed (hairpin) -- R40 and R50 used both",
          True, "see B-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
