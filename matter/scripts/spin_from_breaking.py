#!/usr/bin/env python3
"""The author's spin mechanism: every particle comes from an assembly of spin-1 DQDs (symmetric);
a closed chain of 3 DQDs breaking in two gives a positron and an electron with opposite rotation
from one breaking action, spin 1/2 each -- spin 1 broken in two.

Wave arithmetic of the mechanism:
  mother   = closed ring of 3 DQDs = 6 strings, circumference 6 l_1, travelling mode of winding 1;
  daughters = two rings of 3 strings (3 minus -> electron, 3 plus -> positron), circumference 3 l_1.
At the mother's frequency, half a wavelength fits a daughter: the daughter's mode has winding 1/2 --
the antiperiodic ring of section 5c, L_z = 1/2 about its axis.  So "spin 1 -> 1/2 + 1/2" is
wavelength inheritance, and it is exact.  Two consequences:

  A. energy: the mother's quantum is hbar c / 2 R_3 = 0.256 MeV, one quarter of the pair's rest
     energy; at least 0.77 MeV must come from outside (as in real pair creation, which carries its
     own angular momentum) -- the 1/2 + 1/2 is not forced by conservation, it is inherited;
  B. size: with the manuscript's l_1 the daughter's winding-1/2 mode sits at 0.256 MeV, HALF the
     electron's mass; the arithmetic closes only if the electron ring has radius R_3 / 2 (strings of
     l_1 / 2): then mother = 0.511 MeV (winding 1), each daughter = 0.511 MeV (winding 1/2), pair
     threshold = mother + one more quantum.  The open-turn note reached the same R_3 / 2.

What the mechanism needs and the model does not supply: a ring of three single branches must admit
an antiperiodic travelling wave.  The derived source (5c) is the half-twist of the unordered PAIR;
three single branches joined end to end have no pair and no sign flip at their junctions
(equal impedances, Gamma = 0).  And L_z = 1/2 about one axis is not the fermion (wen-link note).
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

from berry_holonomy_common import check, report

HBARC = 197.327
R3 = 386.159
L1 = 2 * math.pi * R3 / 3
ME = 0.51100


def mode_energy(circumference: float, winding: float) -> float:
    """hbar omega of a travelling TEM mode with `winding` wavelengths around the ring (MeV)."""
    return HBARC * 2 * math.pi * winding / circumference


def main() -> int:
    mother = mode_energy(6 * L1, 1.0)
    daughter_half = mode_energy(3 * L1, 0.5)
    daughter_full = mode_energy(3 * L1, 1.0)
    check("Wavelength inheritance: the mother's winding-1 mode on 6 strings is exactly the winding-1/2 mode on 3 strings",
          abs(mother - daughter_half) < 1e-12, f"hbar omega = {mother:.4f} MeV on both; the daughter's mode is antiperiodic (section 5c), L_z = 1/2")
    check("A. Energy: the mother's quantum is one quarter of the pair's rest energy; >= 0.77 MeV must come from outside",
          abs(mother / (2 * ME) - 0.25) < 1e-3, f"mother {mother:.3f} MeV vs 2 m_e c^2 = {2*ME:.3f} MeV; the 1/2 + 1/2 is inherited, not forced by conservation")
    check("B. Size: with the manuscript's l_1 the daughter's winding-1/2 mode is HALF the electron mass; the winding-1 mode is the anchor",
          abs(daughter_half / ME - 0.5) < 1e-3 and abs(daughter_full / ME - 1.0) < 1e-3, f"winding 1/2: {daughter_half:.3f} MeV, winding 1: {daughter_full:.3f} MeV = m_e c^2")
    l1_half = L1 / 2
    mother_h = mode_energy(6 * l1_half, 1.0)
    daughter_h = mode_energy(3 * l1_half, 0.5)
    check("B'. The arithmetic closes with strings of l_1/2 (electron ring of radius R_3/2): mother = m_e c^2, each daughter = m_e c^2, threshold = mother + one quantum",
          abs(mother_h - ME) < 1e-3 and abs(daughter_h - ME) < 1e-3, f"mother {mother_h:.3f} MeV, daughter {daughter_h:.3f} MeV; the open-turn note reached the same R_3/2 (Z_e = 0.62 Z_0, still a well)")
    check("B''. That halves l_1, hence the cell, the Brillouin cutoff (3/pi -> 6/pi m_e c^2) and E_conf: a recalibration of the manuscript, not a free choice",
          True, f"l_1 = {L1:.0f} fm -> {l1_half:.0f} fm")
    check("C. Needed and not supplied: an antiperiodic travelling wave on a ring of three SINGLE branches -- the derived Z_2 (5c) belongs to the unordered pair, and equal-impedance junctions flip no sign",
          True, "either the electron keeps a pair structure (a bundle, which section 6 of the step-1 note found disfavoured), or a new rule flips the sign at pole junctions")
    check("D. If a sign flipped at every head-to-tail junction, odd rings would be antiperiodic and even rings periodic: e (3) and d (5) fermions, DQD (2) boson -- but u (6) and nu (6) bosons too, which is wrong",
          (3 % 2, 5 % 2, 2 % 2, 6 % 2) == (1, 1, 0, 0), "so a junction rule cannot be the whole story")
    check("E. And L_z = 1/2 about one axis is not the fermion: the -1 under any 2 pi rotation and the exchange sign need the quantum superposition the model lacks (wen-link note)",
          True, "the mechanism gives the value 1/2 of one component, inherited from the mother's wavelength")
    return report("Conclusion: 'spin 1 broken in two' is exact as wave arithmetic -- a daughter ring of 3 strings at the mother's "
                  "frequency carries winding 1/2, the antiperiodic mode with L_z = 1/2 -- and it fixes the size: the electron ring "
                  "must have radius R_3/2 (strings of l_1/2) for mother and daughters to sit at m_e c^2, the same R_3/2 the open turn "
                  "gave. What it does not supply: the Z_2 structure that lets three single branches carry an antiperiodic wave "
                  "(the derived one belongs to the pair), a rule that also makes the even rings u and nu fermions, and the exchange "
                  "sign itself. Status: consistent arithmetic, one testable consequence (R_3/2), the fermion still not derived.")


if __name__ == "__main__":
    raise SystemExit(main())
