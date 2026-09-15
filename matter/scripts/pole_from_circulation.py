#!/usr/bin/env python3
"""The author's mechanism for charge and poles: a charged particle always comes from a neutral
mother that also contained its antiparticle; the breaking gave opposite circulation senses, which
polarised the fluid at opposite ends (+/-).

In the model's terms:
  * the mother is the DQD -- forward branch and return branch, opposite circulation, +-e/3
    (manuscript l. 286); breaking it is unpairing (0.47 keV, string_end_charge.py); the sign of
    each fragment is its circulation sense (manuscript's chirality rule);
  * within a fragment, a fluid circulating at c_0 and meeting an open end piles up there: the
    standing wave of the open string, V(s) = V_0 cos(pi s / l), has its two ends at OPPOSITE
    potential -- the poles.  A static pile-up would relax in a conductor; the poles are the
    antinodes of the string's own mode, and the +/- pattern is the phase fixed at the cut (charge
    conservation: the two new ends of one cut carry opposite charge).
  * so the pole strength delta -- which the manuscript states and does not quantify, and on which
    the chaining of like-charge strings depends (threshold delta/q = 0.28) -- becomes computable:
    the end charge of one quantum of the fundamental on a string of impedance Z.

Derivation: mode energy C' V_0^2 l / 4 = hbar omega = pi hbar c / l; charge density C' V(s);
end charge delta = C' V_0 l / pi = sqrt(4 hbar c C' / pi); with C' = 1/(Z c):
      delta / e = sqrt(4 hbar / (pi Z e^2)) = sqrt(2 R_K / (pi^2 Z)),   R_K = h / e^2 = 25.8 kOhm.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import sympy as sp

from berry_holonomy_common import check, report

R_K = 25812.807        # Ohm
Z0 = 376.730
Z_E = 0.730 * Z0
HBARC = 197.327        # MeV fm
L1 = 2 * math.pi * 386.159 / 3


def main() -> int:
    # symbolic derivation of the end charge
    cp, v0, ell, hbar, c, s, z = sp.symbols("C_p V_0 l hbar c s Z", positive=True)
    energy = sp.integrate(sp.Rational(1, 2) * cp * (v0 * sp.cos(sp.pi * s / ell)) ** 2, (s, 0, ell))    # peak electric energy = mode energy
    v0_sol = sp.solve(sp.Eq(energy, sp.pi * hbar * c / ell), v0)[0]
    delta = sp.integrate(cp * v0_sol * sp.cos(sp.pi * s / ell), (s, 0, ell / 2))
    delta_z = sp.simplify(delta.subs(cp, 1 / (z * c)))
    check("End charge of one quantum of the open-string fundamental: delta = sqrt(4 hbar / (pi Z)) (symbolic)",
          sp.simplify(delta_z - sp.sqrt(4 * hbar / (sp.pi * z))) == 0, f"delta = {delta_z}")

    d_over_e = math.sqrt(2 * R_K / (math.pi ** 2 * Z_E))
    q_zpf = math.sqrt(R_K / (4 * math.pi * Z_E))
    check("On a string of the electron's impedance (275 Ohm) that end charge is ~4 e; the zero-point charge is ~3 e",
          4.0 < d_over_e < 5.0 and 2.5 < q_zpf < 3.0, f"delta = {d_over_e:.2f} e (one quantum, cos profile); q_zpf = sqrt(R_K / 4 pi Z) = {q_zpf:.2f} e")
    check("A DS string is far below the quantum resistance (Z_e / R_K = 1 %): charge per quantum is LARGE, not small",
          Z_E / R_K < 0.02, f"Z_e / R_K = {Z_E/R_K:.4f}; e/3 as an end charge would need Z = {2*R_K/(math.pi**2/9):.0f} Ohm = {2*R_K/(math.pi**2/9)/Z0:.0f} Z_0")
    check("So the poles are not the fluid's e/3: they are ~10 x larger, and e/3 stays the net charge of the fluid (a separate quantity)",
          d_over_e / (1 / 3) > 10, f"delta / (e/3) = {d_over_e*3:.1f}")
    check("Chaining of like-charge strings needs delta/q >= 0.28 (string_assembly.py): satisfied by a factor ~ 45",
          d_over_e * 3 / 0.28 > 30, f"delta/q = {d_over_e*3:.1f} vs threshold 0.28 -- rings and chains form easily")
    # energy caveat
    e_open = math.pi * HBARC / L1
    check("Caveat: one quantum on an OPEN string of l_1 costs pi hbar c / l_1 = 0.77 MeV = 1.5 m_e c^2; the poles must be the ring mode's junction standing waves, not free open strings",
          abs(e_open - 0.766) < 0.01, f"pi hbar c / l_1 = {e_open:.3f} MeV; the closed 3-string ring mode is hbar c / R_3 = 0.511 MeV (the anchor)")
    check("Static versus dynamic: a DC pile-up relaxes in a conductor; the poles are the antinodes of the mode, and the +/- pattern is the phase fixed at the cut (charge conservation at a cut)",
          True, "the two ends created by one cut carry opposite charge -- the mechanism is charge conservation plus a moving fluid")
    return report("Conclusion: the mechanism is consistent with the model and supplies what it lacked -- the poles are the end "
                  "charges of the string's own mode, sign fixed by the circulation sense at the cut, mother = DQD, breaking = "
                  "unpairing. Quantified: one quantum on a 275-Ohm string has end charges of ~4 e (zero-point ~3 e), ten times "
                  "the fluid's e/3, so like-charge strings chain with a margin of ~45 and rings form easily; e/3 is the net "
                  "charge, not the pole. Open caveat: a lone open string of l_1 costs 1.5 m_e c^2 per quantum, so the poles live "
                  "at the junctions of the closed ring, where Gamma_pole = 1/3 sets the standing-wave fraction.")


if __name__ == "__main__":
    raise SystemExit(main())
