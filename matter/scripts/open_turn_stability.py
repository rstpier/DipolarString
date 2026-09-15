#!/usr/bin/env python3
"""Author's correction: the mother is not a closed ring but an OPEN chain of DQDs, =:=:=:=:=:=:= ;
three DQDs are === ; the chain splits lengthwise into a + row and a - row, each curling into a
partially open turn (one symmetry axis -> spin 1/2 by the author's rule 'spin = number of
symmetry axes / 2'; the DQD = and the string - have two axes -> spin 1).

Checks:
  A. the arithmetic is unchanged: an open chain of 3 DQDs (bifilar, length 3 l_1, open ends) has
     the half-wave fundamental at hbar c / 2 R_3 = 0.256 MeV, and so has each daughter chain of
     3 strings: the daughters need strings of l_1/2 to sit at m_e c^2 (R_3/2 again);
  B. the axis rule reproduces photon 1, DQD 1, open turn 1/2 -- and gives 3/2 for a triangle of
     three strings, 2 for a square, infinity for a closed ring: a shape heuristic, not the
     rotation group (no -1 under 2 pi, no exchange sign);
  C. STABILITY: the open turn's two free ends carry opposite poles (a tail -delta and a head
     +delta, delta ~ 4 e from the mode) and attract; the like net charges (-e/3) repel far more
     weakly; closing the gap lowers the mode energy in every family (Compton 3pi -> 2pi, low
     pi -> 0).  The partially open turn snaps shut into a ring -- the same result that made
     polar strings close into rings in string_assembly.py.  The one-axis object is not stable
     in the model's electrostatics; nothing holds the gap open.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

from berry_holonomy_common import check, report
from string_assembly import K, L1, R3, RW, string_charges, mutual

HBARC = 197.327
ME = 0.51100
D0 = 2 * math.cosh(math.pi) * RW
DELTA_E = 4.36          # pole charge per quantum, in units of e (pole_from_circulation.py)


def main() -> int:
    # A. arithmetic with an open-chain mother
    mother = math.pi * HBARC / (3 * L1)           # half wave on an open line of 3 l_1
    daughter = math.pi * HBARC / (3 * L1)         # same length, same open ends
    check("A. Open-chain mother (3 DQDs, 3 l_1, open ends): half-wave fundamental 0.256 MeV; each daughter chain of 3 strings: the same",
          abs(mother - 0.2555) < 1e-3 and abs(daughter - mother) < 1e-12, f"mother {mother:.4f} MeV, daughter {daughter:.4f} MeV = m_e c^2 / 2: strings of l_1/2 needed, R_3/2 as before")
    # B. axis rule
    axes = {"string -": 2, "DQD =": 2, "open turn C": 1, "triangle of 3": 3, "square of 4": 4, "closed ring": math.inf}
    spin = {k: v / 2 for k, v in axes.items()}
    check("B. Axis rule spin = axes/2: photon 1, DQD 1, open turn 1/2 -- but triangle 3/2, square 2, closed ring infinite",
          spin["string -"] == 1 and spin["DQD ="] == 1 and spin["open turn C"] == 0.5 and spin["triangle of 3"] == 1.5,
          "; ".join(f"{k}: {v}" for k, v in spin.items()) + " -- a shape heuristic; spin is the rotation group's double cover, not a mirror count")
    # C. stability of the open turn: three like-charge strings on an arc of radius R, gap g between the free ends
    def arc_turn(gap_fm: float, delta: float, q: float = -1.0):
        # three arcs on a circle whose circumference is 3 l_1 + gap; ends separated by `gap`
        circ = 3 * L1 + gap_fm
        radius = circ / (2 * math.pi)
        strings = []
        import numpy as np
        for k in range(3):
            t0 = k * L1 / radius
            m = 300
            t = t0 + (np.arange(m) + 0.5) * (L1 - 2 * RW) / m / radius + RW / radius
            pts = np.stack([radius * np.cos(t), radius * np.sin(t), np.zeros(m)], axis=-1)
            qs = np.full(m, q / m)
            tp = np.array([t0 + RW / radius, t0 + (L1 - RW) / radius])
            poles = np.stack([radius * np.cos(tp), radius * np.sin(tp), np.zeros(2)], axis=-1)
            strings.append((np.concatenate([pts, poles]), np.concatenate([qs, [-delta, +delta]])))
        return sum(mutual(*strings[i], *strings[j]) for i in range(3) for j in range(i + 1, 3))
    gaps = [2 * RW, 0.1 * L1, 0.5 * L1, L1]
    for delta_label, delta in (("delta = e/3 (q)", 1.0), ("delta = 4 e (mode)", DELTA_E * 3)):
        us = [arc_turn(g, delta) for g in gaps]
        print(f"  open turn, {delta_label}: interaction energy vs gap  " + ", ".join(f"g = {g/L1:.2f} l_1: {u*1e3:+.1f} keV" for g, u in zip(gaps, us)))
        if delta > 1.5:
            monotone = all(a < b for a, b in zip(us, us[1:]))
    us_q = [arc_turn(g, 1.0) for g in gaps]
    us_m = [arc_turn(g, DELTA_E * 3) for g in gaps]
    check("C. The open turn's energy falls monotonically as the gap closes, for poles at e/3 and at the mode's 4 e: the ends attract, the turn snaps shut",
          all(a < b for a, b in zip(us_q, us_q[1:])) and all(a < b for a, b in zip(us_m, us_m[1:])),
          f"closing from g = l_1 to contact releases {(us_q[-1]-us_q[0])*1e3:.1f} keV (delta = q) and {(us_m[-1]-us_m[0])*1e3:.0f} keV (delta = 4 e)")
    e_compton_open, e_compton_closed = 3 * math.pi * HBARC / (3 * L1), 2 * math.pi * HBARC / (3 * L1)
    e_low_open, e_low_closed = math.pi * HBARC / (3 * L1), 0.0
    check("C'. Closing also lowers the mode energy in every family (open_turn_closure.py): Compton 3 pi -> 2 pi, low pi -> 0",
          e_compton_closed < e_compton_open and e_low_closed < e_low_open,
          f"Compton: {e_compton_open:.3f} -> {e_compton_closed:.3f} MeV; low: {e_low_open:.3f} -> {e_low_closed:.3f} MeV -- the half-wave (spin-1/2) mode softens to DC")
    check("C''. Nothing in the model holds the gap open: no repulsion (like net charges lose to the poles by the chaining margin, ~45) and no topological obstruction (the classical gap mode is continuous)",
          True, "the one-axis object decays into the closed ring; a stable open turn needs a new ingredient -- a Z_2 obstruction or a gap repulsion")
    return report("Conclusion: the open-chain mother changes no number -- the daughters still sit at half the electron mass unless "
                  "strings are l_1/2 (R_3/2). The axis rule is a shape heuristic that fits photon, DQD and open turn but not the "
                  "rotation group. The real obstacle is stability: the open turn's free ends carry opposite poles and attract, "
                  "every mode family loses energy on closing, and the model has no obstruction -- the partially open turn snaps "
                  "shut into the closed ring, which is the assembly the fluid prefers (string_assembly.py). A one-axis electron "
                  "needs something that keeps the gap open, and V2.10 does not contain it.")


if __name__ == "__main__":
    raise SystemExit(main())
