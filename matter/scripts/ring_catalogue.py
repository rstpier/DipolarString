#!/usr/bin/env python3
"""Rings of n strings carrying +-e/3 each, joined end to end by their poles: what charges, which
arrangements, and what the mass readings give.  Author's proposal: quarks as n = 5 (charge 1/3) and
n = 6 (charge 2/3).

  A. Charge rule: q = (n_plus - n_minus) e/3 with n = n_plus + n_minus, so the PARITY of n fixes
     the charge class: odd n -> odd thirds (+-1/3, +-1, +-5/3 ...), even n -> even thirds
     (0, +-2/3, +-4/3, ...).  Leptons and d-type quarks are odd rings, u-type quarks and neutral
     states even rings.  The proposal (d: 5, u: 6) obeys it; so would (d: 3, u: 2) or (d: 7, u: 4).
  B. Stability: mixed rings are bound at every pole junction (arcs on the circle R_n = n l_1 / 2 pi,
     poles delta = q, tubes in contact), for every arrangement of the signs.
  C. Arrangements: distinct sign orderings around the ring up to rotation and reflection --
     n = 6 (4+, 2-) has 3, n = 5 (3-, 2+) has 2, n = 3 (2-, 1+) has 1.  Three for u, not for d:
     the count does not reproduce colour for both.
  D. Masses: the inductance scaling (n/3)^2 Z_e gives u (n = 6) = 4.0 m_e = 2.0 MeV (PDG 2.2) and
     d (n = 5) = 2.8 m_e = 1.4 MeV (PDG 4.7); the mode-energy reading makes larger rings lighter
     than the electron.  Current-quark masses are scheme-dependent and a weak target either way.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import itertools
import math
import numpy as np

from berry_holonomy_common import check, report
from string_assembly import K, L1, M, R3, RW, mutual

ME_MEV = 0.51100
PDG = {"u": 2.16, "d": 4.67}      # MeV, MS-bar at 2 GeV


def charges_for(n: int):
    return sorted({(2 * k - n) for k in range(n + 1)})          # in units of e/3: n_plus - n_minus


def distinct_arrangements(n: int, n_plus: int):
    seen = set()
    for pos in itertools.combinations(range(n), n_plus):
        pattern = tuple(1 if i in pos else -1 for i in range(n))
        canon = min(min(pattern[k:] + pattern[:k], tuple(reversed(pattern[k:] + pattern[:k]))) for k in range(n))
        seen.add(canon)
    return sorted(seen)


def ring_energy(pattern, delta: float = 1.0):
    n = len(pattern)
    radius = n * L1 / (2 * math.pi)
    span = 2 * math.pi / n
    gap_ang = 2 * RW / radius
    strings = []
    for k, sgn in enumerate(pattern):
        t0 = k * span
        m = M
        t = t0 + gap_ang / 2 + (np.arange(m) + 0.5) * (span - gap_ang) / m
        pts = np.stack([radius * np.cos(t), radius * np.sin(t), np.zeros(m)], axis=-1)
        qs = np.full(m, sgn / m)
        tp = np.array([t0 + gap_ang / 2 + RW / radius, t0 + span - gap_ang / 2 - RW / radius])
        poles = np.stack([radius * np.cos(tp), radius * np.sin(tp), np.zeros(2)], axis=-1)
        strings.append((np.concatenate([pts, poles]), np.concatenate([qs, [-delta, +delta]])))
    total = 0.0
    junctions = []
    for i in range(n):
        for j in range(i + 1, n):
            u = mutual(*strings[i], *strings[j])
            total += u
            if j == i + 1 or (i == 0 and j == n - 1):
                junctions.append(u)
    return total, junctions


def main() -> int:
    # A. charge rule and parity
    table = {n: charges_for(n) for n in range(2, 9)}
    print("  charges available to a ring of n strings (units e/3):")
    for n, cs in table.items():
        print(f"    n = {n}: " + ", ".join(f"{c:+d}" for c in cs))
    check("A. Charge of a ring = (n_plus - n_minus) e/3: the parity of n fixes the charge class (odd n -> odd thirds, even n -> even thirds)",
          all(all(c % 2 == n % 2 for c in cs) for n, cs in table.items()), "a 1/3 charge needs an odd ring, a 2/3 charge an even ring, a neutral state an even ring")
    check("A'. The proposal d = n 5 (charge -1/3: 2+, 3-) and u = n 6 (charge +2/3: 4+, 2-) obeys the rule; so do (d 3, u 2) and (d 7, u 4)",
          -1 in table[5] and 2 in table[6] and -1 in table[3] and 2 in table[2], "the rule fixes parity, not n: the choice 5/6 needs another criterion")

    # B. stability of mixed rings
    results = {}
    for name, pattern_list in (("e (3-)", [(-1, -1, -1)]), ("d (5: 2+ 3-)", distinct_arrangements(5, 2)),
                               ("u (6: 4+ 2-)", distinct_arrangements(6, 4)), ("nu? (6: 3+ 3-)", distinct_arrangements(6, 3))):
        results[name] = [(p, *ring_energy(p)) for p in pattern_list]
    all_bound = all(all(j < 0 for j in junc) for v in results.values() for _, _, junc in v)
    print("  ring energies at delta = q (keV): total interaction, and the range of junction energies")
    for name, v in results.items():
        for p, tot, junc in v:
            print(f"    {name:16s} {''.join('+' if s > 0 else '-' for s in p):8s}  total {tot*1e3:+7.1f}   junctions {min(junc)*1e3:+.2f} .. {max(junc)*1e3:+.2f}")
    check("B. Every mixed ring is bound at every pole junction, whatever the order of the signs (delta = q)", all_bound,
          "end-to-end pole contact dominates; the sign order changes the total by the like/unlike neighbour terms")

    # C. arrangements
    n_arr = {k: len(distinct_arrangements(*v)) for k, v in (("u 6 (4+2-)", (6, 4)), ("d 5 (2+3-)", (5, 2)), ("d 3 (1+2-)", (3, 1)), ("nu 6 (3+3-)", (6, 3)))}
    check("C. Distinct sign orderings (up to rotation and reflection): u 6 has 3, d 5 has 2, d 3 has 1 -- three states for u only, so 'colour = arrangement' fails for d",
          n_arr["u 6 (4+2-)"] == 3 and n_arr["d 5 (2+3-)"] == 2 and n_arr["d 3 (1+2-)"] == 1, ", ".join(f"{k}: {v}" for k, v in n_arr.items()))
    lowest = {name: min(v, key=lambda t: t[1]) for name, v in results.items()}
    alternating = all(sum(1 for a, b in zip(p, p[1:] + p[:1]) if a != b) == max(sum(1 for a, b in zip(q, q[1:] + q[:1]) if a != b) for q, _, _ in v)
                      for name, v in results.items() for p, _, _ in [lowest[name]])
    check("C'. The lowest-energy arrangement ALTERNATES the signs as much as possible (the most unlike neighbours) -- computed, not assumed",
          alternating, "; ".join(f"{name}: {''.join('+' if s > 0 else '-' for s in p)} ({tot*1e3:+.1f} keV)" for name, (p, tot, _) in lowest.items()))

    # D. masses
    m_ind = {"u": (6 / 3) ** 2 * ME_MEV, "d": (5 / 3) ** 2 * ME_MEV}
    m_mode = {"u": 3 / 6 * ME_MEV, "d": 3 / 5 * ME_MEV}
    check("D. Inductance scaling (n/3)^2 Z_e: u (6) = 2.0 MeV vs PDG 2.2; d (5) = 1.4 MeV vs PDG 4.7 -- one hit, one miss by 3",
          abs(m_ind["u"] / PDG["u"] - 0.95) < 0.05 and abs(m_ind["d"] / PDG["d"] - 0.30) < 0.05,
          f"u: {m_ind['u']:.2f} MeV ({m_ind['u']/PDG['u']:.2f} x PDG), d: {m_ind['d']:.2f} MeV ({m_ind['d']/PDG['d']:.2f} x PDG); the scaling itself is the calibrated mass law, and the 5/6 ordering makes d lighter than u -- the wrong order")
    check("D'. Mode-energy reading (m ~ 1/n): u = 0.26 MeV, d = 0.31 MeV -- larger rings are lighter than the electron, excluded",
          m_mode["u"] < ME_MEV and m_mode["d"] < ME_MEV, f"u: {m_mode['u']:.2f} MeV, d: {m_mode['d']:.2f} MeV")
    check("D''. Current-quark masses are scheme-dependent (MS-bar, 2 GeV) and dominated by binding in hadrons: a weak target for any string count",
          True, "the parity rule (A) is the robust content of the proposal; n itself is not fixed by it")
    return report("Conclusion: a ring of n strings at e/3 each has charge (n+ - n-)/3, so the parity of n fixes the charge "
                  "class -- odd rings for 1/3 and 1 (d, e), even rings for 2/3 and 0 (u, neutral). The proposal d = 5, u = 6 "
                  "obeys it, as would 3/2 or 7/4: parity is derived, the value of n is not. Mixed rings are bound at every "
                  "junction; the arrangement count gives three states for u (6) but two for d (5) and one for d (3), so it is "
                  "not colour. The inductance scaling gives u right and d off by 3 with d lighter than u; the mode reading "
                  "makes both lighter than the electron.")


if __name__ == "__main__":
    raise SystemExit(main())
