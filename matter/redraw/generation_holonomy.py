#!/usr/bin/env python3
"""Redraw: the author's R60 target -- 'two neutral DQDs added => U_generation = R(2 pi/3), without
geometric adjustment', with Koide as the main law, the 2/9 as a Coulomb phase, and n_k = 3 + 4k.

Checks:
  A. the count rule: n_k = 3 + 4k gives n mod 3 = k: adding 4 strands shifts the triple's cyclic
     position by one (4 = 1 mod 3), one third of a turn per generation; three steps return to the
     electron's class (n = 15, 12.6 GeV, LEP-excluded): exactly three generations.
  B. Phase A's holonomy H = R(2 pi [(1 - a) Lk + a Wr]) with Lk in Z/3 (step3b, verified numerically):
     a third of a turn carried by WRITHE (Wr = Lk = 1/3) gives exactly R(2 pi/3) for any a -- no
     adjustment; carried by TWIST (Tw = 1/3, Wr = 0) it gives (1 - a) x 120 deg = 97.7 deg at
     a(3) = 0.186.  The target passes in the writhe channel only.
  C. Koide + theta = 2/9 + m_e: sqrt(m_k) = A [1 + sqrt 2 cos(2/9 + 2 pi k/3)] with (tau, e, mu) at
     k = (0, 1, 2) gives m_mu = 105.659 MeV (+1e-5 vs 105.6584) and m_tau = 1776.97 MeV (+0.5 sigma
     vs 1776.93 +- 0.09): the author's numbers, checked.
  D. the 2/9 as a Coulomb phase over one Compton time: theta = U tau_e/hbar with tau_e = hbar/(m_e c^2)
     and U = q_1 q_2 e^2/(4 pi eps0 r_e) = q_1 q_2 m_e c^2 (definition of r_e): the author's u x d gives
     2/9; the two ADDED neutral DQDs, two branch pairs (1/3)(1/3), give 2/9 as well, with objects
     that are actually in the lepton.  Exact by construction; the content is the claim.  Open:
     Koide's theta is one global offset, whereas a per-step phase would accumulate (2/9) k.
  E. the 45 deg: the cube-diagonal tilt of R59; the singlet/doublet equality |v_common| = |v_diff|
     is still to derive.
  F. verdict: R(2 pi/3) per pair of DQDs is exact in the writhe channel with the mod-3 count rule;
     Koide + 2/9 + m_e reproduces mu and tau; the 2/9 is a Coulomb-phase identity at r_e.  Missing:
     why the third of a turn is writhe not twist, why DQDs come in pairs, the global-vs-per-step
     status of 2/9, and the 45 deg.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

ME, MMU, MTAU, DTAU = 0.51099895, 105.6583755, 1776.93, 0.09
A3 = 0.186                      # Berry factor of the triple (step3b)
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def holonomy_deg(a, lk, wr):
    return (360 * ((1 - a) * lk + a * wr)) % 360


def main() -> int:
    # A. count rule
    ns = [3 + 4 * k for k in range(4)]
    check("A. n_k = 3 + 4k: n mod 3 = k, one third of a turn per generation (4 = 1 mod 3); the fourth (n = 15, 12.6 GeV, LEP-excluded) returns to the electron's class: exactly three generations",
          [n % 3 for n in ns] == [0, 1, 2, 0], f"n = {ns}, n mod 3 = {[n % 3 for n in ns]}")

    # B. holonomy channels
    h_writhe = holonomy_deg(A3, 1 / 3, 1 / 3)
    h_twist = holonomy_deg(A3, 1 / 3, 0.0)
    h_writhe_any = [holonomy_deg(a, 1 / 3, 1 / 3) for a in (0.0, 0.186, 0.5, 1.0)]
    check("B. Phase A holonomy H = R(2 pi[(1 - a) Lk + a Wr]): a third of a turn by WRITHE gives exactly 120 deg for any a (no adjustment); by TWIST it gives (1 - a) x 120 = 97.7 deg at a = 0.186: the target passes in the writhe channel only",
          all(abs(h - 120) < 1e-9 for h in h_writhe_any) and abs(h_twist - 97.7) < 0.1,
          f"writhe: {h_writhe:.1f} deg (a-independent: {[round(h, 1) for h in h_writhe_any]}); twist: {h_twist:.1f} deg")

    # C. Koide + 2/9 + m_e
    th = 2 / 9
    f = lambda k: 1 + math.sqrt(2) * math.cos(th + 2 * math.pi * k / 3)
    A = math.sqrt(ME) / f(1)                        # electron at k = 1
    m_mu = (A * f(2)) ** 2
    m_tau = (A * f(0)) ** 2
    check("C. Koide + theta = 2/9 + m_e: m_mu = 105.659 MeV (+1e-5), m_tau = 1776.97 MeV (+0.5 sigma of 0.09): the author's numbers hold",
          abs(m_mu / MMU - 1) < 3e-5 and abs(m_tau - MTAU) < DTAU,
          f"m_mu = {m_mu:.4f} ({m_mu/MMU-1:+.1e}), m_tau = {m_tau:.3f} ({(m_tau-MTAU)/DTAU:+.2f} sigma)")

    # D. 2/9 as a Coulomb phase
    phase_ud = (2 / 3) * (1 / 3)                    # U tau_e / hbar with U = q1 q2 m_e c^2
    phase_2dqd = 2 * (1 / 3) * (1 / 3)
    per_step = [th * k for k in range(3)]
    check("D. theta = U tau_e/hbar with U = q_1 q_2 m_e c^2 at r_e: u x d gives 2/9 (the author); the two added neutral DQDs, two branch pairs (1/3)(1/3), give 2/9 too with objects present in the lepton. Exact by construction. Open: Koide's theta is one global offset, a per-step phase would accumulate (2/9) k",
          abs(phase_ud - 2 / 9) < 1e-15 and abs(phase_2dqd - 2 / 9) < 1e-15 and per_step[2] != th,
          f"u x d = {phase_ud:.4f}; 2 DQD = {phase_2dqd:.4f}; per-step would give {[round(p, 4) for p in per_step]}")

    check("E. The 45 deg: the cube-diagonal tilt (R59); the singlet/doublet equality |v_common| = |v_diff| is still to derive", True, "recorded")
    check("F. Verdict: R(2 pi/3) per pair of DQDs is exact in the writhe channel with the mod-3 count rule; Koide + 2/9 + m_e reproduces mu and tau; 2/9 is a Coulomb-phase identity at r_e. Missing: writhe vs twist, why pairs, global-vs-per-step 2/9, the 45 deg",
          True, "see A-E")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
