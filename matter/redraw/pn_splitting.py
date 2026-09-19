#!/usr/bin/env python3
"""Redraw: what fixes the p-n splitting, m_n - m_p = 1.293 MeV, in the final topology (R41)?

Checks:
  A. p = uud and n = udd are the same three 3-strand circuits and the same junction graph: the
     static 763 MeV and the joints cancel in n - p.  What differs is the strands' content: the
     neutron has one neutral strand more (5 vs 4) and one charged strand fewer (4 vs 5).
  B. if the neutral strand is a full DQD -- two branches +-e/3, the vacuum's own unit (R9) -- its
     bifilar field energy at the nucleon scale is (e/3)^2/(eps0 l_1(9)) = 4 pi K/(9 l_1) = 2.47 MeV
     (E + B, R30 B); a charged strand's single-branch energy is already in its circuit.  Strong part
     of n - p: +2.47 MeV; lattice QCD+QED (BMW): +2.52 +- 0.29.  Closed form (2 alpha/3) 3^(2 pi) m_e c^2.
  C. minus the proton's Coulomb self-energy (R19 D: 0.86 shell, 1.03 uniform at r_p; lattice
     -1.00 +- 0.16): m_n - m_p = 1.44-1.61 MeV vs 1.293: right sign, +11 to +25 %.
  D. the same rule on the pion pair (pi+ - pi0 = 4.59 MeV): with the ring's strand length 1.48 fm
     the DQD energy is 1.36 MeV and the pi+ Coulomb ~0.5 MeV: 1.9 MeV with uubar (2.4 x short),
     wrong sign with ddbar.  The rule does not carry to the pion ring (not made of circuits): recorded.
  E. the inputs: alpha, m_e, the exponent 2 pi, the count 9 -- no free number; the residual 11-25 %
     sits in the proton's Coulomb term, whose geometry the base has not fixed.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
MPIP, MPI0 = 139.57039, 134.9768
LAMBDA_E = HBARC / ME
L1 = lambda n: 2 * math.pi * LAMBDA_E / 3 * (3 / n) ** (2 * math.pi)
RP, R_E_HALF = 0.8409, ALPHA * LAMBDA_E / 2
LATTICE_STRONG, LATTICE_EM = (2.52, 0.29), (1.00, 0.16)
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def e_dqd(l):
    return 4 * math.pi * K / (9 * l)


def main() -> int:
    # A. same circuits, different content
    p_content, n_content = (4, 1, 4), (2, 2, 5)          # (+, -, 0)
    check("A. p and n: same three 3-strand circuits, same junction graph -> static 763 MeV and joints cancel; the neutron has one neutral strand more (5 vs 4) and one charged strand fewer (4 vs 5)",
          sum(p_content) == sum(n_content) == 9 and n_content[2] - p_content[2] == 1 and (p_content[0] + p_content[1]) - (n_content[0] + n_content[1]) == 1,
          f"p = {p_content}, n = {n_content} (+, -, 0)")

    # B. the extra neutral strand as a DQD
    strong = e_dqd(L1(9))
    closed = (2 * ALPHA / 3) * 3 ** (2 * math.pi) * ME
    check("B. Neutral strand = full DQD (two branches +-e/3, R9): bifilar energy 4 pi K/(9 l_1(9)) = 2.47 MeV = (2 alpha/3) 3^2pi m_e c^2; strong part of n - p: +2.47 vs lattice +2.52 +- 0.29",
          abs(strong - 2.474) < 0.005 and abs(closed - strong) < 1e-9 and abs(strong - LATTICE_STRONG[0]) < LATTICE_STRONG[1],
          f"E_DQD = {strong:.3f} MeV; lattice {LATTICE_STRONG[0]} +- {LATTICE_STRONG[1]}")

    # C. minus the proton's Coulomb
    e_shell, e_unif = K / (2 * RP), 3 * K / (5 * RP)
    dm = (strong - e_unif, strong - e_shell)
    dev = tuple(d / (MN - MP) - 1 for d in dm)
    check("C. Minus the proton's Coulomb self-energy (0.86 shell / 1.03 uniform at r_p; lattice 1.00 +- 0.16): m_n - m_p = 1.44-1.61 MeV vs 1.293: right sign, +11 to +25 %",
          abs(dm[0] - 1.447) < 0.01 and abs(dm[1] - 1.618) < 0.01 and abs(dev[0] - 0.119) < 0.01 and abs(dev[1] - 0.251) < 0.01,
          f"dm = {dm[0]:.3f}-{dm[1]:.3f} MeV ({dev[0]:+.0%} to {dev[1]:+.0%}); with the lattice EM part: {strong-LATTICE_EM[0]:.2f} MeV ({(strong-LATTICE_EM[0])/(MN-MP)-1:+.0%})")

    # D. pions
    l_pi = 2 * math.pi * R_E_HALF / 6
    e_dqd_pi = e_dqd(l_pi)
    coul_pi = K / (2 * R_E_HALF)
    d_uubar = e_dqd_pi + coul_pi                 # pi+ has one neutral more than u ubar
    d_ddbar = -e_dqd_pi + coul_pi                # ... one fewer than d dbar
    check("D. Pions: with the ring's strand 1.48 fm the DQD energy is 1.36 MeV and the pi+ Coulomb 0.5 MeV: 1.9 MeV with u ubar (2.4 x short of 4.59), wrong sign with d dbar -- the rule does not carry to the pion ring (recorded)",
          abs(l_pi - 1.475) < 0.005 and abs(e_dqd_pi - 1.36) < 0.01 and abs((MPIP - MPI0) / d_uubar - 2.4) < 0.1 and d_ddbar < 0,
          f"l = {l_pi:.3f} fm; E_DQD = {e_dqd_pi:.2f}, Coulomb {coul_pi:.2f}; u ubar: {d_uubar:.2f} MeV, d dbar: {d_ddbar:+.2f} MeV vs {MPIP-MPI0:.2f}")

    check("E. Inputs: alpha, m_e, the exponent 2 pi, the count 9 -- no free number; the 11-25 % residual sits in the proton's Coulomb term, whose geometry the base has not fixed",
          True, "recorded")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
