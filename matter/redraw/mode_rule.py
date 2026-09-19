#!/usr/bin/env python3
"""Redraw: test the rule 'charged strands carry the mode' (R28), and say what the mode is.

The mode (BASE.md R10-R13): the fluid runs at c along the string and reflects at the two poles;
the standing half-wave between the poles is the mode, E_mode = hbar c pi / l_arc.  For the electron
(3/4 turn, R = 4 lambda-bar/3) the arc is 2 pi lambda-bar = 2426 fm and E_mode = m_e c^2/2 = 255.5 keV:
the circulating half of the mass (g = 2), the other half static at the poles.  Per strand: 809 fm,
the base's l_1 = 2 pi R_3/3.

Checks:
  A. the electron's mode: hbar c pi / (2 pi lambda-bar) = m_e c^2 / 2 exactly; l_1 = 809 fm.
  B. scale law from the ladder: l_1(n) = l_1(3) (3/n)^(2 pi): mu 3.94 fm, tau 0.23 fm; l_1(9) =
     0.813 fm = r_p (-3 %) (recorded coincidence).
  C. leptons and neutrino under the rule: e, mu, tau unchanged (three charged strands, arc 3 l_1(n));
     nu has no charged arc, no mode, m = 0 (< 0.45 eV, KATRIN 2025).
  D. nucleons' mode part under the rule at the nucleon's count n = 9: proton arc 5 l_1 = 4.06 fm
     -> 153 MeV; neutron arc 4 l_1 = 3.25 fm -> 191 MeV.  The base's own mode part m/|g|
     (g_bookkeeping): 168 MeV (p), 246 MeV (n): -9 %, -22 %, and the ordering n > p is right.
     The ~800 MeV static centre is outside the rule.
  E. pions: at n = 6 the count-set scale gives l_1 = 10.4 fm and a 20 MeV mode for pi+ (3 charged
     strands); the closed ring's own mode at r_e/2 (R16) is 140 MeV: closure sets its own scale.
  F. the rule as a MASS law (each charged strand its own mode, m = 2 x sum): at n = 3, u = 0.34 MeV,
     d = 0.17 MeV, wrong order (m_d - m_u = +2.3 MeV needed); at n = 9, p = 848 MeV (-10 %), n = 678
     (-28 %), wrong order.  As a mode law it holds to 10-20 %; as a mass law it fails.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
MP, MN = 938.27209, 939.56542
MPIP = 139.57039
RP = 0.8409
G_P, G_N = 2 * 2.79284734, 2 * 1.91304273
ALPHA = 1 / 137.035999
P = 2 * math.pi
LAMBDA_E = HBARC / ME
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def l1(n):
    return 2 * math.pi * LAMBDA_E / 3 * (3 / n) ** P


def mode(arc):
    return HBARC * math.pi / arc


def main() -> int:
    # A. electron mode
    arc_e = 0.75 * 2 * math.pi * (4 * LAMBDA_E / 3)
    e_mode = mode(arc_e)
    check("A. The mode = standing half-wave of the fluid between the two poles; electron arc (3/4 turn at 4 lambda-bar/3) = 2 pi lambda-bar = 2426 fm, E_mode = hbar c pi / arc = m_e c^2 / 2 exactly (the circulating half, g = 2); per strand l_1 = 809 fm",
          abs(arc_e - 2 * math.pi * LAMBDA_E) < 1e-9 and abs(e_mode / (ME / 2) - 1) < 1e-12 and abs(l1(3) - 808.8) < 0.5,
          f"arc = {arc_e:.0f} fm, E_mode = {e_mode*1e3:.1f} keV = {e_mode/ME:.3f} m_e c^2, l_1 = {l1(3):.1f} fm")

    # B. scale law
    check("B. Scale law from the ladder l_1(n) = l_1(3)(3/n)^2pi: mu 3.94 fm, tau 0.23 fm; l_1(9) = 0.813 fm = r_p (-3 %), recorded",
          abs(l1(7) - 3.94) < 0.01 and abs(l1(11) - 0.230) < 0.002 and abs(l1(9) / RP - 0.966) < 0.003,
          f"l_1(7) = {l1(7):.2f} fm, l_1(11) = {l1(11):.3f} fm, l_1(9) = {l1(9):.3f} fm = {l1(9)/RP:.3f} r_p")

    # C. leptons and neutrino
    m_lep = {n: 2 * mode(3 * l1(n)) for n in (3, 7, 11)}
    check("C. Under the rule e, mu, tau are unchanged (three charged strands, arc 3 l_1(n), m = 2 E_mode): mu -0.8 %, tau +1.0 %; the neutrino has no charged arc, hence no mode: m = 0 (< 0.45 eV, KATRIN 2025)",
          abs(m_lep[3] / ME - 1) < 1e-12 and abs(m_lep[7] / MMU - 1) < 0.011 and abs(m_lep[11] / MTAU - 1) < 0.011,
          f"e {m_lep[3]:.4f}, mu {m_lep[7]:.1f} ({m_lep[7]/MMU-1:+.1%}), tau {m_lep[11]:.0f} ({m_lep[11]/MTAU-1:+.1%}) MeV; nu: no arc")

    # D. nucleons' mode part
    mode_p, mode_n = mode(5 * l1(9)), mode(4 * l1(9))
    base_p, base_n = MP / G_P, MN / G_N
    check("D. Nucleons at n = 9: proton mode on 5 charged strands = 153 MeV, neutron on 4 = 191 MeV; the base's own mode part m/|g| is 168 (p) and 246 (n) MeV: -9 %, -22 %, ordering n > p right; the ~800 MeV static centre is outside the rule",
          abs(mode_p - 152.6) < 1 and abs(mode_n - 190.7) < 1 and abs(mode_p / base_p - 0.909) < 0.005 and abs(mode_n / base_n - 0.777) < 0.005 and mode_n > mode_p and base_n > base_p,
          f"rule: p {mode_p:.1f}, n {mode_n:.1f} MeV; base m/|g|: p {base_p:.1f}, n {base_n:.1f} MeV; static: p {MP-base_p:.0f}, n {MN-base_n:.0f} MeV")

    # E. pions
    mode_pi = mode(3 * l1(6))
    ring_pi = 2 * ME / ALPHA
    check("E. Pions at n = 6: count-set scale l_1 = 10.4 fm, pi+ mode on 3 charged strands = 20 MeV vs the closed ring's 140 MeV (R16): closure sets its own scale, the count law does not apply to rings",
          abs(l1(6) - 10.38) < 0.05 and abs(mode_pi - 19.9) < 0.3 and abs(ring_pi - 140.05) < 0.1,
          f"l_1(6) = {l1(6):.2f} fm, pi+ mode = {mode_pi:.1f} MeV, ring = {ring_pi:.1f} MeV, measured {MPIP:.1f}")

    # F. as a mass law
    # each charged strand its own mode at scale n, normalised so that three strands give m_e at n = 3
    per_strand = lambda n: ME / 3 * (n / 3) ** P
    mu_q, md_q = 2 * per_strand(3), per_strand(3)
    p_mass, n_mass = 5 * per_strand(9), 4 * per_strand(9)
    check("F. As a MASS law (each charged strand its own mode): at n = 3, u = 0.34, d = 0.17 MeV, wrong order (m_d - m_u must be +2.3 MeV); at n = 9, p = 848 MeV (-10 %), n = 678 (-28 %), wrong order: a mode law to 10-20 %, not a mass law",
          abs(mu_q - 0.341) < 0.002 and abs(md_q - 0.170) < 0.002 and abs(p_mass / MP - 0.904) < 0.005 and abs(n_mass / MN - 0.722) < 0.005,
          f"u {mu_q:.3f}, d {md_q:.3f} MeV; p {p_mass:.0f} ({p_mass/MP-1:+.0%}), n {n_mass:.0f} ({n_mass/MN-1:+.0%}) MeV")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
