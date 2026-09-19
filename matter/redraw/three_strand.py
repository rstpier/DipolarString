#!/usr/bin/env python3
"""Redraw: test the neutral string -- rebuild the quarks with three strands, each strand carrying
-1/3, 0 or +1/3 (Bilson-Thompson's helons), and see what the base gains and loses.

Checks:
  A. charge spectrum: three strands in {-1/3, 0, +1/3} give exactly {0, +-1/3, +-2/3, +-1}, the
     observed charges of the elementary fermions, and nothing else (the chain rule allowed +-4/3 at
     n = 4 and +-5/3 at n = 5, never observed).
  B. assignment: e- = (-,-,-), nu = (0,0,0), u = (+,+,0), d = (-,0,0), antiparticles by sign flip;
     colour = position of the odd strand: u and d have 3 arrangements, e and nu have 1.
  C. hadrons: p = uud = (4+, 1-, 4x0), n = udd = (2+, 2-, 5x0), both 9 strands (same count, as p ~ n
     requires); pi+ = (3+, 3x0), pi- = (3-, 3x0), pi0 = (2+, 2-, 2x0), all 6 strands; Delta++ =
     (6+, 3x0) = +2; a baryon is colourless when its three odd strands sit at three different
     positions: 6 of 27 arrangements.
  D. ladder: mu = (-,-,-) + 4 neutral = 7 strands, tau = + 8 = 11: the generation step is four
     neutral strands, two neutral DQDs (R3's mother).  The neutrino (0,0,0) at n = 3 would weigh m_e
     under m = m_e (n/3)^(2 pi); KATRIN (2025): < 0.45 eV.  Only the rule 'charged strands carry the mode,
     the total count sets the scale' keeps nu = 0 with e, mu, tau unchanged (POSTULATED).
  E. p - n: same strand count, so R24's 'one more string' is void; m_n - m_p = (m_d - m_u) - 1.0 (EM)
     gives m_d - m_u = 2.3 MeV (lattice 2.5): the d, with ONE charged strand, must outweigh the u
     with two: mass is not a strand count.
  F. neutrino size under the base's spin condition S = R E/c = hbar/2: R >= 219 nm for m < 0.45 eV.
  G. what the neutral string voids: the parity rule (u even, d odd, p odd, n even), the identity
     'p carries the pi+ ring's strings' (R23), and the R24 count argument.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction as F

HBARC = 197.3269804
ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
MP, MN = 938.27209, 939.56542
NU_LIMIT_EV = 0.45      # KATRIN 2025 (90 % CL)
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def content(*quarks):
    c = Counter()
    for q in quarks:
        for s in q:
            c[s] += 1
    return c[F(1, 3)], c[F(-1, 3)], c[F(0)]


def fmt(t):
    return f"({t[0]}+, {t[1]}-, {t[2]}x0)"


def main() -> int:
    plus, minus, zero = F(1, 3), F(-1, 3), F(0)
    # A. spectrum
    spec3 = {sum(c) for c in itertools.product((minus, zero, plus), repeat=3)}
    chain4 = {sum(c) for c in itertools.product((minus, plus), repeat=4)}
    chain5 = {sum(c) for c in itertools.product((minus, plus), repeat=5)}
    check("A. Three strands in {-1/3, 0, +1/3} give exactly {0, +-1/3, +-2/3, +-1}: the observed fermion charges and nothing else; the chain rule allowed +-4/3 (n = 4) and +-5/3 (n = 5), never observed",
          spec3 == {F(k, 3) for k in range(-3, 4)} and F(4, 3) in chain4 and F(5, 3) in chain5,
          "3-strand: " + ", ".join(str(q) for q in sorted(spec3)) + f"; chains: 4/3 in n=4 {F(4,3) in chain4}, 5/3 in n=5 {F(5,3) in chain5}")

    # B. assignment and colour
    e, nu, u, d = (minus,) * 3, (zero,) * 3, (plus, plus, zero), (minus, zero, zero)
    arrangements = {name: len(set(itertools.permutations(q))) for name, q in (("e", e), ("nu", nu), ("u", u), ("d", d))}
    charges = {name: sum(q) for name, q in (("e", e), ("nu", nu), ("u", u), ("d", d))}
    check("B. e- = (-,-,-), nu = (0,0,0), u = (+,+,0), d = (-,0,0): charges -1, 0, +2/3, -1/3; colour = position of the odd strand: u and d have 3 arrangements, e and nu have 1",
          charges == {"e": F(-1), "nu": F(0), "u": F(2, 3), "d": F(-1, 3)} and arrangements == {"e": 1, "nu": 1, "u": 3, "d": 3},
          "; ".join(f"{k}: q = {charges[k]}, {arrangements[k]} colour(s)" for k in charges))

    # C. hadrons
    ubar, dbar = tuple(-s for s in u), tuple(-s for s in d)
    p, n = content(u, u, d), content(u, d, d)
    pip, pim, pi0 = content(u, dbar), content(d, ubar), content(u, ubar)
    dpp = content(u, u, u)
    colourless = sum(1 for pos in itertools.product(range(3), repeat=3) if len(set(pos)) == 3)
    check("C. p = (4+, 1-, 4x0), n = (2+, 2-, 5x0), both 9 strands; pi+ = (3+, 3x0), pi- = (3-, 3x0), pi0 = (2+, 2-, 2x0), 6 strands; Delta++ = (6+, 3x0) = +2; a baryon is colourless in 6 of 27 arrangements (odd strands at three different positions)",
          p == (4, 1, 4) and n == (2, 2, 5) and sum(p) == sum(n) == 9 and pip == (3, 0, 3) and pim == (0, 3, 3) and pi0 == (2, 2, 2) and dpp == (6, 0, 3) and colourless == 6,
          f"p = {fmt(p)}, n = {fmt(n)}, pi+ = {fmt(pip)}, pi- = {fmt(pim)}, pi0 = {fmt(pi0)}, Delta++ = {fmt(dpp)}; colourless {colourless}/27")

    # D. ladder and the neutrino
    P = 2 * math.pi
    mu_pred, tau_pred = ME * (7 / 3) ** P, ME * (11 / 3) ** P
    nu_naive = ME * (3 / 3) ** P
    check("D. mu = e + 4 neutral (7), tau = e + 8 neutral (11): the generation step is two neutral DQDs; the ladder keeps mu -0.8 %, tau +1.0 % if the total count sets the scale, but then nu (0,0,0) at n = 3 weighs m_e vs < 0.45 eV (KATRIN 2025): only 'charged strands carry the mode' saves it (POSTULATED)",
          abs(mu_pred / MMU - 1) < 0.011 and abs(tau_pred / MTAU - 1) < 0.011 and nu_naive / (NU_LIMIT_EV * 1e-6) > 1e6,
          f"mu {mu_pred/MMU-1:+.1%}, tau {tau_pred/MTAU-1:+.1%}; naive nu = {nu_naive:.3f} MeV = {nu_naive/(NU_LIMIT_EV*1e-6):.0e} x the KATRIN limit")

    # E. p - n
    md_minus_mu = (MN - MP) + 1.00
    check("E. Same strand count for p and n: R24's count argument is void; m_d - m_u = (m_n - m_p) + 1.0 MeV (EM) = 2.3 MeV (lattice 2.5): the d with one charged strand must outweigh the u with two -- mass is not a strand count",
          abs(md_minus_mu - 2.29) < 0.01 and sum(1 for s in d if s != 0) < sum(1 for s in u if s != 0),
          f"m_d - m_u = {md_minus_mu:.2f} MeV; charged strands: u 2, d 1")

    # F. neutrino size
    R_nu = HBARC / (2 * NU_LIMIT_EV * 1e-6) * 1e-15 * 1e9   # fm -> nm
    check("F. Base spin condition S = R E/c = hbar/2 for a neutrino lighter than 0.45 eV (KATRIN 2025): R >= 219 nm -- a light object is huge, consistent, no prediction",
          abs(R_nu - 219) < 1, f"R >= {R_nu:.0f} nm")

    # G. what is voided
    parity_rule_holds = (sum(u) == F(2, 3) and len(u) % 2 == 0)      # old rule: +2/3 needs even n
    p_is_pip = (p[0], p[1]) == (pip[0], pip[1])
    check("G. Voided by the neutral string: the parity rule (u = 3 strands with +2/3), the identity p = pi+ ring strings (R23), the R24 count argument",
          not parity_rule_holds and not p_is_pip, f"u has 3 strands and +2/3; p charged content {(p[0], p[1])} vs pi+ {(pip[0], pip[1])}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
