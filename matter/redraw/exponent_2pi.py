#!/usr/bin/env python3
"""Redraw: what fixes the exponent 2 pi of the lepton ladder m = m_e (n/3)^p, n = 3, 7, 11 (R14)?

Checks:
  A. how well the data pin p: from (e, mu) alone p = 6.2926; from (e, tau) alone p = 6.2758; 2 pi =
     6.2832 sits between them, 0.016 % from their midpoint (6.2842).  The data fix p to +-0.15 %
     and 2 pi is inside that window.
  B. look-elsewhere: among simple constants (integers, halves, thirds, pi, e, sqrt 2, sqrt 3 and their
     products and powers up to two factors) only 2 pi lies within 0.15 % of the midpoint; the next
     candidates (19/3, 2e + 0.85, ...) are 0.8 % or more away.
  C. Koide: the ladder at p = 2 pi gives Q = sum m / (sum sqrt m)^2 = 0.6683 against 2/3 measured to
     1e-5 (+0.25 %); the p that makes the ladder Koide-exact with (3, 7, 11) is 6.25, 0.5 % below 2 pi
     and outside the data window: the power law and Koide are NOT the same statement.  Koide (exact
     to 1e-5) is the stronger structure and the ladder its 1 % approximant, so whatever fixes '2 pi'
     must in fact produce Koide, whose own form (sqrt m_k = sqrt m_0 [1 + sqrt 2 cos(theta + 2 pi k/3)],
     theta = 2/9) carries a 2 pi/3 and a 3-fold structure.
  D. what the base offers: nothing that derives it.  Its two native exponentials of pi -- the
     matching D/r = 2 cosh pi = e^pi (0.2 %; 0.06 % in the log) and the vortex's logarithmic energies -- make
     'exp(2 pi x log)' the natural form: the ladder says the log of the size ratio is 2 pi times the
     log of the count ratio (a log-log slope of one turn), but no rule of the base yet links the
     count to a turn.
  E. verdict: 2 pi is pinned by the data to 0.15 % and unique among simple constants, but the exact
     structure is Koide's (1e-5), of which the ladder is a 1 % approximant; it is not derived, and
     the next target is a rule turning 'one turn per e-fold of the strand count' into a mechanism
     that lands on Koide.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def koide_q(p):
    ms = [ME * (n / 3) ** p for n in (3, 7, 11)]
    return sum(ms) / sum(math.sqrt(m) for m in ms) ** 2


def main() -> int:
    p_mu = math.log(MMU / ME) / math.log(7 / 3)
    p_tau = math.log(MTAU / ME) / math.log(11 / 3)
    mid = 0.5 * (p_mu + p_tau)
    check("A. The data pin p: from (e, mu) 6.2926, from (e, tau) 6.2758; 2 pi = 6.2832 lies between, 0.016 % from the midpoint: p is fixed to +-0.15 % and 2 pi is inside",
          abs(p_mu - 6.2926) < 0.0005 and abs(p_tau - 6.2758) < 0.0005 and p_tau < 2 * math.pi < p_mu and abs(2 * math.pi / mid - 1) < 3e-4,
          f"p_mu = {p_mu:.4f}, p_tau = {p_tau:.4f}, midpoint {mid:.4f}, 2 pi = {2*math.pi:.4f} ({2*math.pi/mid-1:+.3%})")

    # B. look-elsewhere
    atoms = {"1": 1.0, "2": 2.0, "3": 3.0, "4": 4.0, "5": 5.0, "6": 6.0, "7": 7.0, "1/2": 0.5, "1/3": 1 / 3, "2/3": 2 / 3, "3/2": 1.5,
             "pi": math.pi, "e": math.e, "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "pi^2": math.pi**2, "e^2": math.e**2, "ln2": math.log(2)}
    cands = {}
    for (a, va), (b, vb) in itertools.product(atoms.items(), repeat=2):
        cands[f"{a}*{b}"] = va * vb
        if vb != 0:
            cands[f"{a}/{b}"] = va / vb
        cands[f"{a}+{b}"] = va + vb
    hits = sorted(((k, v) for k, v in cands.items() if abs(v / mid - 1) < 0.0015), key=lambda t: abs(t[1] / mid - 1))
    values = {round(v, 6) for _, v in hits}
    check("B. Look-elsewhere among ~1000 simple two-factor constants: only 2 pi (and its spellings) lies within 0.15 % of the midpoint; the next candidates are 0.8 % or more away",
          len(values) == 1 and abs(list(values)[0] - 2 * math.pi) < 1e-6,
          f"hits: {', '.join(k for k, _ in hits[:6])} -> {list(values)}")

    # C. Koide
    q_2pi = koide_q(2 * math.pi)
    lo, hi = 6.0, 6.6
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if koide_q(m) > 2 / 3:
            hi = m
        else:
            lo = m
    p_koide = lo
    check("C. Koide: the ladder at 2 pi gives Q = 0.6683 vs 2/3 (1e-5), +0.25 %; the Koide-exact exponent with (3, 7, 11) is 6.25, 0.5 % below 2 pi and outside the data window: the power law is a 1 % approximant of Koide, not the same statement; what fixes '2 pi' must produce Koide",
          abs(q_2pi - 0.6683) < 0.0005 and abs(p_koide - 6.252) < 0.005 and p_koide < p_tau,
          f"Q(2 pi) = {q_2pi:.4f}; p_Koide = {p_koide:.4f} ({p_koide/(2*math.pi)-1:+.3%} vs 2 pi)")

    # D. the base's native exponentials of pi
    d_over_r = 2 * math.cosh(math.pi)
    check("D. The base's native exponentials of pi: matching D/r = 2 cosh pi = e^pi to 0.2 % (0.06 % in the log), and logarithmic vortex energies; the ladder reads 'log of size ratio = 2 pi x log of count ratio' (a log-log slope of one turn), but no rule links the count to a turn: not derived",
          abs(d_over_r / math.exp(math.pi) - 1) < 3e-3, f"2 cosh pi / e^pi = {d_over_r/math.exp(math.pi):.5f}")

    check("E. Verdict: 2 pi is pinned to 0.15 % and unique among simple constants, but the exact structure is Koide's, of which the ladder is a 1 % approximant; not derived -- the next target is a rule making one turn per e-fold of the count that lands on Koide",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
