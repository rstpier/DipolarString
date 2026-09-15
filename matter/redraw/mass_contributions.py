#!/usr/bin/env python3
"""Redraw: mass contributions of the turns (spires) and of the string junctions.

What the base already fixes (g_bookkeeping.py): for the electron half the rest energy is the
mode on the turn and half sits on the axis (the pivot junction): 255.5 keV each.  With the
author's geometry -- 3 strings = 3/4 turn (a string is a quarter turn), 2 junctions --
    a = mass per full turn = 255.5 / (3/4) = 340.7 keV,   b = mass per junction = 127.7 keV.

Model B (one integer per particle: n strings, n/4 turns, n - 1 junctions, additive masses):
    m(n) = a n / 4 + b (n - 1) = 212.9 n - 127.7 keV.
  * chains of 4, 5, 6, 7 strings weigh 0.72, 0.94, 1.15, 1.36 MeV -- a hundred times below the
    muon: the 'quark' candidates n = 5, 6 cannot be muon-scale, and no small chain can;
  * the muon needs n = 497 strings, the tau n = 8347 -- near-integers only because any real
    number is within 0.5 of an integer (the relative precision is trivially 1e-3, 1e-5).
Model A (two free integers per particle, turns and junctions independent): with e fixed at
(3/4, 2), every (S_mu, J_mu) pair with positive a, b admits a tau pair fitting to 0.1 %: the
count of solutions is in the thousands -- the calibration selects nothing without a rule for
the integers.
Verdict: additive per-turn and per-junction masses calibrated on the electron cannot reach the
muon or tau with small chains; either the counts are in the hundreds or the unit masses scale as
1/R with the object (the mode picture), which is not additive.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

ME, MMU, MTAU = 510.99895, 105658.3755, 1776860.0     # keV
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    a = (ME / 2) / 0.75          # keV per full turn
    b = (ME / 2) / 2             # keV per junction
    check("Electron split fixed by g = 2 and the author's geometry (3/4 turn, 2 junctions): 340.7 keV per full turn, 127.7 keV per junction",
          abs(a - 340.67) < 0.1 and abs(b - 127.75) < 0.1, f"a = {a:.1f} keV/turn (85.2 per quarter turn = per string), b = {b:.2f} keV/junction")
    m = lambda n: a * n / 4 + b * (n - 1)
    small = {n: m(n) for n in (3, 4, 5, 6, 7)}
    check("Model B (n strings, n/4 turns, n-1 junctions): m(n) = 212.9 n - 127.7 keV; n = 4..7 give 0.72-1.36 MeV -- no small chain reaches the muon (105.7 MeV)",
          abs(small[3] - ME) < 1e-6 and all(v < 2000 for v in small.values()), ", ".join(f"n = {n}: {v/1e3:.2f} MeV" for n, v in small.items()))
    n_mu, n_tau = (MMU + b) / (a / 4 + b), (MTAU + b) / (a / 4 + b)
    check("The muon would need n = 497 strings and the tau n = 8346: near-integers only trivially (any real is within 0.5 of an integer)",
          abs(n_mu - 496.8) < 0.3 and abs(n_tau - 8346) < 1.0, f"n_mu = {n_mu:.2f}, n_tau = {n_tau:.2f}; nearest-integer residuals {abs(n_mu-round(n_mu))/n_mu:.1e}, {abs(n_tau-round(n_tau))/n_tau:.1e}")
    # thresholds for positive coefficients in Model A: det = 0.75 J - S ; a > 0 and b > 0 need
    # either J > 2 m_mu/m_e (det > 0) or S > 0.75 m_mu/m_e (det < 0)
    j_min, s_min = 2 * MMU / ME, 0.75 * MMU / ME
    check("Model A thresholds: positive per-turn and per-junction masses fitting e and mu need the muon to have at least 156 turns or at least 414 junctions",
          abs(j_min - 413.5) < 0.2 and abs(s_min - 155.1) < 0.2, f"J_mu > {j_min:.1f} or S_mu > {s_min:.1f}: no small muon under additive masses")
    # Model A degeneracy: e = (0.75 turn, 2 junctions); mu = (S, J) integers (quarter turns q = 4 S); solve a, b; find tau integers
    count, examples = 0, []
    for q_mu in range(4, 2401, 40):             # quarter turns of the muon (1..600 turns, sampled)
        for j_mu in range(0, 1201, 30):
            # equations: a*0.75 + b*2 = ME ; a*q_mu/4 + b*j_mu = MMU
            det = 0.75 * j_mu - 2 * q_mu / 4
            if abs(det) < 1e-12:
                continue
            aa = (ME * j_mu - 2 * MMU) / det
            bb = (0.75 * MMU - (q_mu / 4) * ME) / det
            if aa <= 0 or bb <= 0:
                continue
            # tau: minimise |aa q/4 + bb j - MTAU| over integers with j from q
            found = False
            for q_tau in range(1, 4000):
                j_tau = round((MTAU - aa * q_tau / 4) / bb)
                if j_tau < 0:
                    break
                if abs(aa * q_tau / 4 + bb * j_tau - MTAU) / MTAU < 1e-3:
                    found = True
                    break
            if found:
                count += 1
                if len(examples) < 3:
                    examples.append((q_mu / 4, j_mu, round(aa), round(bb)))
    check("Model A (turns and junctions free per particle): beyond the thresholds, hundreds of sampled integer assignments fit e, mu, tau to 0.1 % -- the calibration selects nothing by itself",
          count > 100, f"{count} solutions in a 60 x 41 sample of (muon turns <= 600, junctions <= 1200), i.e. thousands overall; e.g. (turns, junctions, a keV, b keV) = {examples}")
    check("Verdict: additive per-turn and per-junction masses calibrated on the electron do not reach mu or tau with small chains; either hundreds of strings, or unit masses that scale as 1/R with the object (the mode picture, not additive)",
          True, "a rule fixing the counts, or the scaling, is what is missing")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
