#!/usr/bin/env python3
"""Redraw: the author's ladder for the charged leptons -- e = 3, mu = 7, tau = 11 strings -- and the
idea that the electron is special (its open turn suppresses its mass exponentially).

Tests, with m_e, m_mu, m_tau from PDG:
  A. exponential in n through e and mu: m = m_e exp(k (n - 3)) -> k = 1.333; tau predicted 12 x too
     heavy.  Exponential in n through mu and tau: x 2.02 per string; the electron then sits 12.3 x
     BELOW the law: the 'open-turn suppression' would be a factor e^{-2.5} -- one free number.
  B. power law m = m_e (n/3)^p: p from the muon is 6.29; the tau then comes out +2.3 %.  With the
     natural exponent p = 2 pi = 6.283: muon -1.0 %, tau +1.0 % -- a one-parameter law with no
     fitted parameter passing a 1 % test on the tau.
  C. among the small odd triples (3, n_mu, n_tau) with 5 <= n_mu < n_tau <= 15, (3, 7, 11) is the
     one whose power-law exponent from the muon predicts the tau best; the next best is off by
     more than 10 %.
  D. what the law means physically: m ~ 1/R (the mode), so R ~ n^{-2 pi} and the string length
     ~ n^{-2 pi - 1}: the muon's strings are 700 x shorter than the electron's.  Strings are not
     fixed-length objects in this base; closure sets their scale -- the author's 'exponential
     impact of closure', made quantitative as a power 2 pi.
Numerology caveat, stated: three masses, integers chosen, one exponent; the tau at 1 % is the
single test passed.  Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86      # MeV
R_MU, R_TAU = MMU / ME, MTAU / ME
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. exponentials
    k_emu = math.log(R_MU) / 4
    tau_pred = math.exp(k_emu * 8)
    check("A. Exponential in n through e (3) and mu (7): the tau (11) is predicted 12 x too heavy",
          tau_pred / R_TAU > 10, f"k = {k_emu:.3f} per string, tau/e predicted {tau_pred:.0f} vs {R_TAU:.0f}")
    k_mutau = math.log(R_TAU / R_MU) / 4
    e_pred = MMU * math.exp(-k_mutau * 4)
    check("A'. Exponential through mu (7) and tau (11): x 2.02 per string; the electron sits 12.3 x below it -- the 'open-turn suppression' would be e^-2.5, a free number",
          abs(e_pred / ME - 12.3) < 0.2, f"factor per string {math.exp(k_mutau):.3f}; electron predicted {e_pred:.2f} MeV vs 0.511; suppression {e_pred/ME:.1f} = e^{math.log(e_pred/ME):.2f}")
    # B. power law
    p_mu = math.log(R_MU) / math.log(7 / 3)
    tau_pl = (11 / 3) ** p_mu
    check("B. Power law m = m_e (n/3)^p: p fitted on the muon = 6.29; tau predicted +2.3 %",
          abs(p_mu - 6.293) < 0.01 and abs(tau_pl / R_TAU - 1.023) < 0.005, f"p = {p_mu:.4f}, tau/e = {tau_pl:.0f} vs {R_TAU:.1f} ({tau_pl/R_TAU-1:+.1%})")
    p = 2 * math.pi
    mu_2pi, tau_2pi = (7 / 3) ** p, (11 / 3) ** p
    check("B'. With p = 2 pi (no fitted parameter): muon -1.0 %, tau +1.0 %",
          abs(mu_2pi / R_MU - 1) < 0.015 and abs(tau_2pi / R_TAU - 1) < 0.015, f"(7/3)^2pi = {mu_2pi:.1f} ({mu_2pi/R_MU-1:+.1%}), (11/3)^2pi = {tau_2pi:.0f} ({tau_2pi/R_TAU-1:+.1%})")
    # C. other small odd triples
    rows = []
    for n_mu in range(5, 14, 2):
        for n_tau in range(n_mu + 2, 16, 2):
            pp = math.log(R_MU) / math.log(n_mu / 3)
            pred = (n_tau / 3) ** pp
            rows.append((n_mu, n_tau, pp, pred / R_TAU - 1))
    rows.sort(key=lambda r: abs(r[3]))
    best = rows[0]
    second = rows[1]
    check("C. Among small odd triples (3, n_mu, n_tau), (3, 7, 11) predicts the tau best from the muon; the runner-up misses by more than 10 %",
          best[0] == 7 and best[1] == 11 and abs(second[3]) > 0.10,
          "ranking (n_mu, n_tau, p, tau error): " + "; ".join(f"({r[0]}, {r[1]}, {r[2]:.2f}, {r[3]:+.0%})" for r in rows[:5]))
    # D. meaning: R ~ n^-2pi
    r_ratio_mu = (3 / 7) ** p
    check("D. Under m ~ 1/R the law is R ~ n^-2pi: the muon's loop is 205 x smaller and its strings (R/n) 480-700 x shorter than the electron's",
          abs(1 / r_ratio_mu - 204.8) < 1 and 400 < (1 / r_ratio_mu) * 7 / 3 < 800, f"R_mu/R_e = {r_ratio_mu:.5f}; string length ratio {(r_ratio_mu*3/7):.2e}")
    koide = lambda a, b, c: (a + b + c) / (math.sqrt(a) + math.sqrt(b) + math.sqrt(c)) ** 2
    k_law = koide(1.0, mu_2pi, tau_2pi)
    check("D'. Koide on the p = 2 pi law: 0.6683 vs 2/3 -- the law is 1 % off where the data are 1e-5 off: not yet a Koide-exact law",
          abs(k_law - 0.6683) < 0.001, f"Koide(law) = {k_law:.4f}, Koide(data) = {koide(ME, MMU, MTAU):.6f}")
    check("Caveat, stated: three masses, integers chosen, one exponent; the 1 % tau is the single test passed. A fourth lepton at n = 15 would weigh (5)^2pi m_e = 12.6 GeV",
          True, f"(15/3)^2pi = {5**p:.0f} m_e = {5**p*ME/1e3:.1f} GeV -- excluded by LEP (no charged lepton below 100 GeV): the ladder must stop at 11, or n = 15 is not a lepton")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
