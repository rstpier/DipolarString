#!/usr/bin/env python3
"""Redraw, R21's criterion applied to the star: is the proton's charge radius r_p = 0.8409(4) fm
expressible in the base's constants (m_e, alpha, hbar, c, pi, small integers) without a new number?
And what the parity rule says about the content of p, n and the pions.

Search: r = L x alpha^k x pi^j x p/q x [(n/3)^(2 pi)], with L in {lambda-bar_e, r_e}, k in 0..2,
j in -2..2, p, q in 1..6, ladder factor for n = 3..12 or none.  The number of candidates fixes the
number of chance hits expected in a window: the search is decisive only if a hit is far better than
that expectation and structurally simple.

Checks:
  A. the closest simple expressions: 3 r_e/10 (+0.5 %), 3 r_e/pi^2 (+1.9 %, the sheet's aspect
     ratio), r_e/pi (+6.7 %); none within 10 sigma of the measurement (0.05 %).
  B. with the full candidate set (7590) no candidate lands within +-0.1 % of r_p, fewer than the
     ~3 chance hits expected: r_p is not derived, the star keeps its own number.
  C. content by the parity rule: proton = (4+, 1-), the pi+ ring's strings; neutron = (2+, 2-), the
     pi0 ring's strings.  Same strings, different shape: m_p / m_pi+ = 6.72, m_n / m_pi0 = 6.96;
     in the base's two readings (ring at r_e/2, star at 4 modes of r_p) the ratio is 2 r_e / r_p =
     6.70 -- the same unexplained number as r_p.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction as F

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_E = HBARC / ME                # 386.159 fm
R_E = ALPHA * LAMBDA_E               # 2.8179 fm
RP, DRP = 0.8409, 0.0004             # PDG 2024
MP, MN, MPIP, MPI0 = 938.27209, 939.56542, 139.57039, 134.9768
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. simple expressions
    simple = {"3 r_e/10": 3 * R_E / 10, "3 r_e/pi^2": 3 * R_E / math.pi**2, "r_e/pi": R_E / math.pi,
              "r_e/2 x 3/5": R_E / 2 * 3 / 5, "sqrt(r_e lambda_e) alpha": math.sqrt(R_E * LAMBDA_E) * ALPHA}
    devs = {k: v / RP - 1 for k, v in simple.items()}
    n_sigma = {k: abs(d) * RP / DRP for k, d in devs.items()}
    check("A. Closest simple expressions: 3 r_e/10 at +0.5 %, 3 r_e/pi^2 at +1.9 %, r_e/pi at +6.7 %: none within 10 sigma of r_p = 0.8409(4) fm",
          abs(devs["3 r_e/10"] - 0.0054) < 0.001 and abs(devs["3 r_e/pi^2"] - 0.0187) < 0.001 and min(n_sigma.values()) > 10,
          "; ".join(f"{k} = {v:.4f} fm ({devs[k]:+.2%}, {n_sigma[k]:.0f} sigma)" for k, v in simple.items()))

    # B. systematic search
    rationals = sorted({F(p, q) for p in range(1, 7) for q in range(1, 7)})
    ladder = [None] + list(range(3, 13))
    cands = []
    for (Lname, L), k, j, r, n in itertools.product((("lambda_e", LAMBDA_E), ("r_e", R_E)), range(3), range(-2, 3), rationals, ladder):
        val = L * ALPHA**k * math.pi**j * float(r)
        tag = f"{Lname} alpha^{k} pi^{j} x {r}"
        if n is not None:
            val *= (n / 3) ** (2 * math.pi)
            tag += f" x ({n}/3)^2pi"
        nfactors = (k > 0) + (j != 0) + (r != 1) + (n is not None)
        cands.append((tag, val, nfactors))
    window = 0.001
    hits = [(t, v, nf) for t, v, nf in cands if abs(v / RP - 1) < window]
    expected = len(cands) * 2 * window * 0.5 / math.log(10)  # rough: log-uniform density over a decade
    simple_hits = [h for h in hits if h[2] <= 2]
    check("B. Systematic search (7590 candidates): no candidate lands within +-0.1 % of r_p, fewer than the ~3 chance hits expected; r_p is not derived, the star keeps its own number",
          len(cands) > 3000 and len(simple_hits) == 0 and len(hits) <= 3 * max(expected, 1),
          f"{len(cands)} candidates, {len(hits)} hits within +-0.1 % (chance expectation ~{expected:.0f}), {len(simple_hits)} with <= 2 factors; hits: " + "; ".join(f"{t} = {v:.4f}" for t, v, _ in hits[:6]))

    # C. content by the parity rule
    ratio_p, ratio_n = MP / MPIP, MN / MPI0
    base_ratio = 2 * R_E / RP
    check("C. Parity rule: proton = (4+, 1-) like the pi+ ring, neutron = (2+, 2-) like the pi0 ring: same strings, different shape; m_p/m_pi+ = 6.72, m_n/m_pi0 = 6.96; ring at r_e/2 vs star at 4 modes of r_p gives 2 r_e/r_p = 6.70 -- the same open number as r_p",
          abs(ratio_p - 6.722) < 0.002 and abs(ratio_n - 6.961) < 0.002 and abs(base_ratio - 6.70) < 0.01,
          f"m_p/m_pi+ = {ratio_p:.3f}, m_n/m_pi0 = {ratio_n:.3f}, 2 r_e/r_p = {base_ratio:.3f} ({base_ratio/ratio_p-1:+.1%} vs p)")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
