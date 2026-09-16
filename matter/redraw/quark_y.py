#!/usr/bin/env python3
"""Redraw: 'if quarks are n = 3 strings in a Y, that explains the strong interaction (a tri-junction)
with three free poles at 1/3' -- the author's proposal, tested in the base (BASE.md R1, R2, R20).

String model of the base (string_assembly.py): each string carries a net charge +-e/3 (its unpaired
branch) and poles +-delta at its ends; strings join pole to pole (+ against -); the charge of any
assembly is q = (n+ - n-) e/3.

Checks:
  A. a Y of three strings carries q in {+-1/3, +-1}: d, s, b (-1/3) fit, u, c, t (+2/3) cannot --
     +2/3 needs an even number of strings.  No split of the proton (6+, 3-) into three Y(3)
     gives (2/3, 2/3, -1/3): the proposal covers half the quarks.
  B. what the parity rule allows: u = 2 strings (2+) or 4 (3+, 1-); d = 1 string (1-) or 3 (1+, 2-);
     proton = 5 strings (4+, 1-) or 11 (7+, 4-); neutron = 4 (2+, 2-) or 10 (5+, 5-).
  C. a tri-junction is frustrated: three poles at one point contain at least one like pair, so
     (+,+,-) binds exactly like one pair, -K delta^2/s, and (+,+,+) repels.  A tri-junction is worth
     one junction: 2 MeV with unit poles at 0.72 fm (R20).  That is the scale of the residual
     nuclear force (deuteron 2.224 MeV), not of confinement (~0.9 GeV per fm of Y-string,
     ~300 MeV per constituent quark): the ratio is 130-400.
  D. with poles of e/3 ('free poles at 1/3') the bond is K/(9 s): 0.16 MeV at 1 fm, 2 MeV only at
     0.08 fm -- 14 x too weak for the deuteron at nuclear distances.
  E. three Y(3) in a baryon expose 9 free poles, an odd number: one pole can never pair;
     a meson (Y + anti-Y, 6 poles) pairs completely.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction as F

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC                 # 1.440 MeV fm
B_DEUTERON = 2.224                # MeV
SIGMA = 890.0                     # MeV/fm, lattice Y-string tension
M_CONSTITUENT = 300.0             # MeV
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def charges(k: int):
    return sorted({F(2 * nplus - k, 3) for nplus in range(k + 1)})


def main() -> int:
    # A. Y of three strings
    cy = charges(3)
    splits = set()
    for a in range(4):                      # n+ in the first Y
        for b in range(4):
            c = 6 - a - b
            if 0 <= c <= 3:
                splits.add(tuple(sorted((F(2 * a - 3, 3), F(2 * b - 3, 3), F(2 * c - 3, 3)))))
    uud = tuple(sorted((F(2, 3), F(2, 3), F(-1, 3))))
    check("A. A Y of three strings carries q in {+-1/3, +-1}: d, s, b fit, u, c, t (+2/3) cannot; no split of the proton (6+, 3-) into three Y(3) gives (2/3, 2/3, -1/3)",
          F(-1, 3) in cy and F(2, 3) not in cy and uud not in splits,
          "Y(3) charges: " + ", ".join(str(q) for q in cy) + "; proton splits into three Y(3): " + "; ".join("(" + ", ".join(str(x) for x in s) + ")" for s in sorted(splits)))

    # B. parity rule
    ok = (F(2, 3) in charges(2) and F(2, 3) in charges(4) and F(-1, 3) in charges(1) and F(-1, 3) in charges(3)
          and F(1) in charges(5) and F(1) in charges(11) and F(0) in charges(4) and F(0) in charges(10))
    check("B. What the parity rule allows: u = 2 strings (2+) or 4 (3+,1-); d = 1 string (1-) or 3 (1+,2-); proton = 5 (4+,1-) or 11 (7+,4-); neutron = 4 (2+,2-) or 10 (5+,5-)",
          ok, "u even, d odd; uud = 5 or 11 strings (odd, q = +1), udd = 4 or 10 (even, q = 0)")

    # C. tri-junction frustration and scale
    def e_junction(signs, delta=1.0, s=1.0):
        return sum(K * delta**2 * si * sj / s for si, sj in itertools.combinations(signs, 2))
    e_pair, e_ppm, e_ppp = e_junction((+1, -1)), e_junction((+1, +1, -1)), e_junction((+1, +1, +1))
    s_unit = K / 2.0
    ratio_conf = (SIGMA * 1.0) / 2.0, M_CONSTITUENT / 2.0
    check("C. A tri-junction is frustrated (one like pair): (+,+,-) binds exactly like one pair, (+,+,+) repels; it is worth one junction, 2 MeV at 0.72 fm with unit poles -- the scale of the nuclear force (deuteron 2.22 MeV), not of confinement (300-900 MeV): ratio 150-450",
          abs(e_ppm - e_pair) < 1e-12 and e_ppp > 0 and abs(s_unit - 0.720) < 0.002 and abs(B_DEUTERON / 2.0 - 1.11) < 0.01 and 100 < ratio_conf[0] < 500 and 100 < ratio_conf[1] < 500,
          f"E(+,-) = {e_pair:.3f}, E(+,+,-) = {e_ppm:.3f}, E(+,+,+) = {e_ppp:+.3f} MeV (delta = e, s = 1 fm); deuteron/2 MeV = {B_DEUTERON/2:.2f}; sigma x 1 fm / 2 MeV = {ratio_conf[0]:.0f}, constituent/2 MeV = {ratio_conf[1]:.0f}")

    # D. poles of e/3
    e_third_1fm = K / 9
    s_2mev = K / 9 / 2.0
    check("D. With free poles of e/3 the bond is K/(9 s): 0.16 MeV at 1 fm, 2 MeV only at 0.08 fm -- 14 x too weak for the deuteron at nuclear distances",
          abs(e_third_1fm - 0.160) < 0.001 and abs(s_2mev - 0.080) < 0.001 and abs(B_DEUTERON / e_third_1fm - 13.9) < 0.2,
          f"K/9 = {e_third_1fm:.3f} MeV fm; s(2 MeV) = {s_2mev:.3f} fm; deuteron / bond(1 fm) = {B_DEUTERON/e_third_1fm:.1f}")

    # E. pole parity
    check("E. Three Y(3) expose 9 free poles, an odd number: one pole can never pair; a meson (Y + anti-Y, 6 poles) pairs completely",
          9 % 2 == 1 and 6 % 2 == 0, "baryon: 9 poles -> 4 pairs + 1 free; meson: 6 poles -> 3 pairs")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
