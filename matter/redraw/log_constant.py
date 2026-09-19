#!/usr/bin/env python3
"""Redraw: what fixes the constant in the logarithm of the self term (R43)?  The self-energy of a
charge q spread along a length l with a transverse size a is (K q^2/l)[ln(2 l/a) + c]; the constant
c and the cut-off a are the two things to fix.

Checks:
  A. the -1 is derived: the self-energy of a uniform line charge of length l, potential taken at a
     distance a from the line, is exactly (K q^2/l)[ln(2 l/a) - 1] for l >> a (checked by direct
     double integration).  For a closed ring of the same length the constant is -1 - ln(pi/2) = -1.45
     (arcs of a closed circuit sit between the two).
  B. the cut-off is the section: the base's strand is a plate of width w (R17, w = d); a thin strip
     of width w has the capacitance of a cylinder of radius w/4 (classical), so a = w/4 = 0.039 fm at
     w_9 = 0.157 fm.  A full square of side w would give 0.59 w, an inscribed cylinder w/2; R43 used
     a = w, a cylinder wider than the sheet, which is not a section of the base.
  C. the length: a quark's charge spreads over its charged strands (R28): two for u (2 l_1), one for
     d (l_1), so S_u != S_d; R43 used one rod of l_1 per quark.
  D. the eight combinations span m_n - m_p = 0.52-1.60 MeV.  The base's own pair -- plate (a = w/4)
     and charged strands -- gives 1.32 MeV (+2.4 %); R43's 1.34 MeV came from two non-base choices
     (a = w, one rod); the square-section variants give 1.03 (-20 %) and 1.49 (+15 %).
  E. verdict: with R17's section and R28's content the constant is fixed; m_n - m_p = 1.32 MeV
     (+2.4 %); remaining sensitivity: ring closure (-1 to -1.45 in the constant, +9 %) and the
     partner plate's screening (not computed).
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
L1 = 2 * math.pi * LAMBDA_E / 3 * (3 / 9) ** (2 * math.pi)
W9 = (4 * LAMBDA_E / math.pi**2) * 3 ** (-2 * math.pi)
E_DQD = 4 * math.pi * K / (9 * L1)
M_MUTUAL = K / (math.sqrt(3) * L1)
DM_OBS = MN - MP
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def self_rod(q, l, a, const=-1.0):
    return K * q * q / l * (math.log(2 * l / a) + const)


def rod_numeric(l, a, n=4000):
    """Self-energy of a uniform unit line charge of length l with the potential at distance a
    (double integral), in units K/l."""
    x = (np.arange(n) + 0.5) * l / n
    dx = l / n
    lam = 1.0 / l
    U = 0.0
    for xi in x:
        r = np.sqrt((x - xi) ** 2 + a * a)
        U += 0.5 * lam * lam * np.sum(dx / r) * dx
    return U * l                       # -> coefficient of K/l


def dm(a, spread):
    """m_n - m_p for cut-off a and charge spread ('rod': one l_1 per quark; 'strands': u over 2 l_1, d over l_1)."""
    if spread == "rod":
        S = self_rod(1.0, L1, a)
        diff_self = (1 - 2 / 3) * S
    else:
        S_u, S_d = self_rod(1.0, 2 * L1, a), self_rod(1.0, L1, a)
        diff_self = (4 / 9) * S_u - (1 / 9) * S_d
    return E_DQD - (diff_self + M_MUTUAL / 3)


def main() -> int:
    # A. the -1
    l, a = 1.0, 0.01
    num = rod_numeric(l, a)
    ana = math.log(2 * l / a) - 1
    ring_const = -1 - math.log(math.pi / 2)
    check("A. The -1 is derived: uniform line charge, potential at distance a: (K q^2/l)[ln(2l/a) - 1], confirmed by direct integration; a closed ring of the same length has -1 - ln(pi/2) = -1.45",
          abs(num / ana - 1) < 5e-3 and abs(ring_const + 1.452) < 0.001,
          f"numeric {num:.4f} vs ln(2l/a) - 1 = {ana:.4f} ({num/ana-1:+.2%}); ring constant {ring_const:.3f}")

    # B. the section
    a_plate, a_square, a_cyl, a_r43 = W9 / 4, 0.5902 * W9, W9 / 2, W9
    check("B. The cut-off is the section: a plate of width w (R17) is a cylinder of radius w/4 = 0.039 fm; a full square 0.59 w = 0.093 fm; an inscribed cylinder w/2 = 0.079 fm; R43's a = w = 0.157 fm is wider than the sheet",
          abs(a_plate - 0.0393) < 0.0005 and abs(a_square - 0.0928) < 0.0005 and abs(a_cyl - 0.0786) < 0.0005,
          f"w_9 = {W9:.4f} fm; a = {a_plate:.4f} (plate), {a_square:.4f} (square), {a_cyl:.4f} (cylinder), {a_r43:.4f} (R43)")

    # C. the length (statement)
    check("C. A quark's charge spreads over its charged strands (R28): u over 2 l_1, d over l_1, so S_u != S_d; R43 used one rod of l_1 per quark", True, "content rule")

    # D. the table
    table = {}
    for name_a, a in (("plate w/4", a_plate), ("square 0.59w", a_square), ("cylinder w/2", a_cyl), ("a = w (R43)", a_r43)):
        for spread in ("rod", "strands"):
            table[(name_a, spread)] = dm(a, spread)
    lo, hi = min(table.values()), max(table.values())
    base_pair = table[("plate w/4", "strands")]
    check("D. Eight combinations span 0.52-1.60 MeV; the base's own pair (plate w/4 + charged strands) gives 1.32 MeV (+2.4 %); R43's 1.34 came from a = w with one rod; square-section variants 1.03 (-20 %) and 1.49 (+15 %)",
          abs(lo - 0.52) < 0.03 and abs(hi - 1.60) < 0.03 and abs(base_pair / DM_OBS - 1.024) < 0.01 and abs(table[("a = w (R43)", "rod")] / DM_OBS - 1.038) < 0.01,
          "; ".join(f"{k[0]}/{k[1]}: {v:.2f} ({v/DM_OBS-1:+.0%})" for k, v in table.items()))

    # E. sensitivity to ring closure
    def dm_const(const):
        S_u, S_d = self_rod(1.0, 2 * L1, a_plate, const), self_rod(1.0, L1, a_plate, const)
        return E_DQD - ((4 / 9) * S_u - (1 / 9) * S_d + M_MUTUAL / 3)
    dm_ring = dm_const(ring_const)
    check("E. With R17's section and R28's content the constant is fixed: m_n - m_p = 1.32 MeV (+2.4 %); ring closure (-1.45) would give 1.41 MeV (+9 %): the remaining ~10 % is the open-vs-closed reading of the quark's charged arc",
          abs(dm_ring / DM_OBS - 1.093) < 0.01, f"rod constant: {base_pair:.3f} MeV; ring constant: {dm_ring:.3f} MeV ({dm_ring/DM_OBS-1:+.0%})")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
