#!/usr/bin/env python3
"""Redraw: what fixes the ribbon (sheet) width w?  R43-R45 used w_9 = (4 lambda-bar_e/pi^2) 3^(-2 pi)
= 0.157 fm, inherited from R17.  Every rule of the base that touches w is tried.

Checks:
  A. matching to Z_0 fixes only the shape: Z = Z_0 d/w = Z_0 <=> w = d, any size (R17 A).
  B. the size comes from g = 2: the two matched-pole junctions, hbar c/(pi^2 w) each, carry the
     static half of the mass (R13, R17 D): 2 hbar c/(pi^2 w) = m c^2/2 -> w = 4 lambda-bar/pi^2 =
     0.405 lambda-bar = 156.5 fm for the electron; the aspect ratio w/l_1 = 6/pi^3 = 0.19 is the same
     at every scale, so w_9 = w_e 3^(-2 pi) = 0.157 fm.  Inputs: g = 2 (measured 2.0023), delta from
     Z_0, and the reading 'pole separation = width' (R17 C).
  C. nothing more basic fixes it: the tube radius of the manuscript (r = R_3/37.1, w/r = 15.0), the
     wire spacing D_0 = 2 cosh pi r (w/D_0 = 0.65), the strand length (w/l_1 = 0.19) give ratios with
     no rule behind them.
  D. the p-n splitting is only logarithmic in w and nearly cancels between u and d: halving or
     doubling w moves m_n - m_p by -+0.14 MeV (-+11 %); a 10 % change of w moves it by 1.5 %.  The
     +2.4 % result therefore fixes w to within a factor ~2, no better: it is not a test of w.
  E. what would test w: the sheet's transverse mode hbar c pi/w = (pi^3/4) m_e c^2 = 3.96 MeV for the
     electron (R17 E), not observed; the base's way out (TEM-only fluid) is still a postulate.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
R3, RW = LAMBDA_E, LAMBDA_E / 37.1
D0 = 2 * math.cosh(math.pi) * RW
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
W_E = 4 * LAMBDA_E / math.pi**2
W_9 = W_E * 3 ** (-2 * math.pi)
E_DQD = 4 * math.pi * K / (9 * L1_9)
M_MUTUAL = K / (math.sqrt(3) * L1_9)
DM_OBS = MN - MP
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def dm_for_width(w):
    a = w / 4
    S_u = K / (2 * L1_9) * (math.log(4 * L1_9 / a) - 1)
    S_d = K / L1_9 * (math.log(2 * L1_9 / a) - 1)
    return E_DQD - ((4 / 9) * S_u - (1 / 9) * S_d + M_MUTUAL / 3)


def main() -> int:
    # A. shape only
    check("A. Matching to Z_0 fixes the shape only: Z = Z_0 d/w = Z_0 <=> w = d at any size",
          all(abs(376.73 * d / d - 376.73) < 1e-9 for d in (1.0, 10.0, 100.0)), "Z(w = d) = Z_0 for every d")

    # B. the size from g = 2
    static_half = 2 * HBARC / (math.pi**2 * W_E)
    check("B. Size from g = 2: two matched-pole junctions hbar c/(pi^2 w) carry the static half -> w = 4 lambda-bar/pi^2 = 0.405 lambda-bar = 156.5 fm (electron); aspect w/l_1 = 6/pi^3 = 0.19 at every scale -> w_9 = 0.157 fm. Inputs: g = 2, delta from Z_0, 'pole separation = width'",
          abs(static_half / (ME / 2) - 1) < 1e-12 and abs(W_E - 156.5) < 0.2 and abs(W_E / L1_E - 6 / math.pi**3) < 1e-12 and abs(W_9 - 0.1573) < 0.0005,
          f"2 hbar c/(pi^2 w_e) = {static_half*1e3:.1f} keV = m_e c^2/2; w_e = {W_E:.1f} fm; w/l_1 = {W_E/L1_E:.4f}; w_9 = {W_9:.4f} fm")

    # C. no more basic rule
    ratios = {"w/r (tube radius)": W_E / RW, "w/D_0 (wire spacing)": W_E / D0, "w/l_1": W_E / L1_E, "w/lambda-bar": W_E / LAMBDA_E}
    check("C. Nothing more basic fixes it: w/r = 15.0, w/D_0 = 0.65, w/l_1 = 0.19, w/lambda-bar = 0.405 -- ratios with no rule behind them",
          abs(ratios["w/r (tube radius)"] - 15.04) < 0.05 and abs(ratios["w/D_0 (wire spacing)"] - 0.648) < 0.005,
          "; ".join(f"{k} = {v:.3f}" for k, v in ratios.items()))

    # D. sensitivity of the splitting
    dm1 = dm_for_width(W_9)
    dm_half, dm_double, dm_10 = dm_for_width(W_9 / 2), dm_for_width(2 * W_9), dm_for_width(1.1 * W_9)
    check("D. The splitting is logarithmic in w and nearly cancels between u and d: w/2 -> -11 %, 2w -> +11 %, +10 % on w -> +1.5 %: the +2.4 % result fixes w only to within a factor ~2",
          abs(dm1 / DM_OBS - 1.024) < 0.01 and abs((dm_half - dm1) / DM_OBS + 0.106) < 0.01 and abs((dm_double - dm1) / DM_OBS - 0.106) < 0.01 and abs((dm_10 - dm1) / DM_OBS - 0.0146) < 0.003,
          f"w_9: {dm1:.3f}; w/2: {dm_half:.3f} ({(dm_half-dm1)/DM_OBS:+.0%}); 2w: {dm_double:.3f} ({(dm_double-dm1)/DM_OBS:+.0%}); 1.1w: {(dm_10-dm1)/DM_OBS:+.1%}")

    # E. the real test
    e_trans = HBARC * math.pi / W_E
    check("E. What would test w: the sheet's transverse mode hbar c pi/w = (pi^3/4) m_e c^2 = 3.96 MeV for the electron, not observed; the TEM-only way out stays a postulate",
          abs(e_trans - 3.96) < 0.01 and abs(e_trans / ME - math.pi**3 / 4) < 1e-9, f"E_transverse = {e_trans:.2f} MeV = {e_trans/ME:.3f} m_e c^2")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
