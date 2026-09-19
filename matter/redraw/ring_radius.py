#!/usr/bin/env python3
"""Redraw: what fixes the nucleon's ring radius at 0.67 fm (R48)?  R = hbar c / (2 E_circ) with
E_circ(N) = (Delta - N)/2 = 147 MeV, i.e. R = 0.671 fm; in the base's own language R = lambda-bar_p/f
with the circulating fraction f = 0.313 (the electron has f = 3/4).

Checks:
  A. the target and its base reading: R = 0.671 fm, f = 0.313, E_circ/m = f/2 = 0.157 (as for the
     electron, 3/8 = f/2).
  B. every length of the base at scale 9 against 0.671 fm: l_1 = 0.813 (x0.826), the circuit radius
     R_c = 3 l_1/2 pi = 0.388 (x1.73 = sqrt 3 to 0.1 %), 2 R_c = 0.776, sqrt 2 R_c = 0.549, lambda-bar_p
     = 0.210 (x3.19, pi to 1.6 %), R25's 0.587 and 0.789, the pion ring 1.409.
  C. the star of three circuit-rings: rings through the centre have their centres at R_c and a
     quark spacing sqrt 3 R_c = 0.672 fm -- the needed radius to 0.1 %; three mutually tangent rings
     have an outer envelope R_c (1 + 2/sqrt 3) = 0.836 fm = r_p to 0.6 % (r_p is an input nowhere),
     centres at 0.448 fm and an rms radius 0.593 fm.
  D. no mechanism selects sqrt 3 R_c: a flow along the quark triangle's edges carries S = R_in E/c
     with the inradius R_in = R_c/2 (E_circ would be 508 MeV); the unpaired circuit's own radius
     gives 254 MeV; neither is 147.  The identification R = sqrt 3 R_c = (3 sqrt 3/2 pi) l_1(9) is a
     coincidence at 0.1 % without a derivation.
  E. if accepted, everything follows from l_1(9) alone: N = 910 MeV (-3 %), Delta = 1204 (-2.3 %),
     mu(Delta++) = 6.4 mu_N, r_p = 0.836 fm (+0.6 %); but mu_p would need 0.87 e circulating at that
     radius, not a clean charge (recorded).
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME, MP, MDELTA = 0.51099895, 938.27209, 1232.0
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
L1 = (2 * math.pi * LAMBDA_E / 3) * (3 / 9) ** (2 * math.pi)
R_C = 3 * L1 / (2 * math.pi)
STATIC = math.pi * HBARC / L1
RP = 0.8409
MU_P = 2.79284734
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    e_circ = (MDELTA - MP) / 2
    R = HBARC / (2 * e_circ)
    f = LAMBDA_P / R
    check("A. Target: R = hbar c/(2 x 147 MeV) = 0.671 fm; in base language R = lambda-bar_p/f with f = 0.313 (electron 3/4), E_circ/m = f/2 = 0.157",
          abs(R - 0.672) < 0.002 and abs(f - 0.313) < 0.002 and abs(e_circ / MP - f / 2) < 1e-12,
          f"R = {R:.4f} fm, f = {f:.4f}, E_circ/m_p = {e_circ/MP:.4f}")

    lengths = {"l_1(9)": L1, "R_c = 3 l_1/2pi": R_C, "2 R_c": 2 * R_C, "sqrt2 R_c": math.sqrt(2) * R_C, "lambda_p": LAMBDA_P,
               "R25 e-ring": MU_P * LAMBDA_P, "R25 R+": 3.749 * LAMBDA_P, "pion ring r_e/2": ALPHA * LAMBDA_E / 2}
    ratios = {k: R / v for k, v in lengths.items()}
    check("B. Base lengths at scale 9 vs 0.671 fm: R_c = 0.388 fm is sqrt 3 times smaller (0.1 %), lambda-bar_p pi times smaller (1.6 %); the rest are unrelated",
          abs(ratios["R_c = 3 l_1/2pi"] / math.sqrt(3) - 1) < 0.002 and abs(ratios["lambda_p"] / math.pi - 1) < 0.02,
          "; ".join(f"{k} = {v:.3f} fm (R/x = {ratios[k]:.3f})" for k, v in lengths.items()))

    spacing = math.sqrt(3) * R_C
    envelope = R_C * (1 + 2 / math.sqrt(3))
    centres = 2 * R_C / math.sqrt(3)
    rms = math.sqrt(centres**2 + R_C**2)
    check("C. Star of three circuit-rings: through the centre, quark spacing sqrt 3 R_c = 0.672 fm (the needed radius to 0.1 %); mutually tangent, outer envelope R_c(1 + 2/sqrt 3) = 0.836 fm = r_p to 0.6 % (r_p is an input nowhere), centres at 0.448 fm, rms 0.593 fm",
          abs(spacing / R - 1) < 0.002 and abs(envelope / RP - 1) < 0.01 and abs(centres - 0.448) < 0.002 and abs(rms - 0.593) < 0.002,
          f"sqrt3 R_c = {spacing:.4f} fm ({spacing/R-1:+.1%}); envelope = {envelope:.4f} fm ({envelope/RP-1:+.1%} vs r_p); centres {centres:.3f}, rms {rms:.3f}")

    e_edges = HBARC / (2 * (R_C / 2))
    e_own = HBARC / (2 * R_C)
    check("D. No mechanism selects sqrt 3 R_c: a flow along the triangle's edges carries S = R_in E/c with R_in = R_c/2 (E_circ = 508 MeV); the unpaired circuit's own radius gives 254 MeV; neither is 147. R = sqrt 3 R_c is a 0.1 % coincidence, not a derivation",
          abs(e_edges - 508) < 2 and abs(e_own - 254) < 1, f"edges: {e_edges:.0f} MeV; own radius: {e_own:.0f} MeV; needed 147")

    mN, mD = STATIC + HBARC / (2 * spacing), STATIC + 3 * HBARC / (2 * spacing)
    mu_dpp = 2 * spacing / LAMBDA_P
    q_p = MU_P * LAMBDA_P / spacing
    check("E. If accepted, all from l_1(9): N = 910 MeV (-3 %), Delta = 1204 (-2.3 %), mu(Delta++) = 6.4 mu_N, r_p = 0.836 fm (+0.6 %); but mu_p would need 0.87 e circulating at that radius (recorded)",
          abs(mN / MP - 0.970) < 0.003 and abs(mD / MDELTA - 0.977) < 0.003 and abs(mu_dpp - 6.39) < 0.05 and abs(q_p - 0.874) < 0.005,
          f"N = {mN:.0f}, Delta = {mD:.0f}, mu(Delta++) = {mu_dpp:.2f}, q_p = {q_p:.3f} e")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
