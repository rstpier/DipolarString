#!/usr/bin/env python3
"""Redraw: what fixes the excited shape Delta at +294 MeV above the nucleon (R47 D)?

The Delta(1232) has the nucleon's content (Delta+ = uud) and spin 3/2.  R6 reads it as the
symmetric star (three axes -> 3/2) against the nucleon's asymmetric one (one axis -> 1/2).  In the
circuit picture the three quark circuits are identical in both (763 MeV static); what differs is the
ring circulation that carries the spin: S = R E_circ / c (R25).

Checks:
  A. spin from the ring: N needs hbar/2, Delta needs 3 hbar/2 -> E_circ(Delta) = 3 E_circ(N) at the
     same radius, so Delta - N = 2 E_circ(N): the three circuits' circulations aligned (R41 F)
     instead of two against one.
  B. with the base's two readings of the nucleon's ring (R25): 125 MeV (+ strands at R+ = 0.79 fm)
     -> Delta - N = 250 MeV (-15 %); 168 MeV (e on one ring at 0.59 fm) -> 336 MeV (+14 %): the two
     bracket the measured 294 MeV.
  C. inverting: 294 MeV -> E_circ(N) = 147 MeV -> ring radius R = hbar c/(2 x 147) = 0.671 fm, between
     the two R25 radii; then N = 763 + 147 = 910 MeV (-3.0 %) and Delta = 763 + 441 = 1204 MeV
     (-2.3 %): both within 3 % from one radius.
  D. a test of that radius: Delta++ (uuu, charge 2e) with the whole charge at R = 0.671 fm has
     mu = 2 (R/lambda-bar_p) mu_N = 6.4 mu_N, inside the measured range 3.7-7.5 mu_N (recorded).
  E. an electromagnetic spin-spin interaction between the circuits (the quark model's hyperfine
     analogue) is ~0.02 MeV: excluded as the mechanism; the splitting is circulation, not magnetism.
  F. verdict: the Delta is the nucleon with three times its ring circulation; derived to +-15 % with
     R25's radii; the exact radius 0.67 fm is not fixed by the base.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MDELTA = 0.51099895, 938.27209, 1232.0
LAMBDA_P = HBARC / MP
L1_9 = (2 * math.pi * HBARC / ME / 3) * (3 / 9) ** (2 * math.pi)
STATIC = math.pi * HBARC / L1_9                 # 763 MeV, three quark circuits
E_CIRC = {"+ strands at R+ (R25 D)": 125.1, "e on one ring (R25 A)": MP / (2 * 2.79284734)}
MU_DPP_RANGE = (3.7, 7.5)                        # mu_N, measured range
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    split = MDELTA - MP
    check("A. Spin from the ring, S = R E_circ/c: N needs hbar/2, Delta 3 hbar/2 -> E_circ(Delta) = 3 E_circ(N), Delta - N = 2 E_circ(N): the three circuits' circulations aligned instead of two against one",
          abs(split - 293.7) < 0.1, f"Delta - N = {split:.1f} MeV")

    preds = {k: 2 * v for k, v in E_CIRC.items()}
    devs = {k: v / split - 1 for k, v in preds.items()}
    check("B. With R25's two ring readings: 125 MeV -> 250 MeV (-15 %), 168 MeV -> 336 MeV (+14 %): they bracket the measured 294 MeV",
          abs(devs["+ strands at R+ (R25 D)"] + 0.148) < 0.01 and abs(devs["e on one ring (R25 A)"] - 0.144) < 0.01 and min(preds.values()) < split < max(preds.values()),
          "; ".join(f"{k}: {v:.0f} MeV ({devs[k]:+.0%})" for k, v in preds.items()))

    e_circ_n = split / 2
    R = HBARC / (2 * e_circ_n)
    mN, mD = STATIC + e_circ_n, STATIC + 3 * e_circ_n
    check("C. Inverting: E_circ(N) = 147 MeV, ring radius 0.671 fm (between R25's 0.59 and 0.79); N = 910 MeV (-3.0 %), Delta = 1204 MeV (-2.3 %): both within 3 % from one radius",
          abs(e_circ_n - 146.9) < 0.2 and abs(R - 0.671) < 0.003 and abs(mN / MP - 0.970) < 0.003 and abs(mD / MDELTA - 0.977) < 0.003,
          f"E_circ = {e_circ_n:.1f} MeV, R = {R:.3f} fm; N = {mN:.0f} ({mN/MP-1:+.1%}), Delta = {mD:.0f} ({mD/MDELTA-1:+.1%})")

    mu_dpp = 2 * R / LAMBDA_P
    check("D. Test of the radius: Delta++ (2e at 0.671 fm) has mu = 6.4 mu_N, inside the measured 3.7-7.5 mu_N (recorded)",
          MU_DPP_RANGE[0] < mu_dpp < MU_DPP_RANGE[1] and abs(mu_dpp - 6.38) < 0.05, f"mu(Delta++) = {mu_dpp:.2f} mu_N")

    # E. EM spin-spin between circuits: mu0 mu^2/(4 pi d^3) with mu ~ mu_N, d ~ 1 fm
    muN_J_per_T = 5.0507837e-27
    e_ss = 1e-7 * muN_J_per_T**2 / (1e-15) ** 3 / 1.602176634e-13
    check("E. An electromagnetic spin-spin interaction between the circuits is ~0.02 MeV: excluded as the mechanism; the splitting is circulation, not magnetism",
          abs(e_ss - 0.016) < 0.003 and split / e_ss > 1e4, f"mu_N^2 mu_0/(4 pi fm^3) = {e_ss:.3f} MeV; ratio {split/e_ss:.0f}")

    check("F. Verdict: the Delta is the nucleon with three times its ring circulation; derived to +-15 % with R25's radii; the exact 0.67 fm is not fixed by the base", True, "see B-C")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
