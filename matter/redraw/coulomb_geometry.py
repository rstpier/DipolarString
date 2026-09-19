#!/usr/bin/env python3
"""Redraw: what fixes the geometry of the proton's Coulomb term (R42 left a +11 to +25 % residual on
m_n - m_p)?  In the final topology the quark charges q_i = (2/3, 2/3, -1/3) for p and (2/3, -1/3, -1/3)
for n give E_C = S sum q_i^2 + M sum_{i<j} q_i q_j with S the self term (one unit charge in the
quark's geometry) and M the mutual term (two unit charges at the quark spacing).

Checks:
  A. what is needed: E_C(p) - E_C(n) = 2.474 - 1.293 = 1.181 MeV; sum q^2 = 1 (p), 2/3 (n);
     sum q_i q_j = 0 (p), -1/3 (n): (S + M)/3 = 1.181 -> S + M = 3.54 MeV.
  B. the mutual term is parameter-free in the base: charge at the poles (R5), arms at 120 deg with
     ends at l_1(9) from the centre, spacing sqrt 3 l_1 = 1.41 fm: M = K/(sqrt 3 l_1) = 1.02 MeV; the
     proton's inter-quark Coulomb vanishes exactly, the neutron's is -0.34 MeV.
  C. the self term needs the charge's size: spread along the strand (rod of length l_1(9) = 0.81 fm)
     with the sheet width scaled from the electron's, w_9 = (4 lambda-bar_e/pi^2) 3^(-2 pi) = 0.157 fm
     (R17), as thickness: S = (K/l_1)(ln(2 l_1/w) - 1) = 2.37 MeV -> S + M = 3.39 MeV ->
     m_n - m_p = 1.34 MeV (+3.8 %).  The log constant and the spacing move it between 1.09 and
     1.34 MeV (-16 to +4 %): the base fixes the geometry to about +-10 %.
  D. the same geometry gives the proton's charge radius: charges at l_1(9) from the centre,
     sqrt<r^2>_p = l_1(9) = 0.813 fm (-3.4 % vs 0.841); the neutron's <r^2> = 0 vs -0.116 fm^2 (its
     negative skin needs the outer loop of R19/R25).
  E. alternatives without the base's thickness: point poles need a = 0.29-0.34 fm, touching spheres
     a = 0.45 fm (d = 0.89 fm), uniform sphere at the measured r_p (+12 %): fits, not the base.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
L1 = lambda n: 2 * math.pi * LAMBDA_E / 3 * (3 / n) ** (2 * math.pi)
L9 = L1(9)
W9 = (4 * LAMBDA_E / math.pi**2) * 3 ** (-2 * math.pi)
E_DQD = 4 * math.pi * K / (9 * L9)          # 2.474 MeV, R42
DM_OBS = MN - MP
RP, R2N = 0.8409, -0.1155
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


Q_P, Q_N = (2 / 3, 2 / 3, -1 / 3), (2 / 3, -1 / 3, -1 / 3)


def sums(q):
    return sum(x * x for x in q), sum(q[i] * q[j] for i in range(3) for j in range(i + 1, 3))


def main() -> int:
    # A. needed
    need = E_DQD - DM_OBS
    s2p, mp_ = sums(Q_P)
    s2n, mn_ = sums(Q_N)
    sm_needed = 3 * need
    check("A. Needed: E_C(p) - E_C(n) = 2.474 - 1.293 = 1.181 MeV; sum q^2 = 1 (p), 2/3 (n); sum q_i q_j = 0 (p), -1/3 (n): S + M = 3.54 MeV",
          abs(need - 1.181) < 0.005 and abs(s2p - 1) < 1e-12 and abs(s2n - 2 / 3) < 1e-12 and abs(mp_) < 1e-12 and abs(mn_ + 1 / 3) < 1e-12 and abs(sm_needed - 3.54) < 0.02,
          f"need {need:.3f} MeV; S + M = {sm_needed:.2f} MeV")

    # B. mutual term
    d = math.sqrt(3) * L9
    M = K / d
    check("B. Mutual term, parameter-free: charge at the poles, arms at 120 deg, spacing sqrt 3 l_1(9) = 1.41 fm: M = 1.02 MeV; the proton's inter-quark Coulomb vanishes exactly, the neutron's is -0.34 MeV",
          abs(d - 1.408) < 0.005 and abs(M - 1.023) < 0.005 and abs(mn_ * M + 0.341) < 0.003,
          f"d = {d:.3f} fm, M = {M:.3f} MeV; p: {mp_*M:.3f}, n: {mn_*M:.3f} MeV")

    # C. self term with the base's thickness
    S = (K / L9) * (math.log(2 * L9 / W9) - 1)
    dm = E_DQD - (S + M) / 3
    variants = {
        "ln(2l/w)-1, d = sqrt3 l": E_DQD - (S + M) / 3,
        "ln(l/w), d = sqrt3 l": E_DQD - ((K / L9) * math.log(L9 / W9) + M) / 3,
        "ln(2l/w)-1, d = l": E_DQD - (S + K / L9) / 3,
    }
    check("C. Self term with the charge along the strand (rod 0.81 fm) and the sheet width scaled from the electron (w_9 = 0.157 fm, R17): S = 2.37 MeV, S + M = 3.39, m_n - m_p = 1.34 MeV (+3.8 %); the log constant and the spacing span 1.09-1.34 MeV (-16 to +4 %): the base fixes the geometry to about +-10 %",
          abs(W9 - 0.1572) < 0.001 and abs(S - 2.37) < 0.02 and abs(dm / DM_OBS - 1.038) < 0.01 and min(variants.values()) > 1.05 and max(variants.values()) < 1.36,
          f"w_9 = {W9:.4f} fm; S = {S:.3f} MeV; " + "; ".join(f"{k}: {v:.3f} ({v/DM_OBS-1:+.0%})" for k, v in variants.items()))

    # D. charge radii from the same geometry
    r2p = sum(q * L9**2 for q in Q_P)
    r2n = sum(q * L9**2 for q in Q_N)
    check("D. Same geometry, charge radii: sqrt<r^2>_p = l_1(9) = 0.813 fm (-3.4 % vs 0.841); <r^2>_n = 0 vs -0.116 fm^2 (the negative skin needs the outer loop of R19/R25)",
          abs(math.sqrt(r2p) / RP - 0.966) < 0.005 and abs(r2n) < 1e-12,
          f"sqrt<r^2>_p = {math.sqrt(r2p):.3f} fm ({math.sqrt(r2p)/RP-1:+.1%}); <r^2>_n = {r2n:.3f} fm^2 vs {R2N}")

    # E. alternatives
    a_shell = K / (2 * (sm_needed - M))
    a_sphere = 3 * K / (5 * (sm_needed - M))
    a_touch = (3 * K / 5 + K / 2) / sm_needed
    e_unif = 3 * K / (5 * RP)
    check("E. Without the base's thickness: point poles need a = 0.29 (shell) - 0.34 (sphere) fm, touching spheres a = 0.45 fm (d = 0.89 fm), uniform sphere at the measured r_p gives +12 %: fits, not the base",
          abs(a_shell - 0.286) < 0.005 and abs(a_sphere - 0.343) < 0.005 and abs(a_touch - 0.447) < 0.005 and abs((E_DQD - e_unif) / DM_OBS - 1.12) < 0.01,
          f"a_shell = {a_shell:.3f}, a_sphere = {a_sphere:.3f}, a_touch = {a_touch:.3f} fm (d = {2*a_touch:.2f}); sphere at r_p: {E_DQD-e_unif:.3f} MeV ({(E_DQD-e_unif)/DM_OBS-1:+.0%})")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
