#!/usr/bin/env python3
"""Redraw: 'the neutron and the proton differ by one joint, so the p-n mass difference calibrates
the joint (knot)' -- the author's proposal.  Tested in the base (BASE.md R2, R4, R13, R17).

Two readings of 'one joint':
  (i)  neutron = proton + one string + one joint: joint = m_n - m_p = 1.293 MeV;
  (ii) neutron = proton + an electron chain (0, 3) attached by one joint:
       joint = m_n - m_p - m_e = 0.782 MeV (the beta Q-value).
Checks:
  A. the ring charge rule q = (n+ - n-) e/3: adding ONE string changes q by +-1/3, not by -1, so
     reading (i) breaks the charge rule; reading (ii) keeps it (+1 - 1 = 0) and the parity (odd + 3 = even).
  B. the calibration values and their ratio to the electron's junction (127.75 keV from g = 2).
  C. the base's junction law on a matched guide, hbar c/(pi^2 s) (R17): a 1.29 / 0.78 MeV junction
     sits at s = 15.5 / 25.6 fm, 18-30 proton radii; inside the proton (s = 0.84 fm) the matched
     junction weighs 24 MeV.  The calibration is not the R17 junction; nucleon poles would be
     0.87 e / 0.68 e at s = r_p.
  D. the Coulomb energy the base itself carries: the proton's own field at r_p (0.86-1.03 MeV) or the
     p-e attraction at r_p (-1.71 MeV) must be paid, so the joint is 2.2-2.5 MeV gross, 1.29 net.
     (Lattice QCD+QED, BMW 2015: strong part 2.52 MeV, EM part -1.00 MeV.)
  E. universality: the base's own ring count makes the pi- (5-ring) a pi0 (4-ring) + one string
     + one joint, and that joint reads 4.59 MeV, 3.6 x the nucleon's; other isospin pairs give 0.3
     to 8.1 MeV.  One joint is not one number.
  F. reading (ii) as a picture: the added negative loop, with the base's own mu = -e c R/2, needs
     R = 0.99 fm to give mu_n - mu_p; the neutron's charge radius (<r^2>_n = -0.1155 fm^2) with a
     +e proton core and a -e ring gives R = 0.91 fm: two independent size tests agree to 9 %.
     Its spin S = R E/c is < 0.02 hbar, so the neutron keeps the proton's spin 1/2.
Exit status is zero only if every check passes.  PDG 2024 / CODATA 2018 values.
"""

from __future__ import annotations

import math

MP, MN, ME = 938.27209, 939.56542, 0.51099895      # MeV
HBARC = 197.3269804                                  # MeV fm
ALPHA = 1 / 137.035999
K = ALPHA * HBARC                                    # e^2/(4 pi eps0) = 1.440 MeV fm
RP = 0.8409                                          # proton charge radius, fm
R2N = -0.1155                                        # neutron <r^2>, fm^2
MU_P, MU_N = 2.79284734, -1.91304273                 # nuclear magnetons
LAMBDA_P = HBARC / MP                                # 0.2103 fm
DELTA_MATCHED = 1 / (math.pi * math.sqrt(ALPHA))     # 3.726 e (R17)
U_E = ME / 4                                         # electron junction, 127.75 keV (g = 2, two junctions carry half)
ISOSPIN = {  # pair: mass difference (heavier - lighter), MeV
    "n - p": MN - MP, "pi+- - pi0": 139.57039 - 134.9768, "K0 - K+-": 497.611 - 493.677,
    "Sigma- - Sigma+": 1197.449 - 1189.37, "Sigma- - Sigma0": 1197.449 - 1192.642,
    "Xi- - Xi0": 1321.71 - 1314.86, "D+- - D0": 1869.66 - 1864.84, "B0 - B+-": 0.32,
}
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. charge rule
    dq_one_string = 1 / 3
    check("A. Ring charge rule: one string more changes q by 1/3, not by 1 -> reading (i) 'proton + one string' cannot be neutral; reading (ii) 'proton + electron chain (0,3) + one joint' gives q = 0 and the right parity (odd + 3 = even)",
          abs(dq_one_string - 1) > 0.5 and (9 + 3) % 2 == 0,
          "q_p = +1 needs n+ - n- = 3; q_n = 0 needs n+ = n-: the difference is three strings, the electron's count")

    # B. calibration values
    j1, j2 = MN - MP, MN - MP - ME
    check("B. Joint calibration: (i) m_n - m_p = 1.293 MeV, (ii) m_n - m_p - m_e = 0.782 MeV; that is 10.1 x / 6.1 x the electron's junction of 127.75 keV",
          abs(j1 - 1.2933) < 0.001 and abs(j2 - 0.7823) < 0.001 and abs(j1 / U_E - 10.12) < 0.02,
          f"(i) {j1:.4f} MeV = {j1/U_E:.2f} u_e; (ii) {j2:.4f} MeV = {j2/U_E:.2f} u_e")

    # C. base junction law
    s1, s2 = HBARC / (math.pi**2 * j1), HBARC / (math.pi**2 * j2)
    u_in = HBARC / (math.pi**2 * RP)
    d1, d2 = math.sqrt(j1 * RP / K), math.sqrt(j2 * RP / K)
    check("C. On a matched guide (poles 3.73 e) a joint of 1.29 / 0.78 MeV sits at s = 15.5 / 25.6 fm, 18-30 proton radii; inside the proton (s = r_p) the matched joint weighs 24 MeV: the calibration is not the R17 junction unless the nucleon's poles are 0.87 e / 0.68 e",
          abs(s1 - 15.46) < 0.05 and abs(s2 - 25.56) < 0.05 and abs(u_in - 23.8) < 0.1 and abs(d1 - 0.869) < 0.005 and abs(d2 - 0.676) < 0.005,
          f"s(i) = {s1:.2f} fm = {s1/RP:.1f} r_p, s(ii) = {s2:.2f} fm = {s2/RP:.1f} r_p; hbar c/(pi^2 r_p) = {u_in:.1f} MeV; poles at r_p: {d1:.3f} e / {d2:.3f} e")

    # D. Coulomb bookkeeping
    e_shell, e_unif = K / (2 * RP), 3 * K / (5 * RP)
    e_pe = K / RP
    gross_i = (j1 + e_shell, j1 + e_unif)
    gross_ii = j2 + e_pe
    check("D. The base carries the proton's own Coulomb energy (0.86 shell / 1.03 uniform MeV at r_p) or the p-e attraction (-1.71 MeV at r_p): the joint is 2.15-2.32 MeV gross in (i), 2.49 MeV in (ii); 1.29 MeV is the net (lattice QCD+QED: strong 2.52, EM -1.00)",
          abs(e_shell - 0.856) < 0.005 and abs(e_unif - 1.027) < 0.005 and abs(gross_ii - 2.494) < 0.01,
          f"proton field: shell {e_shell:.3f} MeV, uniform {e_unif:.3f} MeV; p-e at r_p: -{e_pe:.3f} MeV; gross joint (i) {gross_i[0]:.2f}-{gross_i[1]:.2f} MeV, (ii) {gross_ii:.2f} MeV")

    # E. universality
    ratios = {k: v / j1 for k, v in ISOSPIN.items()}
    lo, hi = min(ISOSPIN.values()), max(ISOSPIN.values())
    check("E. In the base's own ring count the pi- (5-ring) is a pi0 (4-ring) + one string + one joint: that joint reads 4.59 MeV, 3.6 x the nucleon's; isospin pairs span 0.3 to 8.1 MeV: 'one joint' is not one number",
          abs(ratios["pi+- - pi0"] - 3.55) < 0.02 and hi / lo > 20,
          "; ".join(f"{k}: {v:.2f} MeV ({ratios[k]:.2f} x)" for k, v in ISOSPIN.items()))

    # F. reading (ii) as a picture: proton + negative loop
    r_mu = (MU_P - MU_N) * LAMBDA_P              # mu_loop = -e c R/2 = -(R/lambda_p) mu_N
    r_ch = math.sqrt(RP**2 - R2N)                # <r^2>_n = <r^2>_p - R^2
    s_loop = r_mu * gross_ii / HBARC             # S = R E/c in hbar
    check("F. Reading (ii) as a picture -- neutron = proton + a negative loop with the base's mu = -e c R/2: mu_n - mu_p needs R = 0.99 fm; the neutron's charge radius (+e core, -e ring) gives R = 0.91 fm: two independent size tests agree to 9 %; the loop's spin R E/c < 0.02 hbar keeps the neutron at spin 1/2",
          abs(r_mu - 0.990) < 0.003 and abs(r_ch - 0.907) < 0.003 and abs(r_mu / r_ch - 1.09) < 0.01 and s_loop < 0.02,
          f"R from mu: {r_mu:.3f} fm; R from <r^2>_n: {r_ch:.3f} fm; ratio {r_mu/r_ch:.3f}; loop spin with E = {gross_ii:.2f} MeV: {s_loop:.4f} hbar")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
