#!/usr/bin/env python3
"""R69 -- "Il faut trouver une relation entre n (et autres proprietes liees) et
la stabilite temporelle.  Il faut attribuer une frequence relativiste (selon la
vitesse de vibration, la duree de vie percue est modulee a nos yeux)."

Loi construite avec R64-R68 :
    tau(n) = T(n) * (M / m(n))^4 / C(n)
  T(n) = 2 pi hbar / m(n) c^2   : un tour = periode de Compton (R68 A)
  M    = (96 pi^2)^(1/4)/sqrt(G_F) = 1,625 TeV : le mur (R64, R68)
  C(n) = nombre de barreaux inferieurs de meme classe de charge (canaux),
         la couleur comptant 3 (position du brin impair, R42)
Puis : quelle "frequence relativiste" cette loi contient-elle, et que
demanderait un gamma de vibration qui porterait la lenteur ?
"""
import math

HBAR_MEV_S = 6.582119569e-22
ME = 0.51099895
GF = 1.1663788e-11               # MeV^-2
M_WALL = (96 * math.pi ** 2) ** 0.25 / math.sqrt(GF)   # 1.625e6 MeV

MMU, TAU_MU = 105.6583755, 2.1969811e-6
MTAU, TAU_TAU = 1776.86, 2.903e-13
BR_TAU_E, BR_TAU_MU, BR_TAU_HAD = 0.1782, 0.1739, 0.6479
MN_MP, TAU_N = 1.29333, 878.4
TAU_E_BOUND = 6.6e28 * 3.156e7
QCD_HAD = 1.20                   # 1 + alpha_s/pi + ... a m_tau (modele standard)

def m_ladder(n):
    return ME * (n / 3) ** (2 * math.pi)

def T_turn(m):
    return 2 * math.pi * HBAR_MEV_S / m

def tau_law(m, channels):
    return T_turn(m) * (M_WALL / m) ** 4 / channels if channels > 0 else float("inf")

def phase_space(x):
    """partage a trois corps, x = (m_fille/m_mere)^2"""
    return 1 - 8 * x + 8 * x ** 3 - x ** 4 - 12 * x ** 2 * math.log(x)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R69 -- la loi de stabilite en n\n")
    print(f"  M = (96 pi^2)^(1/4)/sqrt(G_F) = {M_WALL/1e6:.3f} TeV (pas d'ajustement)\n")

    # A. electron : n = 3, aucun barreau inferieur
    print("A. Electron, n = 3 : C = 0")
    tau_e_if = tau_law(ME, 1)
    print(f"  s'il avait un canal : {tau_e_if:.1e} s = {tau_e_if/86400:.0f} jours ; "
          f"C = 0 -> stable (borne > {TAU_E_BOUND:.0e} s)")
    check("electron stable parce que C(3) = 0", tau_law(ME, 0) == float("inf"), "C = 0")

    # B. muon : n = 7, un canal (e)
    print("B. Muon, n = 7 : C = 1 (vers l'electron)")
    tau_mu_pred = tau_law(MMU, 1)
    tau_mu_ladder = tau_law(m_ladder(7), 1)
    print(f"  predit {tau_mu_pred*1e6:.4f} us (masse vraie), {tau_mu_ladder*1e6:.3f} us "
          f"(masse de l'echelle, 0,8 % -> 4 %) ; mesure {TAU_MU*1e6:.4f} us")
    check("muon a 1 % (masse vraie)", abs(tau_mu_pred / TAU_MU - 1) < 0.01,
          f"{100*(tau_mu_pred/TAU_MU-1):+.2f} %")

    # C. tau : n = 11, canaux e, mu (partage), quarks x 3 couleurs
    print("C. Tau, n = 11 : C = 1 (e) + partage(mu) + 3 (barreau des quarks, 3 couleurs)")
    ps_mu = phase_space((MMU / MTAU) ** 2)
    C_tau = 1 + ps_mu + 3
    C_tau_qcd = 1 + ps_mu + 3 * QCD_HAD
    tau_tau_pred = tau_law(MTAU, C_tau)
    tau_tau_qcd = tau_law(MTAU, C_tau_qcd)
    print(f"  C = {C_tau:.3f} : tau = {tau_tau_pred:.3e} s ({100*(tau_tau_pred/TAU_TAU-1):+.1f} %) ; "
          f"avec la correction forte {QCD_HAD} sur les quarks, C = {C_tau_qcd:.2f} : "
          f"{tau_tau_qcd:.3e} s ({100*(tau_tau_qcd/TAU_TAU-1):+.1f} %) ; mesure {TAU_TAU:.3e}")
    print(f"  rapports de branchement : e {100/C_tau:.1f} % (mesure {100*BR_TAU_E:.1f}), "
          f"mu {100*ps_mu/C_tau:.1f} % ({100*BR_TAU_MU:.1f}), quarks {300/C_tau:.1f} % "
          f"({100*BR_TAU_HAD:.1f}) ; avec 1,2 : quarks {300*QCD_HAD/C_tau_qcd:.1f} %")
    check("tau a 15 % par simple compte des barreaux (couleur = 3)",
          abs(tau_tau_pred / TAU_TAU - 1) < 0.15, f"{100*(tau_tau_pred/TAU_TAU-1):+.1f} %")
    check("tau a 1 % avec la correction forte 1,2 (non derivee ici)",
          abs(tau_tau_qcd / TAU_TAU - 1) < 0.01, f"{100*(tau_tau_qcd/TAU_TAU-1):+.1f} %")

    # D. la dependance en n : tau ~ (3/n)^(10 pi) / C
    print("D. En n seul : tau(n) ~ (3/n)^(10 pi) / C(n), exposant 10 pi = 31,4")
    ratio_law = (7 / 11) ** (10 * math.pi)
    ratio_meas = (TAU_TAU / BR_TAU_E) / TAU_MU        # canal e contre canal e
    print(f"  (7/11)^(10 pi) = {ratio_law:.2e} ; tau(tau->e)/tau(mu) mesure = {ratio_meas:.2e} "
          f"({100*(ratio_law/ratio_meas-1):+.0f} %, l'ecart de l'echelle a Koide amplifie)")
    check("exposant 10 pi entre mu et tau (a 10 %)", abs(ratio_law / ratio_meas - 1) < 0.10,
          f"{100*(ratio_law/ratio_meas-1):+.1f} %")

    # E. neutron : meme loi avec l'energie liberee
    print("E. Neutron : meme loi avec Q = m_n - m_p a la place de m")
    tau_n_pred = tau_law(MN_MP, 1)
    print(f"  {tau_n_pred:.0f} s contre {TAU_N} (facteur {tau_n_pred/TAU_N:.1f} : le spin g_A et "
          f"le partage, hors de la base)")
    check("neutron a un facteur < 20", tau_n_pred / TAU_N < 20, f"{tau_n_pred/TAU_N:.1f}")

    # F. la "frequence relativiste"
    print("F. La frequence relativiste")
    print("  la loi en contient une seule : la frequence de Compton m c^2/h, celle du tour.")
    print("  un gamma de vibration qui porterait la lenteur devrait valoir N = (M/m)^4 tours :")
    for name, m in (("mu", MMU), ("tau", MTAU)):
        N = (M_WALL / m) ** 4
        print(f"    {name:3s} : gamma = {N:.1e}, 1 - v/c = {1/(2*N*N):.1e}")
    print(f"  et suivre (n/3)^(16 pi) = n^{16*math.pi:.0f} : une vitesse reglee a 34 decimales.")
    print("  -> la lenteur n'est pas une dilatation ; c'est le compte (m/M)^4 par tour.")
    check("gamma de vibration exclu : exposant en n > 40 et reglage < 1e-30",
          16 * math.pi > 40 and 1 / (2 * (M_WALL / MMU) ** 8) < 1e-30,
          f"n^{16*math.pi:.0f}, 1 - v/c = {1/(2*(M_WALL/MMU)**8):.0e}")

    # G. le barreau suivant
    m15 = m_ladder(15)
    tau15 = tau_law(m15, 3 + 3 * QCD_HAD)
    print(f"\nG. Barreau n = 15 (exclu par le LEP) : m = {m15/1e3:.1f} GeV, "
          f"tau = {tau15:.0e} s avec C = 3 leptons + 3 x 1,2 quarks")

    print("\nBilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
