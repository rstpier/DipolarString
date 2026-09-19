#!/usr/bin/env python3
"""R64 -- "La duree de vie est la respiration du Z de la particule : elle devrait
etre adaptee (Z = Z0) et, si elle respire, elle se disperse."

Lecture testee : une particule est un circuit ferme de vortex (R47, R53) adapte a
Z0.  Si son impedance respire autour de Z0, une fraction de l'energie fuit a
chaque tour, comme a une jonction desadaptee : fuite/tour = |(Z - Z0)/(Z + Z0)|^2
= (dZ/2Z0)^2.  La duree de vie est alors tau = T / fuite, T = trajet / c.

On confronte cette lecture aux durees de vie mesurees : muon, tau, neutron,
pions, Delta.  Aucun parametre libre : les trajets sont ceux du redessin.
"""
import math

HBARC = 197.3269804          # MeV fm
C_FM = 2.99792458e23         # fm / s
HBAR_MEV_S = 6.582119569e-22 # MeV s
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME      # 386.16 fm
L1_3 = 2 * math.pi * LAMBDA_BAR / 3   # 808.8 fm

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

# --- donnees (PDG 2024) -----------------------------------------------------
MMU, TAU_MU = 105.6583755, 2.1969811e-6
MTAU, TAU_TAU, BR_TAU_E = 1776.86, 2.903e-13, 0.1782
MN_MP, TAU_N = 1.29333, 878.4
MPI_C, TAU_PI_C = 139.570, 2.6033e-8
MPI_0, TAU_PI_0 = 134.977, 8.43e-17
GAMMA_DELTA = 117.0          # MeV
TAU_DELTA = HBAR_MEV_S / GAMMA_DELTA

# --- trajets du redessin ------------------------------------------------------
PATH = {
    "mu":    3 * l1(7),                 # circuit a trois brins, R42
    "tau":   3 * l1(11),
    "n":     2 * math.pi * 0.95,        # boucle negative accrochee au proton, R42
    "pi+-":  2 * math.pi * (ALPHA * LAMBDA_BAR / 2),   # anneau a r_e/2
    "pi0":   2 * math.pi * (ALPHA * LAMBDA_BAR / 2),
    "Delta": 2 * math.pi * 0.671,       # anneau du nucleon, R49
}
TAU = {"mu": TAU_MU, "tau": TAU_TAU, "n": TAU_N, "pi+-": TAU_PI_C,
       "pi0": TAU_PI_0, "Delta": TAU_DELTA}

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R64 -- la duree de vie comme respiration de Z\n")

    # A. table : periode, nombre de tours, fuite par tour, desadaptation requise
    print("A. Fuite par tour requise (tau = T / fuite ; dZ/Z0 = 2 sqrt(fuite))")
    print(f"  {'objet':6s} {'trajet fm':>10s} {'T s':>10s} {'tours':>10s} "
          f"{'fuite/tour':>11s} {'dZ/Z0':>9s}")
    leak = {}
    for k in ("Delta", "pi0", "pi+-", "tau", "mu", "n"):
        T = PATH[k] / C_FM
        f = T / TAU[k]
        leak[k] = f
        dz = 2 * math.sqrt(f) if f < 1 else float("nan")
        print(f"  {k:6s} {PATH[k]:10.3f} {T:10.2e} {TAU[k]/T:10.2e} {f:11.2e} {dz:9.1e}")
    leak_tau_e = leak["tau"] * BR_TAU_E
    print(f"  tau (canal e seul, BR {BR_TAU_E}) : fuite/tour {leak_tau_e:.2e}, "
          f"dZ/Z0 {2*math.sqrt(leak_tau_e):.1e}")
    print(f"  -> 26 ordres de grandeur entre le Delta et le neutron pour une seule "
          f"'respiration'.\n")
    check("Delta : se defait en moins d'un tour (pas une fuite par respiration)",
          leak["Delta"] > 1, f"fuite/tour = {leak['Delta']:.2f} > 1")

    # B. respiration semblable a elle-meme : meme fuite/tour pour mu et tau
    print("B. Respiration semblable a elle-meme (tout suit l1(n)) : tau ~ 1/m")
    tau_tau_e_meas = TAU_TAU / BR_TAU_E
    tau_tau_e_pred = TAU_MU * MMU / MTAU
    ratio = tau_tau_e_pred / tau_tau_e_meas
    print(f"  tau(tau->e nu nu) predit {tau_tau_e_pred:.2e} s, mesure "
          f"{tau_tau_e_meas:.2e} s : rapport {ratio:.1e}\n")
    check("respiration auto-semblable EXCLUE (tau -> e faux de > 1e4)",
          ratio > 1e4, f"rapport {ratio:.1e}")

    # C. ce que les donnees imposent : loi en m^5
    print("C. Ce que mu et tau imposent")
    g_ratio = (1 / tau_tau_e_meas) / (1 / TAU_MU)
    m5 = (MTAU / MMU) ** 5
    print(f"  Gamma(tau->e)/Gamma(mu) = {g_ratio:.4e} ; (m_tau/m_mu)^5 = {m5:.4e} "
          f"({100*(g_ratio/m5-1):+.1f} %)")
    dz_ratio = math.sqrt(leak_tau_e / leak["mu"])
    m2 = (MTAU / MMU) ** 2
    print(f"  dZ/Z0 (tau, canal e) / dZ/Z0 (mu) = {dz_ratio:.1f} ; (m_tau/m_mu)^2 = {m2:.1f}")
    print(f"  -> avec T ~ l1 ~ 1/m, la fuite/tour doit croitre comme m^4, "
          f"la desadaptation comme m^2.\n")
    check("loi en m^5 entre mu et tau (a 1 %)", abs(g_ratio / m5 - 1) < 0.01,
          f"{100*(g_ratio/m5-1):+.2f} %")
    check("desadaptation dZ/Z0 proportionnelle a m^2 (a 2 %)",
          abs(dz_ratio / m2 - 1) < 0.02, f"{dz_ratio:.1f} vs {m2:.1f}")

    # D. l'echelle fixe que cela exige
    print("D. L'echelle fixe requise : fuite/tour = (m/M)^4")
    M = MMU / leak["mu"] ** 0.25
    L = HBARC / M
    w9 = (4 * LAMBDA_BAR / math.pi ** 2) * (l1(9) / L1_3)   # w_e a l'echelle 9, R17/R54
    smallest = {"w_9 (ruban du quark)": w9, "l1(9)": l1(9), "r_e": ALPHA * LAMBDA_BAR}
    print(f"  M = {M/1e6:.2f} TeV, L = hbar c / M = {L:.1e} fm")
    for name, val in smallest.items():
        print(f"  {name:22s} = {val:8.3f} fm  ({val/L:.0f} x L)")
    n_star = 3 * (L1_3 / L) ** (1 / (2 * math.pi))
    print(f"  sur l'echelle l1(n) il faudrait n = {n_star:.1f} (objet sans role)")
    me_a3 = ME / ALPHA ** 3
    print(f"  m_e/alpha^3 = {me_a3/1e6:.2f} TeV ({100*(me_a3/M-1):+.0f} %, prefacteur "
          f"inconnu : coincidence, pas un mecanisme)\n")
    check("la base n'a aucune longueur a l'echelle requise (plus petite > 100 L)",
          min(smallest.values()) / L > 100,
          f"L = {L:.1e} fm, plus petite longueur de base {min(smallest.values()):.3f} fm")

    # E. neutron : taille de la boucle ou energie liberee ?
    print("E. Le neutron avec la meme loi")
    gamma_mu = 1 / TAU_MU
    # (i) version 'taille' : fuite = (L / trajet)^4
    leak_size = (L / PATH["n"]) ** 4
    tau_n_size = (PATH["n"] / C_FM) / leak_size
    # (ii) version 'energie' : Gamma_n = Gamma_mu (Q/m_mu)^5
    tau_n_energy = 1 / (gamma_mu * (MN_MP / MMU) ** 5)
    print(f"  (i)  par la taille de la boucle (6 fm)     : tau_n = {tau_n_size:.1e} s "
          f"(mesure 878 s, rapport {TAU_N/tau_n_size:.0e})")
    print(f"  (ii) par l'energie liberee Q = {MN_MP} MeV : tau_n = {tau_n_energy:.0f} s "
          f"(rapport {tau_n_energy/TAU_N:.1f} ; le modele standard met ce facteur dans "
          f"g_A et l'espace des phases)\n")
    check("neutron : la version 'taille' est fausse de > 1e5",
          TAU_N / tau_n_size > 1e5, f"rapport {TAU_N/tau_n_size:.0e}")
    check("neutron : la version 'energie' tombe a un facteur < 20",
          tau_n_energy / TAU_N < 20, f"facteur {tau_n_energy/TAU_N:.1f}")

    # F. pions
    print("F. Les pions avec la meme loi (energie)")
    tau_pic_pred = 1 / (gamma_mu * (MPI_C / MMU) ** 5)
    tau_pi0_pred = 1 / (gamma_mu * (MPI_0 / MMU) ** 5)
    print(f"  pi+- : predit {tau_pic_pred:.1e} s, mesure {TAU_PI_C:.1e} "
          f"(facteur {tau_pic_pred/TAU_PI_C:.0f}) ; l'anneau du pion n'est pas un circuit "
          f"de lepton (R42)")
    print(f"  pi0  : predit {tau_pi0_pred:.1e} s, mesure {TAU_PI_0:.1e} "
          f"(facteur {tau_pi0_pred/TAU_PI_0:.0e}) ; fuite electromagnetique, autre "
          f"mecanisme\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
