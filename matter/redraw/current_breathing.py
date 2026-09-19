#!/usr/bin/env python3
"""R66 -- "Respiration plutot sur le courant qui circule, qui fait osciller
l'impedance intrinseque."

Lecture testee : l'impedance du circuit n'est pas une constante du milieu mais
depend du courant qui y circule, Z(I).  Sur l'echelle l1(n) le courant d'un
circuit est I = q c / (3 l1(n)), proportionnel a la masse.  L'ecart Z(I) - Z0
est la desadaptation, et la fuite par tour vaut (dZ/2Z0)^2 (R64, R65).

Questions : quelle puissance de I l'ecart doit-il suivre ?  quel courant fixe
cela demande-t-il ?  que devient l'electron ?  le neutron ?  tau -> e contre
tau -> mu ?
"""
import math

HBARC = 197.3269804
C_FM = 2.99792458e23
C_SI = 2.99792458e8
E_SI = 1.602176634e-19
HBAR_MEV_S = 6.582119569e-22
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
Q_VORTEX = 1 / (2 * math.sqrt(ALPHA))        # 5.85 e, quantum de circulation (R28)

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

def current_A(path_fm):
    """courant d'un quantum q circulant a c sur un trajet ferme (A)"""
    return Q_VORTEX * E_SI * C_SI / (path_fm * 1e-15)

MMU, TAU_MU = 105.6583755, 2.1969811e-6
MTAU, TAU_TAU, BR_TAU_E, BR_TAU_MU = 1776.86, 2.903e-13, 0.1782, 0.1739
MN_MP, TAU_N = 1.29333, 878.4
TAU_E_BOUND = 6.6e28 * 3.156e7

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R66 -- Z(I) : l'impedance suit le courant\n")
    path = {"e": 3 * L1_3, "mu": 3 * l1(7), "tau": 3 * l1(11), "q9": 3 * l1(9),
            "n_loop": 2 * math.pi * 0.95}
    I = {k: current_A(v) for k, v in path.items()}
    T = {k: v / C_FM for k, v in path.items()}
    print("  courants des circuits du redessin (q = 5,85 e a c) :")
    for k in ("e", "mu", "tau", "q9", "n_loop"):
        print(f"    {k:7s} trajet {path[k]:9.3f} fm  I = {I[k]:.2e} A")
    print()

    # A. exposant : dZ/Z0 = (I/I_c)^p ; mu et tau (canal e) fixent p
    leak_mu = T["mu"] / TAU_MU
    leak_tau_e = T["tau"] / (TAU_TAU / BR_TAU_E)
    dz_mu, dz_tau = 2 * math.sqrt(leak_mu), 2 * math.sqrt(leak_tau_e)
    p = math.log(dz_tau / dz_mu) / math.log(I["tau"] / I["mu"])
    print("A. Exposant de la reponse Z(I)")
    print(f"  dZ/Z0 : muon {dz_mu:.2e}, tau {dz_tau:.2e} ; I_tau/I_mu = {I['tau']/I['mu']:.2f}")
    print(f"  p = ln(dZ_tau/dZ_mu) / ln(I_tau/I_mu) = {p:.3f}")
    tau_lin = TAU_MU * (MMU / MTAU) ** 3
    print(f"  Z lineaire en I (p = 1) donnerait tau(tau->e) = {tau_lin:.1e} s "
          f"contre {TAU_TAU/BR_TAU_E:.1e} : faux de {tau_lin/(TAU_TAU/BR_TAU_E):.0f}")
    print("  -> Z(I) = Z0 [1 + (I/I_c)^2] : la premiere correction paire, celle qu'impose un "
          "milieu\n     sans sens de circulation privilegie (R3) ; fuite (I/I_c)^4/4, "
          "taux ~ m^5.\n")
    check("exposant p = 2 (a 2 %)", abs(p - 2) < 0.02, f"p = {p:.3f}")
    check("Z lineaire en I exclu (tau -> e faux de > 100)",
          tau_lin / (TAU_TAU / BR_TAU_E) > 100, f"{tau_lin/(TAU_TAU/BR_TAU_E):.0f}")

    # B. le courant critique
    I_c = I["mu"] / math.sqrt(dz_mu)
    M_equiv = MMU * I_c / I["mu"]
    print("B. Le courant fixe requis")
    print(f"  I_c = I_mu / sqrt(dZ_mu/Z0) = {I_c:.2e} A  (equivalent M = {M_equiv/1e6:.2f} TeV, R64)")
    for k in ("e", "q9", "n_loop"):
        print(f"    I_c / I_{k:6s} = {I_c/I[k]:.1e}")
    print("  -> aucun courant de la base n'approche I_c ; c'est l'echelle du TeV ecrite en amperes.\n")
    check("I_c est > 100 fois le plus grand courant de la base",
          I_c / max(I.values()) > 100, f"I_c/I_max = {I_c/max(I.values()):.0f}")

    # C. l'electron sous la meme regle
    leak_e = (I["e"] / I_c) ** 4 / 4
    tau_e_rule = T["e"] / leak_e
    print("C. L'electron avec la meme regle")
    print(f"  fuite/tour (I_e/I_c)^4/4 = {leak_e:.1e}, tau_e = {tau_e_rule:.1e} s = "
          f"{tau_e_rule/86400:.0f} jours (borne > {TAU_E_BOUND:.0e} s)")
    print("  -> la desadaptation seule ne suffit pas : il faut un etat plus bas ou tomber ;\n"
          "     l'electron est le barreau n = 3, il n'a pas de paire de DQD a lacher.\n")
    check("regle 'desadaptation => fuite' seule tuerait l'electron (< 1e7 s) : "
          "regle d'etat fondamental requise", tau_e_rule < 1e7, f"{tau_e_rule:.1e} s")

    # D. le neutron : courant du parent ou energie liberee ?
    leak_n_parent = (I["n_loop"] / I_c) ** 4 / 4
    tau_n_parent = T["n_loop"] / leak_n_parent
    tau_n_energy = TAU_MU * (MMU / MN_MP) ** 5
    print("D. Le neutron")
    print(f"  courant du parent (boucle a 0,95 fm, I = {I['n_loop']:.1e} A) : "
          f"tau_n = {tau_n_parent:.1e} s (mesure 878 : faux de {TAU_N/tau_n_parent:.0e})")
    print(f"  energie liberee Q = {MN_MP} MeV a la puissance 5 : tau_n = {tau_n_energy:.0f} s "
          f"(facteur {tau_n_energy/TAU_N:.1f})")
    print("  -> le courant qui compte est celui de ce qui part (Q), pas celui du parent.\n")
    check("neutron : le courant du parent est faux de > 1e5",
          TAU_N / tau_n_parent > 1e5, f"{TAU_N/tau_n_parent:.0e}")

    # E. tau -> mu contre tau -> e : lacher une paire ou deux paires de DQD ?
    r_meas = BR_TAU_MU / BR_TAU_E
    r_energy5 = ((MTAU - MMU) / (MTAU - ME)) ** 5
    x = (MMU / MTAU) ** 2
    r_phase = 1 - 8 * x + 8 * x ** 3 - x ** 4 - 12 * x ** 2 * math.log(x)
    print("E. tau -> mu (une paire de DQD lachee, 11 -> 7) contre tau -> e (deux paires, 11 -> 3)")
    print(f"  mesure Gamma(tau->mu)/Gamma(tau->e) = {r_meas:.3f}")
    print(f"  lacher paire par paire : e << mu (deux fuites successives) : exclu")
    print(f"  energie liberee a la puissance 5 : {r_energy5:.3f} ({100*(r_energy5/r_meas-1):+.0f} %)")
    print(f"  courant du parent seul, corrige de l'espace des phases a trois corps : "
          f"{r_phase:.4f} ({100*(r_phase/r_meas-1):+.1f} %)")
    print("  -> un seul effondrement vers n'importe quel barreau inferieur, au taux du parent ;\n"
          "     le partage de l'energie entre les trois corps (espace des phases) manque a la base.\n")
    check("tau -> mu / tau -> e : paire par paire exclu, parent + espace des phases a 1 %",
          abs(r_phase / r_meas - 1) < 0.01 and abs(r_energy5 / r_meas - 1) > 0.1,
          f"phases {r_phase:.4f}, energie^5 {r_energy5:.3f}, mesure {r_meas:.3f}")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
