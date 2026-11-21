#!/usr/bin/env python3
"""R123 -- Adoption de l'axiome d'amplitudes (A3 de la note du 9 octobre ; A9 dans la numerotation du manuscrit) :
les consequences immediates et ce que les donnees en disent.

Decision de l'auteur (21 novembre 2026) : « On adopte A3 ». Sous A3, l'etat du milieu est une amplitude complexe sur
ses configurations : superposition, intrication et regle de Born sont dans le cadre, et les particules sont les
quanta des champs du milieu. Le programme change de question : non plus « reconstruire les particules », mais « que
le milieu ajoute-t-il a la theorie quantique des champs ? ». Trois choses, chacune avec sa borne :

  A. Une coupure : le photon est le quantum de l'onde collective du reseau, sa dispersion est celle du reseau. Trame
     cubique a pas l, premier voisin : omega = (2c/l) sin(k l / 2) = c k (1 - k^2 l^2 / 24 + ...), soit
     E^2 = p^2 c^2 (1 - (E/E_QG2)^2) avec E_QG2 = sqrt(12) hbar c / l. Pour l = l_1 = 808,7 fm : E_QG2 = 0,85 MeV (le
     manuscrit donne 0,69 MeV par sa convention). LHAASO, GRB 221009A (fevrier 2026) : E_QG2 > 6e-8 E_Pl = 7,3e11 GeV.
     Ecart 1e15 ; le pas du reseau porteur des photons doit etre l < sqrt(12) hbar c / E_QG2 = 9e-28 m = 1e-15 l_1.
     La sortie « photon = soliton N = 1 hors du reseau » n'existe plus sous A3. l_1 ne peut pas etre le pas du milieu.
  B. Un spineur : l'electron est le quantum d'un champ de spin 1/2 du milieu. Sur la trame cubique de la base (noeud
     SCN de Johns), toute excitation de spin 1/2 verrouillee en helicite va a 1/3 de la vitesse de lien contre 1/2 pour
     le photon (R106-R107, identite universelle R110) : v_e / c = 2/3, exclu a 1e-11 (LEP). Le reseau qui porte Weyl et
     Maxwell a une vitesse (automate de Weyl, R111) exige un noeud non reciproque de symetrie T_h avec un quart d'onde
     physique (R112), hors des contraintes C1-C5 et contre A5 : adopter A3 oblige a remplacer le noeud et A5.
  C. Une energie de point zero : un reseau quantifie a pas l porte rho_zp = sum_k hbar omega_k / 2 ~ pi^2 hbar c / (8 l^4)
     (deux polarisations, dispersion lineaire jusqu'au bord de zone). Pour l_1 : 9e22 J/m^3, contre rho_Lambda =
     5,3e-10 J/m^3 : 1,7e32. L'axiome A8 du manuscrit (energie native non gravitante) devient obligatoire, et c'est
     le probleme de la constante cosmologique sous un autre nom (note du 9 octobre : 1e55-1e120 selon la coupure).
  Ce qu'A3 donne sans le milieu : superposition, intrication sans canal (CHSH 2 sqrt 2 > 2 ; conforme aux notes des
  4-5 octobre et du 1er novembre), regle de Born (dans l'axiome, ou Gleason). Ce que le milieu doit encore fournir :
  un repere privilegie dont les parametres alpha_1, alpha_2 ne sont pas calcules (probleme ouvert n° 6 du manuscrit,
  desormais le premier calcul du programme) ; un champ de spin 1/2 (B) ; la raison pour laquelle le point zero ne
  gravite pas (C).
  Verdict : A3 ADOPTE (decision). Le milieu quantifie tel qu'ecrit -- trame cubique a l_1 -- est EXCLU comme support
  des photons (A) et de l'electron (B) ; l_1 survit comme longueur de Compton de l'electron (2 pi lambda_bar_e / 3,
  propriete du quantum), pas comme pas du milieu. Ce qui reste de la CD sous A3 : ses rapports (secteur hadronique,
  inchange), ses etats limites gravitationnels (inchanges), et un programme : un milieu dont le pas est < 1e-27 m,
  un noeud qui porte un spineur, et des alpha_1, alpha_2 a calculer.
"""
import math

HBARC_MEVFM = 197.3269804
HBARC_J_M = 3.16152677e-26
LAM_E = HBARC_MEVFM / 0.51099895          # fm
L1 = 2 * math.pi * LAM_E / 3              # fm
E_PL_GEV = 1.220890e19
E_QG2_LHAASO = 6e-8 * E_PL_GEV            # GeV
H0 = 67.4e3 / 3.0857e22                   # s^-1
G, C = 6.6743e-11, 2.99792458e8
OMEGA_L = 0.69

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R123 -- adoption de l'axiome d'amplitudes : consequences immediates\n")

    # A. coupure et LIV
    print("A. La coupure du reseau contre LHAASO")
    ks = [0.05, 0.1, 0.2]
    for x in ks:  # x = k l
        w = 2 * math.sin(x / 2); lin = x * (1 - x * x / 24)
        print(f"      k l = {x:.2f} : omega l / c = {w:.6f} ; c k (1 - k^2 l^2/24) = {lin:.6f} ; ecart {abs(w/lin-1):.1e}")
    e_qg2_l1 = math.sqrt(12) * HBARC_MEVFM / L1      # MeV
    ratio = E_QG2_LHAASO * 1e3 / e_qg2_l1
    l_max = math.sqrt(12) * HBARC_J_M / (E_QG2_LHAASO * 1.602176634e-10)   # m
    print(f"    E_QG2(l_1) = sqrt(12) hbar c / l_1 = {e_qg2_l1:.3f} MeV (manuscrit : 0,69 MeV) ; LHAASO : E_QG2 > {E_QG2_LHAASO:.2e} GeV ; ecart {ratio:.1e}")
    print(f"    pas maximal d'un reseau porteur des photons : l < {l_max:.1e} m = {l_max/(L1*1e-15):.1e} l_1 (l_Planck = 1,6e-35 m)\n")
    check("A. sous A3 le photon a la dispersion du reseau : l_1 exclu comme pas par LHAASO (ecart > 1e14), l < 1e-27 m",
          ratio > 1e14 and l_max < 1e-27, f"ecart {ratio:.1e}, l < {l_max:.1e} m")

    # B. le spineur
    print("B. L'electron comme quantum d'un champ du milieu")
    v_ratio = (1 / 3) / (1 / 2)
    print(f"    trame cubique SCN : v(spin 1/2) / v(photon) = (1/3)/(1/2) = {v_ratio:.4f} ; LEP : |v_e/c - 1| < 1e-11 ; ecart {abs(1-v_ratio):.2f} >> 1e-11 (R106-R107, R110)")
    print("    automate de Weyl (R111) : une vitesse pour Weyl et Maxwell, noeud T_h non reciproque + quart d'onde physique (R112) : hors C1-C5, contre A5\n")
    check("B. la trame de la base ne porte pas l'electron comme quantum (2/3 c) ; le noeud et A5 doivent etre remplaces", abs(1 - v_ratio) > 1e-3, f"{v_ratio:.4f}")

    # C. point zero
    print("C. L'energie de point zero du reseau quantifie")
    rho_zp = math.pi ** 2 * HBARC_J_M / (8 * (L1 * 1e-15) ** 4)
    rho_crit = 3 * H0 ** 2 / (8 * math.pi * G) * C ** 2
    rho_L = OMEGA_L * rho_crit
    print(f"    rho_zp(l_1) = pi^2 hbar c / (8 l_1^4) = {rho_zp:.2e} J/m^3 ; rho_Lambda = {rho_L:.2e} J/m^3 ; rapport {rho_zp/rho_L:.1e} -> A8 obligatoire\n")
    check("C. point zero du reseau a l_1 = 1e32 rho_Lambda : A8 (energie native non gravitante) devient obligatoire", rho_zp / rho_L > 1e30, f"{rho_zp/rho_L:.1e}")

    # D. ce qu'A3 donne et ce qui reste
    print("D. Ce qu'A3 donne et ce que le milieu doit encore fournir")
    chsh = 2 * math.sqrt(2)
    print(f"    superposition -> CHSH max {chsh:.3f} > 2 sans canal (conforme aux notes d'octobre-novembre) ; regle de Born dans l'axiome")
    print("    a fournir : alpha_1, alpha_2 du repere du milieu (premier calcul) ; un champ de spin 1/2 (B) ; pourquoi le point zero ne gravite pas (C)\n")
    check("D. A3 fournit Bell sans canal (2 sqrt 2 > 2) ; le milieu reste redevable de alpha_1, alpha_2, du spineur et du point zero", chsh > 2, f"{chsh:.3f}")

    print("Verdict : A3 ADOPTE. Le milieu quantifie tel qu'ecrit (trame cubique a l_1) est EXCLU comme support des photons (LHAASO, 1e15)")
    print("et de l'electron (2/3 c) ; l_1 survit comme longueur de Compton, pas comme pas du milieu. Restent : les rapports du secteur")
    print("hadronique, les etats limites gravitationnels, et un programme : pas < 1e-27 m, noeud a spineur, alpha_1 et alpha_2 a calculer.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
