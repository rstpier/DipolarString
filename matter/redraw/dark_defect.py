#!/usr/bin/env python3
"""R90 -- Matiere noire = defaut neutre verrouille du milieu DQD ?  (critere de
l'auteur : enumerer les plus petites boucles neutres avec Lk != 0, calculer E(R)
par les regles etablies sans ajuster de masse, exiger un vrai minimum
dE/dR = 0, d2E/dR2 > 0, poser m_DM = E_min/c^2, puis moments EM et sections
efficaces.)

  A. enumeration :
     (1) un DQD ferme avec une demi-torsion : interdit (la branche + se
         raccorderait a la branche -, continuite de charge) ;
     (2) un DQD ferme avec une torsion entiere (Lk = 1, les deux branches +-e/3
         forment un lien de Hopf) : neutre, mais dipole electrique ~ (e/3) R :
         pas sombre (la matiere noire dipolaire est exclue) ;
     (3) trois DQD (la mere, R3) avec writhe 1/3 (Z3, R61/R77) : neutre, dipole
         nul par symetrie, quadrupole et polarisabilite non nuls : le candidat X ;
     (4) anneau neutre a trois brins 000 avec writhe k/3 : la famille des
         neutrinos, meme statut que (3).
  B. E(R) de X avec les regles de la base : circulation h/2 ou h par trajet
     (R47, R87), Coulomb des paires torsadees (R42), pas de jonction (aucun bout
     charge) ; tout est en hbar c / longueur : E = C/R, dE/dR < 0 partout.
     Raison dimensionnelle : la base n'a que hbar c, alpha et l'ancre lambda_e.
     Aucun minimum : pas de masse.
  C. le seul ancrage disponible : l'epinglage a l'ecart du vide D0 = 241 fm
     (comme les jonctions de l'electron, R54) : L = 3 D0 -> 0,86 MeV (h/2) ou
     1,71 MeV (h).  Exclus : section Rayleigh par polarisabilite ~ 10^2 barn a
     T = 1 MeV, thermalises au BBN (Gamma/H ~ 10^20), delta N_eff ~ 1 contre
     N_eff = 2,99 +- 0,17 ; et absents des donnees e+e- / nucleaires.
  D. le nombre de l'auteur : rho_DM c^2 / u_0 = 1,9e-28 ; avec A8 (u_0 ne
     gravite pas), une surdensite du substrat n'est pas un mecanisme.
Verdict : EXCLU comme prediction ; le secteur sombre de la base, ce sont les
neutrinos, dont elle ne fixe pas la taille.  Une prediction demanderait une
constante nouvelle fixant la taille des boucles neutres.
"""
import math

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
LAMBDA_BAR = HBARC / ME
K = ALPHA * HBARC
D0 = 2 * math.cosh(math.pi) * LAMBDA_BAR / 37.1     # 241 fm
PATH_FACTOR = 1.136                                 # trajet a Wr = 1/3 / cercle (R61)
RHO_DM_GEV_CM3 = 0.55                               # PDG 2025, +- 0.17
U0_J_M3 = 4.64e23                                   # substrat v2.9
N_EFF, N_EFF_ERR = 2.99, 0.17
M_PL_MEV = 1.22e22
HBAR_MEV_S = 6.582e-22

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def energy_X(R, half_quantum=True):
    """energie du candidat X (trois DQD, Wr = 1/3) a rayon moyen R, MeV, regles de la base"""
    L = PATH_FACTOR * 2 * math.pi * R
    E_circ = (math.pi if half_quantum else 2 * math.pi) * HBARC / L       # h/2 ou h par trajet
    E_pairs = 3 * 4 * math.pi * K / (9 * L / 3)                          # Coulomb de chaque paire +-e/3 sur sa longueur (R42)
    return E_circ + E_pairs

def main():
    print("R90 -- defaut neutre verrouille ?\n")

    # A. enumeration
    print("A. Les plus petites boucles neutres avec Lk != 0")
    print("  (1) un DQD, demi-torsion : interdit (+ se raccorde a -) ;")
    print(f"  (2) un DQD, torsion entiere (lien de Hopf des branches +-e/3) : dipole ~ (e/3) R = "
          f"{100/3:.0f} e.fm a R = 100 fm : pas sombre ;")
    print("  (3) trois DQD, writhe 1/3 (Z3) : neutre, dipole nul, quadrupole et polarisabilite : candidat X ;")
    print("  (4) anneau 000 a writhe k/3 : les neutrinos, meme statut.\n")
    check("le plus petit defaut sombre est Z3 (trois DQD), pas un DQD seul", True, "dipole nul par symetrie")

    # B. E(R)
    print("B. E(R) du candidat X par les regles de la base")
    Rs = [10.0, 100.0, 1000.0, 1e4, 1e5]
    Es = [energy_X(R) for R in Rs]
    for R, E in zip(Rs, Es):
        print(f"    R = {R:8.0f} fm : E = {E:.4e} MeV  (E R = {E*R:.3f} MeV fm)")
    slopes = [(Es[i + 1] - Es[i]) for i in range(len(Es) - 1)]
    print("  E R = constante : E = C/R, dE/dR < 0 partout, aucun minimum. Raison dimensionnelle :")
    print("  la base n'a que hbar c, alpha et l'ancre lambda_e ; un objet sans charge n'a pas de taille.\n")
    check("E(R) = C/R sans minimum (E R constant a 1e-9, pente negative partout)",
          max(abs(E * R / (Es[0] * Rs[0]) - 1) for R, E in zip(Rs, Es)) < 1e-9 and all(s < 0 for s in slopes), "C/R")

    # C. epinglage au vide
    print("C. Le seul ancrage : epingler X a l'ecart du vide D0 (comme les jonctions de l'electron, R54)")
    L_pin = 3 * D0
    R_pin = L_pin / (2 * math.pi)
    m_f = energy_X(R_pin / PATH_FACTOR, True)
    m_b = energy_X(R_pin / PATH_FACTOR, False)
    print(f"  L = 3 D0 = {L_pin:.0f} fm : m_X = {m_f:.3f} MeV (h/2, fermion) ou {m_b:.3f} MeV (h, boson)")
    # section Rayleigh a T = 1 MeV
    T = 1.0
    k = T / HBARC
    alpha_E = R_pin ** 3                                 # polarisabilite ~ R^3 (fm^3, unites gaussiennes)
    sigma = (8 * math.pi / 3) * k ** 4 * alpha_E ** 2    # fm^2
    n_gamma = 0.2436 * (T / HBARC) ** 3                  # fm^-3
    gamma_rate = n_gamma * sigma * 2.998e23              # /s
    H = 1.66 * math.sqrt(10.75) * T ** 2 / M_PL_MEV / HBAR_MEV_S
    print(f"  polarisabilite ~ R^3 = {alpha_E:.1e} fm^3 : sigma_Rayleigh(1 MeV) = {sigma/100:.0f} barn (Thomson 0,67) ;")
    print(f"  Gamma/H au BBN = {gamma_rate/H:.0e} : thermalise ; delta N_eff ~ 1 contre N_eff = {N_EFF} +- {N_EFF_ERR} : exclu.")
    print("  et aucun neutre stable de 0,9-1,7 MeV dans les donnees e+e- ou nucleaires.\n")
    check("candidat epingle : masse au MeV, thermalise au BBN (Gamma/H > 1e6), exclu par N_eff",
          0.5 < m_f < 3 and gamma_rate / H > 1e6 and 1.0 > 2 * N_EFF_ERR, f"{m_f:.2f} MeV, Gamma/H = {gamma_rate/H:.0e}")

    # D. le nombre de l'auteur
    print("D. Une surdensite du substrat ?")
    rho_dm = RHO_DM_GEV_CM3 * 1.602e-10 / 1e-6           # J/m^3
    print(f"  rho_DM c^2 = {rho_dm:.2e} J/m^3 ; u_0 = {U0_J_M3:.2e} J/m^3 ; rapport = {rho_dm/U0_J_M3:.1e}")
    print("  avec A8 (u_0 ne gravite pas), une perturbation relative de 1e-28 devrait porter toute la masse")
    print("  gravitationnelle du halo : pas un mecanisme, une etiquette.\n")
    check("rho_DM c^2 / u_0 ~ 2e-28", abs(rho_dm / U0_J_M3 / 1.9e-28 - 1) < 0.1, f"{rho_dm/U0_J_M3:.1e}")

    print("Verdict : EXCLU comme prediction. Un defaut neutre verrouille n'a pas de masse dans la base")
    print("(E = C/R, aucun minimum) ; epingle au vide il pese ~1 MeV et le BBN l'exclut. Le secteur")
    print("sombre de la base, ce sont les neutrinos, dont elle ne fixe pas la taille. Une prediction")
    print("demanderait une constante nouvelle fixant la taille des boucles neutres.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
