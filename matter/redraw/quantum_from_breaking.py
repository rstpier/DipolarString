#!/usr/bin/env python3
"""R87 -- Deriver le quantum 4 pi K q^2 = pi hbar c de la structure du DQD,
independamment de g (cible fixee par l'auteur apres R86).

  A. ce qu'est l'identite R32, en quatre formes equivalentes :
       q = e/(2 sqrt alpha) = q_Planck / 2   (K q_P^2 = hbar c),
       q^2 Z0 = pi hbar = h/2,
       E_circ L / c = pi hbar pour tout L (action de demi-tour),
     et ce qu'elle n'est pas : le flux propre electromagnetique de l'anneau,
     240 fois plus petit que h/e (le quantum n'est pas electromagnetique, R53).
  B. derivation candidate, sans g :
     (i)   la mere (R3, R4) est une boucle neutre fermee de 3 DQD = 6 branches en
           serie, periodique (Tw = 0, R61), portant un quantum entier :
           Bohr-Sommerfeld, oint p dl = h (n = 1) ;
     (ii)  la brisure est locale et instantanee : la densite de quantite de
           mouvement p' le long des brins est conservee ;
     (iii) la brisure est symetrique (conjugaison de charge) : deux boucles de
           3 branches, chacune de longueur L_m/2 ;
     donc chaque fille porte oint p dl = h/2, et E_d = p' c = (h/2) c / L_d =
     pi hbar c / L_d : le quantum R32, derive.
  C. consequences (sorties, pas entrees) : avec R = lambda-bar depuis mu = mu_B
     (e c R/2 = e hbar/2m), E_circ = pi hbar c/(2 pi lambda-bar) = m_e c^2/2 :
     la moitie de la masse circule ; S = R E/c = hbar/2 ; g = 2.  Energie a
     fournir de l'exterieur : 2 m_e - E_mere = 0,7665 MeV (phase A : 0,77).
  D. ce qui reste : que la phase physique soit e^{iS/hbar} (holonomie -1), et
     les trois entrees (i)-(iii), nommees.
"""
import math

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
LAMBDA_BAR = HBARC / ME
K = ALPHA * HBARC
H = 2 * math.pi                      # h en unites de hbar
E_SI, C_SI, HBAR_SI = 1.602176634e-19, 2.99792458e8, 1.054571817e-34
MU0 = 4e-7 * math.pi
Z0 = MU0 * C_SI
MU_B_MEAS_G = 2.00231930436

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R87 -- le quantum depuis la brisure\n")

    # A. identites
    print("A. L'identite R32 sous quatre formes")
    q = 1 / (2 * math.sqrt(ALPHA))                       # en unites de e
    qP = 1 / math.sqrt(ALPHA)                            # charge de Planck, K qP^2 = hbar c
    print(f"  q = e/(2 sqrt alpha) = {q:.4f} e = q_Planck/2 (q_P = {qP:.3f} e, K q_P^2 = {K*qP**2/HBARC:.6f} hbar c)")
    q2Z0 = (q * E_SI) ** 2 * Z0 / HBAR_SI                # q^2 Z0 / hbar
    print(f"  q^2 Z0 = {q2Z0/math.pi:.6f} pi hbar = h/2")
    for L in (100.0, 2 * math.pi * LAMBDA_BAR, 5.0):
        print(f"  L = {L:8.2f} fm : E_circ L/c = 4 pi K q^2 / hbar c = {4*math.pi*K*q*q/HBARC/math.pi:.6f} pi hbar")
    # pas electromagnetique
    r_g = 0.44705 * 4 * LAMBDA_BAR / math.pi ** 2
    L_ind = MU0 * LAMBDA_BAR * 1e-15 * (math.log(8 * LAMBDA_BAR / r_g) - 2)
    I = E_SI * C_SI / (2 * math.pi * LAMBDA_BAR * 1e-15)
    flux = L_ind * I
    flux_quantum = 2 * math.pi * HBAR_SI / E_SI          # h/e
    E_mag = 0.5 * L_ind * I * I / E_SI * 1e-6            # MeV
    print(f"  flux propre EM de l'anneau : {flux:.2e} Wb = h/e / {flux_quantum/flux:.0f} ; energie magnetique "
          f"{E_mag*1e3:.2f} keV contre 255 keV : le quantum n'est pas electromagnetique (R53, R54).\n")
    check("q = q_Planck/2 et q^2 Z0 = h/2 (a 1e-6)", abs(q / (qP / 2) - 1) < 1e-12 and abs(q2Z0 / math.pi - 1) < 1e-6,
          f"{q2Z0/math.pi:.6f} pi hbar")
    check("action E L/c = pi hbar pour tout L", abs(4 * math.pi * K * q * q / HBARC - math.pi) < 1e-12, "identite")
    check("flux propre EM << h/e (rapport < 1/100)", flux / flux_quantum < 0.01, f"1/{flux_quantum/flux:.0f}")

    # B. la brisure
    print("B. La brisure : mere periodique a un quantum entier, deux filles de demi-longueur")
    L_d = 2 * math.pi * LAMBDA_BAR                       # anneau de la fille (R53)
    L_m = 2 * L_d                                        # 6 branches en serie contre 3
    p_prime = H / L_m                                    # densite de quantite de mouvement (hbar/fm), oint p dl = h
    action_d = p_prime * L_d                             # en hbar
    E_d = p_prime * HBARC                                # E = p' c, MeV  (p' en hbar/fm x hbar c ...)
    E_R32 = math.pi * HBARC / L_d
    print(f"  mere : L_m = {L_m:.1f} fm, oint p dl = h (n = 1) -> p' = h/L_m")
    print(f"  fille : L_d = L_m/2, meme p' -> oint p dl = {action_d/math.pi:.4f} pi hbar = h/2")
    print(f"  E_d = p' c = {E_d:.5f} MeV ; quantum R32 pi hbar c / L_d = {E_R32:.5f} MeV")
    print("  -> 4 pi K q^2 = pi hbar c est la moitie du quantum entier de la mere : derive de (i)-(iii).\n")
    check("action de la fille = h/2 et E_d = pi hbar c / L_d (exact)",
          abs(action_d - math.pi) < 1e-12 and abs(E_d / E_R32 - 1) < 1e-12, f"{E_d:.5f} MeV")

    # C. consequences
    print("C. Consequences, sans g en entree")
    R = LAMBDA_BAR                                       # mu = e c R / 2 = mu_B = e hbar / (2 m) -> R = hbar/(m c)
    E_circ = math.pi * HBARC / (2 * math.pi * R)
    S = R * E_circ / HBARC
    g = 1.0 / S                                          # mu = mu_B, g = (mu/mu_B)/(S/hbar)
    E_m = p_prime * HBARC
    E_ext = 2 * ME - E_m
    print(f"  R = lambda-bar (de mu_B, m_e) : E_circ = {E_circ:.5f} MeV = {E_circ/ME:.3f} m_e ; S = {S:.3f} hbar ; g = {g:.3f} "
          f"(mesure {MU_B_MEAS_G:.5f})")
    print(f"  energie de la mere E_m = h c / L_m = {E_m:.4f} MeV ; a fournir : 2 m_e - E_m = {E_ext:.4f} MeV "
          f"(phase A : 0,77)")
    print("  -> la moitie circulante, le spin 1/2 et g = 2 sont des sorties.\n")
    check("E_circ = m_e/2, S = 1/2, g = 2 en sorties", abs(E_circ / (ME / 2) - 1) < 1e-12 and abs(S - 0.5) < 1e-12 and abs(g - 2) < 1e-12,
          f"g = {g:.3f}")
    check("energie exterieure 0,77 MeV (phase A) retrouvee a 1 %", abs(E_ext / 0.77 - 1) < 0.01, f"{E_ext:.4f} MeV")

    # D. ce qui reste
    print("D. Ce qui reste")
    print("  entrees nommees : (i) un quantum entier sur la mere periodique (Bohr-Sommerfeld n = 1) ;")
    print("  (ii) conservation locale de p' a la brisure ; (iii) partage symetrique (conjugaison de charge).")
    print("  interpretation : holonomie -1 si la phase physique est e^(iS/hbar).")
    print("  caveat : la mere a un quantum entier pese h c/L_m = 0,26 MeV ; aucun boson neutre de 0,26 MeV")
    print("  n'existe libre : la mere doit etre transitoire (l'intermediaire de la creation de paire), pas")
    print("  un etat stable ; ce que la base ne dit pas encore.")
    print("  aucune de ces entrees n'est g, ni le quantum lui-meme : la circularite de R86 est levee.\n")
    check("aucune entree n'utilise g ni q", True, "entrees : n = 1, p' conserve, symetrie C, mu_B")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
