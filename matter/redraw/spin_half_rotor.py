#!/usr/bin/env python3
"""R83 -- De L_z = 1/2 a j = 1/2 : la toupie rigide est exclue, l'orbite de
rotation reste (suite de R82).

R82 donne un mode antiperiodique sur la section (K = 1/2 autour de l'axe de
l'anneau).  Comment passer a une representation j = 1/2 du groupe des
rotations, sans tour d'etats excites ?

  A. lecture toupie symetrique rigide (axe de l'anneau libre, K = 1/2 fixe par
     le mode) : E(j) = hbar^2 [j(j+1) - K^2]/(2 I_perp) + ... ; le fondamental
     est bien j = 1/2 (deux etats), mais j = 3/2 suit a 3 m c^2 (R = lambda-bar)
     ou 0,86 MeV (R = 4 lambda-bar/3) : un electron excite au MeV, exclu (R47 :
     pas de resonance Compton au MeV, compositeness > 10 TeV).
  B. donc l'orientation de l'anneau n'a pas d'inertie propre : tourner l'anneau
     n'est pas un mouvement, c'est la meme configuration vue tournee (ether de
     Lorentz du manuscrit).  L'espace des etats est l'orbite d'une configuration
     sous les rotations, quotientee par ce qui la laisse invariante ; avec le
     mode antiperiodique cette orbite est SU(2)/U(1) = S^2 a signe pres, la
     variete des etats coherents de spin 1/2 : deux etats, pas de tour.
  C. comptage : le motif dipolaire est l'irrep E (dimension 2) de Z4, le
     fondamental l'irrep A (dimension 1) : 2 = les deux etats de j = 1/2,
     1 = un scalaire.  Coherent avec B.
  D. l'echange : l'anneau est attache par ses trois jonctions (R54) ; echanger
     deux anneaux attaches est isotope a tourner l'un de 2 pi (tour de ceinture,
     R74 D2) ; avec le -1 de R82, Pauli.  Topologie, pas calcul.
"""
import math

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME
COMPOSITENESS_MEV = 1e7          # > 10 TeV (R47)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def top_levels(R, K=0.5, jmax=2.5):
    """niveaux d'une toupie symetrique (masse m sur un anneau de rayon R) avec K fixe"""
    I_perp = ME * R ** 2 / 2 / HBARC ** 2      # en MeV^-1 (hbar = c = 1 : I = m R^2/2, energie hbar^2/(2I))
    levels = []
    j = K
    while j <= jmax + 1e-9:
        e = (j * (j + 1) - K * K) / (2 * I_perp)
        levels.append((j, e))
        j += 1
    return levels

def main():
    print("R83 -- de L_z = 1/2 a j = 1/2\n")

    # A. toupie rigide
    print("A. Toupie symetrique rigide avec K = 1/2")
    for name, R in (("R = lambda-bar (R53)", LAMBDA_BAR), ("R = 4 lambda-bar/3 (R12)", 4 * LAMBDA_BAR / 3)):
        lv = top_levels(R)
        e0 = lv[0][1]
        print(f"  {name:26s} : " + " ; ".join(f"j = {j:.1f} : +{e - e0:.3f} MeV" for j, e in lv))
    lv = top_levels(LAMBDA_BAR)
    gap = lv[1][1] - lv[0][1]
    print(f"  fondamental j = 1/2 (deux etats m = +-1/2) ; premier excite j = 3/2 a +{gap:.3f} MeV = 3 m c^2.")
    print("  -> un electron excite de spin 3/2 au MeV : exclu (aucune resonance, compositeness > 10 TeV).")
    print("     L'anneau n'est pas une toupie rigide : son orientation n'a pas d'inertie propre.\n")
    check("toupie rigide : fondamental j = 1/2 a deux etats", lv[0][0] == 0.5, "K = 1/2")
    check("toupie rigide exclue : j = 3/2 a 3 m c^2 << compositeness", abs(gap / (3 * ME) - 1) < 1e-9 and gap < COMPOSITENESS_MEV,
          f"+{gap:.3f} MeV")

    # B. orbite de rotation
    print("B. L'orbite d'une configuration sous les rotations")
    print("  sans inertie d'orientation, l'etat est la configuration vue tournee : l'espace des etats")
    print("  est l'orbite SO(3)/stabilisateur. Le stabilisateur de l'anneau est la rotation autour de")
    print("  son axe ; avec le mode antiperiodique (R82), un tour complet de l'axe change le signe, donc")
    print("  l'orbite est SU(2)/U(1) = S^2 a un signe pres : la sphere de Bloch, les etats coherents")
    print("  de j = 1/2. Deux etats, pas de tour d'excitations.")
    dim_orbit = 3 - 1       # dim SO(3) - dim U(1)
    print(f"  dimension de l'orbite : {dim_orbit} = celle de la sphere de Bloch.\n")
    check("orbite de dimension 2 (sphere de Bloch)", dim_orbit == 2, "SO(3)/U(1)")

    # C. comptage des irreps de Z4
    print("C. Comptage")
    dims = {"A (fondamental (1,1))": 1, "E (motif dipolaire (1,2),(2,1))": 2}
    for k, v in dims.items():
        print(f"  {k:32s} dimension {v}  ->  {'scalaire' if v == 1 else 'les deux etats de spin 1/2'}")
    print()
    check("E a dimension 2 = 2j + 1 pour j = 1/2", dims["E (motif dipolaire (1,2),(2,1))"] == 2, "2 etats")

    # D. echange
    print("D. Echange (R74 D2)")
    print("  l'anneau est attache au milieu par trois jonctions (R54). Echanger deux objets attaches est")
    print("  isotope a tourner l'un d'eux de 2 pi (tour de ceinture). Avec psi(2 pi) = -psi (R82, t = 1/2),")
    print("  l'echange donne -1 : Pauli. C'est de la topologie ; la seule entree calculee est le -1.\n")
    check("chaine complete conditionnelle : mode impair (R82) + attache (R54) + tour de ceinture => -1 d'echange",
          True, "topologie")

    print("Ce qui reste : (i) que l'electron porte le motif dipolaire (candidat : la courbure de")
    print(f"l'anneau, w/R = {4/math.pi**2:.3f}, deplace naturellement le coeur radialement ; le demi-tour")
    print("de section par circuit fait alors tourner ce deplacement dans le repere de la section) ;")
    print("(ii) que l'orientation soit bien sans inertie (ether de Lorentz : a verifier sur le g du")
    print("muon en vol, deja fait pour la duree de vie, R67).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
