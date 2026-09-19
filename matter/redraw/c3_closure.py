#!/usr/bin/env python3
"""R76 -- Fermeture C3 des generations a partir des paires de DQD (fusion de
R60, R61 et R70, sur commentaire de l'auteur).

Le test n'est plus "pourquoi 3, 7, 11 ?" mais :
    U_pair^3 = I,   U_pair != I,   U_pair^2 != I,
avec U_pair l'operation physique "ajouter deux DQD" (4 brins, 4 = 1 mod 3 :
un tiers de tour dans le canal writhe, R61 ; une paire est le plus petit
ajout neutre qui garde le spin 1/2, R70).  Puis : n = 15 est-il une quatrieme
particule ou le retour dans le secteur de n = 3 ?

  A. U_pair comme decalage cyclique des trois positions : U^3 = I, U, U^2 != I.
  B. un operateur de RACINE de masse Z3-symetrique (circulant hermitien
     C = aI + bU + b*U^2, valeurs propres lambda_k = sqrt(m_k)) a exactement
     trois valeurs propres, a + 2|b| cos(phi + 2 pi k/3) : la forme de
     Koide-Brannen sort de la symetrie C3 seule ; a^2 = 2|b|^2 (Koide, 45 deg)
     et phi = 2/9 sont POSES, pas derives ; et 'ajouter deux DQD => U' est une
     DEFINITION de U, pas une derivation.  Avec ces entrees : m_mu, m_tau a 0,01 %.
  C. la quatrieme application : U^4 = U, pas de quatrieme valeur propre ; sur
     l'echelle (n/3)^(2 pi), n = 15 serait un lepton charge de 12,6 GeV vivant
     1e-17 s (R69), exclu par le LEP (> 100,8 GeV) : l'echelle comme loi de
     masse est exclue par l'absence de quatrieme generation, le circulant non.
  D. trois secteurs : N_nu = 2,9963 +- 0,0074 (largeur invisible du Z, PDG 2024).
"""
import math
import numpy as np

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
LEP_LEPTON_BOUND = 100.8e3        # MeV, lepton charge lourd
N_NU, N_NU_ERR = 2.9963, 0.0074      # PDG 2024 (largeur invisible du Z)
TWO_PI = 2 * math.pi

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R76 -- fermeture C3 des generations\n")

    # A. U_pair
    print("A. U_pair = decalage cyclique des trois positions (4 brins = 1 mod 3)")
    U = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
    I = np.eye(3)
    powers = {k: np.linalg.matrix_power(U, k) for k in (1, 2, 3, 4)}
    print(f"  U != I : {not np.allclose(powers[1], I)} ; U^2 != I : {not np.allclose(powers[2], I)} ; "
          f"U^3 = I : {np.allclose(powers[3], I)} ; U^4 = U : {np.allclose(powers[4], U)}")
    R = np.array([[math.cos(TWO_PI/3), -math.sin(TWO_PI/3)], [math.sin(TWO_PI/3), math.cos(TWO_PI/3)]])
    print(f"  canal writhe (R61) : R(2 pi/3)^3 = I : {np.allclose(np.linalg.matrix_power(R, 3), np.eye(2))}\n")
    check("U_pair^3 = I, U_pair != I, U_pair^2 != I",
          np.allclose(powers[3], I) and not np.allclose(powers[1], I) and not np.allclose(powers[2], I),
          "ordre 3 exactement")

    # B. le circulant de masse
    print("B. Operateur de racine de masse Z3-symetrique : C = a I + b U + conj(b) U^2, lambda_k = sqrt(m_k)")
    print("   (entrees posees : a^2 = 2|b|^2 et phi = 2/9 ; la symetrie C3 ne donne que la forme)")
    phi = 2 / 9
    # Koide : a = sqrt2 |b| ; a fixe par m_e via la plus petite valeur propre
    f = np.array([1 + math.sqrt(2) * math.cos(phi + TWO_PI * k / 3) for k in range(3)])
    a = math.sqrt(ME) / f.min()
    b = (a / math.sqrt(2)) * np.exp(1j * phi)
    C = a * I + b * U + np.conj(b) * np.linalg.matrix_power(U, 2)
    ev = np.sort(np.linalg.eigvalsh(C))
    masses = ev ** 2
    print(f"  hermitien : {np.allclose(C, C.conj().T)} ; valeurs propres sqrt(m) : {ev.round(4)}")
    print(f"  masses : e {masses[0]:.6f}, mu {masses[1]:.4f} ({100*(masses[1]/MMU-1):+.4f} %), "
          f"tau {masses[2]:.3f} ({100*(masses[2]/MTAU-1):+.4f} %)")
    # normes singulet / doublet du vecteur sqrt(m)
    v = ev
    s = np.ones(3) / math.sqrt(3)
    v_sing = np.dot(v, s) * s
    v_doub = v - v_sing
    ang = math.degrees(math.atan2(np.linalg.norm(v_doub), np.linalg.norm(v_sing)))
    print(f"  |singulet| = {np.linalg.norm(v_sing):.4f}, |doublet| = {np.linalg.norm(v_doub):.4f}, "
          f"angle a (1,1,1) = {ang:.3f} deg ; a^2 = {a*a:.4f}, 2|b|^2 = {2*abs(b)**2:.4f}")
    print("  -> Koide (45 deg) <=> a^2 = |b|^2 + |b|^2 : le terme propre au carre egale la somme des")
    print("     deux sauts au carre ; c'est la regle que la dynamique doit produire.\n")
    check("le circulant a exactement 3 valeurs propres et redonne mu, tau a 0,01 %",
          len(ev) == 3 and abs(masses[1] / MMU - 1) < 1e-4 and abs(masses[2] / MTAU - 1) < 1e-4,
          f"mu {100*(masses[1]/MMU-1):+.4f} %, tau {100*(masses[2]/MTAU-1):+.4f} %")
    check("Koide <=> a = sqrt2 |b| <=> angle 45 deg", abs(ang - 45) < 1e-6 and abs(a * a / (2 * abs(b) ** 2) - 1) < 1e-9,
          f"{ang:.6f} deg")

    # C. la quatrieme application
    print("C. Quatrieme application")
    m15 = ME * (15 / 3) ** TWO_PI
    tau15 = 1e-17
    print(f"  circulant : U^4 = U, le quatrieme etat est le premier ; pas de valeur propre nouvelle.")
    print(f"  echelle : n = 15 donnerait un lepton charge de {m15/1e3:.1f} GeV, tau ~ {tau15:.0e} s (R69),")
    print(f"  visible au LEP ; borne : > {LEP_LEPTON_BOUND/1e3:.1f} GeV. L'echelle comme loi de masse est exclue")
    print("  par l'absence de quatrieme generation ; le circulant la predit absente.\n")
    check("l'echelle predit un 4e lepton sous la borne du LEP (exclue comme loi de masse)",
          m15 < LEP_LEPTON_BOUND, f"{m15/1e3:.1f} GeV < {LEP_LEPTON_BOUND/1e3:.1f}")
    check("U^4 = U : pas de quatrieme etat propre", np.allclose(powers[4], U), "ordre 3")

    # D. trois secteurs
    print("D. Compte des secteurs")
    print(f"  largeur du Z : N_nu = {N_NU} +- {N_NU_ERR} ; le circulant en a 3.\n")
    check("trois secteurs (N_nu = 3 a 3 sigma)", abs(N_NU - 3) < 3 * N_NU_ERR, f"{N_NU} +- {N_NU_ERR}")

    # E. ce qui manque
    print("E. Ce qui manque (le point dur du commentaire)")
    print("  U^3 = I sur le spectre ne dit pas encore que six DQD a holonomie nulle ne se lient pas :")
    print("  lecture candidate, six DQD fermes sur eux-memes sans tiers de tour sont un morceau de vide")
    print("  (R3), ils ne portent rien et se detachent a cout nul ; deux DQD portent un tiers de tour")
    print("  qui les verrouille (la desintegration est le defaire, (m/M)^4 par tour, R69).")
    print("  Cela demande une energie de liaison par holonomie : non calculee.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
