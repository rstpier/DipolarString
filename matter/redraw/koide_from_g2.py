#!/usr/bin/env python3
"""R81 -- Le 45 deg de Koide comme conservation : X_singulet = X_doublet avec un
quantum par representation irreductible (relecture de l'auteur : "chercher une
quantite X conservee lors de la formation du fermion ... si cela sort de la
rupture du DQD mere ou du partage circulation/jonctions deja utilise pour
g = 2, Koide commencerait a etre derive").

Resultat : X est l'energie, et le partage est celui de g = 2.
  A. pour un circulant (diagonale constante), ||C_iso||^2 = somme des poids sur
     site, ||C_dev||^2 = somme des poids de saut : le canal singulet est ce qui
     reste sur une position (les jonctions), le canal doublet ce qui passe
     d'une position a l'autre (la circulation).  Exact.
  B. tr C^2 = m_e + m_mu + m_tau : la masse de la famille est le poids total du
     reseau, sur site + sauts.  Exact.
  C. g = 2 pour chaque lepton (R53, R54 : moitie statique aux jonctions, moitie
     en circulation) somme sur la famille : E_stat = E_circ = (somme m)/2.
  D. appariement : statique <-> sur site, circulant <-> saut.  Alors C donne
     a^2 = 2|b|^2 sans sqrt2 a la main : c'est Koide.
  E. l'appariement est familial, pas par particule (m_e/2 << a^2) ; il ne dit
     rien de phi = 2/9 ; et pour une fraction statique f quelconque il donne
     Q = 1/(3f).
Ce qui reste pose : que l'operateur de racine de masse soit la matrice
d'amplitudes du reseau (jonctions sur la diagonale, circulation hors diagonale).
"""
import math
import numpy as np

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
G_E, G_MU = 2.00231930436, 2.00233184123          # mesures ; g_tau = 2 a 2 %
U = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
I3 = np.eye(3)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def circulant(a, b):
    return a * I3 + b * U + np.conj(b) * U @ U

def main():
    print("R81 -- Koide depuis g = 2\n")
    SM = ME + MMU + MTAU

    # A. canaux = sur site / saut
    print("A. Pour un circulant : ||iso||^2 = poids sur site, ||dev||^2 = poids de saut")
    a, b = 1.7, 0.9 * np.exp(1j * 0.37)
    C = circulant(a, b)
    iso = (np.trace(C) / 3) * I3
    dev = C - iso
    onsite = sum(abs(C[i, i]) ** 2 for i in range(3))
    hop = sum(abs(C[i, j]) ** 2 for i in range(3) for j in range(3) if i != j)
    print(f"  exemple a = {a}, |b| = {abs(b):.2f} : ||iso||^2 = {np.linalg.norm(iso)**2:.4f} = sur site {onsite:.4f} ;"
          f" ||dev||^2 = {np.linalg.norm(dev)**2:.4f} = sauts {hop:.4f}")
    print("  -> singulet = ce qui reste sur une position ; doublet = ce qui passe entre positions.\n")
    check("||iso||^2 = poids sur site et ||dev||^2 = poids de saut (circulant)",
          abs(np.linalg.norm(iso) ** 2 - onsite) < 1e-12 and abs(np.linalg.norm(dev) ** 2 - hop) < 1e-12,
          "identite")

    # B. tr C^2 = somme des masses
    print("B. tr C^2 = somme des m_k : la masse de la famille est le poids total du reseau")
    ev = np.linalg.eigvalsh(C)
    print(f"  exemple : tr C^2 = {np.trace(C @ C).real:.4f} = sum lambda^2 = {np.sum(ev**2):.4f}\n")
    check("tr C^2 = sum lambda_k^2", abs(np.trace(C @ C).real - np.sum(ev ** 2)) < 1e-10, "identite")

    # C. g = 2 somme sur la famille
    print("C. g = 2 : moitie statique, moitie circulante, pour chaque lepton (R53, R54)")
    E_stat = sum(m / 2 for m in (ME, MMU, MTAU))
    E_circ = SM - E_stat
    print(f"  g_e = {G_E:.6f}, g_mu = {G_MU:.6f} : E_stat,k = E_circ,k = m_k/2 a 0,1 %")
    print(f"  famille : E_stat = {E_stat:.2f} MeV = E_circ = {E_circ:.2f} MeV = (sum m)/2\n")
    check("g = 2 a 0,2 % pour e et mu", abs(G_E / 2 - 1) < 2e-3 and abs(G_MU / 2 - 1) < 2e-3, "mesure")

    # D. appariement => Koide
    print("D. Appariement : statique <-> sur site (jonctions), circulant <-> saut (circulation)")
    a2 = E_stat / 3                      # poids sur site par position = a^2
    b2 = E_circ / 6                      # poids par saut = |b|^2
    print(f"  a^2 = E_stat/3 = {a2:.3f} MeV ; |b|^2 = E_circ/6 = {b2:.3f} MeV ; a^2 / (2|b|^2) = {a2/(2*b2):.6f}")
    # comparer avec le a de Koide-Brannen ajuste (R76)
    phi = 2 / 9
    f = np.array([1 + math.sqrt(2) * math.cos(phi + 2 * math.pi * k / 3) for k in range(3)])
    a_koide = math.sqrt(ME) / f.min()
    print(f"  a de Koide + 2/9 + m_e (R76) : a^2 = {a_koide**2:.3f} MeV ({100*(a2/a_koide**2-1):+.4f} %)")
    Q = SM / (math.sqrt(ME) + math.sqrt(MMU) + math.sqrt(MTAU)) ** 2
    print(f"  -> a^2 = 2|b|^2 exactement par g = 2 ; Q mesure = {Q:.6f} = 2/3 a {100*(Q/(2/3)-1):+.4f} %.")
    print("     Le sqrt2 n'est pas pose : il est le 'moitie-moitie' de g = 2 lu dans la base des positions.\n")
    check("appariement + g = 2 => a^2 = 2|b|^2 (a^2 de g = 2 = a^2 de Koide a 1e-4)",
          abs(a2 / a_koide ** 2 - 1) < 1e-4, f"{100*(a2/a_koide**2-1):+.4f} %")

    # E. portee
    print("E. Portee de l'appariement")
    print(f"  par particule il est faux : m_e/2 = {ME/2:.3f} MeV contre a^2 = {a2:.1f} MeV (rapport {ME/2/a2:.1e}) ;")
    print("  il vaut pour la somme sur la famille : les trois jonctions de chaque lepton, sommees sur")
    print(f"  e, mu, tau, font a^2 par position ({(ME+MMU+MTAU)/6:.2f} MeV).")
    for fstat in (0.5, 0.813, 1/3):
        print(f"  fraction statique f = {fstat:.3f} : Q = 1/(3f) = {1/(3*fstat):.3f}"
              + ("  (Koide)" if abs(fstat - 0.5) < 1e-9 else "  (nucleon R72 : 81 % statique)" if fstat > 0.8 else ""))
    print("  phi = 2/9 n'est pas touche (les poids ne dependent pas de la phase).\n")
    check("l'appariement est familial (m_e/2 / a^2 < 1e-2) et Q = 1/(3f) redonne 2/3 a f = 1/2",
          ME / 2 / a2 < 1e-2 and abs(1 / (3 * 0.5) - 2 / 3) < 1e-12, f"{ME/2/a2:.1e}")

    print("Ce qui reste pose : l'operateur de racine de masse est la matrice d'amplitudes du reseau")
    print("a trois positions (jonction = noeud, circulation = lien), avec la masse = poids total.")
    print("R80 lisait les canaux comme deux circuits : c'etait la mauvaise lecture ; les canaux sont")
    print("noeud contre lien, et le doublet EST la circulation de l'anneau.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
