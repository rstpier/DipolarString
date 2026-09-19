#!/usr/bin/env python3
"""R94a -- Deriver le vrai analyseur SU(2) du doublet E, sans matrice de Pauli
(test fixe par l'auteur : prendre le doublet E de R82/R84, tourner physiquement
la geometrie de l'analyseur de theta, voir si les deux puissances sortent en
cos^2(theta/2), sin^2(theta/2)).

Le champ E de l'anneau d'axe z : au point d'azimut phi, plan de section
engendre par r(phi) et z ; deplacement du coeur d(phi) = cos(phi/2) r(phi)
+ sin(phi/2) z (demi-tour de section par circuit, R82/R84 : d(2 pi) = -d(0)).
Un anneau d'axe n = R z porte le champ R d(phi) aux points R r(phi).
Trois constructions geometriques d'analyseur oriente selon a = z, l'anneau
incline de theta (rotation autour de x), puissances = carres de projections
ou de recouvrements REELS :
  (i)   projection du deplacement sur a, moyennee le long de l'anneau ;
  (ii)  recouvrement des champs vectoriels apres rotation rigide,
        <d . R d> ;
  (iii) projection dans le plan de section au point de couplage sur l'axe a
        projete, moyennee.
Cible de Born : P(+) = cos^2(theta/2) : 1 a theta = 0, 1/2 a 90 deg, 0 a 180 deg.
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def field(phi):
    r_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    z = np.array([0.0, 0.0, 1.0])
    return np.cos(phi / 2)[:, None] * r_hat + np.sin(phi / 2)[:, None] * z

def rot_x(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def main():
    print("R94a -- l'analyseur SU(2) depuis la geometrie reelle du doublet E\n")
    phi = (np.arange(4000) + 0.5) * 4 * math.pi / 4000          # le champ est 4 pi-periodique
    d0 = field(phi)
    a_hat = np.array([0.0, 0.0, 1.0])
    thetas = [0.0, math.pi / 3, math.pi / 2, 2 * math.pi / 3, math.pi]
    born = [math.cos(t / 2) ** 2 for t in thetas]

    # (i) projection du deplacement sur l'axe de l'analyseur
    print("(i)  projection du deplacement sur a, moyenne le long de l'anneau : <(R d . a)^2> / <|d|^2>")
    P1 = []
    for t in thetas:
        d = d0 @ rot_x(t).T
        P1.append(float(np.mean((d @ a_hat) ** 2) / np.mean(np.sum(d * d, axis=1))))
    for t, p, b in zip(thetas, P1, born):
        print(f"    theta = {math.degrees(t):5.0f} deg : P = {p:.4f}   (Born {b:.4f})")
    print("  -> 1/2 a theta = 0 au lieu de 1 : pas Born.\n")
    check("(i) echoue : P(theta = 0) = 1/2, pas 1", abs(P1[0] - 0.5) < 1e-6, f"{P1[0]:.4f}")

    # (ii) recouvrement des champs apres rotation rigide
    print("(ii) recouvrement des champs vectoriels : <d . R d> / <|d|^2>")
    P2 = []
    for t in thetas:
        d = d0 @ rot_x(t).T
        P2.append(float(np.mean(np.sum(d0 * d, axis=1)) / np.mean(np.sum(d0 * d0, axis=1))))
    for t, p in zip(thetas, P2):
        print(f"    theta = {math.degrees(t):5.0f} deg : recouvrement = {p:+.4f}  (analytique (1 + 3 cos theta)/4 = {(1+3*math.cos(t))/4:+.4f} ; "
              f"cos(theta/2) = {math.cos(t/2):.4f})")
    print("  -> (1 + 3 cos theta)/4, negatif au-dela de 109 deg : ni cos(theta/2) ni une puissance.\n")
    check("(ii) echoue : recouvrement = (1 + 3 cos theta)/4, pas cos(theta/2)",
          all(abs(p - (1 + 3 * math.cos(t)) / 4) < 1e-6 for t, p in zip(thetas, P2)), "analytique")

    # (iii) projection dans le plan de section au point de couplage
    print("(iii) projection dans le plan de section, sur l'axe a projete, moyennee sur le point de couplage")
    P3 = []
    for t in thetas:
        R = rot_x(t)
        d = d0 @ R.T
        r_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1) @ R.T
        n_hat = R @ a_hat
        # axe a projete dans le plan de section (r_hat, n_hat) au point phi
        a_r = r_hat @ a_hat
        a_n = float(n_hat @ a_hat)
        norm = np.sqrt(a_r ** 2 + a_n ** 2)
        a_proj = (a_r[:, None] * r_hat + a_n * n_hat) / norm[:, None]
        P3.append(float(np.mean((np.sum(d * a_proj, axis=1)) ** 2) / np.mean(np.sum(d * d, axis=1))))
    for t, p, b in zip(thetas, P3, born):
        print(f"    theta = {math.degrees(t):5.0f} deg : P = {p:.4f}   (Born {b:.4f})")
    print("  -> 1/2 a theta = 0 : le point de couplage voit sin^2(phi/2), pas l'angle de l'axe.\n")
    check("(iii) echoue : P(theta = 0) = 1/2", abs(P3[0] - 0.5) < 1e-6, f"{P3[0]:.4f}")

    print("Ce que l'echec dit : aucune projection ou recouvrement REEL du deplacement du coeur ne")
    print("porte le demi-angle de l'axe ; le demi-angle de R93 A etait impose. Le cos(theta/2) est le")
    print("recouvrement de deux etats du fibre c_1 = 1 (R85), qui exige la structure COMPLEXE :")
    print("le doublet E combine a la phase de circulation (A = d_x + i d_y co-tournant, R53) et")
    print("transporte avec la connexion du fibre. Un analyseur qui ne voit que le deplacement reel ne")
    print("peut pas donner Born ; il doit etre sensible a la phase de circulation. Non dessine.\n")
    check("aucune construction reelle ne donne cos^2(theta/2) a theta = 0 (toutes donnent 1/2)",
          all(abs(p - 0.5) < 1e-6 for p in (P1[0], P3[0])), "verrou nomme : la phase complexe")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
