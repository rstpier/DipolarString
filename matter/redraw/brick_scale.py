#!/usr/bin/env python3
"""R71 -- Point 1 de coherence : une brique, une seule taille ?

Le DQD du vide a un ecart D0 = 241 fm et un ruban w_e = 156 fm ; les brins du
quark font 0,16 fm.  Regle candidate : une ligne adaptee (Z = Z0) n'a pas de
taille propre, parce que l'impedance et l'energie par longueur d'une ligne
bifilaire ne dependent que des rapports de forme (d/w).  Les tailles
appartiennent aux circuits (ancre m_e + echelle), pas a la brique.

  A. methode des moments 2D : C' de deux conducteurs carres (w, d) et
     (100 w, 100 d) -> identiques ; Z et l'energie par longueur aussi.
  B. w/l1 = 6/pi^3 a tout barreau (par construction, R17/R54) et D0 = 241 fm
     (g = 2, R54) : tout pend a l'ancre lambda_e, rien a la brique.
  C. combien de cellules du vide occupent l'electron et le quark.
"""
import math
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
W_E = 4 * LAMBDA_BAR / math.pi ** 2
D0 = 2 * math.cosh(math.pi) * LAMBDA_BAR / 37.1   # 241 fm, ecart du DQD (manuscrit)

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def square_perimeter(cx, w, n_side):
    """segments (centres, longueur) sur le perimetre d'un carre de cote w centre en cx"""
    h = w / n_side
    pts = []
    for i in range(n_side):
        s = -w / 2 + (i + 0.5) * h
        pts += [(cx + s, -w / 2), (cx + s, w / 2), (cx - w / 2, s), (cx + w / 2, s)]
    return np.array(pts), h

def capacitance_per_length(w, d, n_side=24):
    """C' (unites 2 pi eps0) de deux carres de cote w, centres a distance d"""
    p1, h = square_perimeter(-d / 2, w, n_side)
    p2, _ = square_perimeter(+d / 2, w, n_side)
    pts = np.vstack([p1, p2])
    N = len(pts)
    dx = pts[:, 0][:, None] - pts[:, 0][None, :]
    dy = pts[:, 1][:, None] - pts[:, 1][None, :]
    r = np.sqrt(dx * dx + dy * dy)
    P = np.where(r > 0, -np.log(np.where(r > 0, r, 1.0)), 1.0 - math.log(h / 2))
    V = np.concatenate([np.full(N // 2, 0.5), np.full(N // 2, -0.5)])
    q = np.linalg.solve(P, V)
    return q[: N // 2].sum() / 1.0          # Q par unite de longueur pour dV = 1

def main():
    print("R71 -- la brique a-t-elle une taille ?\n")

    # A. invariance d'echelle de la ligne bifilaire
    print("A. Ligne bifilaire carree (w, d = 2w) a deux echelles")
    c1 = capacitance_per_length(1.0, 2.0)
    c2 = capacitance_per_length(100.0, 200.0)
    c3 = capacitance_per_length(1.0, 3.0)
    print(f"  C'(w = 1, d = 2)   = {c1:.6f} (2 pi eps0)")
    print(f"  C'(w = 100, d = 200) = {c2:.6f}   rapport {c2/c1:.9f}")
    print(f"  C'(w = 1, d = 3)   = {c3:.6f} : seul le rapport d/w compte")
    print("  Z = sqrt(L'/C') = 1/(c C') et l'energie par longueur lambda^2/(2C') suivent :")
    print("  -> une ligne adaptee n'a pas de taille ; la comprimer 1000 fois ne coute rien.\n")
    check("C' invariante d'echelle (a 1e-9)", abs(c2 / c1 - 1) < 1e-9, f"rapport {c2/c1:.9f}")
    check("C' depend de d/w (d = 3w differe de d = 2w)", abs(c3 / c1 - 1) > 0.05,
          f"{c3/c1:.3f}")

    # B. ce qui fixe les tailles : l'ancre et l'echelle
    print("B. Les tailles du redessin")
    ratios = [W_E * l1(n) / L1_3 / l1(n) for n in (3, 7, 9, 11)]
    print(f"  w/l1 = {ratios[0]:.5f} = 6/pi^3 = {6/math.pi**3:.5f} a tout barreau (par construction)")
    print(f"  ecart du vide D0 = 2 cosh(pi) lambda/37,1 = {D0:.0f} fm (rayon de tube lambda/38,1 "
          f"retrouve depuis g = 2, R54) ; w_e = {W_E:.1f} fm")
    print("  -> les deux pendent a l'ancre lambda_e (m_e) ; la brique n'apporte aucun nombre.\n")
    check("w/l1 identique a tous les barreaux", max(ratios) - min(ratios) < 1e-12,
          f"{ratios[0]:.5f}")

    # C. cellules du vide occupees
    print("C. Combien de cellules du vide ?")
    d_e = 2 * LAMBDA_BAR
    d_q = l1(9)
    print(f"  electron : diametre 2 lambda = {d_e:.0f} fm = {d_e/D0:.1f} cellules D0")
    print(f"  quark    : l1(9) = {d_q:.2f} fm = 1/{D0/d_q:.0f} cellule ; ruban w_9 = {W_E*d_q/L1_3:.3f} fm")
    print("  -> le quark est trois branches comprimees 300 fois dans une cellule du vide,")
    print("     ce que A autorise sans cout ; l'electron enjambe trois cellules.\n")
    check("electron > 1 cellule, quark < 1/100 cellule", d_e / D0 > 1 and d_q / D0 < 0.01,
          f"{d_e/D0:.1f} et 1/{D0/d_q:.0f}")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
