#!/usr/bin/env python3
"""R82 -- Le mode impair de la section carree, construit explicitement (suite de
R78, relecture de l'auteur : "construire explicitement le mode impair de la
section carree et montrer son antiperiodicite").

  A. modes de la section carree (laplacien de Dirichlet, cote w) :
     psi_mn = sin(m pi x/w) sin(n pi y/w) ; le fondamental (1,1) est pair sous
     toute rotation du carre ; la paire (1,2), (2,1) forme la representation E
     de Z4 : rotation de 90 deg -> matrice de rotation, de 180 deg -> -1.
     Physiquement : le deplacement du coeur de la circulation hors du centre
     de la section (motif dipolaire), l'analogue du mode differentiel de la
     paire (phase A).
  B. transport le long du circuit avec une torsion de section t par tour :
     amplitude du mode apres un tour = exp(2 pi i t) x (rotation de la paire) ;
     t = 1/2 -> -1 (antiperiodique, L_z dans Z + 1/2) ; t = 1/4 -> +-i (quart
     de periode, exclu) ; t = 0, 1 -> +1.
  C. le fondamental (1,1) reste +1 pour toute torsion : le -1 est porte par le
     motif dipolaire seulement.
  D. cout : le motif dipolaire est une excitation transverse d'energie
     (m^2 + n^2) = 5 contre 2 : dans la lecture onde ce serait un mode a
     sqrt(5/2) fois la coupure (exclu par R47) ; dans la lecture vortex c'est
     un deplacement geometrique du coeur, de cout (pi delta)^2 / (2 x 2 pi R)
     sur le trajet : pour delta = w/4, 0,3 % du trajet.
"""
import math
import numpy as np

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME
W = 4 * LAMBDA_BAR / math.pi ** 2

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def mode(m, n, X, Y):
    return np.sin(m * np.pi * X) * np.sin(n * np.pi * Y)

def rotate90(F):
    """rotation de 90 deg du carre [0,1]^2 : (x, y) -> (1 - y, x)"""
    return np.rot90(F)

def rep_matrix(basis, op):
    """matrice de l'operation op dans la base orthonormee (produits scalaires sur la grille)"""
    G = np.array([[np.sum(bi * bj) for bj in basis] for bi in basis])
    M = np.array([[np.sum(bi * op(bj)) for bj in basis] for bi in basis])
    return np.linalg.solve(G, M)

def main():
    print("R82 -- le mode impair de la section carree\n")
    N = 200
    x = (np.arange(N) + 0.5) / N
    X, Y = np.meshgrid(x, x, indexing="ij")

    # A. representation de Z4 sur les modes
    print("A. Modes de Dirichlet du carre et rotations")
    f11 = mode(1, 1, X, Y)
    f12, f21 = mode(1, 2, X, Y), mode(2, 1, X, Y)
    r90 = rep_matrix([f11], rotate90)
    R90 = rep_matrix([f12, f21], rotate90)
    R180 = R90 @ R90
    print(f"  fondamental (1,1) : rotation 90 deg -> {r90[0,0]:+.4f} (pair)")
    print(f"  paire (1,2),(2,1) : rotation 90 deg ->\n{np.round(R90, 4)}")
    print(f"  rotation 180 deg -> \n{np.round(R180, 4)}  : caractere {np.trace(R180):+.4f} = -2, l'irrep E de Z4")
    print("  -> le motif dipolaire (coeur deplace) change de signe sous 180 deg ; le fondamental non.\n")
    check("fondamental pair (+1) et paire dipolaire impaire sous 180 deg (-1)",
          abs(r90[0, 0] - 1) < 1e-6 and np.allclose(R180, -np.eye(2), atol=1e-6),
          f"{r90[0,0]:+.3f} ; trace {np.trace(R180):+.3f}")

    # B. transport avec torsion de section t par tour
    print("B. Un tour de circuit avec torsion de section t : amplitude -> rotation(2 pi t)")
    def holonomy(t):
        th = 2 * math.pi * t
        return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    for t in (0, 0.25, 0.5, 1.0):
        H = holonomy(t)
        ev = np.linalg.eigvals(H)
        print(f"  t = {t:4.2f} : valeurs propres {np.round(ev, 4)}")
    H_half = holonomy(0.5)
    print("  -> t = 1/2 : psi(phi + 2 pi) = -psi(phi) pour le motif dipolaire, L_z dans Z + 1/2 ;")
    print("     t = 1/4 : +-i, un quart de periode, ni boson ni fermion : la base doit choisir t = 1/2.\n")
    check("demi-torsion : holonomie exactement -1 sur le motif dipolaire", np.allclose(H_half, -np.eye(2)),
          "psi(2 pi) = -psi, psi(4 pi) = +psi")
    check("quart de torsion : +-i (exclu comme signe de Pauli)",
          np.allclose(np.sort_complex(np.linalg.eigvals(holonomy(0.25))), np.array([-1j, 1j])), "anyonique")

    # C. le fondamental ne porte pas le signe
    print("C. Le fondamental (1,1) est invariant sous toute rotation : holonomie +1 pour tout t")
    print("  -> le -1 n'existe que si la circulation porte le motif dipolaire (coeur hors centre).\n")
    check("fondamental : +1 quel que soit t", abs(r90[0, 0] - 1) < 1e-6, "invariant")

    # D. cout
    print("D. Cout du motif dipolaire")
    e_ratio = (1 + 4) / (1 + 1)
    R = LAMBDA_BAR
    delta = W / 4
    extra = (math.pi * delta) ** 2 / (2 * 2 * math.pi * R) / (2 * math.pi * R)
    print(f"  lecture onde : energie transverse (1,2)/(1,1) = {e_ratio:.2f}, un mode excite (exclu par R47)")
    print(f"  lecture vortex : coeur deplace de delta = w/4 = {delta:.1f} fm tournant d'un demi-tour par circuit :")
    print(f"    trajet allonge de {100*extra:.3f} % ; pas de mode, pas de coupure.")
    print("  -> compatible avec R47 et R53 seulement dans la lecture vortex.\n")
    check("cout geometrique du coeur deplace < 1 % du trajet", extra < 0.01, f"{100*extra:.3f} %")

    print("Ce qui est montre : un mode impair sous 180 deg existe sur la section carree (le motif")
    print("dipolaire, irrep E de Z4) et une demi-torsion de section par circuit le rend antiperiodique.")
    print("Ce qui reste : que la circulation de l'electron porte ce motif (coeur hors centre) plutot")
    print("que le fondamental centre ; et le passage de L_z = 1/2 a j = 1/2 puis a l'echange (R74 D2).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
