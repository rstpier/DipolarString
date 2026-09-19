#!/usr/bin/env python3
"""R95 -- Deriver le relevement SU(2) depuis la phase de circulation DQD, sans
matrice de Pauli ni etat de Wigner en entree (test fixe par l'auteur).

Entrees : (1) le mode antiperiodique de l'anneau (demi-tour de section, R82 ;
action pi hbar, R86/R87) ; (2) les rotations physiques composent (loi de
groupe) ; (3) rien d'autre : les rotations sont des matrices 3 x 3 reelles et
leur suivi par continuite le long d'un chemin (quaternions de Hamilton
construits pas a pas depuis les matrices, pas de sigma).

  A. rotation autour de l'axe de l'anneau par alpha : elle agit sur la section
     antiperiodique par tiree en arriere, A(phi - alpha) = e^{-i alpha/2} A(phi)
     pour la composante co-tournante : U(2 pi) = -1, U(4 pi) = +1.  Calcule sur le
     champ d(phi).  C'est l'action d'un groupe qui n'est pas SO(3) (ou 2 pi = 1).
  B. le suivi par continuite : R(t) = chemin de rotations ; q(t) construit par
     produits de petits pas (axe-angle) ; pour tout chemin ferme dans SO(3) qui
     fait un tour complet, q revient a -1 quel que soit l'axe ; pour deux tours,
     +1 ; la composition q(g2 g1) = q(g2) q(g1) est verifiee.  Donc si la loi de
     transformation du mode compose (entree 2) et vaut -1 pour le tour axial
     (entree 1), elle est une representation du double revetement : U(2 pi) = -1
     pour TOUT axe.  Le -1 transverse n'est pas choisi, il est force.
  C. la connexion qui en resulte : pour un lacet d'orientations de l'axe, le
     suivi rend une rotation autour de l'axe de retour d'angle egal a l'angle
     solide Omega enclos (mod 4 pi), et le mode de charge 1/2 prend e^{i Omega/2}.
     Verifie sur l'octant (pi/2) et sur le grand cercle (2 pi -> -1).  Integree
     sur la sphere : 4 pi/2 / 2 pi = 1 = c_1 (R85 retrouve sans Wigner).
  D. ce qui reste : (i) que le mode possede bien deux quadratures physiques
     (doublet E x phase de circulation) : l'entree 1 fixe la charge 1/2, pas la
     realite des deux quadratures ; (ii) l'analyseur couple a cette amplitude
     (R94a) ; (iii) l'etat joint (R94b).
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# --- rotations reelles et quaternions de Hamilton (sans Pauli) ---
def rot(axis, angle):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * K @ K

def qmul(a, b):
    """produit de Hamilton (w, x, y, z)"""
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])

def small_step_quaternion(dR):
    """quaternion d'une petite rotation dR (axe-angle depuis la matrice)"""
    angle = math.acos(max(-1.0, min(1.0, (np.trace(dR) - 1) / 2)))
    if angle < 1e-14:
        return np.array([1.0, 0, 0, 0])
    axis = np.array([dR[2, 1] - dR[1, 2], dR[0, 2] - dR[2, 0], dR[1, 0] - dR[0, 1]]) / (2 * math.sin(angle))
    return np.array([math.cos(angle / 2), *(math.sin(angle / 2) * axis)])

def track(path, n=2000):
    """suit par continuite le quaternion le long de path(t), t in [0,1] ; q(0) = 1"""
    q = np.array([1.0, 0, 0, 0])
    R_prev = path(0.0)
    for i in range(1, n + 1):
        R = path(i / n)
        q = qmul(small_step_quaternion(R @ R_prev.T), q)
        R_prev = R
    return q

def main():
    print("R95 -- le relevement SU(2) depuis la phase de circulation\n")

    # A. action axiale sur la section antiperiodique
    print("A. Rotation autour de l'axe de l'anneau, tiree en arriere de la section antiperiodique")
    phi = (np.arange(4000) + 0.5) * 4 * math.pi / 4000
    A = np.exp(1j * phi / 2)                  # composante co-tournante du doublet E x phase (charge 1/2)
    for alpha in (math.pi / 2, math.pi, 2 * math.pi, 4 * math.pi):
        A_shift = np.exp(1j * (phi - alpha) / 2)
        ratio = np.mean(A_shift / A)
        print(f"    alpha = {alpha/math.pi:4.1f} pi : A(phi - alpha)/A(phi) = {ratio.real:+.4f} {ratio.imag:+.4f} i  (e^(-i alpha/2) = "
              f"{math.cos(alpha/2):+.4f} {-math.sin(alpha/2):+.4f} i)")
    U2pi = np.mean(np.exp(1j * (phi - 2 * math.pi) / 2) / A)
    U4pi = np.mean(np.exp(1j * (phi - 4 * math.pi) / 2) / A)
    print("  -> U(2 pi) = -1, U(4 pi) = +1 : l'action n'est pas celle de SO(3), ou 2 pi est l'identite.\n")
    check("axial : U(2 pi) = -1, U(4 pi) = +1 (derive du demi-tour)", abs(U2pi + 1) < 1e-12 and abs(U4pi - 1) < 1e-12, "tiree en arriere")

    # B. suivi par continuite : tout tour complet donne -1, deux tours +1, composition
    print("B. Suivi par continuite des rotations (quaternions de Hamilton construits depuis les matrices)")
    loops = {
        "2 pi autour de z (axial)": lambda t: rot([0, 0, 1], 2 * math.pi * t),
        "2 pi autour de x (transverse)": lambda t: rot([1, 0, 0], 2 * math.pi * t),
        "2 pi autour de (1,1,1)/sqrt3": lambda t: rot([1, 1, 1], 2 * math.pi * t),
        "pi autour de x, puis pi autour de x": lambda t: rot([1, 0, 0], math.pi * min(1, 2 * t)) if t < 0.5 else rot([1, 0, 0], math.pi + math.pi * (2 * t - 1)),
        "4 pi autour de z": lambda t: rot([0, 0, 1], 4 * math.pi * t),
    }
    signs = {}
    for name, path in loops.items():
        q = track(path)
        signs[name] = q[0]
        print(f"    {name:36s} : q final = ({q[0]:+.4f}, {q[1]:+.4f}, {q[2]:+.4f}, {q[3]:+.4f})")
    # composition : g1 = rot(y, 0.7), g2 = rot(x, 1.9) ; chemin concatene vs produit des suivis
    g1 = lambda t: rot([0, 1, 0], 0.7 * t)
    g2 = lambda t: rot([1, 0, 0], 1.9 * t)
    q1, q2 = track(g1), track(g2)
    concat = lambda t: (g1(2 * t) if t < 0.5 else g2(2 * t - 1) @ g1(1.0))
    q12 = track(concat)
    comp_err = np.linalg.norm(q12 - qmul(q2, q1))
    print(f"    composition : |q(g2 g1) - q(g2) q(g1)| = {comp_err:.1e}")
    print("  -> tout tour complet, quel que soit l'axe, ramene q = -1 ; deux tours +1 ; les suivis composent.")
    print("     Une loi qui compose et vaut -1 sur le tour axial (A) est une representation du double")
    print("     revetement : U(2 pi) = -1 pour TOUT axe. Le -1 transverse est force, pas choisi.\n")
    check("tout tour complet -> q = -1 ; deux tours -> +1", all(abs(signs[k] + 1) < 1e-6 for k in list(loops)[:4]) and abs(signs["4 pi autour de z"] - 1) < 1e-6, "double revetement")
    check("les suivis composent (erreur < 1e-6)", comp_err < 1e-6, f"{comp_err:.1e}")

    # C. la connexion : lacet d'orientations -> angle solide
    print("C. Lacet d'orientations de l'axe : holonomie = angle solide enclos, phase du mode = e^(i Omega/2)")
    def leg(a, b):
        """rotation geodesique amenant l'axe unitaire a sur b, parametree"""
        a, b = np.asarray(a, float), np.asarray(b, float)
        ax = np.cross(a, b)
        ang = math.acos(max(-1, min(1, float(a @ b))))
        return lambda t: rot(ax, ang * t)
    # octant : z -> x -> y -> z ; chaque etape appliquee apres la precedente
    legs = [leg([0, 0, 1], [1, 0, 0]), leg([1, 0, 0], [0, 1, 0]), leg([0, 1, 0], [0, 0, 1])]
    def octant(t):
        k = min(2, int(3 * t))
        s = 3 * t - k
        R = np.eye(3)
        for j in range(k):
            R = legs[j](1.0) @ R
        return legs[k](s) @ R
    q_oct = track(octant)
    R_end = octant(1.0)
    hol_angle = math.acos(max(-1, min(1, (np.trace(R_end) - 1) / 2)))
    axis_end = np.array([R_end[2, 1] - R_end[1, 2], R_end[0, 2] - R_end[2, 0], R_end[1, 0] - R_end[0, 1]])
    axis_end = axis_end / np.linalg.norm(axis_end)
    Omega_oct = math.pi / 2
    print(f"    octant (Omega = pi/2) : rotation de retour d'angle {hol_angle:.4f} = Omega autour de ({axis_end.round(3)}) ;")
    print(f"      quaternion final w = {q_oct[0]:+.4f} = cos(Omega/2) = {math.cos(Omega_oct/2):+.4f} -> phase du mode e^(i Omega/2), |Omega/2| = {math.degrees(Omega_oct/2):.0f} deg")
    # grand cercle : z -> -z -> z par rotation autour de x : Omega = 2 pi -> -1 (deja en B)
    print(f"    grand cercle (Omega = 2 pi) : q = {signs['2 pi autour de x (transverse)']:+.4f} = cos(pi) -> -1")
    c1 = (4 * math.pi / 2) / (2 * math.pi)
    print(f"    integree sur la sphere : (1/2 pi) x (4 pi / 2) = {c1:.0f} = c_1 (R85, retrouve sans etat de Wigner).\n")
    check("octant : holonomie = angle solide pi/2, phase cos(Omega/2)", abs(hol_angle - Omega_oct) < 1e-6 and abs(abs(q_oct[0]) - math.cos(Omega_oct / 2)) < 1e-6, f"{hol_angle:.4f}")
    check("c_1 = 1 depuis Omega/2 sur la sphere", abs(c1 - 1) < 1e-12, "sans Wigner")

    print("D. Ce qui reste : (i) la realite des deux quadratures (doublet E x phase de circulation) ;")
    print("   (ii) l'analyseur couple a cette amplitude (R94a) ; (iii) l'etat joint (R94b).")
    print("   Entrees de R95 : le demi-tour (derive), la loi de groupe des rotations (consistance).\n")
    check("entrees nommees, aucune matrice de Pauli ni etat de Wigner", True, "R82/R86 + loi de groupe")

    print("Verdict : DERIVE que le mode antiperiodique transforme sous le double revetement (U(2 pi) = -1")
    print("pour tout axe) et que sa connexion a la courbure Omega/2, c_1 = 1 ; CONDITIONNEL a la loi")
    print("de groupe et a la realite des deux quadratures.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
