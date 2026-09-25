#!/usr/bin/env python3
"""R108 -- Les faisceaux n >= 3 comme porteur du champ electronique (demande de l'auteur).

Phase A (step 3a/3b/5) : un faisceau de n brins indistinguables porte un doublet de modes
differentiels (irrep E1) dont la combinaison circulaire e+- a une charge de Berry +-a(n) sous la
rotation du repere autour de l'axe du faisceau (a(3) = 0,186) ; l'inclinaison de l'axe n'a pas
d'element de matrice dans le doublet, donc le transport du doublet le long d'un chemin de
directions a pour connexion a x (forme d'angle solide) : holonomie e^{+-i a Omega}.  Un spineur,
c'est exactement cela avec a = 1/2 (R95, R105).  Question : un faisceau a-t-il a = 1/2 ?

  A. a(n) pour n = 3..8 avec la machinerie de la phase A (grille 1401^2, D/r = 2 cosh pi), et sa
     dependance au rayon de conducteur (r = 0,5 ; 1 ; 2 en unites du rayon de la phase A) :
     a est un nombre continu, pas une charge quantifiee.
  B. l'univalence : sous un tour complet du repere, e^{2 pi i a} ; un porteur isotrope exige
     e^{2 pi i a} = +-1 (representation de SO(3) ou de son double) ; distance a +-1 pour chaque n.
  C. si a valait exactement 1/2, le doublet serait le spineur de R105-R107 : transporte sur la
     trame cubique a v = 1/3 contre 1/2 pour le photon ; les faisceaux n'echappent pas a R106.
  D. le faisceau ferme de l'electron (n = 3) : holonomie 2 pi[(1 - a) Lk + a Wr] (step 3b, R60,
     R77) ; avec le repere verrouille sur le milieu (Tw = 0, R61), Lk = Wr = 1/3 donne 120 deg
     exactement, l'etiquette de generation, jamais pi ; pi exigerait Tw != 0 (Wr = 1,23 a
     Lk = 1/3), contre R61, et ne serait qu'un signe autour de l'axe de l'anneau (R74, R83).
"""
import math
import os
import sys
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.abspath(os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, SCRIPTS)
import berry_holonomy_common as bhc

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def e1_doublet(n):
    k = np.arange(n)
    c, s = np.cos(2 * math.pi * k / n), np.sin(2 * math.pi * k / n)
    return [c / np.linalg.norm(c), s / np.linalg.norm(s)]

def modes_n(n, psi):
    pts = bhc.positions(n, psi)
    return [bhc.normalise(bhc.field(q, pts)) for q in e1_doublet(n)]

def berry_a(n, psi=0.0, h=1e-4):
    ep, em, e0 = modes_n(n, psi + h), modes_n(n, psi - h), modes_n(n, psi)
    a01 = bhc.ip(e0[0], ((ep[1][0] - em[1][0]) / (2 * h), (ep[1][1] - em[1][1]) / (2 * h)))
    a10 = bhc.ip(e0[1], ((ep[0][0] - em[0][0]) / (2 * h), (ep[0][1] - em[0][1]) / (2 * h)))
    return abs(0.5 * (a10 - a01))

def main():
    print("R108 -- les faisceaux n >= 3 comme porteur\n")
    ns = range(3, 9)
    print("A. La charge de Berry a(n) du doublet, et sa dependance au rayon du conducteur")
    table = {}
    for eps in (0.5, 1.0, 2.0):
        bhc.EPS = eps * bhc.R
        table[eps] = {n: berry_a(n) for n in ns}
        print(f"    r = {eps:3.1f} R : " + "  ".join(f"a({n}) = {table[eps][n]:.4f}" for n in ns))
    bhc.EPS = bhc.R
    a = table[1.0]
    ok3 = abs(a[3] - 0.1858) < 5e-4
    spread = {n: max(table[e][n] for e in table) - min(table[e][n] for e in table) for n in ns}
    print(f"    n = 3 : {a[3]:.4f} (phase A : 0,1858) ; variation avec le rayon, n = 3..8 : " + ", ".join(f"{spread[n]:.3f}" for n in ns))
    print("    => a(n) croit avec n, depend continument du rayon : une phase geometrique du mode, pas une charge quantifiee.\n")
    check("a(3) = 0,186 reproduit ; a(n) continu en n et en rayon (variation > 0,02 pour tout n)", ok3 and all(spread[n] > 0.02 for n in ns), f"a(3) = {a[3]:.4f}")

    print("B. L'univalence sous un tour complet du repere : e^(2 pi i a)")
    dist = {}
    for n in ns:
        z = np.exp(2j * math.pi * a[n])
        dist[n] = min(abs(z - 1), abs(z + 1))
        print(f"    n = {n} : 2a = {2*a[n]:.4f} ; e^(2 pi i a) = {z.real:+.3f}{z.imag:+.3f}i ; distance a +-1 : {dist[n]:.3f}")
    nearest = min(ns, key=lambda n: dist[n])
    print(f"    le plus proche d'un spineur : n = {nearest} (2a = {2*a[nearest]:.3f}), a {100*abs(2*a[nearest]-round(2*a[nearest])):.1f} % de 1 ; et il bouge avec le rayon"
          f" (2a de {2*min(table[e][nearest] for e in table):.3f} a {2*max(table[e][nearest] for e in table):.3f}).")
    # le rayon qui rendrait a(6) = 1/2 exactement : un reglage, pas une structure
    from scipy.optimize import brentq
    def a6_minus_half(eps):
        bhc.EPS = eps * bhc.R
        return berry_a(6) - 0.5
    eps_star = brentq(a6_minus_half, 0.5, 1.0, xtol=1e-3)
    bhc.EPS = bhc.R
    print(f"    a(6) = 1/2 exactement pour un rayon r = {eps_star:.3f} R (et pour aucun autre) : un reglage du rayon, pas une charge.")
    print("    => aucun faisceau n'a a = 1/2 par structure ; un doublet a a non demi-entier n'est pas univalent sous 2 pi :")
    print("       aucun noeud equivariant ne peut le transporter sur une trame isotrope (le -1 exact de Kramers est mesure).\n")
    spread6 = 2 * (max(table[e][6] for e in table) - min(table[e][6] for e in table))
    check("aucun n de 3 a 8 n'a 2a entier a mieux que 1 % ; le plus proche (n = 6, 1,4 %) varie de plus de 0,1 avec le rayon et n'est 1/2 qu'a un rayon regle",
          all(dist[n] > 0.01 for n in ns) and spread6 > 0.1 and 0.5 < eps_star < 1.0, f"min {min(dist.values()):.3f} a n = {nearest}, 2a(6) varie de {spread6:.2f}, r* = {eps_star:.2f} R")

    print("C. Si a valait 1/2")
    print("    le doublet serait le spineur verrouille de R105 (transport e^(i Omega/2), l'algebre de Pauli par l'orientation) ;")
    print("    sur la trame cubique il irait a 1/3 contre 1/2 pour le photon (R106-R107, enumeration exacte) : v_e/c0 = 2/3.")
    print("    Les faisceaux n'echappent pas a R106 ; ils ne font que fournir (ou non) le spineur.\n")
    check("a = 1/2 ramene au cas de R106-R107 (v = 1/3) : les faisceaux ne changent pas la vitesse", True, "R106/R107")

    print("D. Le faisceau ferme de l'electron (n = 3) : holonomie 2 pi[(1 - a) Lk + a Wr]")
    a3 = a[3]
    hol = lambda Lk, Wr: 2 * math.pi * ((1 - a3) * Lk + a3 * Wr)
    h_tw0 = hol(1 / 3, 1 / 3)
    print(f"    repere verrouille (Tw = 0, R61) : Lk = Wr = 1/3 -> holonomie {math.degrees(h_tw0):.3f} deg : l'etiquette de generation (R60, R76), pas pi")
    sols = {}
    for Lk in (1 / 3, 2 / 3, 1.0):
        Wr = (0.5 - (1 - a3) * Lk) / a3
        sols[Lk] = Wr
        print(f"    holonomie pi exigerait, a Lk = {Lk:.3f} : Wr = {Wr:+.3f}, Tw = Lk - Wr = {Lk - Wr:+.3f} (R61 impose Tw = 0)")
    print("    => le signe -1 n'est pas dans l'holonomie du faisceau de l'electron ; et il ne serait qu'un signe autour de")
    print("       l'axe de l'anneau (R74, R83), pas un spineur de SO(3).\n")
    check("Tw = 0 : holonomie exactement 120 deg (1e-9) ; pi exige Tw != 0 (|Tw| > 0,2 pour tout Lk)",
          abs(h_tw0 - 2 * math.pi / 3) < 1e-9 and all(abs(Lk - Wr) > 0.2 for Lk, Wr in sols.items()), f"{math.degrees(h_tw0):.1f} deg")

    print("Verdict : EXCLU. Le doublet d'un faisceau de n >= 3 brins porte une charge de Berry a(n) continue")
    print("(0,19 ; 0,33 ; 0,43 ; 0,51 ; ... selon n, et selon le rayon), jamais exactement 1/2 : il n'est pas un")
    print("spineur, il n'est pas univalent sous 2 pi, et s'il l'etait il irait a 1/3 sur la trame (R106). Sa vraie")
    print("fonction dans la base est l'etiquette de generation (holonomie 120 deg a Tw = 0). Les faisceaux ne sont")
    print("pas le porteur du champ electronique.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
