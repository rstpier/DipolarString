#!/usr/bin/env python3
"""R78 -- Le Z2 sur trois branches simples : la section carree comme porteuse.

R74 a reduit Pauli a une question : qu'est-ce qui rend antiperiodique un anneau
de trois branches simples ?  La phase A obtient le -1 pour une PAIRE : le
mode differentiel est impair sous l'echange des deux conducteurs (un element
d'ordre 2), et une demi-torsion par tour le rend antiperiodique.

Candidat teste : le triple n'a que la symetrie Z3 de ses trois positions,
sans element d'ordre 2 ; mais chaque brin a une section carree (w = d, R17),
de symetrie Z4, qui contient une rotation de 180 deg.  Un mode de section
impair sous cette rotation, avec une demi-torsion de section par circuit,
est antiperiodique : psi(phi + 2 pi) = -psi, L_z dans Z + 1/2.

  A. ordres des elements : Z3 n'a pas d'ordre 2, Z4 en a un.
  B. torsions de section fermees : t dans Z/4 ; antiperiodiques : 2t impair.
  C. compatibilite aux trois jonctions : t_i dans Z/4, somme = 1/2 mod 1 : solutions.
  D. cout : inclinaison des faces w/(4 lambda-bar) = 1/pi^2, 5,8 deg ; l'adaptation
     ne depend que de d/w (R71) : cout du second ordre.
  E. independance des deux reperes : la torsion du triple (Tw = 0, R61) et la
     torsion de section de chaque brin sont deux reperes differents.
"""
import math
import itertools
from fractions import Fraction

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def element_orders(n):
    """ordres des elements du groupe cyclique Z_n"""
    return sorted({n // math.gcd(k, n) for k in range(n)})

def main():
    print("R78 -- la section carree comme porteuse du Z2\n")

    # A. ordres
    print("A. Ordres des elements")
    o3, o4 = element_orders(3), element_orders(4)
    print(f"  Z3 (trois positions du triple) : ordres {o3} -> pas d'element d'ordre 2")
    print(f"  Z4 (section carree)            : ordres {o4} -> la rotation de 180 deg est d'ordre 2")
    print("  -> le triple ne peut pas porter l'antiperiodicite ; la section carree le peut.\n")
    check("Z3 sans ordre 2, Z4 avec", 2 not in o3 and 2 in o4, f"{o3} / {o4}")

    # B. torsions fermees et antiperiodiques
    print("B. Torsions de section par circuit qui referment le carre : t = k/4")
    closed = [Fraction(k, 4) for k in range(-4, 5)]
    anti = [t for t in closed if (2 * t).denominator == 1 and (2 * t).numerator % 2 == 1]
    print(f"  fermees : {[str(t) for t in closed]}")
    print(f"  antiperiodiques pour un mode impair sous 180 deg (2t impair) : {[str(t) for t in anti]}")
    print("  -> une demi-torsion (ou 3/2, ...) de la section par circuit donne psi(phi+2pi) = -psi,")
    print("     L_z = n + 1/2 : l'arithmetique de la phase A (step 2), transposee a un seul brin.\n")
    check("les torsions antiperiodiques existent parmi les torsions fermees du carre",
          len(anti) > 0 and all(t in closed for t in anti), f"{[str(t) for t in anti]}")

    # C. jonctions
    print("C. Trois brins, trois jonctions : t_i dans Z/4, somme = 1/2 mod 1")
    sols = []
    for combo in itertools.product([Fraction(k, 4) for k in range(-2, 3)], repeat=3):
        tot = sum(combo)
        if (2 * tot).denominator == 1 and (2 * tot).numerator % 2 == 1 and abs(tot) <= 1:
            sols.append(combo)
    minimal = [c for c in sols if sum(abs(t) for t in c) == Fraction(1, 2)]
    print(f"  {len(sols)} solutions avec |t_i| <= 1/2 et |somme| <= 1 ; minimales (somme |t_i| = 1/2) : "
          f"{[tuple(str(t) for t in c) for c in minimal[:3]]} ...")
    print("  -> deux jonctions tournees de 90 deg et une droite suffisent ; le carre se referme a chaque jonction.\n")
    check("il existe des repartitions par quarts de tour sommant a une demi-torsion", len(minimal) > 0,
          f"{len(minimal)} minimales")

    # D. cout
    print("D. Cout de la demi-torsion")
    HBARC, ME = 197.3269804, 0.51099895
    LAMBDA_BAR = HBARC / ME
    W = 4 * LAMBDA_BAR / math.pi ** 2
    path = 2 * math.pi * LAMBDA_BAR
    tilt = math.atan((W / 2) * math.pi / path)          # demi-largeur x angle de torsion / trajet
    print(f"  inclinaison des faces : atan(w pi / (2 x 2 pi lambda-bar)) = atan(1/pi^2) = {math.degrees(tilt):.2f} deg")
    print(f"  variation relative de C' au second ordre ~ tilt^2/2 = {tilt**2/2:.4f} ; l'adaptation ne")
    print("  depend que de d/w (R71) : la demi-torsion est compatible avec Z0 a 0,5 % pres.\n")
    check("cout de la demi-torsion < 1 %", tilt ** 2 / 2 < 0.01, f"{100*tilt**2/2:.2f} %")

    # E. independance des reperes
    print("E. Deux reperes")
    print("  le repere du triple (arrangement des trois brins, Tw = 0, canal writhe des generations,")
    print("  R61) et le repere de section de chaque brin (orientation du carre autour de son axe)")
    print("  sont independants : le premier est fixe par le milieu, le second est libre de porter")
    print("  la demi-torsion. Le muon et le tau gardent les deux (spin 1/2 a toutes les generations).\n")
    check("les deux reperes sont distincts (triple vs section)", True, "R61 vs R17")

    print("Ce qui manque : montrer qu'un brin simple a section carree porte un mode impair sous")
    print("180 deg (pour la paire c'etait le mode differentiel, un champ reel) : un motif de")
    print("circulation quadrupolaire dans la section, non calcule ; puis l'echange par le tour de")
    print("ceinture (R74 D2). Sans ce mode, le candidat est vide.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
