#!/usr/bin/env python3
"""R89 -- Une prediction qui mord : la taille de l'anneau de l'electron dans
l'hydrogene ("Vois-tu une prediction ?").

Un electron etendu deplace le niveau 1S comme un noyau etendu :
    dE = (2 pi/3) alpha hbar c <r^2> |psi(0)|^2.
Dans la theorie de Dirac, le terme de Darwin (mesure, partie de la structure
fine) vaut (pi/2) alpha hbar c lambda-bar^2 |psi(0)|^2.  Si l'electron EST
l'anneau, la taille de l'anneau doit reproduire le terme de Darwin :
    <r^2> = (3/4) lambda-bar^2,  r_rms = (sqrt3/2) lambda-bar = 334,4 fm.
Sinon le 1S bouge de dizaines de GHz, alors qu'il est mesure au kHz.

  A. l'equivalence exacte Darwin <=> <r^2> = 3/4 lambda-bar^2.
  B. anneau fin a R = lambda-bar (R53) : <r^2> = 1,00 (+33 %).
  C. distribution de R84 (charge = courant sur une ligne adaptee, section
     carree, concentree vers l'interieur) : <r^2> = 0,883 (+18 %), 31 GHz sur 1S.
  D. charge aux trois jonctions (R5, R25) a r_J = D/sqrt3 : <r^2> = 0,12 (-84 %).
  E. sensibilite : le 1S est connu a ~10 kHz ; l'ecart de C est exclu par 10^6.
Caveat : |psi(0)|^2 est celui de la mecanique quantique ; l'atome de la base est
une orbite de Bohr (R75).  La base doit produire son propre calcul atomique ;
l'ordre de grandeur (meV contre kHz) dit qu'il devra fixer <r^2> a 3/4 a mieux
que 1e-5.
"""
import math
import numpy as np
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("cd", os.path.join(HERE, "core_displacement.py"))
cd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cd)

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
LAMBDA_BAR = HBARC / ME
W_E = 4 * LAMBDA_BAR / math.pi ** 2
D_J = 6 * LAMBDA_BAR / math.pi ** 2
EV_TO_GHZ = 2.417989e5
LAMB_1S_PRECISION_KHZ = 10.0

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def r84_distribution(n_side=40):
    pts, h = cd.perimeter_filaments(LAMBDA_BAR, W_E, n_side)
    N = len(pts)
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M[i, j] = cd.self_ind(pts[i, 0], h) if i == j else cd.mutual(pts[i, 0], pts[i, 1], pts[j, 0], pts[j, 1])
    A = np.zeros((N + 1, N + 1))
    A[:N, :N] = M
    A[:N, N] = -1.0
    A[N, :N] = 1.0
    sol = np.linalg.solve(A, np.r_[np.zeros(N), 1.0])
    I = sol[:N]
    return pts, I

def main():
    print("R89 -- l'anneau de l'electron dans l'hydrogene\n")
    E_D = 0.5 * ALPHA ** 4 * ME * 1e6             # terme de Darwin 1S, eV
    target = 0.75

    # A. equivalence
    print("A. Taille finie contre Darwin")
    print("  (2 pi/3) alpha hbar c <r^2> |psi(0)|^2 = (pi/2) alpha hbar c lambda-bar^2 |psi(0)|^2  <=>  <r^2> = 3/4 lambda-bar^2")
    print(f"  Darwin 1S = alpha^4 m c^2/2 = {E_D*1e3:.3f} meV = {E_D*EV_TO_GHZ:.1f} GHz ; r_rms requis = {math.sqrt(target)*LAMBDA_BAR:.1f} fm\n")
    check("Darwin <=> <r^2> = 3/4 lambda-bar^2 (algebre exacte)", abs((math.pi / 2) / (2 * math.pi / 3) - 0.75) < 1e-12, "3/4")

    # B, C, D : trois lectures de la charge
    readings = {}
    readings["anneau fin R = lambda-bar (R53)"] = 1.0
    pts, I = r84_distribution()
    readings["distribution R84 (courant = charge, section carree)"] = float(np.sum(I * (pts[:, 0] ** 2 + pts[:, 1] ** 2))) / LAMBDA_BAR ** 2
    r_J = D_J / math.sqrt(3)
    readings["trois charges aux jonctions, r_J = D/sqrt3 (R5, R25)"] = (r_J / LAMBDA_BAR) ** 2
    print("B-D. <r^2> / lambda-bar^2 selon la lecture de la charge, et decalage du 1S par rapport a Darwin")
    for name, r2 in readings.items():
        shift = E_D * (r2 / target - 1)
        print(f"  {name:55s} <r^2> = {r2:.3f} ({100*(r2/target-1):+.0f} %) -> {shift*EV_TO_GHZ:+.1f} GHz sur le 1S")
    r2_84 = readings["distribution R84 (courant = charge, section carree)"]
    print()
    check("anneau fin : +33 % (exclu)", abs(readings["anneau fin R = lambda-bar (R53)"] / target - 1 - 1 / 3) < 1e-9, "1.00 vs 0.75")
    check("distribution R84 : <r^2> = 0,88, +18 % (exclu au niveau du GHz)", abs(r2_84 - 0.883) < 0.005 and r2_84 / target - 1 > 0.1,
          f"{r2_84:.3f}")
    check("charges aux jonctions : -84 % (exclu)", readings["trois charges aux jonctions, r_J = D/sqrt3 (R5, R25)"] < 0.2, f"{r_J:.0f} fm")

    # E. sensibilite
    print("E. Sensibilite")
    shift_84_khz = E_D * (r2_84 / target - 1) * EV_TO_GHZ * 1e6
    print(f"  ecart de R84 : {shift_84_khz/1e6:.1f} GHz ; precision du 1S : ~{LAMB_1S_PRECISION_KHZ:.0f} kHz -> exclu par {shift_84_khz/LAMB_1S_PRECISION_KHZ:.0e}")
    tol = LAMB_1S_PRECISION_KHZ / (E_D * EV_TO_GHZ * 1e6)
    print(f"  pour tenir : <r^2>/(3/4 lambda-bar^2) = 1 a {tol:.0e} pres.")
    print("  caveat : |psi(0)|^2 est celui de la mecanique quantique ; l'atome de la base est une orbite de")
    print("  Bohr (R75) ; la base doit produire son propre calcul atomique. L'ordre de grandeur reste :")
    print("  une taille de 10^2 fm se voit au meV, l'hydrogene est mesure au neV.\n")
    check("l'ecart de R84 est exclu par > 1e5", shift_84_khz / LAMB_1S_PRECISION_KHZ > 1e5, f"{shift_84_khz/LAMB_1S_PRECISION_KHZ:.0e}")

    print("Verdict : la structure fait une prediction sans nombre libre, r_rms = (sqrt3/2) lambda-bar =")
    print("334 fm si l'anneau est le terme de Darwin ; aucune lecture actuelle de la charge ne la donne")
    print("(363 fm pour R84, +8,5 % en rayon). C'est le premier test qui mord ; il exige un calcul")
    print("atomique propre a la base, et il tranchera.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
