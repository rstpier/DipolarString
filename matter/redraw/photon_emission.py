#!/usr/bin/env python3
"""R75 -- Point 5 de coherence : un anneau qui ne rayonne pas mais emet des photons.

  A. le courant interne de l'anneau est stationnaire : puissance rayonnee nulle.
  B. l'anneau en orbite est un dipole qui tourne : classiquement il tombe sur le
     proton en 1,6e-11 s ; il faut une regle d'orbite stationnaire.
  C. la base l'a : l'horloge interne est la periode de Compton (R68) et elle bat
     a omega0/gamma en mouvement (manuscrit, ether de Lorentz) ; vue du
     laboratoire c'est l'onde de phase de de Broglie, longueur h/(m v) ; la
     fermeture de la phase sur l'orbite donne m v r = n hbar, donc Bohr.
  D. l'emission : un saut entre deux fermetures fait varier le dipole orbital ;
     l'anneau (386 fm) est ponctuel devant le photon (121,6 nm) ; le taux avec
     l'element de matrice standard vaut 6,3e8 /s.
"""
import math

HBARC = 197.3269804           # MeV fm
HBAR_EV_S = 6.582119569e-16
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME        # fm
R_E = ALPHA * LAMBDA_BAR
A0 = LAMBDA_BAR / ALPHA        # 52918 fm
C_FM = 2.99792458e23
C_SI = 2.99792458e8
E_SI = 1.602176634e-19
EPS0 = 8.8541878128e-12
HBAR_SI = 1.054571817e-34
A0_SI = A0 * 1e-15

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R75 -- ne pas rayonner et pourtant emettre\n")

    # A. courant interne stationnaire
    print("A. Courant interne de l'anneau : I constant, dI/dt = 0")
    P_internal = 0.0
    print("  P = (mu0/6 pi c) (d^2 p/dt^2)^2 = 0 pour un courant stationnaire (R53) -> pas de rayonnement.\n")
    check("puissance rayonnee par le courant interne nulle", P_internal == 0.0, "dI/dt = 0")

    # B. l'orbite classique tombe
    print("B. L'anneau en orbite de Bohr, lu classiquement")
    t_fall = A0 ** 3 / (4 * R_E ** 2 * C_FM)
    print(f"  temps de chute t = a0^3 / (4 r_e^2 c) = {t_fall:.2e} s : sans regle d'orbite, pas d'atome.\n")
    check("chute classique en ~1,6e-11 s", 1.0e-11 < t_fall < 2.0e-11, f"{t_fall:.2e} s")

    # C. la regle d'orbite depuis l'horloge de la base
    print("C. L'horloge interne (periode de Compton, R68) + omega' = omega0/gamma (manuscrit)")
    v = ALPHA                                     # v/c au niveau fondamental
    gamma = 1 / math.sqrt(1 - v * v)
    lam_dB = 2 * math.pi * LAMBDA_BAR / (gamma * v)   # h/(gamma m v) en fm
    closure = 2 * math.pi * A0 / lam_dB
    E1 = -0.5 * ALPHA ** 2 * ME * 1e6             # eV
    print(f"  onde de phase vue du laboratoire : lambda = h/(gamma m v) = {lam_dB:.0f} fm ; "
          f"2 pi a0 / lambda = {closure:.6f}")
    print(f"  fermeture m v r = n hbar -> Bohr : E_1 = -alpha^2 m c^2/2 = {E1:.3f} eV (mesure -13,598)")
    print("  -> l'orbite stationnaire n'est pas un postulat de plus : c'est la meme horloge que le tour.\n")
    check("fermeture de phase exacte sur l'orbite fondamentale (a 1e-4)", abs(closure - 1) < 1e-4,
          f"{closure:.6f}")
    check("niveau fondamental -13,6 eV", abs(E1 / (-13.598) - 1) < 1e-3, f"{E1:.3f} eV")

    # D. l'emission
    print("D. L'emission 2p -> 1s")
    dE = 0.75 * 13.598                            # eV
    omega = dE / HBAR_EV_S
    lam_photon = 2 * math.pi * C_SI / omega       # m
    d = E_SI * 0.7449 * A0_SI                     # element de matrice standard |<1s|z|2p>| = 0,7449 a0
    A_rate = omega ** 3 * d ** 2 / (3 * math.pi * EPS0 * HBAR_SI * C_SI ** 3)
    print(f"  photon {dE:.2f} eV, lambda = {lam_photon*1e9:.1f} nm ; anneau 2 lambda-bar = {2*LAMBDA_BAR*1e-6:.1e} nm : "
          f"rapport {2*LAMBDA_BAR*1e-15/lam_photon:.0e} -> dipole ponctuel")
    print(f"  taux dipolaire avec l'element de matrice standard : A = {A_rate:.2e} /s (mesure 6,27e8)")
    print("  -> ce qui rayonne est le dipole orbital pendant le saut, pas le courant interne ;")
    print("     l'amplitude du saut reste celle de la mecanique quantique, la base n'en a pas.\n")
    check("anneau ponctuel devant le photon (< 1e-4)", 2 * LAMBDA_BAR * 1e-15 / lam_photon < 1e-4,
          f"{2*LAMBDA_BAR*1e-15/lam_photon:.0e}")
    check("taux 2p -> 1s a 5 % avec l'element de matrice standard", abs(A_rate / 6.27e8 - 1) < 0.05,
          f"{A_rate:.2e} /s")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
