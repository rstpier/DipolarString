#!/usr/bin/env python3
"""R101 -- L'hydrogene repris avec le facteur de forme de R100.

R100 : la charge du modele suit l'anneau, F(q) = <j_0(q rho)>, <r^2> = 0,12 a 1,0 lambda-bar^2.
Que voit l'hydrogene ?

  A. l'atome de la base (orbite de Bohr-Sommerfeld, R75, R99) : le perihelie de l'orbite (n, k)
     vaut r_min = n a0 (n - sqrt(n^2 - k^2)) ; le plus petit sur tous les etats est a0/2 = 26 459 fm,
     soit 68 fois R.  Theoreme de la coquille : le potentiel de toute distribution a symetrie
     spherique (moyenne d'orientation, R99 C) est exactement -K/r hors d'elle.  L'atome de la base
     ne voit pas l'anneau : decalage nul, quelle que soit la lecture.  L'hydrogene n'est pas le
     test dans la base ; R100 (diffusion) l'est.
  B. avec la fonction d'onde 1S importee (mecanique quantique, le caveat de R89) : le decalage
     exact dE = int |psi|^2 [V_rho - V_point] d^3r pour les trois lectures, compare a
     (2 pi/3) K <r^2> |psi(0)|^2 (developpement en <r^2>) et a ce que l'experience laisse au-dessus
     de Dirac sur le 1S : le deplacement de Lamb, 8,17 GHz (que la base n'a pas non plus).
  C. la cible de R89 retiree : le terme de Darwin est la normalisation m/E d'un point (R100 C),
     pas une taille ; la borne de l'hydrogene, si tout le Lamb etait une taille : <r^2> <= 0,035
     lambda-bar^2, r_rms <= 72 fm ; la borne de diffusion (R100) : 1e-3 fm.  L'anneau a 363-386 fm
     est hors des deux.
"""
import math
import os
import importlib.util
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
ff = load("ring_form_factor")

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
K = ALPHA * HBARC                       # MeV fm
LAMBDA_BAR = HBARC / ME
A0 = LAMBDA_BAR / ALPHA                 # 52 918 fm
EV_TO_GHZ = 2.417989e5
LAMB_1S_GHZ = 8.1729                    # deplacement de Lamb du 1S (mesure), l'exces sur Dirac

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def psi2_1s(r):
    return math.exp(-2 * r / A0) / (math.pi * A0 ** 3)

def V_shell(r, rho):                   # potentiel d'une coquille de charge unite au rayon rho, en unites de K
    return -1.0 / max(r, rho)

def main():
    print("R101 -- l'hydrogene avec le facteur de forme de l'anneau\n")
    rd = ff.readings()

    # A. l'atome de la base
    print("A. L'atome de la base : perihelie des orbites de Sommerfeld et theoreme de la coquille")
    rmin = {}
    for n in range(1, 8):
        for k in range(1, n + 1):
            rmin[(n, k)] = n * A0 * (n - math.sqrt(n * n - k * k))
    (nm, km), r_closest = min(rmin.items(), key=lambda kv: kv[1])
    r_lim = A0 / 2
    print(f"    r_min(1,1) = {rmin[(1,1)]:.0f} fm (cercle) ; r_min(2,1) = {rmin[(2,1)]:.0f} fm ; plus petit calcule (n = {nm}, k = {km}) : {r_closest:.0f} fm ;")
    print(f"    limite k = 1, n -> inf : a0/2 = {r_lim:.0f} fm = {r_lim/LAMBDA_BAR:.1f} R")
    okA = True
    for R in rd:
        Vd = sum(w * V_shell(r_lim, rho) for w, rho in zip(R.w, R.rho))
        dV = Vd - (-1.0 / r_lim)
        okA &= abs(dV) < 1e-15
        print(f"    {R.name:44s} V(a0/2) - V_point = {dV:.1e} K/fm : nul (coquille)")
    print("  => dans l'atome de la base, aucun etat ne penetre l'anneau : decalage exactement nul. L'hydrogene")
    print("     de la base est aveugle a la taille de la charge ; ce n'est pas lui qui tranche, c'est R100.\n")
    check("atome de la base : r_min >= a0/2 = 68 R pour tout (n, k) ; decalage nul par la coquille", okA and r_closest >= r_lim, f"{r_closest:.0f} fm")

    # B. avec |psi_1S|^2 importe
    print("B. Avec la fonction d'onde 1S importee (le caveat de R89)")
    okB = True
    shifts = {}
    for R in rd:
        dE = 0.0
        for w, rho in zip(R.w, R.rho):
            if rho <= 0:
                continue
            val, _ = quad(lambda r: 4 * math.pi * r * r * psi2_1s(r) * (1 / r - 1 / rho), 0, rho)
            dE += w * K * val
        dE_ghz = dE * 1e6 * EV_TO_GHZ
        dE_exp = (2 * math.pi / 3) * K * R.r2() * psi2_1s(0) * 1e6 * EV_TO_GHZ
        shifts[R.name] = dE_ghz
        okB &= abs(dE_ghz / dE_exp - 1) < 0.03
        print(f"    {R.name:44s} dE(1S) exact = {dE_ghz:+.1f} GHz ; developpement <r^2> = {dE_exp:+.1f} GHz ; / Lamb = {dE_ghz/LAMB_1S_GHZ:.1f}")
    print(f"  ce que l'experience laisse au-dessus de Dirac sur le 1S : le deplacement de Lamb, {LAMB_1S_GHZ:.2f} GHz (QED, absent de la base).")
    print("  meme en lui attribuant tout le Lamb, l'anneau (1) et (2) deplacent le 1S 25 a 28 fois trop ; les jonctions (3) 4 fois.\n")
    check("1S importe : decalage exact = developpement en <r^2> a 3 % ; anneau > 20 x le Lamb entier",
          okB and shifts[rd[0].name] / LAMB_1S_GHZ > 20 and shifts[rd[1].name] / LAMB_1S_GHZ > 20,
          f"{shifts[rd[0].name]:.0f}, {shifts[rd[1].name]:.0f} GHz")

    # C. la cible de R89 retiree ; bornes
    print("C. La cible de R89 retiree")
    dE_per_r2 = (2 * math.pi / 3) * K * LAMBDA_BAR ** 2 * psi2_1s(0) * 1e6 * EV_TO_GHZ    # GHz par lambda-bar^2
    r2_max = LAMB_1S_GHZ / dE_per_r2
    r_max = math.sqrt(r2_max) * LAMBDA_BAR
    print(f"    decalage par unite de <r^2> : {dE_per_r2:.1f} GHz / lambda-bar^2 ; si tout le Lamb etait une taille : <r^2> <= {r2_max:.3f} lambda-bar^2,")
    print(f"    r_rms <= {r_max:.0f} fm ; borne de diffusion (R100, LEP) : {ff.R_E_BOUND_FM:.0e} fm. L'anneau : 363 a 386 fm.")
    print("    Le terme de Darwin (3/4 lambda-bar^2) est la normalisation m/E d'un point (R100 C) : la cible 334 fm")
    print("    de R89 n'etait pas un rayon de charge ; l'ecart de 15 % a 18 % n'etait pas la question.\n")
    check("borne de l'hydrogene (tout le Lamb) : r_rms <= 72 fm ; l'anneau a 363-386 fm est hors de la borne par > 5 en rayon",
          abs(r_max - 72) < 3 and 363 / r_max > 5, f"{r_max:.0f} fm")

    print("Verdict : dans l'atome de la base, l'hydrogene ne voit pas l'anneau (perihelie a 68 R, coquille) ;")
    print("avec la fonction d'onde importee, l'anneau deplace le 1S de +200 a +230 GHz contre 8 GHz d'exces")
    print("mesure sur Dirac. L'hydrogene ne tranche qu'en important la mecanique quantique ; R100 tranche sans")
    print("elle, et dans le meme sens : la charge electromagnetique de l'electron n'est pas etalee sur 386 fm.")
    print("La cible de R89 (334 fm) est retiree.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
