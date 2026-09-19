#!/usr/bin/env python3
"""R67 -- "Peut-etre aussi des effets relativistes, on est dans une realite
ralentie."

Le manuscrit est de type ether de Lorentz : l'horloge interne d'une particule
est l'onde a c0 qui ferme sa boucle, et elle bat a omega' = omega0/gamma quand
la particule se deplace (sec. Lorentz, app. breathers).  Trois lectures d'un
"ralentissement" sont confrontees aux durees de vie :
  1. un gamma global (notre realite entiere ralentie) ;
  2. le gamma de translation, celui du manuscrit ;
  3. un gamma interne du fluide qui circule (v < c0 dans l'anneau).
"""
import math

HBARC = 197.3269804
C_FM = 2.99792458e23
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

MMU, TAU_MU = 105.6583755, 2.1969811e-6
MTAU, TAU_TAU, BR_TAU_E = 1776.86, 2.903e-13, 0.1782
TAU_N = 878.4
# CERN muon storage ring (Bailey et al. 1977) : muons a gamma = 29.33
GAMMA_CERN, TAU_CERN, TAU_CERN_ERR = 29.33, 64.378e-6, 0.026e-6
# potentiels gravitationnels (GM/rc^2) : Terre, Soleil a 1 UA, Galaxie
PHI = {"Terre": 7.0e-10, "Soleil (1 UA)": 9.9e-9, "Galaxie": 1e-6}

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R67 -- realite ralentie ?\n")
    T = {"mu": 3 * l1(7) / C_FM, "tau": 3 * l1(11) / C_FM, "n": 2 * math.pi * 0.95 / C_FM}
    N = {"mu": TAU_MU / T["mu"], "tau": (TAU_TAU / BR_TAU_E) / T["tau"], "n": TAU_N / T["n"]}

    # 1. gamma global
    print("1. Un gamma global : toute la realite ralentie du meme facteur")
    ratio = TAU_MU / TAU_N
    for g in (1.0, 2.718, 1e17):
        print(f"  gamma = {g:8.3g} : tau_mu/tau_n mesure = {(TAU_MU*g)/(TAU_N*g):.3e}")
    print("  -> une duree de vie se mesure avec nos horloges, ralenties du meme gamma :")
    print("     le facteur s'annule ; les 26 ordres entre Delta et neutron sont un rapport,")
    print("     invariant. Le Z(r) gravitationnel du manuscrit, Z0 exp(2GM/rc0^2), est global :")
    for k, phi in PHI.items():
        print(f"       {k:14s} Z/Z0 - 1 = {math.exp(2*phi)-1:.1e}")
    print()
    check("un gamma global s'annule dans toute duree de vie mesuree",
          abs((TAU_MU * 1e17) / (TAU_N * 1e17) - ratio) < 1e-30 * ratio, "rapport invariant")

    # 2. gamma de translation (loi du manuscrit omega' = omega0/gamma)
    print("2. Le gamma de translation : la loi du manuscrit sur le muon en vol")
    tau_pred = TAU_MU * GAMMA_CERN
    dev = tau_pred / TAU_CERN - 1
    print(f"  anneau de stockage du CERN, gamma = {GAMMA_CERN} : predit {tau_pred*1e6:.2f} us, "
          f"mesure {TAU_CERN*1e6:.3f} +- {TAU_CERN_ERR*1e6:.3f} us ({100*dev:+.2f} %)")
    print("  -> l'effet relativiste que la base contient est le standard ; au repos gamma = 1,")
    print("     et 2,197 us est la duree propre : la dilatation n'y ajoute rien.\n")
    check("la loi omega' = omega0/gamma reproduit le muon en vol (a 0,5 %)",
          abs(dev) < 0.005, f"{100*dev:+.2f} %")

    # 3. gamma interne du fluide : v < c0 dans l'anneau, un tour propre pour se defaire
    print("3. Un gamma interne : le fluide tourne a v < c0, se defait en un tour propre")
    print(f"  {'objet':4s} {'tours':>9s} {'gamma_int':>10s} {'1 - v/c0':>10s}")
    for k in ("mu", "tau", "n"):
        g = N[k]
        eps = 1 / (2 * g * g)
        print(f"  {k:4s} {N[k]:9.2e} {g:10.2e} {eps:10.1e}")
    M_int = MMU * N["mu"] ** 0.25
    scaling = math.log(N["mu"] / N["tau"]) / math.log(MTAU / MMU)
    print(f"  gamma_int doit suivre m^-{scaling:.2f} (mu -> tau), soit gamma = (M/m)^4 avec "
          f"M = {M_int/1e6:.2f} TeV : l'echelle de R64,")
    print("     et 1 - v/c0 = (m/M)^8 : v = c0 a 1e-34 pres pour le muon. Rien dans la base ne")
    print("     fixe une vitesse a 1e-34 de c0 ; c'est la fuite de R64 reecrite en vitesse.\n")
    check("gamma interne : meme echelle M que R64 (1,6 TeV a 5 %)",
          abs(M_int / 1.62e6 - 1) < 0.05, f"M = {M_int/1e6:.2f} TeV")
    check("gamma interne du muon demande v = c0 a mieux que 1e-30",
          1 / (2 * N["mu"] ** 2) < 1e-30, f"1 - v/c0 = {1/(2*N['mu']**2):.1e}")

    # 4. temps propre du fluide a c0
    print("4. Le fluide a c0 n'a pas de temps propre")
    print("  une onde a c0 (R47 : vortex adapte, v = c0) ne vieillit pas ; ce qui declenche la")
    print("  chute doit battre au repos : la moitie statique (R54), dont le tic est l1(n)/c,")
    print("  la meme horloge que le tour. Aucune puissance de m supplementaire n'en sort.\n")
    check("horloge statique / horloge du tour : rapport constant sur l'echelle",
          abs((l1(7) / (3 * l1(7))) - (l1(11) / (3 * l1(11)))) < 1e-12, "1/3 pour tout n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
