#!/usr/bin/env python3
"""R102 -- D'ou vient alors mu_B ?  Le moment de Bohr comme boucle de charge, contre la diffusion.

R100 a arrete la charge etalee sur l'anneau.  Mais dans la base, un moment magnetique est un
courant de charge qui tourne (R53 : mu = q c R/2 = mu_B pour q R = e lambda-bar).  Peut-on garder
mu_B par une boucle de charge plus petite ?

  A. la borne de Lorentz, exacte : mu = q v R/2 ; mu_B = e c lambda-bar/2 ; pour la charge e a
     v <= c, R >= lambda-bar = 386 fm.  R100 : F(q) = 1 a 1 % jusqu'a 1 GeV/c exige R < 0,05 fm ;
     LEP (q ~ 100 GeV/c) : R < 5e-4 fm.  Une boucle de charge e qui porte mu_B avec cette taille
     tourne a v/c >= 8000.  Une charge q > e n'existe pas dans l'electron (charge totale e).
  B. le facteur de forme magnetique de la boucle elle-meme : pour q perpendiculaire a l'axe,
     F_M = 2 J_1(qR)/(qR) ; la partie de rang 1 (la seule que porte un spin 1/2) est la moyenne sur
     les directions de q : 0,72 a 1 MeV/c, ~1e-6 a 1 GeV/c.  La section e-p de Rosenbluth a Q^2 = 1 GeV^2
     utilise le tenseur leptonique d'un point (G_M = 1) et est verifiee au %.
  C. la reference : Dirac ponctuel, repere de Breit, u(p')^dag alpha u(p) = i (sigma x q) exactement
     a tout Q^2 : G_M = 1, le moment mu_B sans boucle, porte par la structure spinorielle (Gordon :
     aimantation (e/2m) psi-bar Sigma psi, etendue = le paquet, pas lambda-bar), et g = 2 avec.
  D. ce que la base doit fournir : un moment mu_B qui n'est pas un courant de charge spatial ; la
     circulation d'energie a lambda-bar (masses, spin) n'est pas touchee si elle est neutre ; le
     moment, lui, n'a plus de source dans la base.
"""
import math
import numpy as np
from scipy.special import j1

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME
R_1GEV = 0.245 * HBARC / 1000.0        # |F - 1| < 1 % a 1 GeV/c : qR/hbar c < 0,245
R_LEP = 0.245 * HBARC / 1e5

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def FM_perp(x):
    return 1.0 if x < 1e-9 else 2 * j1(x) / x

def FM_rank1(q_mev, R_fm):
    """moyenne sur les directions de q (partie de rang 1 du courant de boucle)."""
    x0 = q_mev * R_fm / HBARC
    th = np.linspace(0, math.pi, 200001)                      # grille fine : l'integrande oscille (x0 ~ 2000 a 1 GeV/c)
    f = np.array([FM_perp(x0 * math.sin(t)) for t in th]) * np.sin(th) / 2
    return float(np.trapz(f, th))

def dirac_spinor(p, m, chi):
    E = math.sqrt(m * m + np.dot(p, p))
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
    sp = p[0] * sx + p[1] * sy + p[2] * sz
    return math.sqrt(E + m) * np.concatenate([chi, sp @ chi / (E + m)]), (sx, sy, sz)

def main():
    print("R102 -- mu_B comme boucle de charge, contre la diffusion\n")

    # A. borne de Lorentz
    print("A. La borne de Lorentz, exacte")
    print(f"  mu = q v R / 2 ; mu_B = e c lambda-bar / 2  =>  (q/e)(v/c) R = lambda-bar = {LAMBDA_BAR:.1f} fm")
    print(f"  charge e, v <= c : R >= lambda-bar. Diffusion (R100) : R < {R_1GEV:.3f} fm (1 GeV/c, 1 %), < {R_LEP:.0e} fm (LEP).")
    vc_1, vc_lep = LAMBDA_BAR / R_1GEV, LAMBDA_BAR / R_LEP
    print(f"  une boucle de charge e portant mu_B a cette taille : v/c = {vc_1:.0f} (1 GeV) a {vc_lep:.0e} (LEP).\n")
    check("mu_B par une boucle de charge e exige R >= lambda-bar ; la diffusion exige R < 0,05 fm : v/c >= 8000",
          vc_1 > 5000, f"v/c = {vc_1:.0f}")

    # B. facteur de forme magnetique de la boucle
    print("B. Le facteur de forme magnetique de la boucle a R = lambda-bar")
    res = {}
    for q in (1.0, 1000.0):
        x = q * LAMBDA_BAR / HBARC
        res[q] = (FM_perp(x), FM_rank1(q, LAMBDA_BAR))
        print(f"  q = {q:6.0f} MeV/c : q perpendiculaire a l'axe F_M = 2J_1(x)/x = {res[q][0]:+.4f} ; rang 1 (moyenne) = {res[q][1]:+.4f} ; |G_M|^2 = {res[q][1]**2:.1e}")
    print("  Rosenbluth e-p a Q^2 = 1 GeV^2 : tenseur leptonique d'un point (G_M = 1), verifie au %.\n")
    check("G_M de la boucle (rang 1) : 0,72 a 1 MeV/c, ~1e-6 a 1 GeV/c, contre 1 au %", abs(res[1.0][1]) < 0.9 and abs(res[1000.0][1]) < 1e-2,
          f"{res[1.0][1]:.2f}, {res[1000.0][1]:.0e}")

    # C. Dirac ponctuel : G_M = 1
    print("C. La reference : Dirac ponctuel, repere de Breit, courant spatial")
    okC = True
    chi = np.array([1, 0], complex)
    for Q2 in (1e-4, 1.0, 1e4, 1e6):
        Q = math.sqrt(Q2)
        qv = np.array([Q, 0, 0])
        u, (sx, sy, sz) = dirac_spinor(-qv / 2, ME, chi)
        up, _ = dirac_spinor(qv / 2, ME, chi)
        zero = np.zeros((2, 2)); alpha = [np.block([[zero, s], [s, zero]]) for s in (sx, sy, sz)]
        J = np.array([np.vdot(up, a @ u) for a in alpha])
        # attendu : J = i (sigma x q) ; q = Q x-hat, <up|sigma_z|up> = 1 -> J_y = i Q, J_x = J_z = 0
        GM = (J[1] / (1j * Q)).real
        okC &= abs(GM - 1) < 1e-9 and abs(J[0]) < 1e-9 and abs(J[2]) < 1e-9
        print(f"    Q^2 = {Q2:8.0e} MeV^2 : J_y / (i Q) = {GM:.12f}, J_x = J_z = 0")
    print("  => u^dag alpha u = i (sigma x q) a tout Q^2 : G_M = 1, mu_B sans boucle, porte par la structure spinorielle")
    print("     (Gordon : aimantation (e/2m) psi-bar Sigma psi, etendue = le paquet) ; g = 2 vient avec.\n")
    check("Dirac ponctuel : G_M = 1 exactement a tout Q^2 (moment sans courant spatial)", okC, "i (sigma x q)")

    # D. ce que la base doit fournir
    print("D. Ce que la base doit fournir")
    print("  la base n'a qu'une source de moment, le courant de charge (R53). mu_B et G_E = G_M = 1 ensemble exigent un")
    print("  moment qui n'est pas un courant spatial : la base n'en a pas. La circulation d'energie a lambda-bar (masses,")
    print("  spin) n'est pas touchee si elle est neutre ; le moment n'a plus de source. Decision de l'auteur.\n")

    print("Verdict : mu_B par une boucle de charge est EXCLU par la diffusion, quelle que soit la taille de la boucle")
    print("(R >= lambda-bar a v <= c contre R < 0,05 fm). L'electron mesure a G_E = G_M = 1 : charge et moment")
    print("ponctuels, taille lambda-bar cinematique (m/E, sigma x q). La base doit fournir un moment sans courant de")
    print("charge, ou adopter la structure de Dirac.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
