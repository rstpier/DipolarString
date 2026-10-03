#!/usr/bin/env python3
"""R116 -- L'exposant 2 pi de l'echelle peut-il sortir de l'adaptation a Z_0 ?

Enonce de l'auteur (3 octobre 2026) : « 2 pi, c'est une circonference r = 1 ; et 9, c'est
n_electron^2 ». Test propose : la longueur de circuit adaptee a Z_0 pour n brins varie-t-elle
comme (3/n)^(2 pi) ? L'echelle de la base : l_1(n) = l_1(3) (3/n)^p, p = 2 pi, epinglee par mu et tau.

  A. Ce que les donnees fixent : p(e, mu) = 6,293, p(e, tau) = 6,276 ; 2 pi = 6,283 est dans la
     fenetre (0,15 %) ; avec p = 2 pi exactement, m_mu/m_e sort a -0,75 % et m_tau/m_e a +0,97 %.
  B. Les brins en serie (la geometrie de la base : brins bout a bout) : n troncons adaptes mis en
     chaine ont Z_in = Z_0 pour tout n et toute longueur ; l'adaptation bifilaire fixe D/r = 2 cosh pi
     (le pi de la base) sans n ni longueur. L'adaptation ne produit aucune echelle en n.
  C. Les brins en parallele (faisceau de n brins sur un cercle, retour coaxial) : Z = (eta/2 pi)
     ln(b/r_eq), r_eq = (n r rho^(n-1))^(1/n) ; adapte a Z_0 : b/r_eq = e^(2 pi) = 535,5. Le 2 pi y est
     un rapport de tailles, pas un exposant de n ; la taille adaptee varie comme n^(p_eff) avec
     |p_eff| <= 1,1 pour n = 3 a 11 et rho/r <= 100, contre 2 pi.
  D. Les rapports adaptes de la base sont transverses : D/r = 2 cosh pi = e^pi (1 + e^(-2 pi)), et
     (2 cosh pi)^2 = e^(2 pi) a 0,4 % par identite ; aucun rapport de LONGUEURS de la base ne vaut
     e^(2 pi) a 5 % (le plus proche, lambda-bar/l_1(9), a -11 %). Le 2 pi de l'echelle est
     numeriquement celui du potentiel de ligne 2D, (lambda/2 pi eps) ln r, que la lecture de l'auteur
     nomme ; un mecanisme devrait produire ln(l_1(3)/l_1(n)) = 2 pi ln(n/3), c'est-a-dire une
     impedance eta ln(n/3) entre deux barreaux ; la base n'en a pas.
  Les deux lectures de l'auteur : 9 = 3 x 3 est deja le compte de la base (trois quarks, trois
  brins) ; « 2 pi = tour du cercle unite » est un nom tant qu'aucun calcul ne produit l'integrale
  sur un tour multipliant ln(n/3).
"""
import math
import numpy as np

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
ETA = 376.730313                 # Z_0 (ohm)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def abcd_section(beta_l, Z=1.0):
    return np.array([[math.cos(beta_l), 1j * Z * math.sin(beta_l)], [1j * math.sin(beta_l) / Z, math.cos(beta_l)]])

def z_in(M, ZL):
    (A, B), (C, D) = M
    return (A * ZL + B) / (C * ZL + D)

def main():
    print("R116 -- l'exposant 2 pi de l'echelle et l'adaptation a Z_0\n")

    # A. la fenetre des donnees
    print("A. Ce que mu et tau fixent")
    p_mu = math.log(MMU / ME) / math.log(7 / 3)
    p_tau = math.log(MTAU / ME) / math.log(11 / 3)
    r_mu = (7 / 3) ** (2 * math.pi) / (MMU / ME) - 1
    r_tau = (11 / 3) ** (2 * math.pi) / (MTAU / ME) - 1
    print(f"    p(e, mu) = {p_mu:.4f}, p(e, tau) = {p_tau:.4f}, 2 pi = {2*math.pi:.4f} ; fenetre {100*(p_mu-p_tau)/(2*math.pi):.2f} %, 2 pi dedans")
    print(f"    avec p = 2 pi exactement : m_mu/m_e {100*r_mu:+.2f} %, m_tau/m_e {100*r_tau:+.2f} % (le 0,15 % porte sur p, pas sur les masses)")
    print(f"    n = 3, 7, 11 sont poses ; un mecanisme doit donner p dans [{p_tau:.3f} ; {p_mu:.3f}].\n")
    check("A. 2 pi dans la fenetre de p (0,3 %) ; p = 2 pi met les masses a -0,75 % et +0,97 %",
          p_tau < 2 * math.pi < p_mu and abs(r_mu + 0.0075) < 0.001 and abs(r_tau - 0.0097) < 0.001, f"{p_tau:.4f} < {2*math.pi:.4f} < {p_mu:.4f}")

    # B. brins en serie
    print("B. Les brins en serie : n troncons adaptes en chaine")
    worst = 0.0
    for n in (3, 7, 9, 11):
        for beta_l in (0.3, 1.0, 2.0, math.pi / 2):
            M = np.eye(2, dtype=complex)
            for _ in range(n):
                M = M @ abcd_section(beta_l, 1.0)
            worst = max(worst, abs(z_in(M, 1.0) - 1.0))
    Dr = 2 * math.cosh(math.pi)
    print(f"    |Z_in/Z_0 - 1| <= {worst:.1e} pour n = 3, 7, 9, 11 et quatre longueurs electriques : la chaine est adaptee quel que soit n")
    print(f"    la section bifilaire adaptee : (eta/pi) arccosh(D/2r) = eta -> D/r = 2 cosh pi = {Dr:.3f} : ni n ni longueur")
    print("  -> l'adaptation ne contient pas n ; elle ne peut pas produire l_1(n) proportionnel a n^(-2 pi).\n")
    check("B. n troncons adaptes en serie : Z_in = Z_0 pour tout n (1e-12) ; D/r = 2 cosh pi sans n : pas d'echelle par l'adaptation",
          worst < 1e-12 and abs(Dr - 23.184) < 0.001, f"{worst:.1e}, D/r = {Dr:.3f}")

    # C. brins en parallele
    print("C. Les brins en parallele : faisceau de n brins (rayon r) sur un cercle de rayon rho, retour coaxial en b")
    print("    Z = (eta/2 pi) ln(b/r_eq), r_eq = (n r rho^(n-1))^(1/n) ; adapte : ln(b/r_eq) = 2 pi, b/r_eq = e^(2 pi) =", f"{math.exp(2*math.pi):.1f}")
    ns = np.array([3, 5, 7, 9, 11], dtype=float)
    peffs = []
    for rho_over_r in (3.0, 10.0, 23.14, 100.0):
        r = 1.0
        rho = rho_over_r
        req = (ns * r * rho ** (ns - 1)) ** (1 / ns)
        b = req * math.exp(2 * math.pi)
        p_eff = np.diff(np.log(b)) / np.diff(np.log(ns))
        peffs.append(np.max(np.abs(p_eff)))
        print(f"    rho/r = {rho_over_r:6.2f} : r_eq/rho = " + ", ".join(f"{x:.3f}" for x in req / rho) + f" ; p_eff(taille adaptee) = " + ", ".join(f"{x:+.3f}" for x in p_eff))
    print(f"    |p_eff| <= {max(peffs):.3f} partout (et decroissant en n), contre 2 pi = 6,283 : la taille adaptee d'un faisceau varie au plus comme n^1.")
    print("  -> ici le 2 pi est un rapport de tailles (b/r_eq = 535), pas un exposant du compte de brins.\n")
    check("C. faisceau : le 2 pi de l'adaptation est un rapport de tailles e^(2 pi) ; l'exposant effectif en n est <= 1,1, pas 2 pi",
          max(peffs) < 1.5 and max(peffs) < 2 * math.pi / 4, f"|p_eff| max {max(peffs):.3f}")

    # D. e^(2 pi) contre les rapports de la base
    print("D. Les rapports adaptes sont transverses ; les rapports de longueurs de la base contre e^(2 pi) = 535,5")
    ALPHA = 1 / 137.035999
    target = math.exp(2 * math.pi)
    ident = Dr ** 2 / target - 1
    print(f"    transverse : D/r = 2 cosh pi = e^pi (1 + e^(-2 pi)) = {Dr:.3f} ; (2 cosh pi)^2 / e^(2 pi) = 1 + 2 e^(-2 pi) + ... = {Dr**2/target:.5f} : identite, pas coincidence")
    ratios = {"3^(2 pi) = l_1(3)/l_1(9)": 3 ** (2 * math.pi), "l_1/w = pi^3/6": math.pi ** 3 / 6,
              "1/alpha = a_0/lambda-bar": 1 / ALPHA, "2 pi/alpha": 2 * math.pi / ALPHA,
              "lambda-bar/l_1(9)": 3 ** (2 * math.pi) * 3 / (2 * math.pi), "l_1(3)/w_e = pi^3/6": math.pi ** 3 / 6}
    near = {k: v / target - 1 for k, v in ratios.items()}
    for k, v in ratios.items():
        print(f"      {k:28s} {v:9.2f}  ({100*near[k]:+7.1f} %)")
    closest = min(near.items(), key=lambda kv: abs(kv[1]))
    print(f"    le plus proche : {closest[0]} a {100*abs(closest[1]):.0f} %. Aucun rapport de longueurs de la base ne vaut e^(2 pi).")
    print("    le 2 pi de l'echelle est numeriquement celui du potentiel de ligne 2D, (lambda/2 pi eps) ln r (la lecture de l'auteur) ;")
    print("    un mecanisme devrait donner ln(l_1(3)/l_1(n)) = 2 pi ln(n/3), soit une impedance eta ln(n/3) entre barreaux : la base n'en a pas.\n")
    check("D. (2 cosh pi)^2 = e^(2 pi) a 0,4 % est une identite ; aucun rapport de longueurs de la base n'est e^(2 pi) a 5 %",
          abs(ident - 2 * math.exp(-2 * math.pi)) < 1e-4 and abs(closest[1]) > 0.05, f"identite {ident:.4f} ; {closest[0]} a {100*abs(closest[1]):.0f} %")

    print("Verdict : NON. L'adaptation a Z_0 fixe des rapports transverses (pi dans arccosh, 2 pi dans ln) et jamais une")
    print("puissance du compte de brins ; en serie n n'y entre pas, en faisceau l'exposant effectif est <= 1,1. L'exposant 2 pi")
    print("reste REPARAMETRE ; « 2 pi = tour du cercle unite » et « 9 = 3^2 » sont des noms, enregistres comme tels.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
