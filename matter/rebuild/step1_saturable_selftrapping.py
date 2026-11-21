#!/usr/bin/env python3
"""Etape 1 de la reconstruction (R122) -- L'auto-piegeage saturable : seuil, echelle, loi de masse, bornes.

Idee de l'auteur (20-21 novembre 2026) : l'energie monte l'impedance Z la ou elle est ; les ondes se courbent vers
l'energie ; avec assez de quanta au meme endroit la bosse se referme (confinement) ; les etats mal fermes fuient
(particules ephemeres). Ici on l'ecrit et on la calcule, sans la base : Z(u) = Z0 (1 + (u/u_s)/(1 + u/u_sat)),
c(u) = c0 / sqrt(1 + ...), un milieu de Kerr saturable. L'enveloppe d'une onde porteuse (k, omega0 = c0 k) obeit a
l'equation de Schrodinger non lineaire saturable ; les etats stationnaires localises sont les lumps.

  Unites. Longueur : lambda_bar = 1/k ; amplitude : f^2 = u/u_s ; saturation : x_s = u_sat/u_s. Equation radiale
  stationnaire sans dimension, f(rho) e^{i mu t} : f'' + (d-1)/rho f' = mu f - g(f^2) f, g(s) = s/(1 + s/x_s).
  Norme N~ = int f^2 d^d rho. Retour aux unites : nombre de quanta N_q = N~ u_s lambda_bar^3 / (hbar omega0)
  = N~ (E_s/E_c)^4, masse M c^2 = N_q E_c = N~ E_s^4 / E_c^3, avec E_c = hbar omega0 l'energie d'un quantum et
  E_s = (u_s (hbar c)^3)^(1/4) l'echelle d'energie de la saturation. A N_q quanta et au seuil :
  M c^2 = N~_min^(1/4) E_s N_q^(3/4) -- la loi de masse du sac (MIT : M ~ B^(1/4) N^(3/4), B = u_s).

  A. d = 1 : soliton exact, N~ = 4 sqrt(mu) : pas de seuil (tout paquet se piege). d = 2, cubique : seuil de Townes,
     N~_c = 11,70 (controle numerique de la methode). d = 3, cubique : N~(mu) decroit, aucun lump stable (effondrement,
     Derrick) ; saturable : N~(mu) a un minimum N~_min(x_s) et une branche stable (dN~/dmu > 0, Vakhitov-Kolokolov).
  B. Le seuil : N~_min = 158 x_s^(-1/2), atteint a mu/x_s = 0,049 (loi d'echelle exacte : f = sqrt(x_s) F ramene l'equation
     a une forme universelle ; verifie pour x_s = 0,3 ; 1 ; 3 ; 10). « Deux quanta ne suffisent pas » (auteur) :
     N_q >= N~_min (E_s/E_c)^4, et pour un lump classique N_q >> 1 il faut E_c << E_s.
  C. Les echelles. u_0 = h c / l_1^4 = 4,6e23 J/m^3 (substrat au repos) -> E_s = 0,39 MeV ; borne HIBEF (R119) :
     u_s > 5,7e26 J/m^3 -> E_s > 2,3 MeV ; Euler-Heisenberg : u_EH = 1,5e29 -> 9,2 MeV ; Schwinger : u_S = eps0 E_S^2/2
     = 7,7e24 J/m^3 ; sac MIT : B^(1/4) = 145 MeV -> 9e33 J/m^3.
  D. Ce que la loi de masse donne (x_s = 1, N~_min^(1/4) = 3,54). Le lump classique le plus leger (N_q ~ 10) pese 20 E_s :
     46 MeV a la borne HIBEF, 90 fois l'electron ; l'electron (N_q = 1, 0,511 MeV) exigerait E_s = 0,14 MeV, soit
     u_s = 9e21 J/m^3, 6e4 fois sous la borne HIBEF : EXCLU. Pour un nucleon de 3 quanta : E_s = 116 MeV (sac MIT :
     B^(1/4) = 145-235 MeV), u_s = 4e33 J/m^3 = 8e9 u_0 : l'echelle du substrat (l_1) n'entre pas ; seconde echelle, posee.
  E. Le vrai vide. La borne HIBEF met u_s au-dessus de u_S (Schwinger) d'un facteur 74 : toute non-linearite assez forte
     pour pieger agit la ou le vide reel cree deja des paires ; la « saturation » du vide QED est la creation de paires,
     c'est-a-dire l'electron lui-meme. Un lump de quanta sans masse qui peuvent sortir est metastable (fuite lente,
     oscillon) : la stabilite exige une charge conservee (topologie, boucle fermee de R120) ou un exterieur interdit
     (confinement, sac).
  Verdict : le mecanisme est juste (lumps stables au-dessus d'un seuil, loi M ~ E_s N^(3/4)) ; il est EXCLU comme origine
  de l'electron (HIBEF + Schwinger) ; il est le modele du sac pour les hadrons, CONDITIONNEL a une echelle E_s ~ 150-200
  MeV posee, sans lien avec l_1. L'etape 1 ne fournit pas d'objet nouveau ; elle dit ce que l'objet doit avoir.
"""
import math

HBARC_MEVFM = 197.3269804
HBARC_J_M = 3.16152677e-26        # J m
HBAR, H, C, EPS0 = 1.054571817e-34, 6.62607015e-34, 2.99792458e8, 8.8541878128e-12
ME_KG, E_CH, ALPHA = 9.1093837015e-31, 1.602176634e-19, 1 / 137.035999
MEV_J = 1.602176634e-13
LAM_E = HBARC_MEVFM / 0.51099895          # fm
L1 = 2 * math.pi * LAM_E / 3              # fm, l_1(3)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def g(s, xs):
    return s if xs is None else s / (1 + s / xs)

def shoot(mu, f0, d, xs, rho_max=25.0, h=0.01):
    """Variable reduite rho = sqrt(mu) r : f'' + (d-1)/rho f' = f - g(f^2) f / mu. RK4 depuis rho ~ 0.
    Retourne (statut, [rho], [f]) ; statut 'cross' (f0 trop grand) ou 'positive' (f0 trop petit : la solution repart)."""
    inv_mu = 1.0 / mu
    def rhs(r, f, p):
        return p, f - g(f * f, xs) * f * inv_mu - (d - 1) / r * p
    r, f, p = 1e-6, f0, 0.0
    rs, fs = [r], [f]
    while r < rho_max:
        k1f, k1p = rhs(r, f, p)
        k2f, k2p = rhs(r + h / 2, f + h / 2 * k1f, p + h / 2 * k1p)
        k3f, k3p = rhs(r + h / 2, f + h / 2 * k2f, p + h / 2 * k2p)
        k4f, k4p = rhs(r + h, f + h * k3f, p + h * k3p)
        f += h / 6 * (k1f + 2 * k2f + 2 * k3f + k4f); p += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p); r += h
        rs.append(r); fs.append(f)
        if f < 0: return "cross", rs, fs
        if f > 1.5 * f0: return "positive", rs, fs
    return "positive", rs, fs

def ground_state(mu, d, xs, iters=50):
    """Bissection sur f(0) entre 'positive' (trop petit) et 'cross' (trop grand) ; norme N~ = int f^2 d^d r en unites physiques."""
    lo = (math.sqrt(mu) if xs is None else math.sqrt(mu / (1 - mu / xs))) * 1.0001
    hi = lo * 50
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        st, _, _ = shoot(mu, mid, d, xs)
        if st == "positive": lo = mid
        else: hi = mid
    f0 = 0.5 * (lo + hi)
    _, rs, fs = shoot(mu, f0, d, xs)
    cut = len(fs) - 1
    for i in range(1, len(fs)):
        if fs[i] > fs[i - 1]: cut = i - 1; break
    w = 2 * math.pi if d == 2 else (4 * math.pi if d == 3 else 2.0)
    norm = 0.0
    for i in range(1, cut + 1):
        a = fs[i - 1] ** 2 * rs[i - 1] ** (d - 1); b = fs[i] ** 2 * rs[i] ** (d - 1)
        norm += 0.5 * (a + b) * (rs[i] - rs[i - 1])
    norm *= w
    A, R = fs[cut], rs[cut]
    norm += w * A ** 2 * R ** (d - 1) / 2.0          # queue exponentielle e^{-2 rho}
    return f0, norm * mu ** (-d / 2)                  # retour a r = rho / sqrt(mu)

def main():
    print("Etape 1 (R122) -- auto-piegeage saturable : seuil, echelle, loi de masse, bornes\n")

    # A. controles : 1D exact, 2D Townes, 3D cubique sans minimum
    print("A. Controles de la methode")
    n1 = 4 * math.sqrt(1.0)
    print(f"    d = 1, cubique, mu = 1 : N~ = 4 sqrt(mu) = {n1:.3f} (exact) ; N~ -> 0 quand mu -> 0 : pas de seuil")
    f0_t, n_townes = ground_state(1.0, 2, None)
    print(f"    d = 2, cubique, mu = 1 (Townes) : f(0) = {f0_t:.4f} (ref. 2,2062), N~ = {n_townes:.3f} (ref. 11,70)")
    n3 = [ground_state(mu, 3, None)[1] for mu in (0.25, 1.0, 4.0)]
    print(f"    d = 3, cubique : N~(mu) = {n3[0]:.2f}, {n3[1]:.2f}, {n3[2]:.2f} pour mu = 0,25 ; 1 ; 4 -> decroit comme mu^-1/2 (N~ sqrt(mu) = {n3[1]:.2f}) : branche instable\n")
    check("A. methode : Townes N~_c = 11,70 a 1 % ; 3D cubique decroissant en mu^-1/2 (pas de lump stable)",
          abs(n_townes / 11.70 - 1) < 0.01 and abs(f0_t / 2.2062 - 1) < 0.01 and abs(n3[0] * 0.5 / n3[1] - 1) < 0.02 and abs(n3[2] * 2 / n3[1] - 1) < 0.02,
          f"N~_Townes = {n_townes:.3f}, f0 = {f0_t:.4f}")

    # B. le seuil saturable en 3D
    print("B. d = 3 saturable, g(s) = s/(1 + s/x_s) : N~(mu) et son minimum")
    nmin = {}
    for xs in (0.3, 1.0, 3.0, 10.0):
        mus = [xs * 10 ** (-3 + 2.954 * i / 21) for i in range(22)]      # 1e-3 x_s ... 0,9 x_s, log
        ns = [ground_state(mu, 3, xs)[1] for mu in mus]
        i = min(range(len(ns)), key=lambda k: ns[k])
        nmin[xs] = (mus[i] / xs, ns[i], ns[-1] > ns[i] and ns[0] > ns[i])
        print(f"    x_s = {xs:5.1f} : N~_min = {ns[i]:8.2f} a mu/x_s = {mus[i]/xs:.4f} (f0^2 = {ground_state(mus[i], 3, xs)[0]**2:.2f}) ; N~(1e-3 x_s) = {ns[0]:.1f} (cubique : {18.90/math.sqrt(mus[0]):.1f}), N~(0,9 x_s) = {ns[-1]:.2e} ; branche stable (dN~/dmu > 0) au-dessus du minimum")
    print("    -> un seuil existe des que la non-linearite sature ; N~_min depend de x_s = u_sat/u_s.\n")
    check("B. 3D saturable : un minimum N~_min existe pour chaque x_s (seuil), encadre de part et d'autre",
          all(v[2] for v in nmin.values()) and 5 < nmin[1.0][1] < 200, f"N~_min(x_s=1) = {nmin[1.0][1]:.2f}")

    # C. les echelles
    print("C. Les echelles d'energie E_s = (u (hbar c)^3)^(1/4)")
    def Es(u): return (u * HBARC_J_M ** 3) ** 0.25 / MEV_J
    def u_of_Es(E_mev): return (E_mev * MEV_J) ** 4 / HBARC_J_M ** 3
    u0 = H * C / (L1 * 1e-15) ** 4
    u_eh = 45 * ME_KG ** 4 * C ** 5 / (8 * ALPHA ** 2 * HBAR ** 3)
    u_hibef = u_eh / math.sqrt(7e4)
    E_S = ME_KG ** 2 * C ** 3 / (E_CH * HBAR)
    u_schw = 0.5 * EPS0 * E_S ** 2
    u_bag = u_of_Es(145.0)
    for lab, u in (("u_0 = hc/l_1^4 (substrat au repos)", u0), ("borne HIBEF 2026 (R119), u_s minimal", u_hibef), ("Schwinger eps0 E_S^2/2", u_schw), ("Euler-Heisenberg u_EH", u_eh), ("sac MIT, B^(1/4) = 145 MeV", u_bag)):
        print(f"      {lab:42s} u = {u:9.2e} J/m^3   E_s = {Es(u):8.3f} MeV")
    print()
    check("C. echelles : E_s(u_0) = 0,39 MeV < m_e ; E_s(HIBEF) = 2,3 MeV > 4 m_e ; u_HIBEF > u_Schwinger",
          0.35 < Es(u0) < 0.42 and 2.0 < Es(u_hibef) < 2.6 and u_hibef > u_schw, f"{Es(u0):.3f}, {Es(u_hibef):.3f} MeV ; u_HIBEF/u_S = {u_hibef/u_schw:.0f}")

    # D. la loi de masse
    print("D. La loi de masse au seuil : M c^2 = N~_min^(1/4) E_s N_q^(3/4)")
    n4 = nmin[1.0][1] ** 0.25
    def M(E_s, Nq): return n4 * E_s * Nq ** 0.75
    m10_hibef = M(Es(u_hibef), 10)
    Es_electron = 0.51099895 / (n4 * 1 ** 0.75)
    Es_nucleon = 938.272 / (n4 * 3 ** 0.75)
    print(f"    x_s = 1 : N~_min^(1/4) = {n4:.2f}")
    print(f"    lump classique de 10 quanta a la borne HIBEF : M = {m10_hibef:.0f} MeV = {m10_hibef/0.511:.0f} m_e")
    print(f"    l'electron comme lump d'un quantum : E_s = {Es_electron:.3f} MeV, u_s = {u_of_Es(Es_electron):.1e} J/m^3 : {u_hibef/u_of_Es(Es_electron):.0f} fois sous la borne HIBEF -> EXCLU")
    print(f"    un nucleon de 3 quanta : E_s = {Es_nucleon:.0f} MeV (sac MIT : B^(1/4) = 145-235 MeV), u_s = {u_of_Es(Es_nucleon):.1e} J/m^3 = {u_of_Es(Es_nucleon)/u0:.0e} u_0 : une seconde echelle, posee\n")
    check("D. loi de masse : electron exclu (E_s requis sous la borne HIBEF) ; nucleon -> E_s ~ 100-250 MeV, l'echelle du sac, >= 10^9 u_0",
          u_of_Es(Es_electron) < u_hibef / 100 and 100 < Es_nucleon < 250 and u_of_Es(Es_nucleon) / u0 > 1e9,
          f"E_s(e) = {Es_electron:.2f} MeV ; E_s(N) = {Es_nucleon:.0f} MeV")

    # E. le vrai vide
    print("E. Le vrai vide : saturation = creation de paires ; stabilite")
    print(f"    u_HIBEF / u_Schwinger = {u_hibef/u_schw:.0f} : toute non-linearite assez forte pour pieger agit la ou le vide reel cree des paires (E > E_S = {E_S:.2e} V/m)")
    print("    un lump de quanta sans masse libres de sortir est metastable (oscillon) ; stable seulement avec une charge conservee (boucle fermee, R120) ou un exterieur interdit (sac).\n")
    check("E. u_s(HIBEF) > u_Schwinger : la saturation du vide reel est la creation de paires, pas un plafond de Z", u_hibef / u_schw > 10, f"{u_hibef/u_schw:.0f}")

    print("Verdict : le mecanisme est juste (seuil, lumps stables, loi M ~ E_s N^(3/4) du sac) ; EXCLU comme origine de l'electron")
    print("(HIBEF + Schwinger) ; pour les hadrons c'est le modele du sac avec une echelle E_s ~ 150-200 MeV POSEE, sans lien avec l_1")
    print("(10^10 u_0). L'etape 1 ne donne pas d'objet nouveau ; elle fixe ce que l'objet doit avoir : une charge conservee ou un")
    print("exterieur interdit, et une echelle qui n'est pas celle du substrat.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
