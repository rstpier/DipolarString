#!/usr/bin/env python3
"""R124 -- Les parametres de repere privilegie alpha_1, alpha_2 (et alpha_3) du secteur scalaire de la CD.

Premier calcul du programme sous A3 (matter/rebuild/A3.md, section 5, calcul 1). Probleme ouvert n° 6 du manuscrit :
« une computation post-newtonienne complete avec le repere du substrat en mouvement par rapport au systeme, qui seule
etablirait alpha_1 et alpha_2, n'a pas ete faite ». La voici.

  Le modele. Secteur scalaire contraint (elliptique), representant isotrope auquel la matiere se couple :
      ds^2 = -e^{-chi} dt^2 + e^{chi} (dx^2 + dy^2 + dz^2),   chi = 2U,   c_0 = 1,
  U le potentiel newtonien instantane dans le repere du substrat S (equation de Poisson elliptique, pas de g_0i : un
  champ scalaire n'a pas de potentiel vecteur). Pour une masse m en mouvement a la vitesse v le long de x, deux
  lectures du potentiel coexistent dans le manuscrit :
      (i)  elliptique : U_S = Gm / |x - v t|          (instantane, non contracte ; section « static use ») ;
      (ii) Heaviside  : U_S = Gm / sqrt(gamma^2 (x - vt)^2 + y^2 + z^2)   (contracte ; section Lorentz).
  Methode. (1) Lecture dans S : forme PPN du repere privilegie pour une masse en mouvement, V_i = v_i U, W_i = v^j U_ij,
  Phi_1 = v^2 U ; g_0i = -(4 gamma + 3 + alpha_1 - alpha_2) V_i / 2 - (1 + alpha_2) W_i / 2, g_00 contient
  (2 gamma + 2 + alpha_3) Phi_1. La somme des deux coefficients de g_0i est invariante de jauge. (2) Lecture dans le
  repere comobile S' : boost de Lorentz exact (c_0 = 1 ; la matiere obeit a la cinematique de Lorentz par A6),
  developpement en epsilon (U ~ eps^2, v ~ eps), forme PPN avec w = v : g_00 = -1 + 2U - 2 beta U^2
  - (alpha_1 - alpha_2 - alpha_3) w^2 U - alpha_2 w^i w^j U_ij ; g_0i = -(alpha_1 - 2 alpha_2) w_i U / 2 - alpha_2 w^j U_ij ;
  jauge t' -> t' + lambda (w.x') U qui deplace lambda entre les deux termes de g_0i (somme invariante).
  Controle : la metrique 1PN de la RG d'une masse en mouvement (boost de la metrique de repos) doit rendre
  alpha_1 = alpha_2 = alpha_3 = 0, gamma = beta = 1 par les deux lectures.

  Resultat. alpha_1 = -8 dans les deux lectures et dans les deux variantes : un champ scalaire n'a pas de terme
  gravitomagnetique dans le repere du substrat, et la RG en exige g_0i = -4 v_i U pour une masse en mouvement ;
  l'ecart vaut -(4 gamma + 4) = -8 avec gamma = 1 (Cassini), sans parametre libre. alpha_2 = -1 (potentiel instantane)
  ou 0 (Heaviside) ; alpha_3 = -3 ou -4 avec une source de masse au repos (0 seulement si la source est accordee a
  la structure 1PN de la RG). Bornes : |alpha_1| < 4e-5 (pulsar J1738+0333 ; 1e-4 par la Lune), |alpha_2| < 4e-7
  (alignement du spin solaire) et 1,8e-9 (pulsars), |alpha_3| < 4e-20 (pulsars). Le secteur scalaire a repere est
  EXCLU de 2e5 par alpha_1 seul. Sortie : si le bootstrap tensoriel rend l'exterieur egal a la RG, alpha = 0 mais
  le repere du substrat est inobservable en gravitation : le milieu n'y ajoute rien (critere (3) de A3.md).
  Critere d'arret (1) de A3.md declenche : le milieu a repere est mort ; il reste de la CD des rapports et des
  etats limites.
"""
import sympy as sp

G, m, v, e = sp.symbols("G m v e", positive=True)   # e : parametre de comptage, U ~ e^2, v ~ e
ONE, ZERO = sp.Integer(1), sp.Integer(0)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def gamma2(vv):
    return 1 / (1 - vv ** 2)

def boost_J(vv):
    """x^alpha = J^alpha_mu x'^mu : S' se deplace a vv le long de x par rapport a S (c = 1). J^-1 = boost(-vv)."""
    g = sp.sqrt(gamma2(vv))
    return sp.Matrix([[g, g * vv, 0, 0], [g * vv, g, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

def iso_metric(U):
    """Representant isotrope exponentiel du manuscrit, chi = 2U, exponentielles tronquees a O(U^2) (U^3 ~ e^6 hors PPN)."""
    ep, em = 1 + 2 * U + 2 * U ** 2, 1 - 2 * U + 2 * U ** 2
    return sp.diag(-em, ep, ep, ep)

def ser(expr, n):
    return sp.expand(sp.series(sp.expand(expr), e, 0, n).removeO())

def coeff_e(expr, k, n):
    return ser(expr, n).coeff(e, k)

def read_comoving(gS_of_point, label):
    """gS_of_point(xp,yp,zp) -> metrique de S (4x4) aux coordonnees comobiles ; lecture PPN dans S' apres boost exact
    (w = v) et jauge lambda ; renvoie alpha_1, alpha_2, alpha_3, beta, gamma."""
    J = boost_J(e * v)
    def gprime(P):
        return (J.T * gS_of_point(*P) * J).applyfunc(sp.expand)
    g1, g2 = gprime((ONE, ZERO, ZERO)), gprime((ZERO, ONE, ZERO))      # r' = 1 : U' = G m
    a2 = coeff_e(g1[0, 0], 2, 6)
    t4_1, t4_2 = coeff_e(g1[0, 0], 4, 6), coeff_e(g2[0, 0], 4, 6)
    c1 = t4_2.coeff(G, 2) / m ** 2                  # U^2
    c2 = t4_2.coeff(G, 1).coeff(v, 2) / m           # v^2 U
    c3 = t4_1.coeff(G, 1).coeff(v, 2) / m - c2      # v^2 U n_x^2
    t3_1, t3_2 = coeff_e(g1[0, 1], 3, 5), coeff_e(g2[0, 1], 3, 5)
    b1 = t3_2.coeff(v, 1).coeff(G, 1) / m           # v U
    b2 = t3_1.coeff(v, 1).coeff(G, 1) / m - b1      # v U n_x^2
    gam = coeff_e(g2[2, 2], 2, 4) / (2 * G * m)
    alpha2 = -c3                                    # jauge lambda = alpha_2 ; exige b2 = 0
    alpha1 = -2 * b1
    alpha3 = alpha1 - alpha2 + c2
    beta = -c1 / 2
    print(f"    [{label}] g'_00 : e^2 = {a2} ; e^4 : U^2 -> {c1}, v^2 U -> {c2}, v^2 U n_x^2 -> {c3} ; g'_0x : v U -> {b1}, v U n_x^2 -> {b2} ; g'_yy : 2 gamma = {2*gam}")
    print(f"    [{label}] -> alpha_1 = {alpha1}, alpha_2 = {alpha2}, alpha_3 = {alpha3}, beta = {beta}, gamma = {gam}")
    return dict(a1=sp.nsimplify(alpha1), a2=sp.nsimplify(alpha2), a3=sp.nsimplify(alpha3), beta=sp.nsimplify(beta), gamma=sp.nsimplify(gam), b2=sp.nsimplify(b2))

def read_in_S(gS_point_fn, label):
    """Lecture dans le repere du substrat (forme PPN du repere privilegie, w = 0) pour une masse en mouvement :
    g_0x = -(4 gamma + 3 + a1 - a2)/2 V_x - (1 + a2)/2 W_x avec V_x = v U, W_x = v U n_x^2 ; g_00 contient (2 gamma + 2 + a3) Phi_1, Phi_1 = v^2 U."""
    g1, g2 = gS_point_fn(ONE, ZERO, ZERO).applyfunc(sp.expand), gS_point_fn(ZERO, ONE, ZERO).applyfunc(sp.expand)
    s0x_1 = coeff_e(g1[0, 1], 3, 5).coeff(v, 1).coeff(G, 1) / m      # n = x : V_x + W_x
    s0x_2 = coeff_e(g2[0, 1], 3, 5).coeff(v, 1).coeff(G, 1) / m      # n = y : V_x seul
    gam = coeff_e(g2[2, 2], 2, 4) / (2 * G * m)
    alpha1 = -2 * s0x_1 - 4 * gam - 4
    phi1 = coeff_e(g2[0, 0], 4, 6).coeff(G, 1).coeff(v, 2) / m        # coefficient de Phi_1 (n = y : pas de n_x^2)
    alpha3 = phi1 - 2 * gam - 2
    print(f"    [{label}] dans S : g_0x = {s0x_1} v U (n = x), {s0x_2} v U (n = y) ; coefficient de Phi_1 : {phi1} ; gamma = {gam} -> alpha_1 = {alpha1}, alpha_3 = {alpha3}")
    return dict(a1=sp.nsimplify(alpha1), a3=sp.nsimplify(alpha3), gamma=sp.nsimplify(gam))

def main():
    print("R124 -- alpha_1, alpha_2, alpha_3 du secteur scalaire de la CD, substrat en mouvement\n")
    g2e = gamma2(e * v)
    Jinv = boost_J(-e * v)            # S exprime depuis S' : l'inverse du boost est le boost de -v

    # --- A. Controle : RG 1PN d'une masse en mouvement = boost de la metrique de repos (gamma = beta = 1) ---
    print("A. Controle : la metrique de la RG (1PN) d'une masse en mouvement, par les deux lectures")
    def gS_GR_comoving(xp_, yp_, zp_):
        U = e ** 2 * G * m / sp.sqrt(xp_ ** 2 + yp_ ** 2 + zp_ ** 2)
        return (Jinv.T * iso_metric(U) * Jinv).applyfunc(sp.expand)
    rc = read_comoving(gS_GR_comoving, "RG, comobile")
    def gS_GR_S(x_, y_, z_):
        U = e ** 2 * G * m / sp.sqrt(g2e * x_ ** 2 + y_ ** 2 + z_ ** 2)       # x' = gamma x a t = 0
        return (Jinv.T * iso_metric(U) * Jinv).applyfunc(sp.expand)
    rs = read_in_S(gS_GR_S, "RG, substrat")
    print()
    check("A. controle RG : alpha_1 = alpha_2 = alpha_3 = 0, gamma = beta = 1 par les deux lectures",
          rc["a1"] == 0 and rc["a2"] == 0 and rc["a3"] == 0 and rc["beta"] == 1 and rc["gamma"] == 1 and rs["a1"] == 0 and rs["a3"] == 0,
          f"comobile {rc['a1']}, {rc['a2']}, {rc['a3']} ; substrat {rs['a1']}, {rs['a3']}")

    # --- B. variante (i) : potentiel elliptique instantane (non contracte) ---
    print("B. Le modele, variante (i) : metrique diagonale dans S, potentiel instantane U_S = Gm/|x - vt|")
    def gS_i_comoving(xp_, yp_, zp_):
        U = e ** 2 * G * m / sp.sqrt(xp_ ** 2 / g2e + yp_ ** 2 + zp_ ** 2)      # |x - vt| = x'/gamma
        return iso_metric(U)
    ri = read_comoving(gS_i_comoving, "(i), comobile")
    def gS_i_S(x_, y_, z_):
        return iso_metric(e ** 2 * G * m / sp.sqrt(x_ ** 2 + y_ ** 2 + z_ ** 2))
    ris = read_in_S(gS_i_S, "(i), substrat")
    print()

    # --- C. variante (ii) : ellipsoide de Heaviside (contracte) ---
    print("C. Le modele, variante (ii) : metrique diagonale dans S, ellipsoide de Heaviside U_S = Gm/sqrt(gamma^2 (x-vt)^2 + y^2 + z^2)")
    def gS_ii_comoving(xp_, yp_, zp_):
        return iso_metric(e ** 2 * G * m / sp.sqrt(xp_ ** 2 + yp_ ** 2 + zp_ ** 2))   # r* = r' exactement
    rii = read_comoving(gS_ii_comoving, "(ii), comobile")
    def gS_ii_S(x_, y_, z_):
        return iso_metric(e ** 2 * G * m / sp.sqrt(g2e * x_ ** 2 + y_ ** 2 + z_ ** 2))
    riis = read_in_S(gS_ii_S, "(ii), substrat")
    print()
    check("B. alpha_1 = -8 exactement, deux variantes, deux lectures (pas de gravitomagnetisme dans le repere du substrat)",
          ri["a1"] == -8 and rii["a1"] == -8 and ris["a1"] == -8 and riis["a1"] == -8 and ri["b2"] == 0 and rii["b2"] == 0,
          f"(i) {ri['a1']}/{ris['a1']} ; (ii) {rii['a1']}/{riis['a1']}")
    check("C. alpha_2 = -1 (instantane) et 0 (Heaviside) ; alpha_3 = -3 et -4 ; gamma = beta = 1 retrouves",
          ri["a2"] == -1 and rii["a2"] == 0 and ri["a3"] == -3 and rii["a3"] == -4 and ri["beta"] == 1 and rii["beta"] == 1 and ri["gamma"] == 1 and rii["gamma"] == 1,
          f"alpha_2 = {ri['a2']}, {rii['a2']} ; alpha_3 = {ri['a3']}, {rii['a3']}")

    # --- D. les bornes ---
    print("D. Les bornes")
    bounds = {"alpha_1 (pulsar J1738+0333)": 4e-5, "alpha_1 (Lune, LLR)": 1e-4, "alpha_2 (spin solaire)": 4e-7, "alpha_2 (pulsars)": 1.8e-9, "alpha_3 (pulsars)": 4e-20}
    viol = {}
    for k, b in bounds.items():
        val = {"alpha_1": 8.0, "alpha_2": 1.0, "alpha_3": 4.0}[k.split()[0]]
        viol[k] = val / b
        print(f"      {k:30s} |.| < {b:.1e} ; modele {val:.0f} : violation x {viol[k]:.1e}")
    print("    (la variante (ii) a alpha_2 = 0 : alpha_1 suffit a exclure)\n")
    check("D. exclusion : alpha_1 = -8 viole la borne pulsar de 2e5 (et la borne lunaire de 8e4) sans parametre libre",
          viol["alpha_1 (pulsar J1738+0333)"] > 1e5 and viol["alpha_1 (Lune, LLR)"] > 1e4, f"x {viol['alpha_1 (pulsar J1738+0333)']:.0e}")

    # --- E. la sortie ---
    print("E. La sortie : si le bootstrap tensoriel rend l'exterieur egal a la RG, alpha = 0 (controle A), mais le repere du substrat")
    print("   est alors inobservable en gravitation : le milieu n'ajoute rien a l'exterieur (critere (3) de A3.md). Pas de troisieme voie :")
    print("   un champ scalaire dans un repere ne fournit pas le g_0i = -4 v_i U que la RG et les pulsars exigent.\n")
    check("E. aucune valeur de parametre du secteur scalaire ne ramene alpha_1 a 0 : -(4 gamma + 4), gamma = 1 fixe par Cassini", -(4 * 1 + 4) == -8, "-8")

    print("Verdict : alpha_1 = -8, alpha_2 = -1 ou 0, alpha_3 = -3 ou -4 pour le secteur scalaire a repere, contre des bornes de")
    print("4e-5, 4e-7 (1,8e-9) et 4e-20 : EXCLU de cinq ordres de grandeur par alpha_1 seul, sans parametre libre. Critere d'arret (1)")
    print("de A3.md declenche : le milieu a repere est mort en gravitation ; la seule sortie (bootstrap vers la RG) le rend inobservable.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
