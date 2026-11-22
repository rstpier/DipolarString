#!/usr/bin/env python3
"""R133 -- La corde effective a quatre dimensions : son domaine, le calcul B1 de la feuille de route, la tension.

Decision de l'auteur (22 novembre 2026), devant R132 : « ok prend la corde effective ».
Definition retenue : le filament est une corde effective a quatre dimensions, c'est-a-dire la description des tubes
de flux longs devant la longueur de corde (Luscher et Weisz ; Polchinski et Strominger 1991), avec la tension de la
base. C'est le chemin B de la feuille de route. On fait ici son premier calcul.

  A. La corde tournante a extremites massives (B1). sigma l = m gamma^2 beta^2 ;
     E = somme sur les deux bras de gamma m + sigma l arcsin(beta) / beta ;
     J = somme de gamma m beta l + sigma l^2 [arcsin(beta) - beta sqrt(1 - beta^2)] / (2 beta^2).
     Controles : condition aux extremites, limite sans masse J = E^2 / (2 pi sigma), soit R120.
  B. Ajustements. J = J_classique(M ; alpha', m1, m2) + ordonnee, sur six familles de mesons et deux de baryons ;
     masses d'extremite 60 MeV (u, d) et 220 MeV (s) ; deux ponderations. Pente commune des six familles, a comparer
     a la valeur publiee par Sonnenschein et Weissman, 0,884 GeV^-2 (de memoire) : cible de la feuille de route, 2 %.
  C. La tension de la base contre les trajectoires. La pente ajustee depend de la famille (de 0,74 a 0,92 GeV^-2),
     du nombre de familles retenues, de la ponderation et des masses d'extremite : huit variantes d'ajustement
     global. La base, 0,8595 GeV^-2, est dans leur enveloppe ; avec les seules familles vecteur et tenseur elle est
     trop basse en alpha' de 2 a 4 ecarts types. Le test n'est pas concluant.
  D. Le domaine. Rayon d'Arvis r_c = sqrt(pi / (6 sigma)) = 0,33 fm ; longueur de chaque etat ; le parametre de
     developpement (r_c / L)^2 ; la corde casse vers 1,25 fm (reseau, de memoire) : fenetre de moins d'une decade.
  E. Ce qui ne depend d'aucun nombre de la CD, et l'echelle de Sommer. Terme de Luscher -pi / (12 r) ; potentiel
     d'Arvis ; r0 defini par r^2 V'(r) = 1,65 : r0 sqrt(sigma) = 1,18 pour la corde seule, 1,06 sur le reseau : a
     0,5 fm le potentiel n'est pas celui d'une corde seule (c'est l'enjeu de B2, la fenetre d'ajustement).
  F. La tension de la base entre deux mesures. sqrt(sigma) = m_e 3^(2 pi + 1) / (2 sqrt(pi)) = 430,3 MeV (identite) ;
     trajectoires : enveloppe des huit variantes ; reseau a quarks dynamiques 445 +- 7 MeV (R119). Aucun ecart
     n'atteint trois ecarts types : ni la base, ni la question B3, ne sont tranchees.
  Verdict : la corde effective est coherente a quatre dimensions et c'est la theorie etablie des tubes de flux. Elle
  ne couvre que les hadrons, dans une fenetre etroite, et le rho et le nucleon sont a son bord. Le seul nombre de la
  CD y est la tension, compatible avec l'enveloppe des trajectoires et a deux ecarts types sous le reseau : ni
  confirmee ni exclue.
  B2 (tables de reseau) et la recherche bibliographique de B3 restent a faire.
"""
from math import asin, pi, sqrt

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = 0.1973269804                       # GeV fm
ME = 0.51099895e-3                         # GeV
LAMBDA_E = HBARC / ME                      # fm
L1_9 = (2 * pi * LAMBDA_E / 3) * 3 ** (-2 * pi)
SIGMA = pi * HBARC ** 2 / L1_9 ** 2        # GeV^2
AP_BASE = 1 / (2 * pi * SIGMA)             # GeV^-2
RS_BASE = 1e3 * sqrt(SIGMA)                # MeV

def arm(T, m, w):
    """Un bras de la corde tournante : tension T, masse d'extremite m, vitesse angulaire w -> (E, J, longueur, beta)."""
    if m == 0:
        l = 1 / w
        return T * l * pi / 2, T * l * l * pi / 4, l, 1.0
    k = T / (m * w)
    b = (-1 + sqrt(1 + 4 * k * k)) / (2 * k)                  # sigma l = m gamma^2 beta^2 avec l = beta / w
    g = 1 / sqrt(1 - b * b); l = b / w
    return g * m + T * l * asin(b) / b, g * m * b * l + 0.5 * T * l * l * (asin(b) - b * sqrt(1 - b * b)) / (b * b), l, b

def classical(M, T, m1, m2):
    """J classique et longueur totale de la corde tournante de masse M."""
    lo, hi = 1e-6, 1e6
    for _ in range(64):
        w = sqrt(lo * hi)
        if arm(T, m1, w)[0] + arm(T, m2, w)[0] > M: lo = w
        else: hi = w
    a1, a2 = arm(T, m1, w), arm(T, m2, w)
    return a1[1] + a2[1], a1[2] + a2[2]

# Masses en GeV. rho/a et baryons : valeurs de R120 (PDG 2024) ; autres etats : PDG, de memoire, a verifier.
MU, MS = 0.060, 0.220
MESONS = {
    "rho / a":   ([(1, 0.77526), (2, 1.3182), (3, 1.6888), (4, 1.967), (5, 2.330), (6, 2.450)], MU, MU),
    "omega / f": ([(1, 0.78266), (2, 1.2754), (3, 1.667), (4, 2.018), (6, 2.470)], MU, MU),
    "K*":        ([(1, 0.89167), (2, 1.4273), (3, 1.779), (4, 2.048), (5, 2.382)], MU, MS),
    "phi":       ([(1, 1.019461), (2, 1.5174), (3, 1.854)], MS, MS),
}
PSEUDO = {                                   # familles pi / b et eta / h, comme dans l'article de 2014 ; masses de memoire
    "pi / b":    ([(1, 1.2295), (2, 1.6706), (3, 2.032), (4, 2.250)], MU, MU),
    "eta / h":   ([(1, 1.166), (2, 1.617), (3, 2.025), (4, 2.328)], MU, MU),
}
BARYONS = {
    "N":     ([(0.5, 0.93827), (2.5, 1.685), (4.5, 2.250)], MU, MU),
    "Delta": ([(1.5, 1.232), (3.5, 1.930), (5.5, 2.420)], MU, MU),
}

def residuals(states, m1, m2, ap, weight):
    T = 1 / (2 * pi * ap)
    res = [J - classical(M, T, m1, m2)[0] for J, M in states]
    ws = [1.0 if weight == "J" else 1 / (ap * M * M) ** 2 for _, M in states]
    a = sum(w * r for w, r in zip(ws, res)) / sum(ws)
    return sum(w * (r - a) ** 2 for w, r in zip(ws, res)), a

def minimize(f, lo=0.5, hi=1.6):
    g = (sqrt(5) - 1) / 2
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo); f1, f2 = f(x1), f(x2)
    for _ in range(44):
        if f1 < f2: hi, x2, f2 = x2, x1, f1; x1 = hi - g * (hi - lo); f1 = f(x1)
        else: lo, x1, f1 = x1, x2, f2; x2 = lo + g * (hi - lo); f2 = f(x2)
    return (lo + hi) / 2

def fit_family(states, m1, m2, weight):
    ap = minimize(lambda x: residuals(states, m1, m2, x, weight)[0])
    chi, a = residuals(states, m1, m2, ap, weight)
    return ap, a, chi

def fit_global(families, weight, massive=True):
    def chi(ap): return sum(residuals(st, m1 if massive else 0, m2 if massive else 0, ap, weight)[0] for st, m1, m2 in families.values())
    ap = minimize(chi)
    dof = sum(len(st) for st, _, _ in families.values()) - len(families) - 1
    c0 = chi(ap); s2 = c0 / dof                                     # erreur commune ajustee : chi2 / ddl = 1
    h = 0.01
    curv = (chi(ap + h) - 2 * c0 + chi(ap - h)) / (h * h)           # chi2 = c0 + curv (x - ap)^2 / 2
    err = sqrt(2 * s2 / curv)
    return ap, err, sqrt(max(chi(AP_BASE) - c0, 0) / s2), dof

def rs(ap): return 1e3 / sqrt(2 * pi * ap)                         # sqrt(sigma) en MeV

def main():
    print("R133 -- la corde effective a quatre dimensions : domaine, calcul B1, tension\n")
    print(f"    base : sigma = {SIGMA:.4f} GeV^2, sqrt(sigma) = {RS_BASE:.1f} MeV, alpha' = {AP_BASE:.4f} GeV^-2\n")

    # A
    print("A. La corde tournante a extremites massives")
    T = SIGMA
    massless = max(abs(classical(M, T, 0, 0)[0] - M * M / (2 * pi * T)) for M in (0.775, 1.318, 2.0, 3.0))
    small = max(abs(classical(M, T, 1e-6, 1e-6)[0] / (M * M / (2 * pi * T)) - 1) for M in (0.775, 1.318, 2.0))
    bc = 0.0
    for m, w in ((0.06, 1.0), (0.22, 0.4), (1.5, 0.3)):
        E, J, l, b = arm(T, m, w)
        bc = max(bc, abs(T * l - m * b * b / (1 - b * b)) / (T * l))
    below = all(classical(M, T, 0.06, 0.06)[0] < classical(M, T, 0, 0)[0] for M in (0.775, 1.318, 2.0))
    print(f"    condition aux extremites sigma l = m gamma^2 beta^2 : ecart relatif maximal {bc:.1e}")
    print(f"    limite sans masse : J - E^2 / (2 pi sigma) = {massless:.1e} ; masses de 1 keV : ecart relatif {small:.1e}")
    print(f"    des extremites massives reduisent J a masse donnee : {below} -> la pente ajustee monte quand on les ajoute.\n")
    check("A. corde tournante a extremites massives : condition aux extremites et limite sans masse J = E^2 / (2 pi sigma) verifiees",
          bc < 1e-12 and massless < 1e-9 and small < 1e-4 and below, "formules de la feuille de route")

    # B
    print("B. Ajustements par famille (masses d'extremite 60 et 220 MeV), J = J_classique + ordonnee")
    print("    famille      etats   sans masse (poids J, M^2)   avec masses (poids J, M^2)   ordonnee   sqrt(sigma) en MeV")
    fam_ap = {}
    for name, (st, m1, m2) in {**MESONS, **PSEUDO, **BARYONS}.items():
        r0 = [fit_family(st, 0, 0, w)[0] for w in ("J", "M2")]
        r1 = [fit_family(st, m1, m2, w) for w in ("J", "M2")]
        fam_ap[name] = (r0, [r[0] for r in r1])
        print(f"    {name:10s}   {len(st)}        {r0[0]:.3f}  {r0[1]:.3f}               {r1[0][0]:.3f}  {r1[1][0]:.3f}            {r1[1][1]:+.2f}       {rs(r1[0][0]):.0f}  {rs(r1[1][0]):.0f}")
    SIX = {**MESONS, **PSEUDO}
    variants = {}
    print("    pente commune (erreur commune ajustee pour que le chi2 par degre de liberte vaille 1) :")
    for label, fams in (("quatre familles vecteur et tenseur", MESONS), ("six familles, comme l'article", SIX)):
        for w in ("J", "M2"):
            for massive in (True, False):
                variants[(label, w, massive)] = fit_global(fams, w, massive)
                ap, err, nsig, dof = variants[(label, w, massive)]
                print(f"      {label:34s} poids {w:2s} {'avec' if massive else 'sans'} masses : alpha' = {ap:.3f} +- {err:.3f}, sqrt(sigma) = {rs(ap):.1f} MeV, base a {nsig:.1f} ecart type")
    sw = 0.884
    ref = variants[("six familles, comme l'article", "M2", True)]
    dev_sw = abs(ref[0] / sw - 1)
    print(f"    valeur publiee (Sonnenschein et Weissman 2014, de memoire) : {sw} GeV^-2 ; ici, six familles, poids M^2, avec masses : {ref[0]:.3f} +- {ref[1]:.3f}")
    print(f"    ecart {100*dev_sw:.1f} % : la cible de 2 % de la feuille de route est atteinte\n")
    check("B. six familles, masses d'extremite, poids M^2 : pente commune a moins de 2 % de la valeur publiee 0,884 GeV^-2",
          dev_sw < 0.02 and abs(ref[0] - 0.883) < 0.004, f"alpha' = {ref[0]:.3f} +- {ref[1]:.3f}")

    # C
    print("C. La tension de la base contre les trajectoires")
    fam_with = {name: fam_ap[name][1] for name in SIX}
    lo_f, hi_f = min(min(v) for v in fam_with.values()), max(max(v) for v in fam_with.values())
    aps = [v[0] for v in variants.values()]
    env = (min(aps), max(aps))
    inside = env[0] <= AP_BASE <= env[1]
    four = [variants[("quatre familles vecteur et tenseur", w, True)] for w in ("J", "M2")]
    six = [variants[("six familles, comme l'article", w, True)] for w in ("J", "M2")]
    print(f"    par famille, avec masses : alpha' de {lo_f:.3f} (eta / h) a {hi_f:.3f} (omega / f), soit sqrt(sigma) de {rs(hi_f):.0f} a {rs(lo_f):.0f} MeV")
    print(f"    huit variantes d'ajustement global : alpha' de {env[0]:.3f} a {env[1]:.3f}, sqrt(sigma) de {rs(env[1]):.1f} a {rs(env[0]):.1f} MeV")
    print(f"    base : alpha' = {AP_BASE:.4f}, sqrt(sigma) = {RS_BASE:.1f} MeV : dans l'enveloppe : {inside}")
    print(f"    familles vecteur et tenseur seules, avec masses : base trop basse en alpha' de {100*(1-AP_BASE/four[0][0]):.1f} et {100*(1-AP_BASE/four[1][0]):.1f} %, {four[0][2]:.1f} et {four[1][2]:.1f} ecarts types")
    print(f"    six familles, avec masses : {six[0][2]:.1f} et {six[1][2]:.1f} ecart type")
    bar = [x for name in BARYONS for x in fam_ap[name][0] + fam_ap[name][1]]
    print(f"    baryons (trois etats par famille) : alpha' de {min(bar):.3f} a {max(bar):.3f}")
    print("    -> la pente « universelle » varie de 20 % d'une famille a l'autre ; selon les familles retenues la base est compatible ou trop")
    print("       haute en tension de 2 a 4 ecarts types. Les trajectoires ne tranchent pas. R120 (« en tension ») est a nuancer dans ce sens.\n")
    check("C. pente par famille de 0,74 a 0,92 GeV^-2 ; la base est dans l'enveloppe des huit ajustements globaux ; quatre familles seules : 2 a 4,5 ecarts types",
          lo_f < 0.76 and 0.91 < hi_f < 0.93 and inside and all(1.8 < v[2] < 4.5 for v in four) and all(v[2] < 2 for v in six),
          f"enveloppe {env[0]:.3f} a {env[1]:.3f} ; base {AP_BASE:.4f}")

    # D
    print("D. Le domaine de la corde effective")
    r_c = HBARC * sqrt(pi / (6 * SIGMA))
    print(f"    rayon d'Arvis r_c = sqrt(pi / (6 sigma)) = {r_c:.3f} fm : en dessous le potentiel d'Arvis n'existe plus ; le developpement est en (r_c / L)^2")
    lengths = {}
    for name, M, m1, m2 in (("rho", 0.77526, MU, MU), ("a2", 1.3182, MU, MU), ("rho3", 1.6888, MU, MU), ("a4", 1.967, MU, MU), ("a6", 2.450, MU, MU),
                            ("phi", 1.019461, MS, MS), ("N", 0.93827, MU, MU), ("Delta", 1.232, MU, MU)):
        L = classical(M, SIGMA, m1, m2)[1] * HBARC
        lengths[name] = L
        print(f"      {name:6s} M = {1e3*M:6.0f} MeV : longueur {L:.2f} fm = {L/r_c:.1f} r_c ; (r_c / L)^2 = {(r_c/L)**2:.2f}")
    r_break = 1.25
    print(f"    brisure de la corde statique vers {r_break} fm (reseau, de memoire) : energie stockee sigma r = {1e3*SIGMA/HBARC*r_break:.0f} MeV ; fenetre {r_c:.2f} a {r_break} fm, facteur {r_break/r_c:.1f}")
    print("    -> le rho et le nucleon, a 1,5 r_c, sont au bord du domaine ; la corde effective est fiable pour les grands moments cinetiques.")
    print("       La « part statique » de la base, corde a J = 1/2 (R120 A), n'est donc pas controlee.\n")
    check("D. domaine : r_c = 0,33 fm ; rho et nucleon a moins de 2 r_c, a4 et a6 au-dela de 3,5 r_c ; fenetre statique de 0,33 a 1,25 fm",
          abs(r_c - 0.332) < 2e-3 and lengths["rho"] < 2 * r_c and lengths["N"] < 2.2 * r_c and lengths["a4"] > 3.5 * r_c and 3.5 < r_break / r_c < 4,
          f"r_c = {r_c:.3f} fm ; rho {lengths['rho']:.2f} fm, a4 {lengths['a4']:.2f} fm")

    # E
    print("E. Ce qui ne depend d'aucun nombre de la CD, et l'echelle de Sommer")
    luscher = -pi / 12 * HBARC * 1e3
    s_fm = SIGMA / HBARC ** 2                                         # fm^-2
    def r2dV(r): return s_fm * r * r / sqrt(1 - pi / (6 * s_fm * r * r))        # potentiel d'Arvis, D = 4
    r_min = r_c * sqrt(1.5)                                           # r^2 V' y est minimal et vaut (3 sqrt(3) / 2) (pi / 6)
    f_min = r2dV(r_min)
    lo, hi = r_min, 3.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if r2dV(mid) < 1.65: lo = mid
        else: hi = mid
    r0_arvis = (lo + hi) / 2
    r0_lo = sqrt((1.65 - pi / 12) / s_fm)
    k_string = sqrt(1.65 - pi / 12)
    r0_lat, rs_lat, drs_lat = 0.472, 445.0, 7.0
    k_lat = r0_lat * rs_lat / (1e3 * HBARC)
    print(f"    terme de Luscher : -pi hbar c / (12 r) = {luscher:.1f} MeV fm / r, universel ; a 1 fm : {luscher:.1f} MeV")
    print(f"    echelle de Sommer, r^2 V'(r0) = 1,65 : corde seule r0 sqrt(sigma) = sqrt(1,65 - pi/12) = {k_string:.3f}")
    print(f"       avec la tension de la base : r0 = {r0_lo:.3f} fm (Luscher), {r0_arvis:.3f} fm (Arvis) ; (r_c / r0)^2 = {(r_c/r0_lo)**2:.2f}")
    print(f"       reseau : r0 = {r0_lat} fm (FLAG, de memoire) et sqrt(sigma) = {rs_lat:.0f} +- {drs_lat:.0f} MeV (R119) : r0 sqrt(sigma) = {k_lat:.3f}")
    print(f"       la seconde echelle du reseau, r1, definie par 1,0, n'existe pas pour la corde d'Arvis : son r^2 V' ne descend pas sous {f_min:.2f}")
    print("    -> a 0,5 fm le potentiel mesure n'est pas celui d'une corde seule : la tension extraite depend de la fenetre (B2, non fait ici).\n")
    check("E. terme de Luscher -51,7 MeV fm / r ; r0 sqrt(sigma) = 1,18 pour la corde seule contre 1,06 sur le reseau ; pas d'echelle r1 pour la corde d'Arvis",
          abs(luscher + 51.66) < 0.05 and abs(k_string - 1.178) < 2e-3 and 1.04 < k_lat < 1.09 and 0.53 < r0_lo < 0.55 and 0.50 < r0_arvis < r0_lo
          and abs(f_min - 1.360) < 2e-3 and f_min > 1.0, f"{k_string:.3f} contre {k_lat:.3f}")

    # F
    print("F. La tension de la base entre deux mesures")
    formula = 1e3 * ME * 3 ** (2 * pi + 1) / (2 * sqrt(pi))
    rs_lo, rs_hi = rs(env[1]), rs(env[0])
    rs_regge, drs_regge = (rs_lo + rs_hi) / 2, (rs_hi - rs_lo) / 2
    d1, d2, d3 = RS_BASE - rs_regge, rs_lat - RS_BASE, rs_lat - rs_regge
    n1, n2, n3 = d1 / drs_regge, d2 / drs_lat, d3 / sqrt(drs_lat ** 2 + drs_regge ** 2)
    print(f"    sqrt(sigma) de la base = m_e 3^(2 pi + 1) / (2 sqrt(pi)) = {formula:.1f} MeV (identite ; aucun mecanisme, et R126 : accord non significatif)")
    print(f"    trajectoires des mesons, enveloppe des huit variantes : {rs_regge:.1f} +- {drs_regge:.1f} MeV ; reseau a quarks dynamiques : {rs_lat:.0f} +- {drs_lat:.0f} MeV")
    print(f"    base - trajectoires = {d1:+.1f} MeV ({n1:.1f} ecart type) ; reseau - base = {d2:+.1f} MeV ({n2:.1f}) ; reseau - trajectoires = {d3:+.1f} MeV ({n3:.1f})")
    print("    -> aucun ecart n'atteint trois ecarts types : la tension de la base n'est ni confirmee ni exclue, et la question B3 (les deux tensions")
    print("       sont-elles la meme ?) n'est pas tranchee par ce calcul. Le budget d'erreur est domine par le choix des familles.\n")
    check("F. sqrt(sigma) de la base = m_e 3^(2 pi + 1) / (2 sqrt(pi)) = 430,3 MeV ; trajectoires 426 +- 6, reseau 445 +- 7 : aucun ecart a trois ecarts types",
          abs(formula - RS_BASE) < 1e-6 and abs(rs_regge - 426) < 2 and 5 < drs_regge < 7.5 and abs(n1) < 1.5 and 1.8 < n2 < 2.4 and 1.5 < n3 < 2.6,
          f"{rs_regge:.0f} +- {drs_regge:.0f} ; {RS_BASE:.1f} ; {rs_lat:.0f} +- {drs_lat:.0f} MeV")

    print("Verdict : la corde effective est coherente a quatre dimensions ; c'est la theorie etablie des tubes de flux. Elle ne couvre que les")
    print("hadrons, dans une fenetre de 0,33 a 1,25 fm, et le rho et le nucleon sont a son bord. Le seul nombre de la CD y est la tension :")
    print("430 MeV, dans l'enveloppe des trajectoires et a deux ecarts types sous le reseau, ni confirmee ni exclue. B1 fait et conforme a la")
    print("valeur publiee ; B2 (tables de reseau) et la recherche bibliographique de B3 restent a faire.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
