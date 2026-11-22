#!/usr/bin/env python3
"""R136 -- Calcul B1 refait a l'identique de l'article de 2014, lu a la source, et ce que la lecture corrige.

Le 22 novembre 2026 l'auteur a autorise le telechargement depuis arXiv. Lus : Sonnenschein et Weissman,
arXiv:1402.5603 (JHEP 08 (2014) 013) ; Bulava et al., arXiv:2403.00754 (Phys. Lett. B 854 (2024) 138754) ; Badalian
et Bakker, arXiv:1901.10280. Les masses ci-dessous sont celles de la table 8 de l'article de 2014 ; la mesure
d'ecart est la sienne, chi2 = (1 / (N - 1)) somme ((M^2 - E^2) / M^2)^2, E etant la masse calculee de l'etat.

  A. Ajustements lineaires (extremites sans masse) des dix trajectoires : pente, ordonnee et chi2 publies, refaits.
  B. Ajustements a extremites massives aux masses publiees : rho/a (65 MeV), pi/b (170), K* (60, 220), phi (200 et
     400), Psi (1500), D (80, 1640), Ds* (400, 1580).
  C. Ajustement universel de l'article : neuf trajectoires, 38 etats, masses fixees a 60, 220 et 1500 MeV, une pente
     et neuf ordonnees. Publie : alpha' = 0,884 GeV^-2, chi2 = 13,13e-4, et les masses calculees de la table 8.
  D. Ce que la lecture corrige. R133 comparait 0,884 a un ajustement de six familles « comme l'article » : l'article
     en ajuste neuf, dont trois avec un quark charme, et garde l'eta. Avec ses etats, les six familles legeres et
     etranges donnent a elles seules 0,884 : le nombre de R133 tient, sa description etait inexacte.
  E. La tension dans l'ajustement universel refait. Pente fixee a celle de la base (0,8595) puis a celle du reseau
     (0,8037), ordonnees libres : hausse du chi2 et ecart en ecarts types (erreur commune ajustee).
  F. Cote reseau, lu a la source. Potentiel V0 + sigma r + gamma / r avec gamma = -0,434(11), contre -pi/12 = -0,262
     pour une corde seule ; sigma t0 = 0,1061(7)(20), sqrt(t0) = 0,1443 fm ; echelles de Sommer de ce potentiel ;
     brisure a 1,211 et 1,336 fm au point physique. Autres determinations citees par l'article : 467(7) et 482(7) MeV
     (2+1+1 saveurs, ajustements de Cornell sous 1 fm) ; jauge pure, sigma t0 de 0,143 a 0,159.
  Verdict : B1 est reproduit a l'identique, a mieux que 1 % sur les pentes et au chi2 publie pres. Dans cet
  ajustement la tension des spectres vaut 424 +- 2,5 MeV ; celle de la base en est a 2,5 ecarts types, celle du
  reseau a 2,9 ; et la « tension du reseau » depend elle-meme de sa definition, de 445 a 482 MeV avec des quarks
  dynamiques.
"""
from math import asin, pi, sqrt

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = 0.1973269804                       # GeV fm
ME = 0.51099895e-3
L1_9 = (2 * pi * HBARC / ME / 3) * 3 ** (-2 * pi)
AP_BASE = 1 / (2 * pi * (pi * HBARC ** 2 / L1_9 ** 2))
AP_LAT = 1 / (2 * pi * 0.445 ** 2)

def arm(T, m, w):
    if m == 0:
        l = 1 / w
        return T * l * pi / 2, T * l * l * pi / 4
    k = T / (m * w)
    b = (-1 + sqrt(1 + 4 * k * k)) / (2 * k)
    g = 1 / sqrt(1 - b * b); l = b / w
    return g * m + T * l * asin(b) / b, g * m * b * l + 0.5 * T * l * l * (asin(b) - b * sqrt(1 - b * b)) / (b * b)

def mass_of(Jc, ap, m1, m2):
    """Masse de la corde tournante de moment cinetique classique Jc > 0."""
    if Jc <= 0: return m1 + m2
    if m1 == 0 and m2 == 0: return sqrt(Jc / ap)
    T = 1 / (2 * pi * ap)
    lo, hi = 1e-6, 1e6
    for _ in range(70):
        w = sqrt(lo * hi)
        if arm(T, m1, w)[1] + arm(T, m2, w)[1] > Jc: lo = w
        else: hi = w
    return arm(T, m1, w)[0] + arm(T, m2, w)[0]

def chi2_terms(states, ap, a, m1, m2):
    return [((M * M - mass_of(J - a, ap, m1, m2) ** 2) / (M * M)) ** 2 for J, M in states]

def golden(f, lo, hi, n=48):
    g = (sqrt(5) - 1) / 2
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo); f1, f2 = f(x1), f(x2)
    for _ in range(n):
        if f1 < f2: hi, x2, f2 = x2, x1, f1; x1 = hi - g * (hi - lo); f1 = f(x1)
        else: lo, x1, f1 = x1, x2, f2; x2 = lo + g * (hi - lo); f2 = f(x2)
    return (lo + hi) / 2

def best_a(states, ap, m1, m2):
    jmin = min(J for J, _ in states)
    f = lambda a: sum(chi2_terms(states, ap, a, m1, m2))
    a = golden(f, jmin - 6.0, jmin - 1e-6, 40)
    return a, f(a)

def fit(states, m1, m2, lo=0.3, hi=1.4):
    ap = golden(lambda x: best_a(states, x, m1, m2)[1], lo, hi, 36)
    a, s = best_a(states, ap, m1, m2)
    return ap, a, s / (len(states) - 1)

# Table 8 de l'article (masses mesurees, MeV) ; variable de trajectoire : J, sauf Psi (L)
TRAJ = {
    "pi/b":  ([(1, 1.229), (2, 1.672), (3, 2.032), (4, 2.250)], "u", "u"),
    "rho/a": ([(1, 0.776), (2, 1.318), (3, 1.689), (4, 1.996), (5, 2.330), (6, 2.450)], "u", "u"),
    "eta/h": ([(0, 0.548), (1, 1.170), (2, 1.617), (3, 2.025), (4, 2.328)], "u", "u"),
    "omega/f": ([(1, 0.783), (2, 1.275), (3, 1.667), (4, 2.018), (5, 2.250), (6, 2.469)], "u", "u"),
    "K*":    ([(1, 0.892), (2, 1.426), (3, 1.776), (4, 2.045), (5, 2.382)], "u", "s"),
    "phi":   ([(1, 1.020), (2, 1.525), (3, 1.854)], "s", "s"),
    "D":     ([(0, 1.865), (1, 2.421), (2, 2.737)], "u", "c"),
    "Ds*":   ([(1, 2.112), (2, 2.572), (3, 2.862)], "s", "c"),
    "Psi":   ([(0, 3.097), (1, 3.494), (2, 3.778)], "c", "c"),
}
LINEAR_PUB = {"pi/b": (0.808, -0.23, 7.99), "rho/a": (0.882, 0.47, 9.90), "eta/h": (0.839, -0.25, 39.63), "omega/f": (0.909, 0.45, 8.85),
              "K*": (0.849, 0.33, 7.15), "phi": (0.814, 0.15, 4.43), "D": (0.480, -1.69, 13.92), "Ds*": (0.522, -1.35, 6.44), "Psi": (0.418, -4.04, 3.41)}
MASSIVE_PUB = [("rho/a", 0.065, 0.065, 0.896, 0.52), ("pi/b", 0.170, 0.170, 0.844, 0.00), ("K*", 0.060, 0.220, 0.885, 0.50),
               ("phi", 0.200, 0.200, 0.882, 0.41), ("phi", 0.400, 0.400, 1.078, 0.82), ("Psi", 1.500, 1.500, 0.979, -0.09),
               ("D", 0.080, 1.640, 1.073, -0.07), ("Ds*", 0.400, 1.580, 1.093, 0.89)]
UNIV_MASS = {"u": 0.060, "s": 0.220, "c": 1.500}
UNIV_PUB_A = {"pi/b": -0.33, "rho/a": 0.52, "eta/h": -0.22, "omega/f": 0.53, "K*": 0.50, "phi": 0.46, "D": -0.19, "Ds*": -0.39, "Psi": -0.06}
# L'ordonnee publiee du Ds* (-0,39) differe de 1,00 exactement de l'ordonnee refaite en J (+0,61) : elle est donnee dans la
# convention du moment orbital L = J - 1 (table 4 de l'article : 0,89 - 1 = -0,11 pour l'ajustement individuel).
UNIV_SHIFT = {"Ds*": 1.0}
TABLE8_CALC = {"pi/b": [1257, 1650, 1965, 2236], "rho/a": [776, 1324, 1701, 2008, 2274, 2511], "eta/h": [545, 1206, 1612, 1933, 2208],
               "omega/f": [768, 1319, 1698, 2006, 2271, 2509], "K*": [892, 1415, 1783, 2084, 2345], "phi": [1019, 1514, 1870],
               "D": [1862, 2408, 2752], "Ds*": [2112, 2563, 2881], "Psi": [3080, 3535, 3824]}

def universal(ap):
    out = {}; tot = 0.0
    for name, (st, q1, q2) in TRAJ.items():
        a, s = best_a(st, ap, UNIV_MASS[q1], UNIV_MASS[q2]); out[name] = a; tot += s
    return tot, out

def main():
    print("R136 -- calcul B1 refait a l'identique de l'article de 2014, lu a la source\n")

    # A
    print("A. Ajustements lineaires (extremites sans masse)")
    print("    trajectoire   publie (alpha', a, chi2 x 1e4)      refait                  ecart sur alpha'")
    dev_lin = {}
    for name, (pa, pb, pc) in LINEAR_PUB.items():
        ap, a, c = fit(TRAJ[name][0], 0, 0)
        dev_lin[name] = (ap / pa - 1, a - pb, c * 1e4 / pc - 1)
        print(f"    {name:9s}     {pa:.3f}  {pb:+.2f}  {pc:6.2f}            {ap:.3f}  {a:+.2f}  {c*1e4:6.2f}        {100*(ap/pa-1):+.2f} %")
    ap, a, c = fit(TRAJ["eta/h"][0][1:], 0, 0)
    dev_eta = (ap / 0.745 - 1, a + 0.01)
    print(f"    eta/h sans l'eta : publie 0.745  -0.01    4.57 ; refait {ap:.3f}  {a:+.2f}  {c*1e4:6.2f}")
    worst_lin = max(abs(v[0]) for v in dev_lin.values())
    print(f"    ecart maximal sur la pente : {100*worst_lin:.2f} % ; sur l'ordonnee : {max(abs(v[1]) for v in dev_lin.values()):.3f}\n")
    check("A. neuf ajustements lineaires publies refaits : pentes a moins de 1 %, ordonnees a moins de 0,03",
          worst_lin < 0.01 and max(abs(v[1]) for v in dev_lin.values()) < 0.03 and abs(dev_eta[0]) < 0.01, f"ecart maximal {100*worst_lin:.2f} %")

    # B
    print("B. Ajustements a extremites massives, aux masses publiees")
    dev_mass = []
    for name, m1, m2, pa, pb in MASSIVE_PUB:
        ap, a, c = fit(TRAJ[name][0], m1, m2)
        dev_mass.append((ap / pa - 1, a - pb))
        print(f"    {name:6s} masses {1e3*m1:5.0f} et {1e3*m2:5.0f} MeV : publie {pa:.3f}  {pb:+.2f} ; refait {ap:.3f}  {a:+.2f} ; ecart {100*(ap/pa-1):+.2f} %")
    worst_mass = max(abs(d[0]) for d in dev_mass)
    print(f"    ecart maximal sur la pente : {100*worst_mass:.2f} % ; sur l'ordonnee : {max(abs(d[1]) for d in dev_mass):.3f}\n")
    check("B. huit ajustements massifs publies refaits : pentes a moins de 2 %, ordonnees a moins de 0,05",
          worst_mass < 0.02 and max(abs(d[1]) for d in dev_mass) < 0.05, f"ecart maximal {100*worst_mass:.2f} %")

    # C
    print("C. Ajustement universel : neuf trajectoires, 38 etats, masses 60, 220 et 1500 MeV")
    nstates = sum(len(st) for st, _, _ in TRAJ.values())
    ap_u = golden(lambda x: universal(x)[0], 0.6, 1.2, 40)
    tot, a_u = universal(ap_u)
    chi_u = tot / (nstates - 1)
    print(f"    pente : publiee 0.884, refaite {ap_u:.4f} ({100*(ap_u/0.884-1):+.2f} %) ; chi2 : publie 13.13e-4, refait {chi_u*1e4:.2f}e-4 ({nstates} etats)")
    a_cmp = {k: a_u[k] - UNIV_SHIFT.get(k, 0.0) for k in TRAJ}
    print("    ordonnees (publiee / refaite) : " + " ; ".join(f"{k} {UNIV_PUB_A[k]:+.2f} / {a_cmp[k]:+.2f}" for k in TRAJ))
    print(f"       (Ds* : refaite {a_u['Ds*']:+.2f} avec J, soit {a_cmp['Ds*']:+.2f} dans la convention du moment orbital employee par l'article pour cette famille)")
    da = max(abs(a_cmp[k] - UNIV_PUB_A[k]) for k in TRAJ)
    dm = []
    for name, (st, q1, q2) in TRAJ.items():
        for (J, M), calc in zip(st, TABLE8_CALC[name]):
            dm.append(1e3 * mass_of(J - a_u[name], ap_u, UNIV_MASS[q1], UNIV_MASS[q2]) - calc)
    print(f"    masses calculees de la table 8 (38 etats) : ecart maximal a la table {max(abs(x) for x in dm):.1f} MeV, moyen {sum(abs(x) for x in dm)/len(dm):.1f} MeV")
    print()
    check("C. ajustement universel refait : alpha' a moins de 0,5 % de 0,884, ordonnees a moins de 0,01, chi2 = 13,13e-4 comme publie, masses de la table 8 a moins de 5 MeV",
          nstates == 38 and abs(ap_u / 0.884 - 1) < 0.005 and da < 0.01 and abs(chi_u * 1e4 - 13.13) < 0.05 and max(abs(x) for x in dm) < 5,
          f"alpha' = {ap_u:.4f}, chi2 = {chi_u*1e4:.2f}e-4")

    # D
    print("D. Ce que la lecture corrige dans R133 et R135")
    six = {k: TRAJ[k] for k in ("pi/b", "rho/a", "eta/h", "omega/f", "K*", "phi")}
    def partial(ap, fam): return sum(best_a(st, ap, UNIV_MASS[q1], UNIV_MASS[q2])[1] for st, q1, q2 in fam.values())
    ap_six = golden(lambda x: partial(x, six), 0.6, 1.2, 40)
    ap_light = golden(lambda x: partial(x, {k: TRAJ[k] for k in ("pi/b", "rho/a", "eta/h", "omega/f")}), 0.6, 1.2, 40)
    ap_charm = golden(lambda x: partial(x, {k: TRAJ[k] for k in ("D", "Ds*", "Psi")}), 0.6, 1.4, 40)
    print(f"    l'ajustement universel porte sur neuf trajectoires, dont trois charmees, et garde l'eta ; R133 en prenait six, sans l'eta ni omega5")
    print(f"    avec les etats de l'article : six familles legeres et etranges {ap_six:.3f} ; quatre familles legeres {ap_light:.3f} ; trois familles charmees {ap_charm:.3f} ; neuf {ap_u:.3f}")
    print("    masse d'extremite du charme : 1500 MeV (R133 et R135 citaient 1490 de memoire) ; 60 et 220 MeV confirmes ; 0,884 confirme")
    print("    -> le nombre de R133 (0,883 +- 0,014) tient : les six familles legeres et etranges donnent a elles seules la pente publiee.")
    print("       Sa description, « six familles, comme l'article », etait inexacte ; la reproduction a l'identique la remplace.\n")
    check("D. pente universelle par sous-ensemble : familles legeres et etranges, familles charmees, ensemble des neuf",
          0.85 < ap_six < 0.92 and 0.85 < ap_light < 0.93 and 0.80 < ap_charm < 1.10, f"{ap_six:.3f} ; {ap_charm:.3f} ; {ap_u:.3f}")

    # E
    print("E. La tension dans l'ajustement universel refait")
    dof = nstates - len(TRAJ) - 1
    s2 = tot / dof
    out = {}
    for label, ap in (("base", AP_BASE), ("reseau", AP_LAT)):
        t, _ = universal(ap)
        out[label] = (t / tot, sqrt(max(t - tot, 0) / s2))
        print(f"    pente fixee a celle de {label:6s} ({ap:.4f}) : chi2 multiplie par {t/tot:.2f}, soit {out[label][1]:.2f} ecarts types ({dof} degres de liberte)")
    h = 0.01
    curv = (universal(ap_u + h)[0] - 2 * tot + universal(ap_u - h)[0]) / (h * h)
    err_u = sqrt(2 * s2 / curv)
    rs = lambda ap: 1e3 / sqrt(2 * pi * ap)
    print(f"    pente universelle refaite : {ap_u:.3f} +- {err_u:.3f} GeV^-2, soit sqrt(sigma) = {rs(ap_u):.1f} +- {rs(ap_u)*err_u/(2*ap_u):.1f} MeV")
    n_lat_full = (445.0 - rs(ap_u)) / sqrt(6.7 ** 2 + (rs(ap_u) * err_u / (2 * ap_u)) ** 2)
    print(f"    avec l'erreur du reseau (+- 6,7 MeV) : reseau - spectres = {445.0 - rs(ap_u):+.1f} MeV, {n_lat_full:.1f} ecarts types")
    print("    -> dans l'ajustement meme de l'article, la base est a 2,5 ecarts types : en tension moderee, pas exclue ; la pente que donnerait")
    print("       la tension du reseau est exclue comme pente des trajectoires.\n")
    check("E. ajustement universel : pente 0,885 +- 0,010, sqrt(sigma) = 424 +- 2,5 MeV ; base a 2,5 ecarts types, reseau a 8,6 (2,9 avec son erreur)",
          2.2 < out["base"][1] < 2.7 and 8 < out["reseau"][1] < 9.2 and abs(rs(ap_u) - 424.1) < 1 and 2.6 < n_lat_full < 3.2,
          f"base {out['base'][1]:.1f} ; reseau {out['reseau'][1]:.1f} ; {n_lat_full:.1f}")

    # F
    print("F. Cote reseau, lu a la source (Bulava et al. 2024)")
    gamma, dgamma = -0.434, 0.011
    sig_t0, sqrt_t0 = 0.1061, 0.1443
    rs_lat = sqrt(sig_t0) / sqrt_t0 * HBARC * 1e3
    r0 = sqrt((1.65 + gamma) / sig_t0) * sqrt_t0
    r1 = sqrt((1.00 + gamma) / sig_t0) * sqrt_t0
    luscher = -pi / 12
    print(f"    potentiel V0 + sigma r + gamma / r : gamma = {gamma}({int(dgamma*1000)}) ; corde seule : -pi/12 = {luscher:.3f} ; rapport {gamma/luscher:.2f}")
    print(f"    sigma t0 = {sig_t0}, sqrt(t0) = {sqrt_t0} fm -> sqrt(sigma) = {rs_lat:.0f} MeV (publie : 445(3)(6))")
    print(f"    echelles de Sommer de ce potentiel : r0 = {r0:.3f} fm (reseau 0,472 a 0,474), r1 = {r1:.3f} fm (reseau 0,3127) ; r0 sqrt(sigma) = {sqrt(1.65 + gamma):.3f}")
    print("    brisure au point physique : 1,211(7)(11) fm (legers) et 1,336(7)(12) fm (etrange)")
    others = {"2+1 saveurs, modele avec brisure": (445, 445), "2+1+1 saveurs, ajustements de Cornell sous 1 fm": (467, 482),
              "jauge pure (sigma t0 de 0,143 a 0,159)": (sqrt(0.143) / sqrt_t0 * HBARC * 1e3, sqrt(0.159) / sqrt_t0 * HBARC * 1e3)}
    for k, (a, b) in others.items(): print(f"      {k:50s} : sqrt(sigma) = {a:.0f} a {b:.0f} MeV")
    print("    -> a distance intermediaire le terme en 1/r est 1,7 fois celui d'une corde ; et la tension du reseau depend de sa definition.\n")
    check("F. reseau : gamma = -0,434 soit 1,66 fois le terme de Luscher ; sqrt(sigma) = 445 MeV ; r0 = 0,489 fm et r1 = 0,333 fm ; 445 a 482 MeV selon la definition",
          abs(gamma / luscher - 1.66) < 0.02 and abs(rs_lat - 445) < 1.5 and abs(r0 - 0.4885) < 0.002 and abs(r1 - 0.3333) < 0.002
          and 510 < others["jauge pure (sigma t0 de 0,143 a 0,159)"][0] < 525, f"r0 = {r0:.3f} fm, r1 = {r1:.3f} fm")

    print("Verdict : le calcul B1 est reproduit a l'identique de l'article, a mieux que 1 % sur les pentes : la cible de la feuille de route est")
    print("atteinte, cette fois contre des nombres lus a la source. La tension des spectres vaut 424 +- 2,5 MeV ; celle de la base en est a")
    print("2,5 ecarts types, celle du reseau a 2,9. Et la tension du reseau depend de sa definition, de 445 a 482 MeV : B3 se pose des deux cotes.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
