#!/usr/bin/env python3
"""R139 -- Calcul B2 : la tension tiree du potentiel statique du reseau avec la forme d'une corde, et l'ecart aux spectres.

Choix fait par delegation le 22 novembre 2026, inscrit au registre avec son critere avant le calcul. Mecanisme : le
reseau mesure l'energie d'un quark et d'un antiquark immobiles en fonction de leur distance ; la pente de cette
energie a grande distance est la tension du fil. Mais la valeur qu'on en tire depend de la forme qu'on ajuste.

Source : Bulava et al., arXiv:2403.00754 (Phys. Lett. B 854 (2024) 138754), table 4, parametres extrapoles au point
physique, en unites de l'echelle sqrt(t0) = 0,1443 fm. L'article ne publie pas les points du potentiel, seulement le
modele ajuste : un etat de corde V^(r) = V0 + sigma r + gamma / r, couple a deux etats de paires de mesons.
B2 est donc fait sur ce modele publie, entre 0,25 et 1,5 fm ou il a ete ajuste aux donnees.

Critere ecrit avant le calcul (BASE, 22 novembre) : ecart entre la tension tiree avec la forme d'une corde,
V0 + sigma r sqrt(1 - pi / (6 sigma r^2)), dans la fenetre de 0,4 a 1,0 fm, et celle des spectres, 423,7 +- 2,5 MeV,
en ecarts types, erreurs ajoutees en quadrature, la dependance a la fenetre comptee dans l'erreur du reseau.
A plus de 3, l'ecart est etabli ; a moins de 2, les deux tensions s'accordent ; entre les deux, non tranche.

  A. Le modele publie, relu : il redonne les deux distances de brisure publiees et la tension de 445 MeV.
  B. La forme d'une corde : sa force ne peut pas descendre sous 1,360 / r^2. Sous 0,43 fm elle ne peut pas egaler
     celle du reseau, quelle que soit la tension.
  C. La tension locale : celle qu'une corde doit avoir pour tirer comme le reseau a la distance r.
  D. Ajustements par fenetre sur l'etat de corde.
  E. Le niveau le plus bas, celui que le fil partage avec les paires de mesons : sa detente, mesuree par le reseau.
  F. Sensibilite aux parametres publies.
  G. Budget d'erreur et critere.
  Verdict : avec la forme d'une corde le reseau donne 455 MeV et non 445 ; l'ecart aux spectres est de 31 MeV,
  3,65 ecarts types selon le critere ecrit, 2,97 avec l'erreur la plus prudente du cote des spectres.
"""
from math import acos, cos, pi, sqrt
import ds_b1_reproduction as b1
import ds_string_relaxation as rel

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = b1.HBARC * 1e3                       # MeV fm
ST0 = 0.1443                                 # sqrt(t0) en fm
A_LAT = 0.0633                               # maille des trois ensembles, en fm
PAR = dict(sig=0.1061, v0=-0.835, gam=-0.434, e1=0.0034, e2=0.1005, gl=0.0256, gs=0.0134)   # table 4, unites de t0
RS_LAT, DRS_LAT = 445.0, 6.7                 # MeV, publie : 445(3)(6)
mev = lambda e: e * HBARC / ST0              # energie en unites de 1 / sqrt(t0) -> MeV
rs = lambda s: sqrt(s) * HBARC / ST0         # tension en unites de 1 / t0 -> sqrt(sigma) en MeV
X = lambda fm: fm / ST0

def bis(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        m = 0.5 * (lo + hi); fm_ = f(m)
        if (fm_ > 0) == (flo > 0): lo, flo = m, fm_
        else: hi = m
    return 0.5 * (lo + hi)

def golden(f, lo, hi, n=80):
    return b1.golden(f, lo, hi, n)

def v_string(x, p=PAR): return p["v0"] + p["sig"] * x + p["gam"] / x
def f_string(x, p=PAR): return p["sig"] - p["gam"] / (x * x)

def levels(x, p=PAR):
    """Les trois niveaux du modele publie : valeurs propres de [[V^, sqrt(2) gl, gs], [sqrt(2) gl, E1, 0], [gs, 0, E2]]."""
    a, b, c = v_string(x, p), p["e1"], p["e2"]; d, e = sqrt(2) * p["gl"], p["gs"]
    q = (a + b + c) / 3
    pp = sqrt(((a - q) ** 2 + (b - q) ** 2 + (c - q) ** 2 + 2 * (d * d + e * e)) / 6)
    det = ((a - q) * (b - q) * (c - q) - d * d * (c - q) - e * e * (b - q)) / pp ** 3
    ph = acos(max(-1.0, min(1.0, det / 2))) / 3
    l1 = q + 2 * pp * cos(ph); l3 = q + 2 * pp * cos(ph + 2 * pi / 3)
    return l3, 3 * q - l1 - l3, l1
def v_ground(x, p=PAR): return levels(x, p)[0]

def arvis(x, s): return s * x * sqrt(1 - pi / (6 * s * x * x))
def f_arvis(x, s): return s / sqrt(1 - pi / (6 * s * x * x))

def local_tension(x, p=PAR):
    """Tension qu'une corde doit avoir pour tirer comme l'etat de corde du reseau a la distance x ; None si impossible."""
    k = pi / (6 * x * x)
    if f_string(x, p) < 1.5 * sqrt(3) * k: return None
    return bis(lambda s: f_arvis(x, s) - f_string(x, p), 1.5 * k * (1 + 1e-12), 5.0)

def window_fit(target, r1, r2, wpow=0.0, step=A_LAT):
    """Moindres carres de V0 + corde sur target, aux distances r1, r1 + step, ..., r2 (fm) ; poids r^wpow."""
    xs = []; r = r1
    while r <= r2 + 1e-9: xs.append(X(r)); r += step
    w = [x ** wpow for x in xs]
    def cost(s, full=False):
        shape = [arvis(x, s) for x in xs]
        v0 = sum(wi * (target(x) - sh) for wi, x, sh in zip(w, xs, shape)) / sum(w)
        res = [target(x) - sh - v0 for x, sh in zip(xs, shape)]
        c = sum(wi * r_ * r_ for wi, r_ in zip(w, res))
        return (c, res) if full else c
    s = golden(cost, pi / (6 * xs[0] ** 2) * (1 + 1e-9), 0.5)
    res = cost(s, True)[1]
    return rs(s), mev(sqrt(sum(r_ * r_ for r_ in res) / len(res))), mev(max(abs(r_) for r_ in res)), len(xs)

def main():
    print("R139 -- calcul B2 : la tension tiree du potentiel statique du reseau avec la forme d'une corde\n")

    # A
    print("A. Le modele publie, relu (Bulava et al. 2024, table 4)")
    xl = bis(lambda x: v_string(x) - PAR["e1"], 5, 12); xs_ = bis(lambda x: v_string(x) - PAR["e2"], 5, 12)
    print(f"    etat de corde V^(r) = V0 + sigma r + gamma / r, sigma t0 = {PAR['sig']}, gamma = {PAR['gam']} : sqrt(sigma) = {rs(PAR['sig']):.1f} MeV (publie 445)")
    print(f"    brisure, ou V^ croise l'energie de deux mesons : {xl:.2f} sqrt(t0) = {xl*ST0:.3f} fm (publie 8,39 et 1,211) ; avec l'etrange {xs_:.2f} = {xs_*ST0:.3f} fm (publie 9,26 et 1,336)")
    lv = levels(X(1.0)); a, b, c = v_string(X(1.0)), PAR["e1"], PAR["e2"]; d, e = sqrt(2) * PAR["gl"], PAR["gs"]
    chk = max(abs((a - l) * (b - l) * (c - l) - d * d * (c - l) - e * e * (b - l)) for l in lv)
    print(f"    les trois niveaux a 1 fm : {mev(lv[0]):.1f}, {mev(lv[1]):.1f}, {mev(lv[2]):.1f} MeV (zero : deux mesons statiques legers) ; residu du determinant {chk:.0e}")
    print("    l'article n'ajuste qu'a une maille, 0,0633 fm, entre 4 et 24 a 27 mailles : de 0,25 a 1,5 ou 1,7 fm. Il ne publie pas les points.\n")
    check("A. modele publie relu : brisure retrouvee a 8,39 et 9,26 sqrt(t0), soit 1,211 et 1,336 fm ; sqrt(sigma) = 445 MeV",
          abs(xl - 8.39) < 0.01 and abs(xs_ - 9.26) < 0.01 and abs(xl * ST0 - 1.211) < 0.002 and abs(rs(PAR["sig"]) - 445) < 1 and chk < 1e-12,
          f"{xl:.3f} et {xs_:.3f} sqrt(t0) ; {rs(PAR['sig']):.1f} MeV")

    # B
    print("B. La forme d'une corde : V0 + sigma r sqrt(1 - pi / (6 sigma r^2))")
    x_far = X(6.0); s = PAR["sig"]
    lus = (arvis(x_far, s) - s * x_far) * x_far
    kmin = 1.5 * sqrt(3) * pi / 6
    xmin = sqrt((kmin + PAR["gam"]) / PAR["sig"])
    print(f"    a grande distance elle redonne le terme universel : (V - sigma r) r = {lus:.4f} a 6 fm, pour -pi/12 = {-pi/12:.4f}")
    print(f"    sa force, sigma / sqrt(1 - pi / (6 sigma r^2)), ne descend jamais sous 3 sqrt(3) pi / (12 r^2) = {kmin:.3f} / r^2, quelle que soit la tension")
    print(f"    la force du reseau vaut sigma + {-PAR['gam']} / r^2 : une corde ne peut l'egaler qu'au-dela de r = {xmin*ST0:.3f} fm")
    print("    -> la fenetre de la feuille de route, 0,4 a 1,0 fm, commence juste sous la limite. Sous 0,43 fm le potentiel n'est pas celui d'une corde.\n")
    check("B. forme d'une corde : terme universel -pi/12 retrouve ; force minimale 1,360 / r^2 ; egalite avec le reseau impossible sous 0,426 fm",
          abs(lus + pi / 12) < 2e-3 and abs(kmin - 1.3603) < 1e-3 and abs(xmin * ST0 - 0.426) < 0.002 and local_tension(X(0.42)) is None and local_tension(X(0.43)) is not None,
          f"{lus:.4f} ; {kmin:.4f} ; {xmin*ST0:.3f} fm")

    # C
    print("C. La tension locale : celle qu'une corde doit avoir pour tirer comme le reseau a la distance r")
    loc = {r: rs(local_tension(X(r))) for r in (0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 4.0)}
    print("    " + " ; ".join(f"{r:.2f} fm : {v:.1f}" for r, v in loc.items()) + " MeV")
    rmax = golden(lambda r: -rs(local_tension(X(r))), 0.45, 1.0, 60)
    print(f"    maximum {rs(local_tension(X(rmax))):.1f} MeV a {rmax:.2f} fm ; elle ne rejoint 445 MeV qu'a l'infini : le terme en 1 / r du reseau vaut {-PAR['gam']/(pi/12):.2f} fois celui d'une corde\n")
    check("C. tension locale : 455 MeV a 0,5 fm, maximum 458 vers 0,6 fm, 452 a 1 fm, 446 a 4 fm ; au-dessus de 445 partout au-dela de 0,46 fm",
          abs(loc[0.5] - 455.1) < 0.5 and 457 < rs(local_tension(X(rmax))) < 459 and abs(loc[1.0] - 452.0) < 0.5 and abs(loc[4.0] - 445.9) < 0.5 and all(v > 445 for r, v in loc.items() if r >= 0.5)
          and rs(local_tension(X(0.46))) > 445 > rs(local_tension(X(0.45))),
          f"{loc[0.5]:.1f} ; {rs(local_tension(X(rmax))):.1f} ; {loc[1.0]:.1f} ; {loc[4.0]:.1f} MeV")

    # D
    print("D. Ajustements par fenetre sur l'etat de corde (une distance par maille ; poids egaux, puis en 1 / r^2)")
    wins = [(0.40, 1.0), (0.45, 1.0), (0.50, 1.0), (0.60, 1.0), (0.70, 1.0), (0.50, 1.2), (0.70, 1.2), (0.80, 1.5)]
    fits = {}
    for r1, r2 in wins:
        t, rms, mx, n = window_fit(v_string, r1, r2); t2 = window_fit(v_string, r1, r2, -2.0)[0]
        fits[(r1, r2)] = (t, rms, mx, t2)
        print(f"    [{r1:.2f} ; {r2:.2f}] fm, {n:2d} points : sqrt(sigma) = {t:.1f} MeV ; ecart a la courbe publiee {rms:.2f} MeV en moyenne, {mx:.2f} au pire ; poids en 1 / r^2 : {t2:.1f}")
    central = fits[(0.40, 1.0)][0]
    fine = window_fit(v_string, 0.40, 1.0, 0.0, 0.01)[0]
    bad = window_fit(v_string, 0.25, 1.0)
    print(f"    fenetre de la feuille de route : {central:.1f} MeV ; avec un pas de 0,01 fm : {fine:.1f} ; en partant de 0,25 fm comme l'article : {bad[0]:.0f} MeV, ecart moyen {bad[1]:.0f} MeV, la forme ne tient plus")
    print("    -> entre 0,45 et 1,0 fm la forme d'une corde, qui a un parametre de moins, suit la courbe publiee a 1 MeV pres. Elle demande")
    print("       455 MeV la ou le terme en 1 / r libre donne 445 : la tension du reseau depend de 10 MeV de la forme ajustee.\n")
    string_vals = [v[0] for k, v in fits.items() if k[1] <= 1.2]
    check("D. etat de corde : 454,8 MeV dans la fenetre 0,4 a 1,0 fm ; de 453 a 456 MeV dans les fenetres sous 1,2 fm ; ecart a la courbe publiee sous 1,3 MeV ; poids sans effet",
          abs(central - 454.8) < 0.3 and 452.5 < min(string_vals) and max(string_vals) < 456.5 and max(v[1] for k, v in fits.items()) < 1.3
          and max(abs(v[0] - v[3]) for v in fits.values()) < 0.7 and abs(fine - central) < 0.5 and bad[1] > 50,
          f"{central:.1f} MeV ; de {min(string_vals):.1f} a {max(string_vals):.1f}")

    # E
    print("E. Le niveau le plus bas : la detente du fil, mesuree par le reseau")
    gap = {r: mev(v_ground(X(r)) - v_string(X(r))) for r in (0.5, 0.7, 1.0, 1.1, 1.2, 1.3)}
    def fratio(r):
        h = 1e-5; x = X(r)
        return (v_ground(x + h) - v_ground(x - h)) / (2 * h) / f_string(x)
    fr = {r: fratio(r) for r in (0.7, 1.0, 1.1, 1.2, 1.3)}
    bb = {}
    for r in (1.0, 1.1, 1.2, 1.3):
        f, df = rel.phi(r / b1.HBARC, rel.BEST); bb[r] = f + df * r / b1.HBARC
    print("    niveau le plus bas moins etat de corde : " + " ; ".join(f"{r:.1f} fm : {v:.1f} MeV" for r, v in gap.items()))
    print("    force du niveau le plus bas, en part de celle de la corde : " + " ; ".join(f"{r:.1f} fm : {v:.2f}" for r, v in fr.items()))
    print("    profil publie utilise en R138, meme grandeur : " + " ; ".join(f"{r:.1f} fm : {v:.2f}" for r, v in bb.items()))
    gfits = {w: window_fit(v_ground, *w) for w in ((0.50, 1.0), (0.60, 1.0), (0.70, 1.0), (0.50, 1.1), (0.50, 1.2))}
    print("    forme d'une corde ajustee sur le niveau le plus bas : " + " ; ".join(f"[{w[0]:.1f} ; {w[1]:.1f}] : {v[0]:.1f} MeV (ecart {v[1]:.1f})" for w, v in gfits.items()))
    print("    -> sur le reseau le fil ne se detend pas lentement jusqu'a 3 fm : son niveau le plus bas tire encore a 96 % a 1,0 fm, puis")
    print("       decroche entre 1,1 et 1,3 fm, quand une paire de mesons coute moins que le fil. Sous 1,0 fm le melange abaisse la tension")
    print("       ajustee de 3 a 4 MeV ; au-dela la forme d'une corde ne decrit plus ce niveau.\n")
    ground_vals = [v[0] for w, v in gfits.items() if w[1] <= 1.0]
    check("E. niveau le plus bas : 11 MeV sous l'etat de corde a 1,0 fm et 45 MeV a 1,2 fm ; force a 96 % a 1,0 fm, 56 % a 1,2 fm, 16 % a 1,3 fm ; fenetres sous 1,0 fm : 450 a 453 MeV",
          abs(gap[1.0] + 11.2) < 0.5 and abs(gap[1.2] + 45.2) < 1.0 and abs(fr[1.0] - 0.956) < 0.005 and abs(fr[1.2] - 0.562) < 0.01 and abs(fr[1.3] - 0.158) < 0.01
          and 450 < min(ground_vals) and max(ground_vals) < 453.5 and gfits[(0.50, 1.2)][1] > 5 and bb[1.2] > 0.9,
          f"{gap[1.0]:.1f} et {gap[1.2]:.1f} MeV ; {fr[1.0]:.3f}, {fr[1.2]:.3f}, {fr[1.3]:.3f} ; {min(ground_vals):.1f} a {max(ground_vals):.1f} MeV")

    # F
    print("F. Sensibilite aux parametres publies")
    sens = {}
    for key, d in (("sig", 0.0007), ("gam", 0.011)):
        hi = dict(PAR); hi[key] += d; lo = dict(PAR); lo[key] -= d
        sens[key] = 0.5 * (window_fit(lambda x: v_string(x, hi), 0.40, 1.0)[0] - window_fit(lambda x: v_string(x, lo), 0.40, 1.0)[0])
    print(f"    sigma t0 +- 0,0007 : {sens['sig']:+.1f} MeV ; gamma +- 0,011 : {sens['gam']:+.1f} MeV (les deux sont correles dans l'article, matrice non reprise ici)")
    drs_lat = DRS_LAT * central / RS_LAT
    print(f"    erreur retenue pour le reseau : celle de l'article, +- {DRS_LAT} MeV sur 445, portee a la meme part de {central:.1f} : +- {drs_lat:.1f} MeV\n")
    check("F. sensibilite : +- 1,5 MeV pour l'erreur statistique de sigma, +- 1,3 MeV pour celle de gamma ; erreur retenue +- 6,8 MeV",
          1.0 < abs(sens["sig"]) < 2.0 and 0.8 < abs(sens["gam"]) < 2.0 and abs(drs_lat - 6.85) < 0.1, f"{sens['sig']:+.2f} et {sens['gam']:+.2f} MeV ; +- {drs_lat:.2f}")

    # G
    print("G. Budget d'erreur et critere ecrit avant le calcul")
    ft = lambda ap: rel.total(1e3 / sqrt(2 * pi * ap), None, "a")[0]
    ap0 = golden(ft, 0.6, 1.2, 30); tot0 = ft(ap0)
    s2 = tot0 / 28; curv = (ft(ap0 + 0.01) - 2 * tot0 + ft(ap0 - 0.01)) / 1e-4
    rs_spec = 1e3 / sqrt(2 * pi * ap0); drs_spec = rs_spec * sqrt(2 * s2 / curv) / (2 * ap0)
    allv = string_vals + ground_vals
    dwin = max(abs(v - central) for v in allv)
    dlat = sqrt(drs_lat ** 2 + dwin ** 2)
    n_crit = (central - rs_spec) / sqrt(dlat ** 2 + drs_spec ** 2)
    env = 0.5 * (1e3 / sqrt(2 * pi * 0.853) - 1e3 / sqrt(2 * pi * 0.903))
    n_prud = (central - rs_spec) / sqrt(dlat ** 2 + drs_spec ** 2 + env ** 2)
    n_pub = (RS_LAT - rs_spec) / sqrt(DRS_LAT ** 2 + drs_spec ** 2)
    rs_base = 1e3 / sqrt(2 * pi * b1.AP_BASE)
    n_base = (central - rs_base) / dlat
    print(f"    spectres (R137, refait ici) : {rs_spec:.1f} +- {drs_spec:.1f} MeV")
    print(f"    reseau, forme d'une corde, fenetre 0,4 a 1,0 fm : {central:.1f} MeV ; fenetres et niveau ajuste : de {min(allv):.1f} a {max(allv):.1f}, soit +- {dwin:.1f} ; avec l'erreur de l'article : +- {dlat:.1f} MeV")
    print(f"    ecart : {central - rs_spec:.1f} MeV, soit {n_crit:.2f} ecarts types. Seuils ecrits : plus de 3, ecart etabli ; moins de 2, accord.")
    print(f"    avec la tension publiee, 445 +- 6,7 : {n_pub:.2f} (R137). Avec l'enveloppe des familles de R133 du cote des spectres, +- {env:.1f} MeV : {n_prud:.2f}")
    print(f"    la valeur de la base, {rs_base:.1f} MeV, est a {n_base:.1f} ecarts types sous cette tension du reseau, et a 2,6 au-dessus des spectres (R137)")
    print("    hors budget : l'article n'a qu'une maille, donc pas de limite du continu ; d'autres ajustements publies donnent 467 a 482 MeV (R136 F),")
    print("    ce qui agrandit l'ecart ; du cote des spectres, la corde tournante est classique.\n")
    check("G. critere : 454,8 +- 8,1 MeV contre 423,7 +- 2,5, ecart de 31 MeV, 3,65 ecarts types, au-dessus du seuil de 3 ; 2,97 avec l'enveloppe des familles ; 2,97 avec la tension publiee",
          abs(rs_spec - 423.7) < 0.1 and abs(drs_spec - 2.5) < 0.1 and 7.6 < dlat < 8.6 and 3.4 < n_crit < 4.0 and n_crit > 3 and 2.8 < n_prud < 3.2 and 2.9 < n_pub < 3.05 and 2.7 < n_base < 3.3,
          f"{central:.1f} +- {dlat:.1f} MeV ; {n_crit:.2f} ; prudent {n_prud:.2f} ; publie {n_pub:.2f}")

    print("Verdict : avec la forme d'une corde, le potentiel publie du reseau donne 455 MeV et non 445 : la tension du reseau depend de")
    print("10 MeV de la forme ajustee. L'ecart aux spectres est de 31 MeV, 3,65 ecarts types : le seuil de 3 ecrit d'avance est franchi.")
    print("Avec l'erreur la plus prudente du cote des spectres il retombe a 2,97, juste sous le seuil : etabli au sens du critere, pas")
    print("au-dela. Aucun nombre de la CD n'entre dans ce calcul.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
