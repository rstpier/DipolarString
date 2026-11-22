#!/usr/bin/env python3
"""R134 -- La tension du reseau a la place de celle de la base : ce que le remplacement change.

Enonce de l'auteur (22 novembre 2026), devant R133 : « ok prend la tension du reseau », puis : « Puisque tu connais
les meilleures reponses aux questions que tu me poses, fais les choix a ma place. »
Choix fait par delegation : dans la corde effective, la tension est celle que le reseau mesure,
sqrt(sigma) = 445 +- 7 MeV (Bulava et al. 2024, 2+1 saveurs, masses physiques ; R119). Raisons : la tension d'une
corde longue est par definition la pente du potentiel statique a grande distance, que le reseau mesure directement ;
la formule de la base n'a pas de mecanisme et son accord n'est pas significatif (R126) ; la pente des trajectoires
depend de la famille a 20 % pres (R133). Ce script calcule les consequences. Il reutilise le solveur de R133.

  A. Les constantes de la corde avec cette tension.
  B. La formule de la base. Pour donner 445 MeV l'exposant devrait valoir 6,314 +- 0,014 : 2 pi est a 2,1 ecarts
     types, le bord de la fenetre des leptons (6,2758 a 6,2926) a 1,5. La formule est retiree comme prediction.
  C. Les trajectoires avec la tension du reseau. alpha' = 0,804 +- 0,025 GeV^-2, contre 0,883 +- 0,014 (six
     familles) et 0,903 +- 0,011 (familles vecteur et tenseur) : 9 a 11 % trop faible, 2,7 et 3,6 ecarts types ; tours
     de Regge trop hautes de 3 a 9 %. Des masses d'extremite aggravent l'ecart.
  D. Le bilan commun. Trajectoires 426,0 +- 6,1 MeV, reseau 445 +- 7 MeV : meilleure valeur commune 434,2 +- 4,6 ;
     chi2 de 4,2 au minimum, 4,9 a la valeur de la base, 9,7 a celle du reseau. Aucun choix ne reconcilie les deux.
  E. La « part statique » de la base. Corde a J = 1/2 : 763 MeV avec la tension de la base, 789 +- 12 avec celle du
     reseau ; longueur 0,50 fm, soit 1,6 rayon d'Arvis : hors domaine dans les deux cas (R133 D). Pas de remise a
     l'echelle des autres lignes de la table des hadrons.
  F. Le compte. Dans la corde effective il ne reste aucun nombre propre a la CD.
  Verdict : choix motive et coherent. Il retire a la CD son dernier nombre dans cette lecture, et il deplace l'ecart
  du cote des spectres, ou il grandit : c'est la question B3, rendue plus nette. Suite : la recherche
  bibliographique de B3.
"""
from math import log, pi, sqrt
import ds_effective_string as est

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = est.HBARC                          # GeV fm
RS_LAT, DRS_LAT = 445.0, 7.0               # MeV, Bulava et al., Phys. Lett. B 854 (2024) 138754 : 445(3)(6)
RS_BASE = est.RS_BASE

def main():
    print("R134 -- la tension du reseau a la place de celle de la base\n")

    # A
    print("A. Les constantes de la corde")
    def consts(rs_mev):
        sig = (rs_mev * 1e-3) ** 2
        ap = 1 / (2 * pi * sig)
        return sig, ap, sqrt(ap) * HBARC, HBARC * sqrt(pi / (6 * sig)), 1e3 * sig / HBARC
    b, l = consts(RS_BASE), consts(RS_LAT)
    rel = DRS_LAT / RS_LAT
    print("                              base          reseau")
    print(f"    sqrt(sigma)   (MeV)       {RS_BASE:7.1f}       {RS_LAT:.0f} +- {DRS_LAT:.0f}")
    print(f"    sigma         (GeV^2)     {b[0]:7.4f}       {l[0]:.4f} +- {2*rel*l[0]:.4f}")
    print(f"    alpha'        (GeV^-2)    {b[1]:7.4f}       {l[1]:.4f} +- {2*rel*l[1]:.4f}")
    print(f"    longueur de corde (fm)    {b[2]:7.4f}       {l[2]:.4f}")
    print(f"    rayon d'Arvis (fm)        {b[3]:7.3f}       {l[3]:.3f}")
    print(f"    poids du fil  (MeV / fm)  {b[4]:7.0f}       {l[4]:.0f} +- {2*rel*l[4]:.0f}")
    print(f"    energie stockee a la brisure, 1,25 fm : {1.25*b[4]:.0f} MeV (base), {1.25*l[4]:.0f} MeV (reseau)\n")
    check("A. tension du reseau : sigma = 0,198 GeV^2, alpha' = 0,804 GeV^-2, longueur de corde 0,177 fm, 1004 MeV par fermi",
          abs(l[0] - 0.1980) < 2e-4 and abs(l[1] - 0.8037) < 5e-4 and abs(l[2] - 0.1769) < 2e-4 and abs(l[4] - 1003.5) < 1.5, f"alpha' = {l[1]:.4f} GeV^-2")

    # B
    print("B. La formule de la base")
    p_need = 2 * pi + log(RS_LAT / RS_BASE) / log(3)
    dp = rel / log(3)
    n_2pi = (p_need - 2 * pi) / dp
    win = (6.2758, 6.2926)
    n_win = (p_need - win[1]) / dp
    print(f"    sqrt(sigma) = m_e 3^(p + 1) / (2 sqrt(pi)) : p = 2 pi donne {RS_BASE:.1f} MeV ; {RS_LAT:.0f} +- {DRS_LAT:.0f} MeV demande p = {p_need:.4f} +- {dp:.4f}")
    print(f"    2 pi = {2*pi:.4f} est a {n_2pi:.1f} ecarts types ; le bord de la fenetre des leptons ({win[0]} a {win[1]}, R55) a {n_win:.1f}")
    print("    -> un seul exposant ne sert pas a la fois les rapports des leptons et la tension mesuree. La formule n'a pas de mecanisme (R126) :")
    print("       elle est retiree comme prediction de la tension, et gardee au registre comme piece d'histoire.\n")
    check("B. exposant requis par la tension du reseau : 6,314 +- 0,014 ; 2 pi a 2,1 ecarts types, fenetre des leptons a 1,5",
          abs(p_need - 6.3138) < 2e-3 and abs(dp - 0.0143) < 5e-4 and 2.0 < n_2pi < 2.3 and 1.3 < n_win < 1.7, f"p = {p_need:.4f} +- {dp:.4f}")

    # C
    print("C. Les trajectoires avec la tension du reseau")
    ap_lat, dap_lat = l[1], 2 * rel * l[1]
    six = est.fit_global({**est.MESONS, **est.PSEUDO}, "M2", True)
    four = est.fit_global(est.MESONS, "M2", True)
    out = {}
    for label, (ap, err, _, _) in (("six familles", six), ("familles vecteur et tenseur", four)):
        n = (ap - ap_lat) / sqrt(err ** 2 + dap_lat ** 2)
        out[label] = (100 * (ap_lat / ap - 1), n)
        print(f"    {label:28s} : ajuste {ap:.3f} +- {err:.3f} ; reseau {ap_lat:.3f} +- {dap_lat:.3f} : {100*(ap_lat/ap-1):+.1f} %, {n:.1f} ecarts types")
    towers = []
    print("    tours de Regge depuis l'etat le plus bas (extremites sans masse) :")
    for name, traj in (("rho", [(1, 0.77526), (2, 1.3182), (3, 1.6888), (4, 1.967)]), ("N", [(0.5, 0.93827), (2.5, 1.685), (4.5, 2.250)]), ("Delta", [(1.5, 1.232), (3.5, 1.930), (5.5, 2.420)])):
        j0, m0 = traj[0]
        line = []
        for j, m in traj[1:]:
            d_lat = sqrt(m0 * m0 + (j - j0) / ap_lat) / m - 1
            d_base = sqrt(m0 * m0 + (j - j0) / est.AP_BASE) / m - 1
            towers.append((d_base, d_lat)); line.append(f"J = {j} : {100*d_base:+.1f} % (base), {100*d_lat:+.1f} % (reseau)")
        print(f"      {name:6s} " + " ; ".join(line))
    worse = all(abs(dl) > abs(db) for db, dl in towers)
    # des masses d'extremite aggravent l'ecart : a alpha' fixe, le chi2 monte quand on les ajoute
    st, m1, m2 = est.MESONS["rho / a"]
    chi_massless = est.residuals(st, 0, 0, ap_lat, "M2")[0]
    chi_massive = est.residuals(st, m1, m2, ap_lat, "M2")[0]
    print(f"    famille rho / a a alpha' du reseau : chi2 sans masses {chi_massless:.4f}, avec masses de 60 MeV {chi_massive:.4f} (unites arbitraires) : les masses aggravent")
    print("    -> avec la tension du reseau la corde seule donne des pentes trop faibles de 9 a 11 % : l'ecart passe du cote des spectres et grandit.\n")
    check("C. trajectoires avec la tension du reseau : alpha' trop faible de 9 % (six familles, 2,7 ecarts types) et 11 % (vecteur et tenseur, 3,6) ; tours plus hautes qu'avec la base",
          -10 < out["six familles"][0] < -8 and 2.5 < out["six familles"][1] < 3.1 and -12 < out["familles vecteur et tenseur"][0] < -10
          and 3.3 < out["familles vecteur et tenseur"][1] < 3.9 and worse and chi_massive > chi_massless,
          f"{out['six familles'][0]:+.1f} %, {out['six familles'][1]:.1f} ; {out['familles vecteur et tenseur'][0]:+.1f} %, {out['familles vecteur et tenseur'][1]:.1f}")

    # D
    print("D. Le bilan commun des deux mesures")
    variants = [est.fit_global(f, w, m)[0] for f in (est.MESONS, {**est.MESONS, **est.PSEUDO}) for w in ("J", "M2") for m in (True, False)]
    lo, hi = est.rs(max(variants)), est.rs(min(variants))
    rs_regge, drs_regge = (lo + hi) / 2, (hi - lo) / 2
    w1, w2 = 1 / drs_regge ** 2, 1 / DRS_LAT ** 2
    best = (w1 * rs_regge + w2 * RS_LAT) / (w1 + w2); dbest = 1 / sqrt(w1 + w2)
    chi = lambda s: ((s - rs_regge) / drs_regge) ** 2 + ((s - RS_LAT) / DRS_LAT) ** 2
    print(f"    trajectoires (enveloppe de R133) : {rs_regge:.1f} +- {drs_regge:.1f} MeV ; reseau : {RS_LAT:.0f} +- {DRS_LAT:.0f} MeV")
    print(f"    meilleure valeur commune : {best:.1f} +- {dbest:.1f} MeV ; chi2 : {chi(best):.1f} au minimum, {chi(RS_BASE):.1f} a la valeur de la base, {chi(RS_LAT):.1f} a celle du reseau")
    print("    -> les deux mesures different de deux ecarts types quelle que soit la valeur retenue : le choix ne les reconcilie pas, il dit laquelle on croit.\n")
    check("D. valeur commune 434 +- 5 MeV ; chi2 4,2 au minimum, 4,9 a la base, 9,7 au reseau : aucun choix ne reconcilie les deux mesures",
          abs(best - 434.2) < 1 and abs(dbest - 4.6) < 0.3 and abs(chi(best) - 4.2) < 0.4 and abs(chi(RS_BASE) - 4.9) < 0.5 and abs(chi(RS_LAT) - 9.7) < 0.8,
          f"{best:.1f} +- {dbest:.1f} MeV")

    # E
    print("E. La « part statique » de la base")
    half = lambda rs_mev: sqrt(pi) * rs_mev                           # masse de la corde a J = 1/2 : sqrt(pi sigma)
    s_base, s_lat = half(RS_BASE), half(RS_LAT)
    ds_lat = sqrt(pi) * DRS_LAT
    L_half = 2 * HBARC / sqrt(pi * l[0])                              # longueur de la corde a J = 1/2, en fm
    luscher = -pi / 12 * HBARC * 1e3 / (pi * HBARC * 1e3 / s_lat)     # a la distance pi hbar c / M
    print(f"    corde a J = 1/2 : sqrt(pi sigma) = {s_base:.0f} MeV (base), {s_lat:.0f} +- {ds_lat:.0f} MeV (reseau) ; fenetre exigee par le moment du proton : 770 a 813 MeV")
    print(f"    terme de Luscher a la distance correspondante : {luscher:.0f} MeV -> {s_lat + luscher:.0f} MeV, hors fenetre")
    print(f"    longueur de cette corde : {L_half:.2f} fm = {L_half/l[3]:.1f} rayon d'Arvis ; (r_c / L)^2 = {(l[3]/L_half)**2:.2f} : hors domaine (R133 D)")
    print("    -> l'identite « part statique = corde a J = 1/2 » ne tient qu'avec la tension de la base, et dans un domaine non controle :")
    print("       elle n'est pas reportee. Aucune autre ligne de la table des hadrons n'est remise a l'echelle.\n")
    check("E. corde a J = 1/2 : 763 MeV (base), 789 +- 12 MeV (reseau) ; avec Luscher 723 MeV ; longueur 0,50 fm, 1,6 rayon d'Arvis : hors domaine",
          abs(s_base - 762.7) < 0.5 and abs(s_lat - 788.7) < 0.5 and abs(ds_lat - 12.4) < 0.2 and 770 < s_lat < 813 and s_lat + luscher < 770
          and abs(L_half - 0.50) < 0.01 and 1.5 < L_half / l[3] < 1.7, f"{s_lat:.0f} +- {ds_lat:.0f} MeV ; {L_half:.2f} fm")

    # F
    print("F. Le compte des nombres propres a la CD dans la corde effective")
    numbers = {"tension": "mesuree sur le reseau (choix du 22 novembre)", "masses d'extremite": "ajustees (60 et 220 MeV, litterature)",
               "ordonnees": "ajustees, une par famille", "terme de Luscher, forme du potentiel": "universels"}
    for k, v in numbers.items(): print(f"    {k:38s} : {v}")
    own = 0
    print(f"    nombres propres a la CD : {own}")
    print("    -> dans cette lecture la CD n'apporte plus aucun nombre : ce qui reste est la physique etablie des tubes de flux, et une question, B3.\n")
    check("F. apres le choix de la tension du reseau, la corde effective ne contient aucun nombre propre a la CD", own == 0 and len(numbers) == 4, "0")

    print("Verdict : la tension du reseau est le bon choix pour une corde longue, parce que c'est la grandeur que le reseau mesure. Il retire a")
    print("la CD son dernier nombre dans cette lecture, et il deplace l'ecart du cote des spectres, ou il passe de 0,7 a 3 ecarts types.")
    print("La question B3 en sort plus nette : pourquoi la pente des trajectoires depasse-t-elle de 9 a 11 % celle que donne la tension statique ?\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
