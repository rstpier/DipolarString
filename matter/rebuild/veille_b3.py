#!/usr/bin/env python3
"""R135 -- Question B3, etat de la litterature au 22 novembre 2026 : les deux tensions sont-elles la meme ?

Choix fait par delegation (R134) : avant tout calcul sur B3, chercher si la reponse est publiee (critere d'arret de la
feuille de route). Recherche partielle, huit requetes. Ce script confronte les nombres trouves a ceux de R133 et R134.

  A. Cote spectres. Badalian et Bakker (2019) : pente de la trajectoire dominante 2 pi sigma = 1,13(1) GeV^2, soit
     sigma = 0,180(2) GeV^2 et sqrt(sigma) = 424 +- 2 MeV. Compatible avec l'enveloppe de R133, 426 +- 6 MeV.
  B. Cote reseau. Bulava et al. (2024), 2+1 saveurs, masses physiques : sqrt(sigma) = 445(3)(6) MeV. Ecart aux
     spectres : 21 MeV, soit 3,0 ecarts types avec l'erreur publiee des spectres, 2,1 avec l'enveloppe de R133.
     En pente : le reseau donne 2 pi sigma = 1,24 +- 0,04 GeV^2, 10 % au-dessus de 1,13.
  C. Brisure de la corde : 1,22 fm (quarks legers) et 1,29 fm (etrange), meme source. R133 avait pris 1,25 fm de
     memoire : confirme.
  D. Echelles de Sommer. r1 = 0,3127(30) fm (FLAG 2021, 2+1 saveurs) ; r0 = 0,474(14) fm (FLAG 2021, 2+1+1) ; mesure
     plus recente r1 = 0,3104(33) fm et r0 / r1 = 1,522(13), soit r0 = 0,472 fm. R133 avait pris 0,472 de memoire :
     confirme. r0 sqrt(sigma) = 1,07 +- 0,04 contre 1,178 pour une corde seule.
  E. Baryons : pente universelle publiee 0,95 GeV^-2 (Sonnenschein et Weissman 2015), contre 0,92 a 0,99 en R133.
     La pente des mesons, 0,884, et les masses d'extremite restent citees de memoire : non confirmees par la recherche.
  F. Mecanismes publies. Aplatissement du potentiel par la creation de paires, qui abaisse les pentes des
     trajectoires (Badalian et Bakker 2019) ; corrections quantiques de la corde tournante a extremites massives,
     qui deplacent l'ordonnee et non la pente (Sonnenschein et Weissman 2018). Aucune source trouvee ne confronte
     directement la tension de 2024 aux pentes de Regge.
  Verdict : B3 n'est pas tranchee dans ce que la recherche a trouve. L'ecart vaut 2 a 3 ecarts types selon le budget
  d'erreur du cote des spectres. Une piste publiee existe et doit etre lue avant tout calcul.
"""
from math import pi, sqrt

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = 197.3269804                        # MeV fm

def main():
    print("R135 -- question B3 : etat de la litterature au 22 novembre 2026\n")

    # A
    slope, dslope = 1.13, 0.01             # GeV^2, Badalian et Bakker, arXiv:1901.10280
    sig = slope / (2 * pi); dsig = dslope / (2 * pi)
    rs_bb = 1e3 * sqrt(sig); drs_bb = rs_bb * dsig / (2 * sig)
    env, denv = 426.0, 6.1                 # R133
    print("A. Cote spectres")
    print(f"    2 pi sigma = {slope}({int(dslope*100)}) GeV^2 -> sigma = {sig:.4f} +- {dsig:.4f} GeV^2, sqrt(sigma) = {rs_bb:.1f} +- {drs_bb:.1f} MeV")
    print(f"    enveloppe de R133 : {env} +- {denv} MeV ; ecart {env - rs_bb:+.1f} MeV, {abs(env - rs_bb)/sqrt(denv**2 + drs_bb**2):.1f} ecart type\n")
    check("A. litterature : sigma = 0,180(2) GeV^2, sqrt(sigma) = 424 +- 2 MeV, compatible avec l'enveloppe de R133",
          abs(sig - 0.1798) < 3e-4 and abs(rs_bb - 424.1) < 0.5 and abs(drs_bb - 1.9) < 0.2 and abs(env - rs_bb) < denv, f"{rs_bb:.1f} +- {drs_bb:.1f} MeV")

    # B
    rs_lat, stat, syst = 445.0, 3.0, 6.0   # Bulava et al., Phys. Lett. B 854 (2024) 138754, arXiv:2403.00754
    drs_lat = sqrt(stat ** 2 + syst ** 2)
    d = rs_lat - rs_bb
    n_pub = d / sqrt(drs_lat ** 2 + drs_bb ** 2)
    n_env = (rs_lat - env) / sqrt(drs_lat ** 2 + denv ** 2)
    slope_lat = 2 * pi * (rs_lat * 1e-3) ** 2; dslope_lat = 2 * slope_lat * drs_lat / rs_lat
    print("B. Cote reseau")
    print(f"    sqrt(sigma) = {rs_lat:.0f}({stat:.0f})({syst:.0f}) MeV, soit +- {drs_lat:.1f}")
    print(f"    reseau - spectres = {d:+.1f} MeV : {n_pub:.1f} ecarts types avec l'erreur publiee des spectres, {n_env:.1f} avec l'enveloppe de R133")
    print(f"    en pente : le reseau donne 2 pi sigma = {slope_lat:.3f} +- {dslope_lat:.3f} GeV^2, soit {100*(slope_lat/slope-1):+.0f} % au-dessus de {slope}\n")
    check("B. reseau 445(3)(6) MeV : ecart aux spectres de 21 MeV, 3,0 ou 2,1 ecarts types selon le budget d'erreur ; pente de 10 % plus forte",
          abs(drs_lat - 6.7) < 0.1 and abs(d - 20.9) < 0.6 and 2.8 < n_pub < 3.2 and 1.9 < n_env < 2.3 and 0.08 < slope_lat / slope - 1 < 0.12, f"{n_pub:.1f} ; {n_env:.1f}")

    # C
    r_light, r_strange, r_used = 1.22, 1.29, 1.25
    print("C. Brisure de la corde")
    print(f"    {r_light} fm (quarks legers), {r_strange} fm (etrange) ; R133 avait pris {r_used} fm de memoire : dans l'intervalle")
    e_break = (rs_lat ** 2 / HBARC) * r_light
    print(f"    energie stockee a {r_light} fm avec la tension du reseau : {e_break:.0f} MeV\n")
    check("C. brisure a 1,22 et 1,29 fm : la valeur de 1,25 fm prise en R133 est confirmee", r_light < r_used < r_strange and 1200 < e_break < 1250, f"{r_light} et {r_strange} fm")

    # D
    r1, dr1 = 0.3127, 0.0030               # FLAG 2021, Nf = 2+1
    r0_flag, dr0_flag = 0.474, 0.014       # FLAG 2021, Nf = 2+1+1 (ETM 14)
    r1_new, ratio, dratio = 0.3104, 1.522, 0.013
    r0_new = r1_new * ratio
    k0 = r0_flag * rs_lat / HBARC; dk0 = k0 * sqrt((dr0_flag / r0_flag) ** 2 + (drs_lat / rs_lat) ** 2)
    k_string = sqrt(1.65 - pi / 12)
    k1 = r1 * rs_lat / HBARC
    print("D. Echelles de Sommer")
    print(f"    r0 = {r0_flag}({int(dr0_flag*1000)}) fm (FLAG 2021) ; r1 x (r0 / r1) = {r1_new} x {ratio} = {r0_new:.4f} fm ; R133 avait pris 0,472 fm de memoire : confirme")
    print(f"    r0 sqrt(sigma) = {k0:.3f} +- {dk0:.3f} ; corde seule : {k_string:.3f}, a {(k_string - k0)/dk0:.1f} ecarts types")
    print(f"    r1 sqrt(sigma) = {k1:.3f} : aucune contrepartie pour la corde d'Arvis (R133 E)\n")
    check("D. r0 = 0,474(14) fm et 0,472 fm par deux voies ; r0 sqrt(sigma) = 1,07 +- 0,04 contre 1,178 pour la corde seule",
          abs(r0_new - 0.4724) < 5e-4 and abs(k0 - 1.069) < 0.004 and 0.03 < dk0 < 0.04 and 2.5 < (k_string - k0) / dk0 < 3.5, f"{k0:.3f} +- {dk0:.3f}")

    # E
    bar_pub = 0.95
    n_r133, d_r133 = (0.964, 0.993), (0.924, 0.929)
    print("E. Baryons")
    print(f"    pente universelle publiee : {bar_pub} GeV^-2 ; R133 : N {n_r133[0]} a {n_r133[1]}, Delta {d_r133[0]} a {d_r133[1]} ; ecarts de {100*(d_r133[0]/bar_pub-1):+.0f} a {100*(n_r133[1]/bar_pub-1):+.0f} %")
    print("    pente des mesons 0,884 et masses d'extremite 60, 220 MeV : citees de memoire, non confirmees par la recherche\n")
    check("E. pente des baryons publiee 0,95 GeV^-2 ; celles de R133 en sont a moins de 5 %", all(abs(x / bar_pub - 1) < 0.05 for x in n_r133 + d_r133), "0,95")

    # F
    print("F. Mecanismes publies et etat de la question")
    print("    aplatissement du potentiel par creation de paires : abaisse les pentes des trajectoires (Badalian et Bakker 2019)")
    print("    corrections quantiques de la corde tournante a extremites massives : deplacent l'ordonnee, pas la pente (Sonnenschein et Weissman 2018)")
    print("    aucune source trouvee ne confronte directement la tension de 2024 aux pentes de Regge")
    open_question = n_env < 3 <= round(n_pub, 1) + 0.05
    print(f"    ecart entre 2 et 3 ecarts types selon le budget d'erreur : {open_question}\n")
    check("F. B3 non tranchee dans ce que la recherche a trouve : ecart de 2,1 a 3,0 ecarts types, une piste publiee a lire", open_question, "2,1 a 3,0")

    print("Verdict : la question B3 est vivante. L'ecart vaut 21 MeV, deux a trois ecarts types selon le budget d'erreur du cote des spectres,")
    print("que R133 a montre domine par le choix des familles. Une piste publiee, l'aplatissement du potentiel, doit etre lue avant tout calcul.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
