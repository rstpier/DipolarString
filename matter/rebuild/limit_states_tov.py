#!/usr/bin/env python3
"""R125 -- Les etats limites sous l'objection R124 : la branche des etoiles de substrat avec l'hydrostatique de la RG.

R124 a exclu le secteur scalaire a repere ; la seule gravitation coherente qui reste au modele est la RG (bootstrap).
Or la branche des etoiles de substrat (manuscrit, eq. star_structure) est integree avec la loi de force scalaire :
    dM/dr = 4 pi r^2 (u - u_0)/c^2,   dP/dr = -(u + P) G M/(r^2 c^2),   P = u/3,   surface P(R) = u_0/3  (A8),
soit, sans dimension (r = a x, M = b m, u = u_0 theta ; a = c^2/(4 sqrt(pi G u_0)), b = c^4/(16 G sqrt(pi G u_0))) :
    dm/dx = x^2 (theta - 1),   dtheta/dx = -theta m / x^2,   compacite r_s/R = m/(2x).
Le manuscrit le dit « newtonien » : il manque les deux facteurs de la RG, (M + 4 pi r^3 P_g/c^2) et 1/(1 - 2GM/rc^2).
Ici on les restaure (equation de Tolman-Oppenheimer-Volkoff) en gardant toutes les conventions du manuscrit (A8 :
masse et pression gravitantes = exces u - u_0, (u - u_0)/3 ; inertie = u + P ; surface a theta = 1) :
    dtheta/dx = -theta [m + x^3 (theta - 1)/3] / (x^2 (1 - m/(2x))).
On integre aussi chaque facteur seul pour voir lequel pese. Le manuscrit lui-meme annoncait que la reformulation
covariante « pourrait deplacer les nombres » ; voici de combien.

  A. Controle : le systeme du manuscrit redonne M_crit = 1,29e6 M_sun a theta_c = 14,5-14,8, r_s/R = 0,82.
  B. TOV-A8 : branche stable (premier maximum) a M_crit = 3,9e5 M_sun (x 0,30), theta_c = 7,6, r_s/R = 0,31 ;
     compacite maximale sur la branche stable 0,31 < 2/3 : AUCUNE configuration stable n'est sombre.
  C. Facteur metrique seul : 5,8e5 M_sun, 0,40 ; facteur de pression seul : 7,3e5 M_sun, 0,53 ; chacun suffit a
     interdire les objets sombres (< 2/3).
  D. Consequences. (i) Les « objets sombres stables » (ombre de trou noir sans horizon, region hachuree de la figure
     de branche, famille 2 sombre) n'existent que sous la loi de force scalaire exclue par R124. (ii) Sous la RG la plus
     compacte des etoiles stables a 1 + z = (1 - 0,31)^(-1/2) = 1,20 : surface brillante, l'argument de luminosite
     (Broderick-Narayan) frappe toute la famille, pas seulement Sgr A*. (iii) La coupure de la fonction de masse passe
     de log M = 6,11 a 5,59, sous la mediane 6,10 de l'ecole des petites masses et a 2 ordres de la masse dynamique
     (7,7). (iv) Les effondrements geles (familles 1 et 3) tenaient a c_eff -> 0 du secteur scalaire : sans lui, la RG
     forme un horizon en temps propre fini ; les deux premiers falsificateurs (pas de sursaut PBH, echos) perdent leur
     mecanisme. (v) La loi M_crit ~ l_1^2 survit (dimensionnelle), pas la noirceur.
  E. Coherence : 0,31 est sous la borne de Buchdahl (8/9) et sous celle des spheres de rayonnement auto-gravitantes
     (2GM/Rc^2 <= 0,49, Sorkin-Wald-Zhang 1981) ; le 0,82 newtonien la depasse : il n'etait pas realisable en RG.
  Verdict : les etats limites tombent avec le secteur scalaire. Ce qui reste est une sphere de rayonnement
  auto-gravitante a fond non gravitant (A8, pose), de masse ~ l_1^2 et lumineuse : pas d'objet sombre sans horizon,
  pas de mecanisme de gel. Critere d'arret d'A3.md : de la CD restent ses rapports hadroniques et ses resultats negatifs.
"""
import math

G, C, H, HBAR, ME = 6.67430e-11, 2.99792458e8, 6.62607015e-34, 1.054571817e-34, 9.1093837015e-31
MSUN = 1.98847e30
L1 = 2 * math.pi * HBAR / (ME * C) / 3
U0 = H * C / L1 ** 4
B = C ** 4 / (16 * G * math.sqrt(math.pi * G * U0))      # unite de masse, kg

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def rhs(model, x, m, th):
    if model == "newton":
        return x * x * (th - 1), -th * m / (x * x)
    den = 1 - m / (2 * x)
    if model == "pressure_only":
        return x * x * (th - 1), -th * (m + x ** 3 * (th - 1) / 3) / (x * x)
    if den <= 1e-9:
        return None
    if model == "metric_only":
        return x * x * (th - 1), -th * m / (x * x * den)
    return x * x * (th - 1), -th * (m + x ** 3 * (th - 1) / 3) / (x * x * den)      # TOV-A8 complet

def integrate(theta_c, model, h=1e-3):
    eps = 1e-4
    m = (theta_c - 1) * eps ** 3 / 3; th = theta_c - theta_c * (theta_c - 1) * eps ** 2 / 6; x = eps
    while x < 1e4:
        step = min(h, 0.02 * x + 1e-4)
        k1 = rhs(model, x, m, th)
        if k1 is None: return None
        k2 = rhs(model, x + step / 2, m + step / 2 * k1[0], th + step / 2 * k1[1])
        if k2 is None: return None
        k3 = rhs(model, x + step / 2, m + step / 2 * k2[0], th + step / 2 * k2[1])
        if k3 is None: return None
        k4 = rhs(model, x + step, m + step * k3[0], th + step * k3[1])
        if k4 is None: return None
        m_new = m + step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        th_new = th + step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if th_new <= 1.0:
            f = (th - 1) / (th - th_new)
            return x + f * step, m + f * (m_new - m)
        x, m, th = x + step, m_new, th_new
    return None

def branch(model, n=121):
    thetas = [1.01 * 10 ** (5 * i / (n - 1)) for i in range(n)]
    rows = []
    for tc in thetas:
        r = integrate(tc, model)
        if r is not None:
            rows.append((tc, r[1], r[0], r[1] / (2 * r[0])))
    i = max(range(len(rows)), key=lambda k: rows[k][1])
    tc, ms, xs, comp = rows[i]
    comp_stable = max(row[3] for row in rows[:i + 1])
    bracketed = 0 < i < len(rows) - 1
    return dict(M=ms * B / MSUN, theta_c=tc, x_s=xs, comp=comp, comp_stable=comp_stable, bracketed=bracketed, rows=rows, i=i)

def main():
    print("R125 -- les etoiles de substrat avec l'hydrostatique de la RG (TOV), conventions A8 du manuscrit\n")
    print(f"  l_1 = {L1:.4e} m ; u_0 = {U0:.3e} J/m^3 ; unite de masse b = {B/MSUN:.4e} M_sun\n")
    res = {k: branch(k) for k in ("newton", "tov_a8", "metric_only", "pressure_only")}
    names = {"newton": "manuscrit (newtonien)", "tov_a8": "TOV-A8 (deux facteurs)", "metric_only": "facteur 1/(1 - 2GM/rc^2) seul", "pressure_only": "facteur M + 4 pi r^3 P_g/c^2 seul"}
    print("A-C. Les branches : premier maximum de M(theta_c)")
    for k in ("newton", "tov_a8", "metric_only", "pressure_only"):
        r = res[k]
        print(f"    {names[k]:36s} M_crit = {r['M']:.3e} M_sun (log {math.log10(r['M']):.2f}) a theta_c = {r['theta_c']:.2f}, x_s = {r['x_s']:.3f}, r_s/R = {r['comp']:.3f} ; compacite max sur la branche stable {r['comp_stable']:.3f} ; maximum encadre : {r['bracketed']}")
    print("    points de la table du manuscrit (theta_c = 1,05 ; 1,9 ; 3,3 ; 5,9 ; 14,5) :")
    for tc in (1.05, 1.9, 3.3, 5.9, 14.5):
        rn, rt = integrate(tc, "newton"), integrate(tc, "tov_a8")
        print(f"      theta_c = {tc:5.2f} : newtonien M = {rn[1]*B/MSUN:.2e}, r_s/R = {rn[1]/(2*rn[0]):.3f} ; TOV-A8 M = {rt[1]*B/MSUN:.2e}, r_s/R = {rt[1]/(2*rt[0]):.3f}")
    print()
    n, t = res["newton"], res["tov_a8"]
    check("A. controle : le systeme du manuscrit redonne M_crit = 1,29e6 M_sun, theta_c = 14,5-15, r_s/R = 0,82",
          abs(n["M"] / 1.29e6 - 1) < 0.02 and 13 < n["theta_c"] < 16 and 0.80 < n["comp"] < 0.83, f"{n['M']:.3e}, {n['theta_c']:.2f}, {n['comp']:.3f}")
    check("B. TOV-A8 : branche stable encadree, M_crit = 3,9e5 M_sun (x 0,30), compacite maximale 0,31 < 2/3 : aucune configuration stable sombre",
          t["bracketed"] and 0.25 < t["M"] / n["M"] < 0.35 and t["comp_stable"] < 2 / 3 and 0.28 < t["comp_stable"] < 0.34, f"{t['M']:.3e}, x {t['M']/n['M']:.2f}, {t['comp_stable']:.3f}")
    check("C. chaque facteur de la RG seul interdit deja les objets sombres (compacites 0,40 et 0,53 < 2/3)",
          res["metric_only"]["comp_stable"] < 2 / 3 and res["pressure_only"]["comp_stable"] < 2 / 3, f"{res['metric_only']['comp_stable']:.3f}, {res['pressure_only']['comp_stable']:.3f}")

    print("D. Consequences")
    z = (1 - t["comp_stable"]) ** -0.5 - 1
    logM = math.log10(t["M"])
    print(f"    (i)  objets sombres stables (r_s/R > 2/3) : aucun sous la RG ; ils n'existaient que sous la loi de force scalaire exclue par R124")
    print(f"    (ii) redshift de surface de l'etoile stable la plus compacte : 1 + z = {1+z:.3f} (manuscrit : 2,3) -> surface brillante, l'argument Broderick-Narayan frappe toute la famille 2")
    print(f"    (iii) coupure de la fonction de masse : log M_crit = {logM:.2f} (manuscrit 6,11) ; mediane LRD petites masses 6,10 ; masse dynamique A2744-QSO1 : 7,7")
    print(f"    (iv) effondrements geles : le gel tenait a c_eff -> 0 du secteur scalaire ; sous la RG un horizon se forme en temps propre fini -> les falsificateurs 1 et 2 perdent leur mecanisme")
    print(f"    (v)  M_crit ~ l_1^2 survit (dimensionnel : b ~ l_1^2), la noirceur non\n")
    check("D. sous la RG : 1 + z = 1,20 (surface brillante), log M_crit = 5,59 < 6,10 (sous la mediane LRD), pas de mecanisme de gel",
          1.15 < 1 + z < 1.25 and 5.5 < logM < 5.7, f"1+z = {1+z:.3f}, log M = {logM:.2f}")

    print("E. Coherence avec les bornes de la RG")
    print(f"    compacite TOV-A8 {t['comp_stable']:.3f} < 0,493 (spheres de rayonnement, Sorkin-Wald-Zhang 1981) < 8/9 (Buchdahl) ; newtonien {n['comp_stable']:.3f} > 0,493 : irrealisable en RG\n")
    check("E. 0,31 sous la borne des spheres de rayonnement (0,49) ; le 0,82 newtonien la depasse", t["comp_stable"] < 0.493 < n["comp_stable"], f"{t['comp_stable']:.3f} / {n['comp_stable']:.3f}")

    print("Verdict : les etats limites tombent avec le secteur scalaire. Sous la seule gravitation coherente qui reste (RG), la branche")
    print("existe mais n'est pas sombre (r_s/R <= 0,31), M_crit = 3,9e5 M_sun, surface a 1 + z = 1,2, pas de gel ; l'ombre sans horizon,")
    print("la famille 2 sombre, l'argument Sgr A* et les falsificateurs 1-3 etaient des consequences de la loi de force scalaire exclue")
    print("par R124. Il reste une sphere de rayonnement a fond non gravitant (A8, pose), lumineuse, de masse ~ l_1^2.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
