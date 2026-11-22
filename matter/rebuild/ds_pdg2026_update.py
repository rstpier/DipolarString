#!/usr/bin/env python3
"""R137 -- Les masses du PDG 2026 lues a la source : controle des valeurs citees de memoire, mise a jour de l'ajustement.

Le 22 novembre 2026 l'auteur a autorise le telechargement depuis le PDG. Lus : les tables de synthese des mesons et
des baryons, edition 2026 (F. Takahashi et al., Particle Data Group, Int. J. Mod. Phys. A 41, 2630011 (2026)), et
les fiches de rho5(2350), f6(2510), K5*(2380), D2(2740)0 et des « autres etats ». Ce script reutilise le calcul de R136.

  A. Controle des masses que R133 citait de memoire : vingt-huit etats, ecart maximal a la table.
  B. Deux corrections : Delta(2420), dont la masse de Breit-Wigner est estimee a 2450 MeV et non 2420 (R120, R133) ;
     K*(892) charge, 891,88 MeV et non 891,67.
  C. L'ajustement universel de l'article de 2014 refait avec les masses de 2026, memes neuf familles, memes etats :
     pente, mesure d'ecart, tension des spectres. Les masses ont bouge de 30 MeV au plus (a4) en douze ans.
  D. La tension de la base et celle du reseau contre cette mise a jour.
  E. Les baryons avec la masse corrigee.
  Verdict : les valeurs citees de memoire sont confirmees a 0,3 MeV pres, sauf deux, corrigees. Avec les masses de
  2026 la tension des spectres ne bouge pas de facon significative.
"""
from math import pi, sqrt
import ds_b1_reproduction as b1

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# (valeur citee de memoire en R133 ou R120, valeur du PDG 2026, incertitude), en MeV
CONTROL = {
    "rho(770)": (775.26, 775.26, 0.23), "a2(1320)": (1318.2, 1318.2, 0.6), "rho3(1690)": (1688.8, 1688.8, 2.1), "a4(1970)": (1967.0, 1967.0, 16.0),
    "rho5(2350)": (2330.0, 2330.0, 35.0), "a6(2450)": (2450.0, 2450.0, 130.0),
    "omega(782)": (782.66, 782.66, 0.13), "f2(1270)": (1275.4, 1275.4, 0.8), "omega3(1670)": (1667.0, 1667.0, 4.0), "f4(2050)": (2018.0, 2018.0, 11.0),
    "f6(2510)": (2470.0, 2470.0, 50.0),
    "K*(892)+-": (891.67, 891.88, 0.23), "K2*(1430)+-": (1427.3, 1427.3, 1.5), "K3*(1780)": (1779.0, 1779.0, 8.0), "K4*(2045)": (2048.0, 2048.0, 9.0),
    "K5*(2380)": (2382.0, 2382.0, 24.0),
    "phi(1020)": (1019.461, 1019.460, 0.016), "f2'(1525)": (1517.4, 1517.3, 2.4), "phi3(1850)": (1854.0, 1854.0, 7.0),
    "b1(1235)": (1229.5, 1229.5, 3.2), "pi2(1670)": (1670.6, 1670.6, 2.9), "b3(2030)": (2032.0, 2032.0, 12.0), "pi4(2250)": (2250.0, 2250.0, 15.0),
    "h1(1170)": (1166.0, 1166.0, 6.0), "eta2(1645)": (1617.0, 1617.0, 5.0), "h3(2025)": (2025.0, 2025.0, 20.0), "eta4(2330)": (2328.0, 2328.0, 38.0),
    "N(1680) 5/2+": (1685.0, 1685.0, 5.0), "N(2220) 9/2+": (2250.0, 2250.0, 50.0), "Delta(1950) 7/2+": (1930.0, 1930.0, 18.0),
    "Delta(2420) 11/2+": (2420.0, 2450.0, 150.0),
}

# L'ajustement universel de l'article, memes etats, masses du PDG 2026 (GeV). D2(2740)0 : 2747 +- 6 MeV, lu dans sa
# fiche (omis de la table de synthese ; l'article de 2014 prenait 2737, la premiere mesure).
TRAJ_2026 = {
    "pi/b":  ([(1, 1.2295), (2, 1.6706), (3, 2.032), (4, 2.250)], "u", "u"),
    "rho/a": ([(1, 0.77526), (2, 1.3182), (3, 1.6888), (4, 1.967), (5, 2.330), (6, 2.450)], "u", "u"),
    "eta/h": ([(0, 0.547862), (1, 1.166), (2, 1.617), (3, 2.025), (4, 2.328)], "u", "u"),
    "omega/f": ([(1, 0.78266), (2, 1.2754), (3, 1.667), (4, 2.018), (5, 2.250), (6, 2.470)], "u", "u"),
    "K*":    ([(1, 0.89188), (2, 1.4273), (3, 1.779), (4, 2.048), (5, 2.382)], "u", "s"),
    "phi":   ([(1, 1.019460), (2, 1.5173), (3, 1.854)], "s", "s"),
    "D":     ([(0, 1.86484), (1, 2.4221), (2, 2.747)], "u", "c"),
    "Ds*":   ([(1, 2.1122), (2, 2.5691), (3, 2.860)], "s", "c"),
    "Psi":   ([(0, 3.096900), (1, 3.51067), (2, 3.7737)], "c", "c"),
}

def universal(traj, ap):
    tot = 0.0; out = {}
    for name, (st, q1, q2) in traj.items():
        a, s = b1.best_a(st, ap, b1.UNIV_MASS[q1], b1.UNIV_MASS[q2]); out[name] = a; tot += s
    return tot, out

def fit_universal(traj):
    ap = b1.golden(lambda x: universal(traj, x)[0], 0.6, 1.2, 40)
    tot, a = universal(traj, ap)
    n = sum(len(st) for st, _, _ in traj.values())
    dof = n - len(traj) - 1
    s2 = tot / dof
    h = 0.01
    curv = (universal(traj, ap + h)[0] - 2 * tot + universal(traj, ap - h)[0]) / (h * h)
    return ap, sqrt(2 * s2 / curv), tot / (n - 1), tot, s2, a

def rs(ap): return 1e3 / sqrt(2 * pi * ap)

def main():
    print("R137 -- masses du PDG 2026 lues a la source : controle et mise a jour\n")

    # A
    print("A. Controle des masses citees de memoire")
    diffs = {k: v[1] - v[0] for k, v in CONTROL.items()}
    exact = [k for k, d in diffs.items() if abs(d) < 0.05]
    moved = {k: d for k, d in diffs.items() if abs(d) >= 0.05}
    print(f"    {len(CONTROL)} etats controles ; {len(exact)} identiques a la table ; differents : " + ", ".join(f"{k} {d:+.2f} MeV" for k, d in moved.items()))
    small = max(abs(d) for k, d in diffs.items() if k not in ("Delta(2420) 11/2+",))
    print(f"    hors Delta(2420), ecart maximal {small:.2f} MeV\n")
    check("A. trente et un etats controles contre le PDG 2026 : vingt-huit identiques, ecart maximal 0,21 MeV hors Delta(2420)",
          len(CONTROL) == 31 and len(exact) == 28 and small < 0.25, f"{len(exact)} sur {len(CONTROL)}")

    # B
    print("B. Deux corrections")
    d_delta = diffs["Delta(2420) 11/2+"]
    print(f"    Delta(2420) 11/2+ : masse de Breit-Wigner estimee de 2300 a 2600 MeV, soit environ 2450 ; R120 et R133 prenaient 2420 : {d_delta:+.0f} MeV")
    print(f"    K*(892) charge, produit en collisions de hadrons : 891,88 +- 0,23 MeV ; R133 prenait 891,67 : {diffs['K*(892)+-']:+.2f} MeV")
    d_f2p = diffs["f2'(1525)"]
    print(f"    troisieme ecart, sans effet : f2'(1525), 1517,3 MeV et non 1517,4 : {d_f2p:+.2f} MeV")
    print("    etats lus dans les fiches et non dans la table de synthese : rho5, a6, f6, K5*, b3, pi4, h3, eta4, omega5\n")
    check("B. corrections : Delta(2420) a 2450 MeV (+30), K*(892) a 891,88 MeV (+0,21)", abs(d_delta - 30) < 1e-9 and abs(diffs["K*(892)+-"] - 0.21) < 1e-9, "+30 MeV ; +0,21 MeV")

    # C
    print("C. L'ajustement universel avec les masses de 2026")
    old = fit_universal(b1.TRAJ)
    new = fit_universal(TRAJ_2026)
    shifts = []
    for name in b1.TRAJ:
        for (j, m0), (_, m1) in zip(b1.TRAJ[name][0], TRAJ_2026[name][0]): shifts.append((abs(m1 - m0) * 1e3, name, j, (m1 - m0) * 1e3))
    shifts.sort(reverse=True)
    print("    masses qui ont le plus bouge depuis l'article : " + ", ".join(f"{n} J = {j} {d:+.0f} MeV" for _, n, j, d in shifts[:5]))
    print(f"    masses de 2014 : alpha' = {old[0]:.4f} +- {old[1]:.4f}, chi2 = {old[2]*1e4:.2f}e-4, sqrt(sigma) = {rs(old[0]):.1f} +- {rs(old[0])*old[1]/(2*old[0]):.1f} MeV")
    print(f"    masses de 2026 : alpha' = {new[0]:.4f} +- {new[1]:.4f}, chi2 = {new[2]*1e4:.2f}e-4, sqrt(sigma) = {rs(new[0]):.1f} +- {rs(new[0])*new[1]/(2*new[0]):.1f} MeV")
    print("    ordonnees 2026 : " + " ; ".join(f"{k} {v - b1.UNIV_SHIFT.get(k, 0):+.2f}" for k, v in new[5].items()))
    dshift = (new[0] - old[0]) / sqrt(old[1] ** 2 + new[1] ** 2)
    print(f"    deplacement de la pente : {new[0]-old[0]:+.4f}, soit {dshift:+.1f} ecart type\n")
    check("C. mise a jour 2026 de l'ajustement universel : pente entre 0,87 et 0,90, deplacee de moins d'un ecart type",
          0.87 < new[0] < 0.90 and abs(dshift) < 1 and shifts[0][0] < 35, f"alpha' = {new[0]:.4f} +- {new[1]:.4f}")

    # D
    print("D. La tension de la base et celle du reseau contre cette mise a jour")
    rs_new, drs_new = rs(new[0]), rs(new[0]) * new[1] / (2 * new[0])
    out = {}
    for label, ap in (("base", b1.AP_BASE), ("reseau", b1.AP_LAT)):
        t, _ = universal(TRAJ_2026, ap)
        out[label] = sqrt(max(t - new[3], 0) / new[4])
        print(f"    pente fixee a celle de {label:6s} : {out[label]:.2f} ecarts types d'ajustement")
    n_lat = (445.0 - rs_new) / sqrt(6.7 ** 2 + drs_new ** 2)
    rs_base = rs(b1.AP_BASE)
    print(f"    tension des spectres : {rs_new:.1f} +- {drs_new:.1f} MeV ; base {rs_base:.1f} MeV : {rs_base - rs_new:+.1f} MeV ; reseau 445 +- 6,7 MeV : {445 - rs_new:+.1f} MeV, {n_lat:.1f} ecarts types")
    print()
    check("D. masses de 2026 : tension des spectres entre 420 et 428 MeV ; base entre 1,5 et 3,5 ecarts types ; reseau entre 2,5 et 3,5 avec son erreur",
          420 < rs_new < 428 and 1.5 < out["base"] < 3.5 and 2.5 < n_lat < 3.5, f"{rs_new:.1f} +- {drs_new:.1f} MeV ; base {out['base']:.1f} ; reseau {n_lat:.1f}")

    # E
    print("E. Les baryons avec la masse corrigee")
    def lin(pts):
        xs = [m * m for _, m in pts]; ys = [j for j, _ in pts]
        xb, yb = sum(xs) / len(xs), sum(ys) / len(ys)
        return sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / sum((x - xb) ** 2 for x in xs)
    d_old = lin([(1.5, 1.232), (3.5, 1.930), (5.5, 2.420)])
    d_new = lin([(1.5, 1.232), (3.5, 1.930), (5.5, 2.450)])
    n_slope = lin([(0.5, 0.93827), (2.5, 1.685), (4.5, 2.250)])
    print(f"    Delta(1232), Delta(1950), Delta(2420) : pente {d_old:.3f} avec 2420 MeV (R120), {d_new:.3f} avec 2450 MeV ; N : {n_slope:.3f} (inchangee)")
    print("    -> la pente du Delta descend de 0,92 a 0,89 : elle rejoint celle des mesons ; seule la famille du nucleon reste au-dessus, a 0,955.\n")
    check("E. pente du Delta : 0,922 avec 2420 MeV, 0,892 avec 2450 MeV ; nucleon 0,955", abs(d_old - 0.922) < 2e-3 and abs(d_new - 0.892) < 3e-3 and abs(n_slope - 0.955) < 2e-3, f"{d_old:.3f} -> {d_new:.3f}")

    print("Verdict : les valeurs citees de memoire sont confirmees par le PDG 2026, a deux corrections pres. Avec les masses de 2026 l'ajustement")
    print("universel ne bouge pas de facon significative : la tension des spectres reste a 424 MeV, la base et le reseau de part et d'autre.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
