#!/usr/bin/env python3
"""R120 -- La corde de Nambu-Goto avec la tension de la base : pente de Regge, tours N/Delta, identites.

Enonce (21 novembre 2026) : importer la « vieille » corde hadronique (Nambu-Goto : energie = tension x longueur,
extremites a c) avec la tension de la base, sigma = pi hbar c / l_1(9)^2 (sqrt(sigma) = 430,3 MeV, sigma = 0,185 GeV^2),
et confronter ce qu'elle impose aux tables : la pente de Regge alpha' = 1/(2 pi sigma), les tours J = alpha' M^2 + alpha_0,
la boucle fermee avec un quantum de circulation (l'extremum que R113 n'avait pas), le terme de Luscher.

  A. Identite. La corde ouverte tournante (extremites sans masse) a M = pi sigma R, J = pi sigma R^2 / 2, donc
     J = alpha' M^2 avec alpha' = 1/(2 pi sigma) ; a J = 1/2 : M = sqrt(pi sigma) = pi hbar c / l_1(9) = 3 E_q = 762,7 MeV,
     la « part statique » de la base (trois circuits de quark). Identite de definitions (sigma = pi hbar c / l_1^2) : la
     regle 1 (un quantum pi hbar c par circuit) est la corde de Nambu-Goto a J = 1/2 ; l'energie d'anneau de la lecture
     « masse + spin » (R72 : 938 - 763 = 176 MeV) est l'ordonnee de Regge alpha_0(N) = 1/2 - alpha' M_N^2 = -0,26.
  B. Pente. alpha'_base = 0,859 GeV^-2 ; ajustements lineaires J(M^2) sur les tables : rho-a2-rho3 0,888 (-3,3 %) ; avec
     a4(1970) 0,915 (-6 %) ; N(939)-N(1685)-N(2250) 0,955 (-10 %) ; Delta(1232)-Delta(1930)-Delta(2420) 0,922 (-7 %).
     La fenetre de l'exposant (6,2758-6,2926) deplace alpha' entre 0,842 et 0,874 : toutes les pentes mesurees sont
     au-dessus. Meme sens que R117/R119 (sqrt(sigma) sous le reseau) : la tension de la base est basse, par le compte 9.
  C. Tours. Avec alpha'_base et l'etat le plus bas de chaque trajectoire : a2 1328 (+0,8 %), rho3 1711 (+1,3 %),
     a4 2023 (+2,8 %) ; N 5/2 1791 (+6 %), N 9/2 2353 (+5 %) ; Delta 7/2 1961 (+1,6 %), Delta 11/2 2484 (+2,6 %).
  D. Boucle fermee avec un quantum de circulation. E(R) = 2 pi sigma R + kappa hbar c / R a un minimum (R113 : aucun sans
     tension) : R* = sqrt(kappa hbar c / 2 pi sigma), E* = 2 sqrt(2 pi sigma kappa hbar c). Demi-quantum de la base
     (kappa = 1/2) : R* = 0,129 fm, E* = 2 sqrt(pi sigma) = 1525 MeV = deux fois la part statique (identite) ; quantum
     entier : 0,183 fm, 2157 MeV. Etat neutre vers 1,5 GeV, non revendique (glueball 0++ du reseau : ~1,7 GeV).
  E. Luscher. Premier terme quantique de la corde en D = 4 : -pi hbar c / (12 r) = -63,6 MeV a r = l_1(9), -8,3 % de 763 :
     la part statique quantique est 699 MeV, hors des 770-813 MeV que mu_p exige (HADRONS, ligne « part statique »).
  Verdict : CONDITIONNEL (sigma de la base) et EN TENSION pour la pente ; IDENTITE pour la part statique et l'anneau ;
  la boucle fermee a desormais un rayon (acquis pour la reconstruction : E = tension x longueur + quantum / R est le
  fonctionnel minimal de l'etape 1) ; le terme de Luscher va dans le mauvais sens pour le nucleon. Pas de decouverte :
  la base avait deja la tension de Regge-Cornell ; R120 montre que sa regle posee est cette corde a J = 1/2.
"""
import math

HBARC = 197.3269804          # MeV fm
ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
MN, MDELTA = 938.27208816, 1232.0
LAM = HBARC / ME
L1_3 = 2 * math.pi * LAM / 3
def l1(p): return L1_3 * 3 ** (-p)
L1 = l1(2 * math.pi)
SIGMA = math.pi * HBARC / L1 ** 2          # MeV/fm
SIGMA_GEV2 = SIGMA * HBARC * 1e-6          # GeV^2
E_Q = math.pi * HBARC / (3 * L1)
def alpha_prime(sig_gev2): return 1 / (2 * math.pi * sig_gev2)   # GeV^-2

# Tables (MeV), masses Breit-Wigner PDG 2024
RHO = [(1, 775.26), (2, 1318.2), (3, 1688.8), (4, 1967.0)]          # rho, a2(1320), rho3(1690), a4(1970)
N_TRAJ = [(0.5, MN), (2.5, 1685.0), (4.5, 2250.0)]                  # N, N(1680) 5/2+, N(2220) 9/2+
D_TRAJ = [(1.5, MDELTA), (3.5, 1930.0), (5.5, 2420.0)]              # Delta, Delta(1950) 7/2+, Delta(2420) 11/2+

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def fit(pts):
    xs = [(m / 1000) ** 2 for _, m in pts]; ys = [j for j, _ in pts]
    xb, yb = sum(xs) / len(xs), sum(ys) / len(ys)
    a = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / sum((x - xb) ** 2 for x in xs)
    return a, yb - a * xb

def main():
    ap = alpha_prime(SIGMA_GEV2)
    print("R120 -- la corde de Nambu-Goto avec la tension de la base\n")
    print(f"  base : l_1(9) = {L1:.4f} fm ; sigma = pi hbar c / l_1^2 = {SIGMA:.1f} MeV/fm = {SIGMA_GEV2:.4f} GeV^2 ; sqrt(sigma) = {1000*math.sqrt(SIGMA_GEV2):.1f} MeV ; alpha' = 1/(2 pi sigma) = {ap:.4f} GeV^-2\n")

    # A. identite : part statique = corde a J = 1/2 ; anneau = ordonnee
    print("A. La corde ouverte tournante : M = pi sigma R, J = pi sigma R^2/2, J = alpha' M^2")
    m_half = math.sqrt(math.pi * SIGMA_GEV2) * 1000
    static = math.pi * HBARC / L1
    a0_N = 0.5 - ap * (MN / 1000) ** 2
    ring = MN - m_half
    print(f"    J = 1/2 : M = sqrt(pi sigma) = {m_half:.2f} MeV ; pi hbar c / l_1(9) = {static:.2f} MeV ; 3 E_q = {3*E_Q:.2f} MeV  -> identite")
    print(f"    nucleon : alpha_0(N) = 1/2 - alpha' M_N^2 = {a0_N:+.3f} ; en masse : M_N - sqrt(pi sigma) = {ring:.1f} MeV = l'energie d'anneau de R72 (176 MeV)")
    print(f"    Delta : alpha_0 = 3/2 - alpha' M_Delta^2 = {1.5 - ap*(MDELTA/1000)**2:+.3f} ; M(J = 3/2) classique = sqrt(3 pi sigma) = {1000*math.sqrt(3*math.pi*SIGMA_GEV2):.0f} MeV ({100*(math.sqrt(3*math.pi*SIGMA_GEV2)*1000/MDELTA-1):+.1f} %)\n")
    check("A. identite : la part statique de la base est la corde de Nambu-Goto a J = 1/2 ; l'anneau (176 MeV) est l'ordonnee alpha_0(N) = -0,26",
          abs(m_half / static - 1) < 1e-12 and abs(3 * E_Q - static) < 1e-9 and abs(ring - 176) < 1.0 and -0.27 < a0_N < -0.25,
          f"sqrt(pi sigma) = {m_half:.2f} = pi hbar c/l_1 = {static:.2f} MeV ; anneau {ring:.1f} MeV ; alpha_0(N) = {a0_N:+.3f}")

    # B. la pente contre les tables
    print("B. La pente de Regge contre les trajectoires mesurees (ajustement lineaire J = alpha' M^2 + alpha_0)")
    fits = {
        "rho - a2 - rho3": fit(RHO[:3]),
        "rho - a2 - rho3 - a4(1970)": fit(RHO),
        "N(939) - N(1685) 5/2 - N(2250) 9/2": fit(N_TRAJ),
        "Delta(1232) - Delta(1930) 7/2 - Delta(2420) 11/2": fit(D_TRAJ),
    }
    for k, (a, a0) in fits.items():
        print(f"      {k:50s} alpha' = {a:.3f} GeV^-2, alpha_0 = {a0:+.2f} ; base {ap:.3f} : {100*(ap/a-1):+.1f} %")
    win = (alpha_prime(SIGMA_GEV2 * 3 ** (2 * (6.2926 - 2 * math.pi))), alpha_prime(SIGMA_GEV2 * 3 ** (2 * (6.2758 - 2 * math.pi))))
    print(f"    fenetre de l'exposant [6,2758 ; 6,2926] : alpha' entre {win[0]:.3f} et {win[1]:.3f} GeV^-2 ; toutes les pentes ajustees sont au-dessus")
    print("    -> la base est basse de 3 % (mesons) a 7-10 % (baryons), meme sens que sqrt(sigma) sous le reseau (R117, R119).\n")
    slopes = [a for a, _ in fits.values()]
    check("B. alpha'_base = 0,859 sous toutes les pentes ajustees : -3 % (rho-a2-rho3), -6 a -10 % (a4, N, Delta) ; la fenetre de l'exposant ne couvre aucune",
          all(ap < a for a in slopes) and abs(ap / fits["rho - a2 - rho3"][0] - 1) < 0.05 and max(abs(ap / a - 1) for a in slopes) < 0.12 and all(win[1] < a for a in slopes),
          f"ecarts {', '.join(f'{100*(ap/a-1):+.1f}' for a in slopes)} % ; fenetre [{win[0]:.3f} ; {win[1]:.3f}]")

    # C. les tours
    print("C. Les tours : alpha'_base et l'etat le plus bas de chaque trajectoire")
    devs = []
    for name, traj in (("rho", RHO), ("N", N_TRAJ), ("Delta", D_TRAJ)):
        j0, m0 = traj[0]
        for j, m in traj[1:]:
            pred = 1000 * math.sqrt((m0 / 1000) ** 2 + (j - j0) / ap)
            devs.append(pred / m - 1)
            print(f"      {name:6s} J = {j:4.1f} : predit {pred:.0f} MeV, table {m:.0f} : {100*(pred/m-1):+.1f} %")
    print("    -> 1-3 % pour les mesons et le Delta, 5-6 % pour le N : des predictions a l'etat le plus bas et a la tension de la base pres.\n")
    check("C. tours de Regge avec alpha'_base : tous les etats predits a moins de 7 % des tables", max(abs(d) for d in devs) < 0.07, f"max {100*max(abs(d) for d in devs):.1f} %")

    # D. la boucle fermee avec un quantum
    print("D. La boucle fermee : E(R) = 2 pi sigma R + kappa hbar c / R (R113 : aucun extremum sans la tension)")
    out = {}
    for kappa, lab in ((0.5, "demi-quantum de la base (pi hbar c par tour)"), (1.0, "quantum entier (2 pi hbar c par tour)")):
        R = math.sqrt(kappa * HBARC / (2 * math.pi * SIGMA))
        E = 2 * math.pi * SIGMA * R + kappa * HBARC / R
        d2 = 2 * kappa * HBARC / R ** 3
        out[kappa] = (R, E, d2)
        print(f"      {lab:46s} R* = {R:.4f} fm, E* = {E:.1f} MeV, d2E/dR2 = {d2:.0f} MeV/fm^2 > 0")
    print(f"    E*(1/2) = 2 sqrt(pi sigma) = {2*m_half:.1f} MeV = deux fois la part statique (identite) ; etat neutre vers 1,5 GeV, non revendique (glueball 0++ reseau ~1,7 GeV).\n")
    check("D. la tension fixe le rayon de la boucle fermee (minimum) : R* = 0,129 fm, E* = 1525 MeV = 2 sqrt(pi sigma) pour le demi-quantum",
          out[0.5][2] > 0 and 0.12 < out[0.5][0] < 0.14 and abs(out[0.5][1] / (2 * m_half) - 1) < 1e-12,
          f"R* = {out[0.5][0]:.4f} fm, E* = {out[0.5][1]:.1f} MeV")

    # E. Luscher
    print("E. Le terme de Luscher en D = 4 : V(r) = sigma r - pi hbar c / (12 r)")
    lus = -math.pi * HBARC / (12 * L1)
    stat_q = static + lus
    print(f"    a r = l_1(9) : sigma r = {static:.1f} MeV, Luscher = {lus:.1f} MeV ({100*lus/static:+.1f} %), part statique quantique {stat_q:.1f} MeV contre 770-813 MeV exiges par mu_p (HADRONS)")
    print("    -> le premier terme quantique de la corde va dans le mauvais sens pour le nucleon ; la base ne le contient pas.\n")
    check("E. Luscher a l_1(9) : -64 MeV (-8 %), la part statique quantique 699 MeV sort de 770-813", -70 < lus < -58 and stat_q < 770, f"{lus:.1f} MeV, {stat_q:.1f} MeV")

    print("Verdict : CONDITIONNEL et EN TENSION pour la pente (la tension de la base, par le compte 9, est basse de 3-10 % en alpha',")
    print("meme sens que sqrt(sigma) sous le reseau) ; IDENTITE pour la part statique (corde a J = 1/2) et l'anneau (ordonnee alpha_0) ;")
    print("ACQUIS pour la reconstruction : avec une tension, la boucle fermee a un rayon, E = tension x longueur + quantum / R est le")
    print("fonctionnel minimal de l'etape 1 ; le terme de Luscher (-64 MeV) va dans le mauvais sens. Pas de decouverte : la regle posee")
    print("« un quantum par circuit » est la corde de Nambu-Goto a J = 1/2, avec la tension de Regge-Cornell que la base avait deja.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
