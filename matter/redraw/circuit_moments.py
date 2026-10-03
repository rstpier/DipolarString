#!/usr/bin/env python3
"""R115 -- Le moment du nucleon porte par les circuits de quark eux-memes, a R_c.

R114 a ferme la lecture « charge circulante sur un anneau separe » (mu_p et r_p incompatibles a
tout rayon). Reste la seule autre structure chargee du nucleon : les trois circuits de quark,
cercles de rayon R_c = 3 l_1(9)/2 pi = 0,388 fm portant chacun un quantum a c (R113 A : S = hbar/2
chacun), avec la charge du quark sur ses brins (R43). Le moment d'un circuit est alors
mu_q = q_q (R_c / lambda-bar_p) w_q mu_N, w_q = poids (sens et fraction) de sa circulation.

  A. sens classiques, deux contre un (paire alignee, impair contre, |w| = 1) : mu_p = 3,08 mu_N
     (+10 %), mu_n = -2,46 (+29 %), rapport -5/4 (-14 %) : la lecture exclue par R72, retrouvee.
  B. poids de SU(6), (2 w_paire, w_impair) = (4/3, -1/3) : rapport -3/2 exact (il ne depend que
     de w_impair/w_paire = -1/2), grandeur mu_p = 1,85 mu_N (-34 %) ; la grandeur exige
     R = 0,587 fm = 1,51 R_c, entre g = 1 (R_c) et g = 2 (2 R_c).
  C. la famille complete : mu_q proportionnel a q_q, somme des poids 2 w_paire + w_impair = 1
     (spin 1/2). Un seul parametre : mu_p exact donne le rapport a -1,28 (-12 %), le rapport exact
     donne mu_p a -34 % ; le minimax est a 9,5 %. Aucun partage a spin 1/2 ne donne les deux.
  D. coincidence enregistree et rejetee : (w_paire, w_impair) = (1, -1/2) donne
     mu_p = (3/2)(R_c/lambda-bar_p) = 2,77 mu_N (-0,9 %), mu_n = -1,85 (-3,6 %), rapport -3/2, mais la
     somme des poids vaut 3/2 : S_z = 3 hbar/4, pas un spin 1/2.
"""
import math

HBARC = 197.3269804
ME, MP = 0.51099895, 938.27209
MU_P, MU_N = 2.79284734, -1.91304273
LAM = HBARC / ME
L1 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)
R_C = 3 * L1 / (2 * math.pi)
E_Q = math.pi * HBARC / (3 * L1)
LAM_P = HBARC / MP
U = R_C / LAM_P                     # moment d'une charge e a c sur un cercle de rayon R_c, en mu_N
Q_U, Q_D = 2 / 3, -1 / 3
R_MEAS = MU_P / MU_N

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def moments(w_pair, w_odd, unit=U):
    """proton = (u, u | d), neutron = (d, d | u) : la paire porte w_pair chacune, l'impair w_odd"""
    mu_p = (2 * w_pair * Q_U + w_odd * Q_D) * unit
    mu_n = (2 * w_pair * Q_D + w_odd * Q_U) * unit
    return mu_p, mu_n

def dev(x, ref):
    return x / ref - 1

def main():
    print("R115 -- le moment du nucleon porte par les circuits de quark a R_c\n")
    print(f"  R_c = {R_C:.4f} fm, E_q = {E_Q:.1f} MeV, lambda-bar_p = {LAM_P:.4f} fm ; unite : e a c sur R_c = {U:.4f} mu_N")
    print(f"  (R_c = hbar c / 2 E_q : un quantum sur un cercle a toute son energie en circulation, g = 1 ; g = 2 serait 2 R_c)")
    print(f"  mesure : mu_p = {MU_P:.3f}, mu_n = {MU_N:.3f} mu_N, rapport {R_MEAS:.3f}\n")

    # A. sens classiques
    print("A. Sens classiques : la paire alignee (+1, +1), l'impair contre (-1)")
    mp, mn = moments(1.0, -1.0)
    print(f"    mu_p = {mp:.3f} mu_N ({100*dev(mp, MU_P):+.1f} %), mu_n = {mn:.3f} ({100*dev(mn, MU_N):+.1f} %), rapport {mp/mn:.3f} ({100*dev(mp/mn, R_MEAS):+.1f} %)")
    print("    = « impair a contre-sens a egalite » de R72 (-1,25), exclue par le rapport ; la grandeur est a 10 %.\n")
    check("A. sens classiques a R_c : mu_p = 3,08 (+10 %), rapport -5/4 (-14 %), la lecture exclue par R72",
          abs(mp / mn + 1.25) < 1e-9 and abs(dev(mp, MU_P) - 0.10) < 0.01, f"{mp:.2f} mu_N, {mp/mn:.3f}")

    # B. SU(6)
    print("B. Poids de SU(6) : (2 w_paire, w_impair) = (4/3, -1/3), somme 1 (spin 1/2)")
    mp6, mn6 = moments(2 / 3, -1 / 3)
    R_need = MU_P * LAM_P / (mp6 / U)       # rayon qui donnerait la grandeur avec ces poids
    print(f"    mu_p = {mp6:.3f} mu_N ({100*dev(mp6, MU_P):+.1f} %), mu_n = {mn6:.3f} ({100*dev(mn6, MU_N):+.1f} %), rapport {mp6/mn6:.3f} ({100*dev(mp6/mn6, R_MEAS):+.1f} %)")
    print(f"    le rapport ne depend que de w_impair/w_paire = -1/2 ; la grandeur exigerait R = {R_need:.4f} fm = {R_need/R_C:.3f} R_c")
    print(f"    (R25 : 0,587 fm ; g = 2 donnerait 2 R_c et mu_p = {2*mp6:.2f}, +{100*dev(2*mp6, MU_P):.0f} %).\n")
    check("B. SU(6) a R_c : rapport -3/2 exact, grandeur -34 % ; la grandeur exige R = 0,587 fm = 1,51 R_c",
          abs(mp6 / mn6 + 1.5) < 1e-9 and abs(dev(mp6, MU_P) + 0.34) < 0.01 and abs(R_need - 0.5874) < 0.001, f"{mp6:.3f} mu_N, R = {R_need:.3f} fm")

    # C. la famille a spin 1/2
    print("C. Toute la famille mu_q = q_q w_q (R_c/lambda-bar_p) mu_N avec 2 w_paire + w_impair = 1 (spin 1/2) : un parametre w_paire")
    best = None
    rows = []
    for i in range(0, 20001):
        w = 0.3 + 0.9 * i / 20000
        wo = 1 - 2 * w
        p, n = moments(w, wo)
        d_p, d_r = dev(p, MU_P), dev(p / n, R_MEAS)
        m = max(abs(d_p), abs(d_r))
        if best is None or m < best[0]:
            best = (m, w, wo, p, n, d_p, d_r)
    w_mu = (MU_P / U + 1 / 3) / 2          # mu_p exact : (2 w - 1/3) U = mu_p
    p1, n1 = moments(w_mu, 1 - 2 * w_mu)
    print(f"    mu_p exact      : w_paire = {w_mu:.3f}, w_impair = {1-2*w_mu:+.3f} -> rapport {p1/n1:.3f} ({100*dev(p1/n1, R_MEAS):+.1f} %)")
    print(f"    rapport -3/2    : w_paire = 2/3, w_impair = -1/3 -> mu_p {100*dev(mp6, MU_P):+.1f} %")
    print(f"    minimax         : w_paire = {best[1]:.3f}, w_impair = {best[2]:+.3f} -> mu_p {100*best[5]:+.1f} %, rapport {100*best[6]:+.1f} %")
    print("  -> a R_c, aucun partage a spin 1/2 ne donne a la fois la grandeur et le rapport ; le mieux est 9,5 % sur les deux.\n")
    check("C. famille complete a spin 1/2 : mu_p exact <-> rapport -12 % ; rapport exact <-> mu_p -34 % ; minimax > 5 %",
          abs(dev(p1 / n1, R_MEAS) + 0.12) < 0.01 and best[0] > 0.05, f"minimax {100*best[0]:.1f} % a w_paire = {best[1]:.2f}")

    # D. la coincidence (1, -1/2)
    print("D. Coincidence : (w_paire, w_impair) = (1, -1/2), la paire entiere, l'impair a moitie contre")
    pc, nc = moments(1.0, -0.5)
    print(f"    mu_p = (3/2)(R_c/lambda-bar_p) = {pc:.3f} mu_N ({100*dev(pc, MU_P):+.1f} %), mu_n = {nc:.3f} ({100*dev(nc, MU_N):+.1f} %), rapport {pc/nc:.3f}")
    print(f"    mais somme des poids 2 w_paire + w_impair = {2*1.0-0.5:.1f} : S_z = {(2*1.0-0.5)/2:.2f} hbar, pas 1/2. Rejetee par le compte ; enregistree.\n")
    check("D. (1, -1/2) : mu_p -0,9 %, mu_n -3,6 %, rapport -3/2, mais S_z = 3/4 : coincidence rejetee par le compte",
          abs(dev(pc, MU_P)) < 0.01 and abs(dev(nc, MU_N)) < 0.04 and abs(pc / nc + 1.5) < 1e-9, f"{pc:.3f}, {nc:.3f} mu_N")

    print("Verdict : NON. A R_c, les circuits de quark donnent le rapport -3/2 avec les poids de SU(6) (regle posee, R72) mais")
    print("une grandeur 34 % trop faible, ou la grandeur avec des poids dont la somme n'est pas un spin 1/2 ; les sens classiques")
    print("donnent -5/4. Ce qui porte le moment du nucleon reste ouvert : ni un anneau charge (R114), ni les circuits a R_c.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
