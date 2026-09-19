#!/usr/bin/env python3
"""R100 -- Le facteur de forme electromagnetique de l'anneau (critere de l'auteur).

Partir uniquement de la charge telle qu'elle existe dans le modele, calculer
    F(q) = (1/e) int rho(r) e^{i q.r} d^3r,   <r^2> = -6 dF/dq^2 |_(q=0),
et ne declarer PASS que si la theorie existante produit naturellement F(q^2) ~ 1
malgre l'anneau de 386 fm, sans mecanisme ajoute pour sauver le resultat.

Les trois lectures de la charge dans la base :
  (1) charge e sur l'anneau fin R = lambda-bar (R53 : c'est e qui circule) ;
  (2) charge = courant de surface de la ligne adaptee, section carree w = 4 lambda-bar/pi^2,
      coeur deplace vers l'interieur (R84) ;
  (3) trois charges e/3 aux jonctions, r_J = D/sqrt3 = 136 fm (R5, R25 ; lecture ancienne).

  A. F(q) et <r^2> pour les trois lectures ; -6 dF/dq^2 = second moment (coherence).
  B. l'orientation : anneau fin, q le long de l'axe : F = 1 exactement (la charge est dans
     le plan perpendiculaire a q) ; q perpendiculaire : J_0(qR).  Mais l'orientation est une
     variable de Bloch de spin 1/2 (R99 C) : aucun moment de rang 2, la densite de charge de
     tout etat est sa moyenne d'orientation, F = <j_0(q rho)> : la sortie « q le long de l'axe »
     est fermee par la base elle-meme.
  C. la reference : l'electron de Dirac ponctuel.  Dans le repere de Breit, u(p')^dag u(p) = 2m
     exactement a tout Q^2 : G_E = 1.  Le terme de Darwin de l'hydrogene, <r^2>_DF = 3/4 lambda-bar^2,
     est le facteur de normalisation m/E = (1 + Q^2/4m^2)^(-1/2), une cinematique, pas une
     distribution de charge : la cible de R89 (334 fm) n'etait pas un rayon de charge.
  D. les sondes : Compton 1 MeV a 90 deg (Klein-Nishina verifie au %), e-p elastique a
     Q^2 = 1 GeV^2 (l'electron y est un point au %), LEP (r_e < 1e-3 fm).
  E. l'echappatoire de l'instantane : au GeV le temps d'interaction (7e-25 s) est bien plus court
     que le tour (8e-21 s) ; la sonde voit une charge ponctuelle quelque part sur l'anneau.
     L'amplitude elastique reste la moyenne <e^{iq.r}> = F(q) ; le reste, 1 - |F|^2, va dans des
     etats internes excites de l'anneau, que l'electron n'a pas (R47) et que la diffusion
     elastique e-p et Bhabha ne montrent pas.  Pas de sortie.
  F. le critere de l'auteur : |F(q) - 1| < 1 % jusqu'a q = 1 GeV/c, sans mecanisme ajoute.
"""
import math
import os
import importlib.util
import numpy as np
from scipy.special import j0

HERE = os.path.dirname(os.path.abspath(__file__))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
hd = load("hydrogen_darwin")

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999      # MeV fm, MeV
LAMBDA_BAR = HBARC / ME
D_J = 6 * LAMBDA_BAR / math.pi ** 2
R_E_BOUND_FM = 1e-3                                             # LEP : < 1e-18 m

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def sph_j0(x):
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x) < 1e-6, 1 - x * x / 6, np.sin(np.where(x == 0, 1, x)) / np.where(x == 0, 1, x))

class Reading:
    """distribution axisymetrique : poids w_i en (r_i, z_i) [fm], somme des poids = 1."""
    def __init__(self, name, r, z, w):
        self.name, self.r, self.z, self.w = name, np.asarray(r, float), np.asarray(z, float), np.asarray(w, float)
        self.rho = np.sqrt(self.r ** 2 + self.z ** 2)
    def r2(self):
        return float(np.sum(self.w * self.rho ** 2))
    def F_sph(self, q_mev):                     # moyenne d'orientation (monopole)
        return float(np.sum(self.w * sph_j0(q_mev * self.rho / HBARC)))
    def F_axis(self, q_mev):                    # q le long de l'axe de l'anneau
        return float(np.sum(self.w * np.cos(q_mev * self.z / HBARC)))
    def F_perp(self, q_mev):                    # q dans le plan de l'anneau
        return float(np.sum(self.w * j0(q_mev * self.r / HBARC) * 1.0))

def readings():
    out = [Reading("(1) anneau fin R = lambda-bar (R53)", [LAMBDA_BAR], [0.0], [1.0])]
    pts, I = hd.r84_distribution()
    out.append(Reading("(2) distribution R84 (courant de surface)", pts[:, 0], pts[:, 1], I))
    r_J = D_J / math.sqrt(3)
    out.append(Reading("(3) trois charges aux jonctions r_J = 136 fm", [r_J] * 3, [0.0] * 3, [1 / 3] * 3))
    return out

def dirac_spinor(p, m, chi):
    E = math.sqrt(m * m + np.dot(p, p))
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
    sp = p[0] * sx + p[1] * sy + p[2] * sz
    return math.sqrt(E + m) * np.concatenate([chi, sp @ chi / (E + m)])

def compton_q(k_mev, theta):
    kp = k_mev / (1 + (k_mev / ME) * (1 - math.cos(theta)))
    return math.sqrt(k_mev ** 2 + kp ** 2 - 2 * k_mev * kp * math.cos(theta))

def main():
    print("R100 -- le facteur de forme electromagnetique de l'anneau\n")
    rd = readings()

    # A. F(q), <r^2>
    print("A. F(q) (moyenne d'orientation) et <r^2> = -6 dF/dq^2")
    qs = [0.1, 0.5, 1.0, 2.0, 10.0, 100.0, 1000.0]
    print("    q [MeV/c]        " + "".join(f"{q:>9.1f}" for q in qs))
    okA = True
    for R in rd:
        Fs = [R.F_sph(q) for q in qs]
        q0 = 1e-3
        r2_from_F = -6 * (R.F_sph(q0) - 1) / (q0 / HBARC) ** 2      # fm^2
        r2_direct = R.r2()
        okA &= abs(r2_from_F / r2_direct - 1) < 1e-4
        print(f"  {R.name:44s}" + "".join(f"{F:>9.4f}" for F in Fs))
        print(f"    <r^2>: -6 dF/dq^2 = {r2_from_F:.1f} fm^2, second moment = {r2_direct:.1f} fm^2 -> r_rms = {math.sqrt(r2_direct):.1f} fm = {math.sqrt(r2_direct)/LAMBDA_BAR:.3f} lambda-bar")
    print()
    check("<r^2> = -6 dF/dq^2 coincide avec le second moment pour les trois lectures (1e-4)", okA, "coherence")

    # B. orientation
    print("B. Orientation : anneau fin, q = 1 MeV/c")
    ring = rd[0]
    q = 1.0
    fa, fp, fs = ring.F_axis(q), ring.F_perp(q), ring.F_sph(q)
    print(f"  q le long de l'axe : F = {fa:.4f} (charge dans le plan perpendiculaire a q : ponctuelle pour cette direction)")
    print(f"  q dans le plan     : F = J_0(qR) = {fp:.4f}")
    print(f"  moyenne (Bloch, spin 1/2 sans rang 2, R99 C) : F = j_0(qR) = {fs:.4f}")
    print("  la seule lecture ou F = 1 exige un axe classique fixe, que R99 a exclu par le 1S : fermee par la base.\n")
    check("l'anneau n'est ponctuel que pour q le long d'un axe classique ; Bloch => F = j_0(qR) < 1",
          abs(fa - 1) < 1e-12 and fs < 0.6, f"axe {fa:.3f}, plan {fp:.3f}, moyenne {fs:.3f}")

    # C. Dirac ponctuel : Breit
    print("C. La reference : Dirac ponctuel dans le repere de Breit")
    okC = True
    for Q2 in (1e-4, 1.0, 1e4, 1e6):                                  # MeV^2
        qv = np.array([0, 0, math.sqrt(Q2)])
        chi = np.array([1, 0], complex)
        u, up = dirac_spinor(-qv / 2, ME, chi), dirac_spinor(qv / 2, ME, chi)
        GE = (np.vdot(up, u) / (2 * ME)).real
        okC &= abs(GE - 1) < 1e-9
        print(f"    Q^2 = {Q2:8.0e} MeV^2 : u(p')^dag u(p) / 2m = {GE:.12f}")
    Q2 = 1e-4
    mE = lambda Q2: 1 / math.sqrt(1 + Q2 / (4 * ME ** 2))
    r2_DF = -6 * (mE(Q2) - 1) / (Q2 / HBARC ** 2) / LAMBDA_BAR ** 2
    print(f"  facteur de normalisation m/E = (1 + Q^2/4m^2)^(-1/2) : -6 d/dQ^2 = {r2_DF:.4f} lambda-bar^2 = le terme de Darwin (3/4)")
    print("  => G_E = 1 exactement ; le 3/4 lambda-bar^2 de l'hydrogene est une cinematique de Dirac, pas une taille.")
    print("  La cible de R89 (r_rms = 334 fm) n'etait pas un rayon de charge : la cible est G_E = 1.\n")
    check("Dirac ponctuel : G_E = 1 a tout Q^2 ; Darwin 3/4 lambda-bar^2 = normalisation m/E, pas une distribution",
          okC and abs(r2_DF - 0.75) < 1e-3, f"{r2_DF:.4f}")

    # D. sondes
    print("D. Les sondes")
    qC = compton_q(1.0, math.pi / 2)
    print(f"  Compton 1 MeV a 90 deg : q = {qC:.3f} MeV/c")
    for R in rd:
        print(f"    {R.name:44s} |F|^2 = {R.F_sph(qC)**2:.3f}  (Klein-Nishina verifie au % : ecart {100*(1-R.F_sph(qC)**2):.0f} %)")
    qP = 1000.0
    print(f"  e-p elastique, Q^2 = 1 GeV^2 (q = {qP:.0f} MeV/c) :")
    for R in rd:
        print(f"    {R.name:44s} |F|^2 = {R.F_sph(qP)**2:.1e}  (l'electron y est un point au %)")
    print(f"  LEP : r_e < {R_E_BOUND_FM:.0e} fm ; l'anneau : r_rms = {math.sqrt(rd[0].r2()):.0f} fm, rapport {math.sqrt(rd[0].r2())/R_E_BOUND_FM:.0e} en rayon\n")
    check("Compton 1 MeV : |F|^2 = 0,18 (anneau) contre 1 au % ; e-p 1 GeV^2 : |F|^2 ~ 1e-7 ; LEP : 4e5 en rayon",
          rd[0].F_sph(qC) ** 2 < 0.5 and rd[0].F_sph(qP) ** 2 < 1e-5, f"{rd[0].F_sph(qC)**2:.2f}, {rd[0].F_sph(qP)**2:.0e}")

    # E. instantane
    T_turn = 2 * math.pi * LAMBDA_BAR / 2.99792458e23          # s (fm / (fm/s))
    t_int = HBARC / qP / 2.99792458e23
    print("E. L'echappatoire de l'instantane")
    print(f"  tour de l'anneau T = {T_turn:.1e} s ; temps d'interaction a 1 GeV/c : {t_int:.1e} s (rapport {T_turn/t_int:.0e})")
    print("  la sonde voit une charge ponctuelle a une position aleatoire sur l'anneau ; l'amplitude elastique est")
    print(f"  la moyenne <e^(iq.r)> = F(q) = j_0(qR) ; le complement 1 - |F|^2 = {1-rd[0].F_sph(qP)**2:.7f} irait dans des etats")
    print("  internes excites de l'anneau : l'electron n'en a pas (R47), et e-p au GeV comme Bhabha au LEP sont")
    print("  elastiques et ponctuels. Pas de sortie sans mecanisme nouveau.\n")

    # F. critere de l'auteur
    print("F. Le critere de l'auteur : |F(q) - 1| < 1 % jusqu'a q = 1 GeV/c, sans mecanisme ajoute")
    worst = {}
    for R in rd:
        qgrid = np.geomspace(0.01, 1000.0, 400)
        dev = max(abs(R.F_sph(q) - 1) for q in qgrid)
        q1 = next(q for q in qgrid if abs(R.F_sph(q) - 1) > 0.01)
        worst[R.name] = (dev, q1)
        print(f"    {R.name:44s} ecart max {dev:.3f} ; F quitte 1 % des q = {q1:.2f} MeV/c")
    print()
    passF = all(dev < 0.01 for dev, _ in worst.values())
    check("critere R100 : F(q^2) ~ 1 (1 %) jusqu'a 1 GeV/c pour la charge du modele, sans mecanisme ajoute",
          passF, "aucune lecture ; F quitte 1 % des q ~ 0,1 MeV/c")

    print("Verdict : la charge du modele suit l'anneau, dans ses trois lectures, et son facteur de forme est")
    print("j_0(qR) : <r^2> = 0,12 a 1,0 lambda-bar^2, F = 0,42 a 1 MeV/c, 5e-4 a 1 GeV/c. L'electron mesure a G_E = 1")
    print("(Compton au %, e-p au %, LEP < 1e-3 fm) et son terme de Darwin est la normalisation m/E d'un point,")
    print("pas une taille. Rien dans la base ne rend la charge ponctuelle : l'orientation de Bloch ferme la")
    print("direction axiale, l'instantane exige des etats excites absents. Le critere de l'auteur ECHOUE :")
    print("cette representation de l'electron (charge e portee par l'anneau de 386 fm) s'arrete ici.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
