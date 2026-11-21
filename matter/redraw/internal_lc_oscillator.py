#!/usr/bin/env python3
"""R121 -- « Le cout d'un deplacement est energie pour masse » : les 3 + 3 dimensions de stockage comme oscillateurs LC.

Enonce de l'auteur (21 novembre 2026), apres la proposition « CD en 10 dimensions : 3 massiques, 3 stockage capacitif,
3 stockage inductif, 1 temps » : le cout d'un deplacement dans une direction est une energie, comme pour une masse.
Pris au mot : chaque direction porte un terme ½ (inertie) x (vitesse)^2. Espace : inertie m. Direction « inductive » :
coordonnee q (charge), inertie L, cout ½ L q_dot^2 = ½ L I^2. Direction « capacitive » : cout q^2 / 2C. Mais ½ L I^2 et
q^2/2C sont l'energie cinetique et l'energie potentielle de la MEME coordonnee (I = q_dot) : les 3 + 3 dimensions sont
3 oscillateurs LC internes attaches au point, de frequence 1/sqrt(LC), pas 6 dimensions.

  A. Identite de Legendre : L = ½ L q_dot^2 - q^2/2C, p = L q_dot, H = p^2/2L + q^2/2C ; la lecture duale (flux Phi)
     echange les roles (½ C Phi_dot^2 + Phi^2/2L). Les nombres du dossier (CONSTRAINTS 1b) : ½ L_e I_e^2 = 2,19 keV avec
     I_e = e c / (2 pi lambda_bar_e) = 19,8 A, hbar / sqrt(L_e C_e) = 101 keV ; d'ou L_e = 1,79e-18 H, C_e = 2,37e-23 F,
     Z_e = sqrt(L_e/C_e) = 275 ohm = 0,73 Z_0 : les trois nombres du manuscrit sont coherents entre eux.
  B. La masse comme energie interne : trois oscillateurs au point zero, (3/2) hbar omega = 151 keV = 0,30 m_e ; un quantum
     par axe, 303 keV ; m_e exigerait n = 3,56 quanta au-dessus du point zero (un oscillateur : 4,56) : aucun niveau.
     Sans point zero, m_e / hbar omega = 5,06 : « cinq quanta a 1 % », une coincidence entre nombres calibres sur
     l'electron lui-meme (R_3 = lambda_bar_e, R_3/r = 37,1), non revendiquee.
  C. Les niveaux excites : m_e + n hbar omega = 612, 713, 814 keV... Un electron excite est exclu jusqu'a > 100 GeV (LEP,
     leptons excites) et l'electron n'a aucune structure interne a 10^-3 fm (G_E = 1, R100) ni aucun degre de liberte
     interne visible dans g - 2 (10^-12) : EXCLU. L'oscillateur doit etre gele dans son fondamental.
  D. Le fondamental gele : si toute la masse est le point zero d'un oscillateur, omega = 2 m_e c^2 / hbar = 1,55e21 rad/s,
     la frequence du Zitterbewegung de Dirac (R103-R104, identite) ; il faudrait L = Z_e/omega = 1,77e-19 H et
     C = 2,34e-24 F, dix fois sous les L_e, C_e du modele : le modele ne le fournit pas sans reajuster L et C.
  Verdict : REDUCTION (6 dimensions de stockage -> 3 oscillateurs), IDENTITE (fondamental = Zitterbewegung, a L et C
  reajustes), EXCLU (niveaux). Le « 10 = 3 + 3 + 3 + 1 » se reduit au point materiel avec oscillateur interne gele :
  un seul nombre (la masse) regle par L et C, rien d'autre d'observable.
"""
import math

HBAR, C, E_CH, EPS0 = 1.054571817e-34, 2.99792458e8, 1.602176634e-19, 8.8541878128e-12
MU0 = 1.25663706212e-6
ME_MEV = 0.51099895
KEV_J = 1.602176634e-16
LAM_E = 3.8615926796e-13          # m
Z0 = math.sqrt(MU0 / EPS0)
E_IND_KEV, E_LC_KEV = 2.19, 101.0  # CONSTRAINTS 1b

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R121 -- les dimensions de stockage comme oscillateurs LC internes\n")

    # A. Legendre + coherence des nombres du dossier
    print("A. Une coordonnee, deux energies ; les nombres du dossier")
    L_test, C_test, q, qd = 2.0, 0.5, 0.3, 0.7
    lag = 0.5 * L_test * qd ** 2 - q ** 2 / (2 * C_test)
    p = L_test * qd
    H = p * qd - lag
    H_expected = p ** 2 / (2 * L_test) + q ** 2 / (2 * C_test)
    print(f"    L = ½Lq'^2 - q^2/2C -> H = p^2/2L + q^2/2C : {H:.6f} = {H_expected:.6f} ; inductif = cinetique, capacitif = potentiel, meme q")
    I_e = E_CH * C / (2 * math.pi * LAM_E)
    L_e = 2 * E_IND_KEV * KEV_J / I_e ** 2
    omega = E_LC_KEV * KEV_J / HBAR
    C_e = 1 / (L_e * omega ** 2)
    Z_e = math.sqrt(L_e / C_e)
    print(f"    I_e = ec/(2 pi lambda_bar_e) = {I_e:.2f} A ; L_e = {L_e:.3e} H ; omega = {omega:.3e} rad/s ; C_e = {C_e:.3e} F ; Z_e = {Z_e:.1f} ohm = {Z_e/Z0:.3f} Z_0 (manuscrit : 0,73 Z_0 = 275 ohm)\n")
    check("A. Legendre : capacitif et inductif sont les deux energies d'une seule coordonnee ; L_e, C_e, Z_e du dossier coherents (275 ohm)",
          abs(H - H_expected) < 1e-12 and abs(Z_e / 275.0 - 1) < 0.01, f"Z_e = {Z_e:.1f} ohm, I_e = {I_e:.2f} A")

    # B. la masse comme energie interne
    print("B. La masse comme energie des oscillateurs internes")
    hw = E_LC_KEV / 1000  # MeV
    zp3 = 1.5 * hw
    n3 = (ME_MEV - zp3) / hw
    n1 = (ME_MEV - 0.5 * hw) / hw
    n0 = ME_MEV / hw
    print(f"    hbar omega = {hw*1000:.0f} keV ; point zero de 3 oscillateurs {zp3*1000:.0f} keV = {zp3/ME_MEV:.2f} m_e ; un quantum par axe {3*hw*1000:.0f} keV")
    print(f"    m_e = (n + 3/2) hbar omega -> n = {n3:.2f} ; (n + 1/2) -> n = {n1:.2f} ; sans point zero n = {n0:.2f} (coincidence « 5 quanta », non revendiquee)\n")
    check("B. aucun niveau d'oscillateur a 101 keV ne tombe sur m_e avec le point zero (n = 3,56 ou 4,56) ; 5,06 sans point zero",
          abs(n3 - round(n3)) > 0.3 and abs(n1 - round(n1)) > 0.3 and abs(n0 - 5) < 0.1, f"n = {n3:.2f}, {n1:.2f}, {n0:.2f}")

    # C. niveaux excites
    print("C. Les niveaux excites de l'electron")
    levels = [ME_MEV + n * hw for n in (1, 2, 3)]
    print(f"    m_e + n hbar omega = {', '.join(f'{1000*l:.0f}' for l in levels)} keV ; bornes : e* > 100 GeV (LEP), G_E = 1 a 10^-3 fm, g - 2 a 10^-12 -> EXCLU\n")
    check("C. electron excite a 612 keV exclu (bornes LEP > 100 GeV ; aucune structure interne)", levels[0] < 1e5, f"{1000*levels[0]:.0f} keV")

    # D. fondamental gele = Zitterbewegung
    print("D. Toute la masse au point zero d'un oscillateur : omega = 2 m_e c^2 / hbar")
    omega_zb = 2 * ME_MEV * 1e3 * KEV_J / HBAR
    L_zb = Z_e / omega_zb; C_zb = 1 / (Z_e * omega_zb)
    print(f"    omega = {omega_zb:.3e} rad/s (Zitterbewegung) ; a Z_e fixe : L = {L_zb:.3e} H, C = {C_zb:.3e} F ; L_e/L = {L_e/L_zb:.1f}, C_e/C = {C_e/C_zb:.1f}\n")
    check("D. fondamental = Zitterbewegung exige L et C dix fois sous les valeurs du modele : identite a L, C reajustes",
          9 < L_e / L_zb < 11 and 9 < C_e / C_zb < 11, f"L_e/L = {L_e/L_zb:.2f}")

    print("Verdict : REDUCTION (6 dimensions de stockage -> 3 oscillateurs LC), IDENTITE (fondamental = Zitterbewegung, L et C reajustes),")
    print("EXCLU (niveaux excites). Le « 10 = 3 + 3 + 3 + 1 » est un point materiel avec un oscillateur interne gele : un nombre, la masse.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
