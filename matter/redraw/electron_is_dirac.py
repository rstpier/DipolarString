#!/usr/bin/env python3
"""R103 -- Decision : l'electron electromagnetique est celui de Dirac.

R100-R102 : charge et moment de l'electron sont ponctuels (G_E = G_M = 1 a < 1e-3 fm), la base n'a
qu'une source de moment (le courant de charge) et elle est exclue.  Deux issues : (a) inventer un
moment sans courant de charge (un mecanisme nouveau, interdit par le critere de R100) ; (b) donner
l'electron electromagnetique a Dirac et ne garder de la base que ce qui ne depend pas de la
position de la charge.  Decision : (b).  Ce script mesure ce que (b) coute et ce qu'il laisse.

  A. les nombres de l'electron dans la base sont des identites de Compton, choisies pour Dirac :
     S = R E_circ/c = hbar/2 et mu = q c R/2 = mu_B avec R = lambda-bar, E_circ = m/2, q = e (R53).
     Le partage m/2 est ce qui fait g = 2 : avec toute la masse en circulation, g = 1 ; avec le
     rayon du Zitterbewegung (lambda-bar/2) et toute la masse, S = hbar/2 mais mu = mu_B/2.
     Les nombres sont donc ajustes a Dirac, pas derives ; et Dirac les a avec G_E = G_M = 1.
  B. Dirac au repos : Zitterbewegung calcule (alpha(t) = e^{iHt} alpha e^{-iHt}, H = beta m) :
     frequence 2 m c^2/hbar, amplitude lambda-bar/2, vitesse c ; action par periode
     m c^2 x pi hbar/(m c^2) = pi hbar : le demi-quantum de R86-R87 (E_circ T = pi hbar) est
     celui du Zitterbewegung.  Le tour de la base (omega = m c^2/hbar, R = lambda-bar) est
     l'horloge de de Broglie ; le Zitterbewegung en est le double.
  C. ce que (b) change : R53 (mu_B) -> IDENTITE par construction ; g = 2 -> entree (Dirac) ;
     R84-R85 (motif de charge) -> sans objet pour la charge, l'orientation reste ; R89 -> retire ;
     R99 A tient ; R54 (moitie statique = trois jonctions) : 3 hbar c/(pi^2 D) = m/2 est une
     identite qui DEFINIT D = 6 lambda-bar/pi^2 ; son seul contenu est la coincidence
     6/pi^2 = 0,608 contre 2 cosh(pi)/37,1 = 0,625 (2,7 %) : COINCIDENCE.
  D. ce qui reste a la base, hors de Dirac : le spectre de charge des brins, la structure des
     generations (Koide, conditionnelle), l'echelle des quarks et le nucleon, la tension forte,
     la force nucleaire, la loi de fuite (identite avec G_F).  L'electron n'y predit plus rien
     au-dela de Dirac.
"""
import math
import numpy as np
from scipy.linalg import expm

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
LAMBDA_BAR = HBARC / ME

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R103 -- decision : l'electron electromagnetique est celui de Dirac\n")

    # A. identites de Compton (unites hbar = c = 1, m = 1 ; R en lambda-bar)
    print("A. Les nombres de l'electron dans la base : identites de Compton ajustees a Dirac")
    cases = [("base R53 : R = lambda-bar, E_circ = m/2, q = e", 1.0, 0.5, 1.0),
             ("toute la masse en circulation a R = lambda-bar", 1.0, 1.0, 1.0),
             ("rayon du Zitterbewegung lambda-bar/2, toute la masse", 0.5, 1.0, 1.0)]
    print(f"    {'lecture':52s} {'S/hbar':>8s} {'mu/mu_B':>8s} {'g':>6s}")
    g_base = None
    for name, R, Ecirc, q in cases:
        S = R * Ecirc                     # R E/c en hbar
        mu = q * R / 2                    # q c R/2 en unites e hbar/(2m) = mu_B  (mu_B = e/2 en ces unites)
        g = (mu / S) / (q / 2)            # mu = g (q/2m) S
        if g_base is None:
            g_base = g
        print(f"    {name:52s} {S:8.3f} {mu:8.3f} {g:6.2f}")
    print("  S = hbar/2 et mu = mu_B ne sortent que du choix (R = lambda-bar, E_circ = m/2, q = e) : g = 2 est le partage")
    print("  moitie-moitie pose en R53, pas une sortie ; Dirac a S = hbar/2, mu = mu_B, g = 2 avec G_E = G_M = 1.\n")
    check("R53 : S = hbar/2, mu = mu_B, g = 2 sont des identites du choix (lambda-bar, m/2, e) ; g = 1 ou mu_B/2 sinon",
          abs(g_base - 2) < 1e-12, "identite par construction")

    # B. Zitterbewegung de Dirac au repos
    print("B. Zitterbewegung de Dirac au repos (unites hbar = c = m = 1)")
    I2, Z2 = np.eye(2), np.zeros((2, 2))
    sz = np.array([[1, 0], [0, -1]], complex)
    beta = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
    alpha_z = np.block([[Z2, sz], [sz, Z2]])
    H = beta                                              # p = 0, m = 1
    psi = np.array([1, 0, 1, 0], complex) / math.sqrt(2)  # superposition +E / -E, meme spin
    ts = np.linspace(0, 4 * math.pi, 4001)
    x = []
    pos = 0.0
    for i, t in enumerate(ts):
        U = expm(-1j * H * t)
        v = np.vdot(U @ psi, alpha_z @ (U @ psi)).real     # <v_z>(t)
        x.append(v)
    v = np.array(x)
    # <v_z>(t) = cos(2 t) pour cet etat ; x(t) = sin(2t)/2 : amplitude 1/2 lambda-bar, frequence 2 m c^2/hbar
    vmax = float(np.max(np.abs(v)))
    xt = np.concatenate([[0], np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(ts))])
    amp = 0.5 * (xt.max() - xt.min())
    # frequence : premier retour de v a sa valeur initiale
    k = next(i for i in range(10, len(v)) if v[i] > 0.999 * vmax and v[i - 1] <= v[i] and v[i] >= v[i + 1])
    period = ts[k]
    omega = 2 * math.pi / period
    action = 1.0 * period                                 # E T avec E = m c^2
    print(f"    |v_z| max = {vmax:.4f} c ; amplitude de x = {amp:.4f} lambda-bar ; omega = {omega:.4f} m c^2/hbar ; action par periode E T = {action/math.pi:.4f} pi hbar")
    print(f"    base (R53, R86) : R = 1 lambda-bar, omega = 1 m c^2/hbar (horloge de de Broglie), E_circ T = (1/2)(2 pi) = 1 pi hbar")
    print("  => meme vitesse c, meme action pi hbar par periode ; rayon et frequence differents d'un facteur 2 (m/2 a lambda-bar")
    print("     contre m a lambda-bar/2). Le demi-quantum de R86-R87 est celui du Zitterbewegung.\n")
    check("Zitterbewegung : v = c, amplitude lambda-bar/2, omega = 2 m c^2/hbar, action pi hbar par periode (= R86)",
          abs(vmax - 1) < 1e-3 and abs(amp - 0.5) < 2e-3 and abs(omega - 2) < 2e-3 and abs(action / math.pi - 1) < 2e-3,
          f"{amp:.3f}, {omega:.3f}, {action/math.pi:.3f} pi")

    # C. R54 : identite + coincidence
    print("C. R54 revisite : la moitie statique en trois jonctions")
    D_id = 6 * LAMBDA_BAR / math.pi ** 2
    E_j = 3 * HBARC / (math.pi ** 2 * D_id)            # MeV
    ratio_id, ratio_ms = 6 / math.pi ** 2, 2 * math.cosh(math.pi) / 37.1
    print(f"    3 hbar c/(pi^2 D) = m/2 <=> D = 6 lambda-bar/pi^2 = {D_id:.1f} fm (E = {E_j:.4f} MeV = m/2 : identite, D est defini par m/2)")
    print(f"    contenu restant : 6/pi^2 = {ratio_id:.4f} contre 2 cosh(pi)/37,1 = {ratio_ms:.4f} (manuscrit) : {100*(ratio_ms/ratio_id-1):+.1f} %, une coincidence.\n")
    check("R54 : identite (D defini par m/2) + coincidence 6/pi^2 vs 2cosh(pi)/37,1 a 2,7 %",
          abs(E_j - ME / 2) < 1e-9 and abs(ratio_ms / ratio_id - 1 - 0.028) < 0.003, f"{100*(ratio_ms/ratio_id-1):+.1f} %")

    # D. decision
    print("D. Decision")
    print("  (b) : l'electron electromagnetique est celui de Dirac (point, G_E = G_M = 1, mu_B et g = 2 par la structure")
    print("  spinorielle, Darwin = m/E). La base garde ce qui ne depend pas de la position de la charge : spectre de")
    print("  charge des brins, structure des generations (conditionnelle), echelle des quarks, nucleon, tension forte,")
    print("  force nucleaire, loi de fuite (identite avec G_F). L'anneau a lambda-bar reste comme circulation d'energie")
    print("  neutre (horloge de de Broglie), sans charge dessus, et ne predit rien au-dela de Dirac pour l'electron.\n")

    print("Verdict : DECISION (b). Les nombres de l'electron de la base sont des identites de Compton ajustees a")
    print("Dirac (S, mu, g), que Dirac possede sans boucle ; leur seule sortie propre, le demi-quantum pi hbar, est")
    print("l'action du Zitterbewegung. Le secteur electron de la base est une reformulation ; le contenu propre de la")
    print("base est ailleurs (brins, generations, hadrons).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
