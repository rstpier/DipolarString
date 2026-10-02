#!/usr/bin/env python3
"""R113 -- Le rayon de l'anneau du nucleon : S = hbar/2 sur la circulation collective le fixe-t-il ?

Le secteur hadronique (R109) a un nombre non fixe : le rayon R de l'anneau qui porte le spin du
nucleon, 0,56 a 0,67 fm selon la lecture (masse, mu_p, Delta - N), qui pese 3 % sur m_p. La
question posee : si l'on exige S = hbar/2 sur la circulation collective des trois circuits de
quark, le rayon sort-il ? Sinon, on l'enregistre comme coincidence (R49) et on s'arrete la.

  A. L'identite. Un quantum de circulation (E = pi hbar c / L) sur un chemin ferme de longueur L
     et d'aire vectorielle A porte S = (2A/L)(E/c) = 2 pi hbar A / L^2 : c'est hbar/2 si et
     seulement si le chemin est un cercle (egalite isoperimetrique), et cela pour TOUT rayon.
     S = hbar/2 ne contient donc aucune information sur R ; l'energie seule le fixe (E = hbar c/2R).
  B. La circulation collective. Le moment cinetique d'un ecoulement ferme autour de n'importe quel
     point est 2A fois sa densite d'impulsion (son impulsion totale est nulle) : les trois circuits
     donnent S = +1/2 + 1/2 - 1/2 = hbar/2 (nucleon) ou 3 hbar/2 (Delta) quelle que soit leur
     distance au centre. La condition est satisfaite pour toute geometrie : elle ne fixe rien.
  C. Le compte. Si les circuits de quark sont des cercles (trois brins bout a bout, comme l'anneau
     de l'electron a l'echelle 3), chacun porte hbar/2 ; un anneau separe (R48, R72) est alors un
     quatrieme circuit, et le total n'est demi-entier que si l'anneau porte un nombre PAIR de
     quanta. Les masses ne fixent que n/R : l'assignation qui garde N et Delta a 3 % est
     (n_N, n_Delta) = (2, 6) a R = 1,34 fm, hors du proton (r_p = 0,84 fm), avec 0,44 e en
     circulation pour mu_p ; la lecture de R48 (1, 3) donne un total entier. Un circuit non
     circulaire (Y bifilaire a l'echelle 9) porte 0,13 hbar (0,04 en sens alterne) : non quantifie.
  D. L'energie. Tous les termes de la base pour l'anneau (quantum, Coulomb de la charge circulante)
     decroissent en 1/R : aucun point stationnaire, aucun rayon d'equilibre.
  E. Les cercles de la geometrie de l'etoile contre les trois lectures : seul sqrt 3 R_c tombe sur
     Delta - N (0,1 %, R49) ; les lectures s'etalent sur 17 % ; coincidence confirmee.
"""
import math
import itertools
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MDELTA = 0.51099895, 938.27209, 1232.0
MU_P = 2.79284734
RP = 0.8409
LAM = HBARC / ME
L1 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)      # l_1(9) = 0.8128 fm
W9 = (4 * LAM / math.pi ** 2) * 3 ** (-2 * math.pi)      # section a l'echelle 9 = 0.1573 fm
R_C = 3 * L1 / (2 * math.pi)                             # rayon du circuit de quark = 0.3881 fm
E_Q = math.pi * HBARC / (3 * L1)                         # 254.2 MeV
E_STATIC = 3 * E_Q                                       # 762.7 MeV
LAM_P = HBARC / MP

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# ----------------------------------------------------------------------------------------------
def polygon(pts):
    """longueur et aire signee (shoelace) d'un chemin ferme donne par ses sommets (N, 2)"""
    x, y = pts[:, 0], pts[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    L = float(np.sum(np.hypot(xn - x, yn - y)))
    A = 0.5 * float(np.sum(x * yn - xn * y))
    return L, A

def S_one_quantum(pts):
    """S/hbar pour un quantum de circulation (E = pi hbar c/L) a la vitesse c sur le chemin"""
    L, A = polygon(pts)
    return 2 * math.pi * A / L ** 2

def circle(R, cx=0.0, cy=0.0, n=400000, sense=+1):
    t = np.linspace(0, 2 * math.pi, n, endpoint=False) * sense
    return np.column_stack([cx + R * np.cos(t), cy + R * np.sin(t)])

def ellipse(a, b, n=400000):
    t = np.linspace(0, 2 * math.pi, n, endpoint=False)
    return np.column_stack([a * np.cos(t), b * np.sin(t)])

def regular_polygon(k, side, n_per=50000):
    Rcirc = side / (2 * math.sin(math.pi / k))
    verts = [(Rcirc * math.cos(2 * math.pi * j / k), Rcirc * math.sin(2 * math.pi * j / k)) for j in range(k)]
    pts = []
    for j in range(k):
        (x0, y0), (x1, y1) = verts[j], verts[(j + 1) % k]
        s = np.linspace(0, 1, n_per, endpoint=False)
        pts.append(np.column_stack([x0 + s * (x1 - x0), y0 + s * (y1 - y0)]))
    return np.vstack(pts)

def trefoil_envelope(Rc, n_per=100000):
    """enveloppe exterieure de trois cercles de rayon Rc passant par l'origine, centres a Rc"""
    pts = []
    for j in range(3):
        phi = 2 * math.pi * j / 3
        cx, cy = Rc * math.cos(phi), Rc * math.sin(phi)
        t = np.linspace(phi - 2 * math.pi / 3, phi + 2 * math.pi / 3, n_per, endpoint=False)
        pts.append(np.column_stack([cx + Rc * np.cos(t), cy + Rc * np.sin(t)]))
    return np.vstack(pts)

def bifilar_Y(arm, w, n_per=50000, alternate=False):
    """circuit en Y bifilaire : trois bras de longueur `arm`, aller sur un conducteur, retour sur
    l'autre, ecart w ; sens de parcours des trois boucles identique ou alterne"""
    pts = []
    for j in range(3):
        phi = 2 * math.pi * j / 3
        u = np.array([math.cos(phi), math.sin(phi)])
        nrm = np.array([-u[1], u[0]]) * (w / 2)
        sgn = -1 if (alternate and j == 1) else +1
        s = np.linspace(0, 1, n_per, endpoint=False)[:, None]
        out = (sgn * nrm) + s * arm * u
        back = (arm * u - sgn * nrm) - s * arm * u
        pts += [out, back]
    return np.vstack(pts)

# ----------------------------------------------------------------------------------------------
def main():
    print("R113 -- le rayon de l'anneau du nucleon : S = hbar/2 sur la circulation collective le fixe-t-il ?\n")
    print(f"  l_1(9) = {L1:.4f} fm ; R_c = 3 l_1/2 pi = {R_C:.4f} fm ; E_q = {E_Q:.1f} MeV ; part statique {E_STATIC:.1f} MeV ; w_9 = {W9:.4f} fm")
    E_mass = MP - E_STATIC
    R_mass = HBARC / (2 * E_mass)
    R_mu = MU_P * LAM_P
    E_mu = HBARC / (2 * R_mu)
    E_delta = (MDELTA - MP) / 2
    R_delta = HBARC / (2 * E_delta)
    readings = {"masse + spin": (E_mass, R_mass), "mu_p (e sur l'anneau)": (E_mu, R_mu), "Delta - N = 2E": (E_delta, R_delta)}
    print("  les trois lectures de l'anneau (R109 B) :")
    for k, (E, R) in readings.items():
        print(f"    {k:24s} E = {E:6.1f} MeV, R = {R:.4f} fm")
    print()

    # A. l'identite
    print("A. S = 2 pi hbar A / L^2 pour un quantum sur un chemin ferme ; hbar/2 ssi cercle, pour tout rayon")
    circ = {f"cercle R = {R:.3f} fm": S_one_quantum(circle(R)) for R in (R_mass, R_mu, R_delta, R_C, 3.0)}
    non = {"ellipse b/a = 0,9": S_one_quantum(ellipse(1.0, 0.9)),
           "ellipse b/a = 0,5": S_one_quantum(ellipse(1.0, 0.5)),
           "triangle equilateral": S_one_quantum(regular_polygon(3, 1.0)),
           "carre": S_one_quantum(regular_polygon(4, 1.0)),
           "enveloppe des trois circuits (trefle)": S_one_quantum(trefoil_envelope(R_C)),
           "Y bifilaire, bras l_1/2, ecart w_9, meme sens": abs(S_one_quantum(bifilar_Y(L1 / 2, W9))),
           "Y bifilaire, sens alterne": abs(S_one_quantum(bifilar_Y(L1 / 2, W9, alternate=True)))}
    for k, v in circ.items():
        print(f"    {k:46s} S/hbar = {v:.9f}")
    for k, v in non.items():
        print(f"    {k:46s} S/hbar = {v:.4f}")
    print("  (Y bifilaire : |S| en valeur absolue, le sens est une convention ; les segments de bout et de jonction, ~w, comptent dans L)")
    print("  -> S = hbar/2 est l'egalite isoperimetrique : vraie pour tout cercle, fausse pour tout autre chemin.")
    print("     Elle ne contient aucune information sur R ; seule l'energie E = hbar c/2R le fixe.\n")
    ok_circ = max(abs(v - 0.5) for v in circ.values()) < 1e-7
    ok_non = all(0 < v < 0.5 for v in non.values())
    check("A. un quantum sur un cercle porte S = hbar/2 a tout rayon (identite, 1e-7) ; tout chemin non circulaire porte moins",
          ok_circ and ok_non, f"cercles {max(abs(v-0.5) for v in circ.values()):.1e} ; non-cercles max {max(non.values()):.3f}")

    # B. la circulation collective
    print("B. Moment cinetique collectif des trois circuits autour du centre du nucleon, selon leur ecartement d")
    pts_n = circle(R_C, n=200000)
    rows = []
    for d in np.linspace(0, 2 * R_C, 9):
        tot = {}
        for label, senses in (("+ + -", (+1, +1, -1)), ("+ + +", (+1, +1, +1))):
            S = 0.0
            for j, sg in enumerate(senses):
                phi = 2 * math.pi * j / 3
                P = circle(R_C, d * math.cos(phi), d * math.sin(phi), n=200000, sense=sg)
                L, A = polygon(P)
                S += 2 * math.pi * A / L ** 2          # aire signee ABSOLUE (origine au centre du nucleon)
            tot[label] = S
        rows.append((d, tot["+ + -"], tot["+ + +"]))
        print(f"    d = {d:.3f} fm ({d/R_C:.2f} R_c) : S(+ + -) = {tot['+ + -']:+.9f} hbar, S(+ + +) = {tot['+ + +']:+.9f} hbar")
    dev = max(max(abs(r[1] - 0.5), abs(r[2] - 1.5)) for r in rows)
    print("  -> l'impulsion totale d'un ecoulement ferme est nulle : son moment cinetique est le meme autour de tout point.")
    print("     S = hbar/2 (deux contre un) et 3 hbar/2 (alignes) sortent a TOUT ecartement : la condition ne fixe rien.\n")
    check("B. S collectif = hbar/2 (+ + -) et 3 hbar/2 (+ + +) independamment de l'ecartement d des circuits (1e-7)",
          dev < 1e-7, f"ecart max {dev:.1e} sur d = 0 a 2 R_c")

    # C. le compte avec un anneau separe
    print("C. Le compte : circuits circulaires (hbar/2 chacun) + un anneau de n quanta (S_anneau = n/2)")
    E1 = HBARC / 2      # energie d'un quantum sur un cercle de 1 fm : hbar c / 2R avec R en fm
    def reachable(n, target):
        for sg in itertools.product((+1, -1), repeat=4):
            if abs(abs(sg[0] * 0.5 + sg[1] * 0.5 + sg[2] * 0.5 + sg[3] * n / 2) - target) < 1e-9:
                return True
        return False
    table = []
    for nN in range(1, 9):
        for nD in range(nN + 1, 13):
            okN, okD = reachable(nN, 0.5), reachable(nD, 1.5)
            if not (okN and okD):
                continue
            # masses ne dependent que de n/R : R depuis Delta - N, puis N ; mu_p -> charge circulante
            R = (nD - nN) * E1 / (MDELTA - MP)
            N = E_STATIC + nN * E1 / R
            D = E_STATIC + nD * E1 / R
            Q = MU_P * LAM_P / R
            table.append((nN, nD, R, N, D, Q))
    print("    assignations a total demi-entier (|S_N| = 1/2, |S_Delta| = 3/2), R fixe par Delta - N = 294 MeV :")
    print("      n_N  n_D   R (fm)   N (MeV)   ecart   Delta (MeV)  ecart   Q_circ pour mu_p")
    best = None
    for nN, nD, R, N, D, Q in table:
        devN, devD = N / MP - 1, D / MDELTA - 1
        flag = " <-- a 3 %" if abs(devN) < 0.03 and abs(devD) < 0.03 else ""
        print(f"      {nN:3d}  {nD:3d}   {R:.4f}   {N:7.1f}  {100*devN:+6.1f} %   {D:7.1f}   {100*devD:+6.1f} %   {Q:.3f} e{flag}")
        if best is None or max(abs(devN), abs(devD)) < best[0]:
            best = (max(abs(devN), abs(devD)), nN, nD, R, Q)
    parity_ok = all(nN % 2 == 0 and nD % 2 == 0 for nN, nD, *_ in table)
    R48 = (E_STATIC + 1 * E1 / R_delta, E_STATIC + 3 * E1 / R_delta)
    print(f"    la lecture de R48, (n_N, n_D) = (1, 3) a R = {R_delta:.3f} fm : N = {R48[0]:.0f}, Delta = {R48[1]:.0f} MeV, mais total ENTIER")
    print(f"    (quatre demi-entiers) : exclue par le compte si les circuits de quark sont des cercles.")
    print(f"    meilleure assignation : (n_N, n_D) = ({best[1]}, {best[2]}), R = {best[3]:.3f} fm = {best[3]/RP:.2f} r_p = {best[3]/(math.sqrt(3)*R_C):.3f} x sqrt3 R_c,")
    print(f"    ecart max {100*best[0]:.1f} %, charge circulante {best[4]:.3f} e pour mu_p (R49 : 0,87 e a 0,67 fm).")
    print(f"    un circuit de quark non circulaire (Y bifilaire, A) porte {non['Y bifilaire, bras l_1/2, ecart w_9, meme sens']:.2f} hbar : total non demi-entier.")
    print("  -> le compte fixe la PARITE de n (paire), les masses fixent n/R : le rayon double (1,34 fm, hors du proton)")
    print("     ou les masses s'ecartent de 10 % ; rien ne fixe R. Les circuits « statiques » de R72 exigent A = 0 :")
    print("     aucun circuit ferme a sens unique de la base ne l'a.\n")
    check("C. avec des circuits circulaires, seul un anneau a nombre pair de quanta garde S demi-entier ; N et Delta a 3 % exigent alors (2, 6) a R = 1,34 fm > r_p",
          parity_ok and best[1] == 2 and best[2] == 6 and best[0] < 0.035 and best[3] > RP,
          f"({best[1]}, {best[2]}), R = {best[3]:.3f} fm, ecart {100*best[0]:.1f} %, Q = {best[4]:.2f} e")

    # D. l'energie
    print("D. L'energie de l'anneau en fonction de R (quantum + Coulomb de la charge circulante, coupure a = w_9/4)")
    a_cut = W9 / 4
    Rgrid = np.linspace(0.15, 3.0, 2000)
    def E_ring(R, n=1, Q=1.0):
        return n * HBARC / (2 * R) + K * Q ** 2 / (2 * math.pi * R) * (np.log(8 * R / a_cut) - 2)
    mono = True
    for n, Q in ((1, 1.0), (2, 1.0), (1, 0.0), (2, 0.44)):
        E = E_ring(Rgrid, n, Q)
        dE = np.diff(E)
        mono &= bool(np.all(dE < 0))
        print(f"    n = {n}, Q = {Q:.2f} e : E({Rgrid[0]:.2f}) = {E[0]:.0f}, E(0,67) = {float(E_ring(0.672, n, Q)):.1f}, E(1,34) = {float(E_ring(1.344, n, Q)):.1f}, E({Rgrid[-1]:.1f}) = {E[-1]:.0f} MeV ; dE/dR < 0 partout : {bool(np.all(dE < 0))}")
    print(f"    Coulomb de e sur l'anneau a 0,67 fm : {float(K/(2*math.pi*0.672)*(math.log(8*0.672/a_cut)-2)):.2f} MeV (1,4 % du quantum)")
    print("  -> tout est homogene de degre -1 en R (le quantum EST le terme de tension) : pas de point stationnaire,")
    print("     pas de rayon d'equilibre. Le rayon est une entree.\n")
    check("D. E(R) strictement decroissante sur [0,15 ; 3] fm pour n = 1, 2 : aucun extremum, aucun rayon d'equilibre", mono, "dE/dR < 0 partout")

    # E. la geometrie contre les lectures
    print("E. Les cercles que la geometrie de l'etoile offre, contre les trois lectures")
    cands = {"R_c/2 (inscrit)": R_C / 2, "R_c (centres, etoile par le centre)": R_C,
             "2 R_c/sqrt3 (centres, etoile tangente)": 2 * R_C / math.sqrt(3), "sqrt3 R_c (ecart entre quarks)": math.sqrt(3) * R_C,
             "2 R_c (bord exterieur)": 2 * R_C, "R_c (1 + 2/sqrt3) (enveloppe tangente)": R_C * (1 + 2 / math.sqrt(3)),
             "l_1(9)": L1, "2 sqrt3 R_c (compte pair)": 2 * math.sqrt(3) * R_C}
    print("      candidat                                   R (fm)   E = hbar c/2R   masse  mu_p   Delta-N")
    hits = {}
    for k, R in cands.items():
        E = HBARC / (2 * R)
        devs = [E / Er - 1 for (Er, _) in readings.values()]
        hits[k] = devs
        print(f"      {k:42s} {R:.4f}   {E:7.1f} MeV   " + "  ".join(f"{100*x:+6.1f} %" for x in devs))
    spread = (E_mass - E_delta) / ((E_mass + E_delta) / 2)
    d_sqrt3 = abs(hits["sqrt3 R_c (ecart entre quarks)"][2])
    print(f"    les trois lectures s'etalent sur {100*spread:.1f} % (147 a 176 MeV) : aucun rayon unique ne les satisfait a mieux que ~9 % ;")
    print(f"    sqrt3 R_c tombe sur Delta - N a {100*d_sqrt3:.2f} % (R49) et sur rien d'autre ; aucun mecanisme ne le selectionne (A, B, D).")
    print("  -> coincidence confirmee, non revendiquee. On s'arrete la.\n")
    check("E. sqrt3 R_c reproduit Delta - N a 0,2 % (R49) ; les lectures s'etalent sur > 10 % ; aucun cercle de la geometrie n'en satisfait deux",
          d_sqrt3 < 0.002 and spread > 0.10 and all(sum(abs(x) < 0.03 for x in devs) <= 1 for devs in hits.values()),
          f"{100*d_sqrt3:.2f} % ; etalement {100*spread:.1f} %")

    print("Verdict : NON. S = hbar/2 est l'egalite isoperimetrique d'un quantum sur un cercle, vraie a tout rayon (A) ;")
    print("la circulation collective des trois circuits donne hbar/2 et 3 hbar/2 a tout ecartement (B) ; aucun terme de la")
    print("base n'a d'extremum en R (D). Le rayon n'est pas fixe ; sqrt3 R_c reste une coincidence (E). En route, le compte :")
    print("si les circuits de quark sont des cercles, l'anneau separe doit porter un nombre pair de quanta, ce qui le met a")
    print("1,34 fm (hors du proton) ou ecarte les masses de 10 % (C) : la lecture (1, 3) de R48 et les circuits « statiques »")
    print("de R72 ne sont compatibles qu'avec des circuits sans aire, que la base n'a pas.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
