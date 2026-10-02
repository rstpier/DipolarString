#!/usr/bin/env python3
"""R114 -- L'anneau a 1,34 fm : recalcul de mu_p (« Mets l'anneau a 1,34 fm et recalcule mu_p »).

R113 C : si les circuits de quark sont des cercles, l'anneau separe doit porter un nombre pair de
quanta, et les masses (N, Delta - N) le mettent a R = 1,344 fm avec (n_N, n_Delta) = (2, 6). Ici,
l'anneau est pose a ce rayon et le moment magnetique est recalcule avec la regle de la base
(R25, R72) : la charge circulante Q tourne a c sur l'anneau, mu = Q (R / lambda-bar_p) mu_N, avec le
partage de SU(6) : Q_p = (2/3)(2/3) + (2/3)(2/3) + (-1/3)(-1/3) = 1, Q_n = -2/3.

  A. mu_p, mu_n, mu(Delta++) a 1,344 fm avec le partage : mu_p = 6,39 mu_N (+129 %), mu_n = -4,26
     (+123 %), Delta++ = 12,8 mu_N (mesure 3,7-7,5) ; le rapport mu_p/mu_n = -3/2 ne voit pas R.
  B. ce que le moment fixe : le produit Q R = mu_p lambda-bar_p = 0,587 fm e, pas R ; a 1,344 fm il
     faut Q = 0,437 e ; aucune fraction de la base ne le donne (4/9 a 1,7 %, sans lecture).
  C. le rayon de charge : un anneau de 0,437 e a 1,344 fm a <r^2> = 0,789 fm^2 > r_p^2 = 0,707 ;
     r_p >= 0,888 fm (+5,6 %) meme avec le reste de la charge au centre ; 0,98 fm (+16 %) avec le
     reste sur les circuits. Exclu par r_p (0,05 %). A 0,672 fm avec 0,874 e : 0,64 fm (-22 %).
  D. mu_p et r_p ensemble fixent (R, Q) = (0,99 fm, 0,59 e) ; avec N et Delta - N le compte n'est
     pas entier (1,8 et 2,9 quanta) ; a n entier le meilleur est (2, 5) a 2,4 %, de total ENTIER ;
     a n pair, (2, 4) met Delta - N a 199 MeV (-32 %). Le compte (R113), mu_p et r_p ne laissent
     aucun anneau.
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME, MP, MDELTA = 0.51099895, 938.27209, 1232.0
MU_P, MU_N = 2.79284734, -1.91304273
MU_DPP = (3.7, 7.5)                       # mu(Delta++) en mu_N, fourchette PDG (R48)
RP = 0.8409
LAM = HBARC / ME
L1 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)
R_C = 3 * L1 / (2 * math.pi)
E_STATIC = math.pi * HBARC / L1
LAM_P = HBARC / MP
R_RING = 2 * math.sqrt(3) * R_C           # 1,344 fm : (2, 6) quanta, R113 C
R_48 = math.sqrt(3) * R_C                 # 0,672 fm : (1, 3) quanta, R48
R2_CIRC = 2 * R_C ** 2                    # <r^2> d'un circuit de l'etoile par le centre (d = R_c)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def mu(Q, R):
    """moment d'une charge Q (en e) tournant a c sur un cercle de rayon R, en magnetons nucleaires"""
    return Q * R / LAM_P

def main():
    print("R114 -- l'anneau a 1,34 fm : recalcul de mu_p\n")
    print(f"  R = 2 sqrt3 R_c = {R_RING:.4f} fm (R113 C), lambda-bar_p = {LAM_P:.4f} fm, R/lambda-bar_p = {R_RING/LAM_P:.3f}")
    E1 = HBARC / (2 * R_RING)
    print(f"  masses a ce rayon : quantum {E1:.1f} MeV ; N = {E_STATIC + 2*E1:.1f} ({100*((E_STATIC+2*E1)/MP-1):+.1f} %), Delta = {E_STATIC + 6*E1:.1f} ({100*((E_STATIC+6*E1)/MDELTA-1):+.1f} %)\n")

    # A. le partage de la base
    print("A. Le moment avec la regle de la base (charge circulante a c, partage 2/3, 2/3, -1/3)")
    w = (2 / 3, 2 / 3, -1 / 3)
    Qp = sum(wi * qi for wi, qi in zip(w, (2 / 3, 2 / 3, -1 / 3)))
    Qn = sum(wi * qi for wi, qi in zip(w, (-1 / 3, -1 / 3, 2 / 3)))
    mup, mun, mudpp = mu(Qp, R_RING), mu(Qn, R_RING), mu(2.0, R_RING)
    print(f"    Q_p = {Qp:.3f} e, Q_n = {Qn:.3f} e")
    print(f"    mu_p = {mup:.2f} mu_N (mesure {MU_P:.3f}, {100*(mup/MU_P-1):+.0f} %) ; mu_n = {mun:.2f} (mesure {MU_N:.3f}, {100*(mun/MU_N-1):+.0f} %)")
    print(f"    mu(Delta++) (2 e) = {mudpp:.1f} mu_N (fourchette {MU_DPP[0]}-{MU_DPP[1]}) ; rapport mu_p/mu_n = {mup/mun:.3f} (mesure {MU_P/MU_N:.3f}, inchange : le rapport ne voit pas R)")
    print(f"    pour memoire, aux rayons des lectures : 0,587 fm -> {mu(1, 0.5874):.3f} ; 0,672 fm -> {mu(1, R_48):.2f} ; 1,344 fm -> {mup:.2f}\n")
    check("A. a 1,344 fm, le partage de la base donne mu_p = 6,39 mu_N (+129 %), mu_n = -4,26, Delta++ hors fourchette",
          abs(mup / MU_P - 2.29) < 0.02 and mudpp > MU_DPP[1], f"mu_p/mesure = {mup/MU_P:.2f}, Delta++ {mudpp:.1f}")

    # B. ce que le moment fixe
    print("B. Ce que le moment fixe : le produit Q R, pas R")
    QR = MU_P * LAM_P
    Qneed = QR / R_RING
    Qneed_n = MU_N * LAM_P / R_RING
    print(f"    Q R = mu_p lambda-bar_p = {QR:.4f} fm e, independant du nombre de quanta (la charge tourne une fois par tour a c)")
    print(f"    a 1,344 fm : Q_p = {Qneed:.3f} e (R49 : 0,874 e a 0,672 fm ; R72 : 1,045 e a 0,562 fm) ; Q_n = {Qneed_n:.3f} e")
    fracs = {"1": 1.0, "2/3": 2 / 3, "1/2": 0.5, "4/9": 4 / 9, "1/3": 1 / 3, "2/9": 2 / 9, "1/(2 sqrt alpha) (fluide)": 1 / (2 * math.sqrt(ALPHA)),
             "alpha^(1/2)": math.sqrt(ALPHA)}
    closest = min(fracs.items(), key=lambda kv: abs(kv[1] / Qneed - 1))
    for k, v in fracs.items():
        print(f"      {k:28s} {v:.4f}  ({100*(v/Qneed-1):+6.1f} %)")
    print(f"    la plus proche : {closest[0]} a {100*abs(closest[1]/Qneed-1):.1f} %, sans lecture dans la base (le partage donne 1).\n")
    check("B. le moment fixe Q R = 0,587 fm e ; a 1,344 fm il faut 0,437 e, qu'aucune fraction de la base ne donne a 1 %",
          abs(Qneed - 0.437) < 0.001 and abs(closest[1] / Qneed - 1) > 0.01, f"Q = {Qneed:.3f} e ; {closest[0]} a {100*abs(closest[1]/Qneed-1):.1f} %")

    # C. le rayon de charge
    print("C. Le rayon de charge avec la charge circulante sur l'anneau")
    def r_charge(R, Q, rest="circuits"):
        r2_rest = R2_CIRC if rest == "circuits" else 0.0
        return math.sqrt(Q * R ** 2 + (1 - Q) * r2_rest)
    rows = []
    for R, Q, label in ((R_RING, Qneed, "1,344 fm, 0,437 e (compte pair)"), (R_48, QR / R_48, "0,672 fm, 0,874 e (R48/R49)"), (0.5874, 1.0, "0,587 fm, 1 e (R25/R72)")):
        rc, r0 = r_charge(R, Q), r_charge(R, Q, rest="centre")
        rows.append((label, rc, r0))
        print(f"    {label:34s} <r^2>_anneau = {Q*R**2:.3f} fm^2 ; r_p = {rc:.3f} fm ({100*(rc/RP-1):+.1f} %) reste sur les circuits, >= {r0:.3f} fm ({100*(r0/RP-1):+.1f} %) reste au centre")
    print(f"    mesure r_p = {RP} fm (0,05 %) ; <r^2> d'un circuit de l'etoile par le centre = 2 R_c^2 = {R2_CIRC:.3f} fm^2")
    print("  -> a 1,344 fm, la seule charge que mu_p exige sur l'anneau depasse deja r_p^2 : exclu quel que soit le reste.")
    print("     a 0,672 fm le meme calcul donne -22 % : l'image « charge circulante sur l'anneau » n'a jamais eu le rayon de charge.\n")
    check("C. l'anneau a 1,344 fm portant les 0,437 e de mu_p a r_p >= 0,888 fm (+5,6 %) meme avec le reste au centre : exclu par r_p",
          rows[0][2] / RP - 1 > 0.05 and Qneed * R_RING ** 2 > RP ** 2, f"r_p >= {rows[0][2]:.3f} fm ; anneau seul {Qneed*R_RING**2:.3f} > {RP**2:.3f} fm^2")

    # D. mu_p et r_p ensemble, puis le compte
    print("D. mu_p et r_p ensemble : (R, Q) ; puis le compte avec N et Delta - N")
    # Q R = QR ; Q R^2 + (1 - Q) R2_CIRC = RP^2  ->  QR R + R2_CIRC - QR R2_CIRC / R = RP^2
    a, b, c = QR, R2_CIRC - RP ** 2, -QR * R2_CIRC
    Rj = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
    Qj = QR / Rj
    E1j = HBARC / (2 * Rj)
    nN_real, dn_real = (MP - E_STATIC) / E1j, (MDELTA - MP) / E1j
    print(f"    R = {Rj:.3f} fm, Q = {Qj:.3f} e (reste sur les circuits) : mu_p et r_p exacts par construction")
    print(f"    quanta demandes par N et Delta - N a ce rayon : n_N = {nN_real:.2f}, n_Delta - n_N = {dn_real:.2f} : non entiers")
    best_any, best_even = None, None
    for nN in range(1, 7):
        for nD in range(nN + 1, 13):
            N, D = E_STATIC + nN * E1j, E_STATIC + nD * E1j
            dev = max(abs(N / MP - 1), abs((nD - nN) * E1j / (MDELTA - MP) - 1))     # sur N et sur l'ecart Delta - N (R48)
            if best_any is None or dev < best_any[0]:
                best_any = (dev, nN, nD, N, D)
            if nN % 2 == 0 and nD % 2 == 0 and (best_even is None or dev < best_even[0]):
                best_even = (dev, nN, nD, N, D)
    for label, b in (("meilleur compte entier", best_any), ("meilleur compte pair (R113 C)", best_even)):
        print(f"    {label:30s} : ({b[1]}, {b[2]}) : N = {b[3]:.0f} ({100*(b[3]/MP-1):+.1f} %), Delta - N = {b[4]-b[3]:.0f} MeV ({100*((b[4]-b[3])/(MDELTA-MP)-1):+.0f} %), Delta = {b[4]:.0f} ({100*(b[4]/MDELTA-1):+.1f} %)")
    print("    (2, 5) : total de spin ENTIER (n_Delta impair) ; (2, 4) : l'ecart Delta - N, la donnee que R48 ajuste, tombe a -32 %.")
    print("  -> quatre nombres (R, Q, n_N, n_Delta) pour quatre donnees (N, Delta - N, mu_p, r_p) : sans le compte, un ajustement a 2,4 %,")
    print("     pas une prediction ; avec le compte pair, Delta - N a -32 %. Le compte, mu_p et r_p ne laissent aucun anneau.\n")
    check("D. mu_p + r_p fixent (0,99 fm, 0,59 e) ; N et Delta - N y demandent des quanta non entiers ; a n pair, Delta - N s'ecarte de > 25 %",
          abs(Rj - 0.99) < 0.01 and abs(nN_real - round(nN_real)) > 0.2 and best_even[0] > 0.25 and best_any[0] < 0.03,
          f"R = {Rj:.3f}, Q = {Qj:.2f}, n_N = {nN_real:.2f}, pair {100*best_even[0]:.0f} %, entier {100*best_any[0]:.1f} %")

    print("Verdict : EXCLU. A 1,34 fm la regle de la base donne mu_p = 6,4 mu_N ; le moment ne fixe que Q R = 0,587 fm e, donc")
    print("0,44 e en circulation, que rien ne donne ; et cette charge a ce rayon depasse a elle seule le rayon de charge du proton.")
    print("Le compte pair (R113), mu_p et r_p sont incompatibles : l'anneau charge a c n'est pas la lecture du moment du nucleon.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
