#!/usr/bin/env python3
"""R132 -- Le maillon infiniment fin : la limite continue de la chaine, ses conditions, ce qu'elle laisse.

Decision de l'auteur (22 novembre 2026), devant R131 : « ok prend le maillon infiniment fin ».
C'est la limite M -> infini, m -> 0 a impulsion P+ = M m fixee : la chaine devient la corde continue.

  A. La limite existe. L'ecart de la chaine a la corde, 1 - sin(x) / x avec x = n pi / 2M, tombe comme
     n^2 pi^2 / (24 M^2) : la brisure de l'invariance de Lorentz de R131 C s'eteint comme 1 / M^2.
  B. Elle n'est invariante de Lorentz qu'en dimension critique. Le premier niveau a D - 2 etats, le nombre d'etats
     d'un vecteur sans masse : sa masse doit etre nulle. Avec l'ordonnee de la chaine (R131 B), a = (D - 2) / 24,
     M^2 = (1 - (D - 2) / 24) / alpha' : nul si et seulement si D = 26. Avec un fermion de feuillet par dimension
     (modes demi-entiers, somme de chaine + pi / 24) : a = (D - 2) / 16, M^2 = (1/2 - (D - 2) / 16) / alpha' : nul si
     et seulement si D = 10. Le coefficient d'anomalie de l'algebre de Lorentz,
     m [1 - (D - 2) / 24] + (1 / m) [(D - 2) / 24 - a], s'annule pour tout m dans ce cas seulement. Le choix de dix
     dimensions (R128) devient une consequence. A quatre dimensions : premier niveau a 1033 MeV, 2 etats pour 3.
     Portee : corde fondamentale, valable a toute longueur ; une corde effective longue reste permise (R127 B).
  C. Les contre-termes. Chaine bosonique : energie propre des maillons 4 M^2 / pi et energie des bouts - M par
     dimension, a compenser exactement ; precision relative requise pi^2 / (48 M^2), soit 1,7e-30 pour 3,5e14
     maillons, infinie a la limite. Un boson et un fermion de memes frequences : tout s'annule sans reglage.
  D. Le dipole continu. N dipoles bout a bout de charge q sur une longueur L : moment total q L pour tout N, charges
     interieures nulles, potentiel identique a celui des deux charges de bout. Le moment de chaque maillon, q L / N,
     tend vers zero : le filament dipolaire est une ligne uniformement polarisee, de polarisation lineique q.
  E. La taille du fil. Ecart quadratique moyen d'un maillon au centre, par dimension : alpha' (ln M + c). Il croit
     comme le logarithme du nombre de maillons et diverge a la limite : un fil infiniment fin est infiniment etale,
     et ce que l'on en voit depend de la resolution. C'est le retrecissement de Regge (pente de diffraction
     B = B0 + 2 alpha' ln s) et l'elargissement logarithmique du tube de flux (Luscher, Munster, Weisz 1981) ; la
     constante c vaut gamma d'Euler + ln(4 / pi). Le seul nombre qui y entre est la tension.
  F. Les longueurs. Il ne reste aucune taille d'element. La seule longueur tiree de la tension est la longueur de
     corde, l_s = l1(9) / (pi sqrt(2)) = 0,183 fm, amplitude des vibrations et non taille d'un maillon.
  Verdict : le maillon infiniment fin est la corde continue (IDENTITE). La limite n'est coherente qu'a dix dimensions
  avec des fermions (DERIVE, classique) : le choix de R128 en decoule. Elle ne laisse aucun element a caracteriser :
  le fil est fixe par sa tension seule. Aucun nombre nouveau.
"""
from fractions import Fraction
from math import log, pi, sin, sqrt, tan

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = 0.1973269804                       # GeV fm
LAMBDA_E = 386.15926796                    # fm
L1_3 = 2 * pi * LAMBDA_E / 3
L1_9 = L1_3 * 3 ** (-2 * pi)
SIGMA = pi * HBARC ** 2 / L1_9 ** 2        # GeV^2
ALPHA_P = 1 / (2 * pi * SIGMA)             # GeV^-2
L_S = sqrt(ALPHA_P) * HBARC                # fm

def main():
    print("R132 -- le maillon infiniment fin : la limite continue de la chaine\n")

    # A
    print("A. La limite existe")
    dev = lambda n, M: 1 - sin(n * pi / (2 * M)) / (n * pi / (2 * M))
    print("    ecart relatif du mode n de la chaine au mode de la corde, 1 - sin(x) / x, x = n pi / 2M")
    print("    M            n = 1        n = 2        n = 5        n^2 pi^2 / (24 M^2) pour n = 5")
    for M in (10, 100, 10 ** 4, 10 ** 6):
        print(f"    {M:<10d}   {dev(1, M):.3e}    {dev(2, M):.3e}    {dev(5, M):.3e}    {25 * pi**2 / (24 * M * M):.3e}")
    scaling = [dev(n, 100) / dev(n, 1000) for n in (1, 2, 5)]
    leading = all(abs(dev(n, 10 ** 4) / (n * n * pi ** 2 / (24 * 1e8)) - 1) < 1e-4 for n in (1, 2, 5))
    print(f"    dix fois plus de maillons : ecart divise par " + ", ".join(f"{s:.1f}" for s in scaling))
    print("    -> la brisure de R131 C s'eteint comme 1 / M^2 : la limite continue existe, c'est la corde.\n")
    check("A. limite continue : ecart chaine - corde = n^2 pi^2 / (24 M^2), divise par 100 quand M est multiplie par 10",
          all(abs(s - 100) < 1 for s in scaling) and leading, "1 / M^2")

    # B
    print("B. La limite n'est invariante de Lorentz qu'en dimension critique")
    rows = []
    for M in (10, 100, 1000):
        integer = 2 * M * sum(sin(n * pi / (2 * M)) for n in range(1, M))
        half = 2 * M * sum(sin((n - 0.5) * pi / (2 * M)) for n in range(1, M + 1))
        rows.append((M, integer - 4 * M * M / pi + M, half - 4 * M * M / pi, abs(half - M / sin(pi / (4 * M)))))
    for M, ci, ch, err in rows:
        print(f"    M = {M:4d} : modes entiers, terme constant {ci:.7f} ; modes demi-entiers, terme constant {ch:.7f}")
    print(f"    limites : -pi / 12 = {-pi/12:.7f} (un boson), +pi / 24 = {pi/24:.7f} (modes demi-entiers) ; en unites de 1 / alpha' : -1/24 et +1/48")
    sums_ok = abs(rows[-1][1] + pi / 12) < 1e-7 and abs(rows[-1][2] - pi / 24) < 1e-7 and all(r[3] < 1e-6 * r[0] ** 2 for r in rows)
    a_bos = lambda D: Fraction(D - 2, 24)
    a_ns = lambda D: Fraction(D - 2, 24) + Fraction(D - 2, 48)
    m2_bos = lambda D: 1 - a_bos(D)
    m2_ns = lambda D: Fraction(1, 2) - a_ns(D)
    sol_bos = [D for D in range(2, 40) if m2_bos(D) == 0]
    sol_ns = [D for D in range(2, 40) if m2_ns(D) == 0]
    print("    premier niveau : D - 2 etats, ceux d'un vecteur sans masse (un vecteur massif en demande D - 1) -> sa masse doit etre nulle")
    print(f"    chaine bosonique : alpha' M^2 = 1 - (D - 2) / 24 ; nul pour D = {sol_bos}")
    print(f"    avec un fermion de feuillet par dimension : a = (D - 2) / 16, alpha' M^2 = 1/2 - (D - 2) / 16 ; nul pour D = {sol_ns}")
    delta = lambda m, D, a: m * (1 - Fraction(D - 2, 24)) + Fraction(1, m) * (Fraction(D - 2, 24) - a)
    closes = {(D, a): all(delta(m, D, a) == 0 for m in range(1, 21)) for D, a in ((26, Fraction(1)), (4, Fraction(1, 12)), (4, Fraction(1)), (10, Fraction(1, 3)))}
    second = all(Fraction(D - 2, 24) - a_bos(D) == 0 for D in range(2, 40))
    print(f"    coefficient d'anomalie m [1 - (D - 2) / 24] + (1 / m) [(D - 2) / 24 - a] nul pour tout m : " + ", ".join(f"(D = {D}, a = {a}) {v}" for (D, a), v in closes.items()))
    print(f"    avec l'ordonnee propre de la chaine le second crochet est nul en toute dimension : {second} ; il reste la seule condition D = 26")
    m1_d4 = sqrt(float(m2_bos(4)) / ALPHA_P) * 1e3
    print(f"    a quatre dimensions : premier niveau a sqrt(11/12 / alpha') = {m1_d4:.0f} MeV, 2 etats pour 3 requis : le resultat de R127 B")
    print("    -> le maillon infiniment fin n'est coherent qu'a 26 dimensions, ou 10 avec des fermions : le choix de R128 en decoule.")
    print("       Portee : corde fondamentale, valable a toute longueur. Une corde effective, longue devant l_s, reste permise a quatre dimensions (R127 B).\n")
    check("B. premier niveau sans masse si et seulement si D = 26 (chaine bosonique) ou D = 10 (avec fermions) ; sommes de chaine -pi/12 et +pi/24 ; D = 4 : 1033 MeV, 2 etats pour 3",
          sums_ok and sol_bos == [26] and sol_ns == [10] and closes[(26, Fraction(1))] and not closes[(4, Fraction(1, 12))] and not closes[(4, Fraction(1))]
          and second and abs(m1_d4 - 1033) < 1, "D = 26 ; D = 10")

    # C
    print("C. Les contre-termes de la limite")
    print("    par dimension, en unites de la tension : boson seul            4 M^2 / pi - M - pi / 12")
    print("                                             boson + fermion (demi-entiers)   - M - pi / 8  (l'energie propre des maillons s'annule)")
    print("                                             boson + fermion (memes frequences)        0")
    M_iso = pi / sqrt(8e-29)
    prec = pi ** 2 / (48 * M_iso ** 2)
    cancel = all(abs((2 * M * sum(sin(n * pi / (2 * M)) for n in range(1, M)) - 2 * M * sum(sin((n - 0.5) * pi / (2 * M)) for n in range(1, M + 1))) - (-M - pi / 8)) < 0.03 / M ** 2 + 1e-8
                 for M in (10, 100, 1000))                                   # egalite a des termes en 1 / M^2 pres
    print(f"    difference boson - fermion verifiee sur la chaine : {cancel}")
    print(f"    chaine bosonique : compensation requise a pi^2 / (48 M^2) pres, soit {prec:.1e} pour {M_iso:.1e} maillons (R131 C), et exacte a la limite")
    print("    -> le maillon infiniment fin demande deux contre-termes ajustes exactement, ou la supersymetrie, qui les rend inutiles.\n")
    check("C. limite bosonique : deux contre-termes (4 M^2 / pi et - M) a ajuster a 1,7e-30 pres pour 3,5e14 maillons ; boson + fermion : l'energie propre s'annule",
          cancel and 1e-30 < prec < 3e-30, f"{prec:.1e}")

    # D
    print("D. Le dipole continu")
    q, L = 1.0, 1.0
    probe = (0.37, 2.3, -1.1)
    def potential(N):
        v = 0.0
        for k in range(N):
            z0, z1 = -L / 2 + k * L / N, -L / 2 + (k + 1) * L / N
            for z, s in ((z1, +q), (z0, -q)):
                v += s / sqrt(probe[0] ** 2 + probe[1] ** 2 + (probe[2] - z) ** 2)
        return v
    ref = potential(1)
    same = all(abs(potential(N) - ref) < 1e-9 for N in (2, 5, 50, 2000))
    moments = [N * q * (L / N) for N in (1, 2, 5, 50, 2000)]
    print(f"    N dipoles bout a bout (charge q, longueur L / N) : moment total {moments} pour N = 1, 2, 5, 50, 2000")
    print(f"    potentiel en un point quelconque identique a celui des deux charges de bout : {same} ; moment d'un maillon q L / N -> 0")
    print("    -> le filament dipolaire est une ligne uniformement polarisee : polarisation par unite de longueur q, charges liees aux deux bouts seulement.\n")
    check("D. chaine de N dipoles : moment total q L et champ des deux charges de bout pour tout N ; moment d'un maillon nul a la limite",
          same and all(abs(m - q * L) < 1e-12 for m in moments), "q L")

    # E
    print("E. La taille d'un fil aux maillons infiniment fins")
    size = lambda M: (pi / (2 * M)) * sum(1 / sin(n * pi / (2 * M)) for n in range(1, M))          # <x^2> / alpha', par dimension
    vals = {M: size(M) for M in (10, 100, 1000, 10 ** 4, 10 ** 5)}
    for M, v in vals.items():
        print(f"    M = {M:6d} : <x^2> = {v:.4f} alpha' ; moins ln M : {v - log(M):.4f}")
    step = vals[10 ** 5] - vals[10 ** 4]
    c = vals[10 ** 5] - log(10 ** 5)
    ls2 = L_S ** 2
    M_iso_ln = log(M_iso) + c
    c_exact = 0.5772156649015329 + log(4 / pi)                    # constante d'Euler + ln(4 / pi)
    print(f"    dix fois plus de maillons : + {step:.4f} alpha' (ln 10 = {log(10):.4f}) -> <x^2> = alpha' (ln M + {c:.4f}) par dimension, divergent a la limite")
    print(f"    la constante est gamma d'Euler + ln(4 / pi) = {c_exact:.4f}")
    print(f"    avec la tension de la base : alpha' = {ALPHA_P:.4f} GeV^-2, soit {ls2:.4f} fm^2 par facteur e sur le nombre de maillons")
    print(f"    pour {M_iso:.1e} maillons : ecart type {sqrt(ls2 * M_iso_ln):.2f} fm par dimension, {sqrt(2 * ls2 * M_iso_ln):.2f} fm dans le plan transverse")
    lo, hi = sqrt(ls2 * (log(1e3) + c)), sqrt(ls2 * (log(1e30) + c))
    print(f"    de mille a 1e30 maillons : {lo:.2f} a {hi:.2f} fm par dimension : l'ordre de grandeur d'un hadron, sans valeur definie")
    print(f"    retrecissement de Regge : pente de diffraction B = B0 + 2 alpha' ln s, coefficient 2 alpha' = {2*ALPHA_P:.2f} GeV^-2 avec la tension de la base")
    print(f"    elargissement du tube de flux avec sa longueur R : w^2 = (D - 2) alpha' ln(R / R0), soit {2*ls2:.4f} fm^2 par facteur e a quatre dimensions")
    print("       (meme nombre que R119 et R120 : la tension ; rien de nouveau a confronter)")
    print("    -> un fil infiniment fin est infiniment etale, lentement ; ce que l'on en voit depend de la resolution. Aucune taille propre.\n")
    check("E. taille du fil : <x^2> = alpha' (ln M + c) par dimension, + ln 10 par decade de maillons, divergente a la limite ; coefficient de Regge 2 alpha' = 1,72 GeV^-2",
          abs(step - log(10)) < 1e-3 and abs(c - c_exact) < 1e-3 and abs(2 * ALPHA_P - 1.719) < 2e-3 and vals[10 ** 5] > vals[10] and 0.4 < lo < 0.6 and 1.4 < hi < 1.7, f"ln M + {c:.3f}")

    # F
    print("F. Les longueurs qui restent")
    ident = abs(L_S - L1_9 / (pi * sqrt(2))) < 1e-12
    print(f"    longueur de Compton de l'electron {LAMBDA_E:.2f} fm : une masse mesuree ; l1(3) = {L1_3:.1f} fm et l1(9) = {L1_9:.4f} fm : construites sur elle avec les entiers poses")
    print(f"    longueur de corde l_s = hbar c sqrt(alpha') = l1(9) / (pi sqrt(2)) = {L_S:.4f} fm : identite {ident}")
    print("    aucune n'est une taille de maillon (R131 C) ; l_s est l'amplitude des vibrations : chaque mode n ajoute alpha' / n a la taille")
    print("    parametres du fil : la tension, et la probabilite de jonction ; le reste (espace, branes) vient de R128 a R130")
    print("    -> la question d'origine, caracteriser l'element du filament, recoit sa reponse dans cette lecture : il n'y a pas d'element.\n")
    check("F. aucune taille d'element : la seule longueur tiree de la tension est l_s = l1(9) / (pi sqrt(2)) = 0,183 fm, amplitude de vibration", ident and abs(L_S - 0.1829) < 1e-3, f"{L_S:.4f} fm")

    print("Verdict : le maillon infiniment fin est la corde continue. La limite existe, mais elle n'est coherente qu'a dix dimensions avec des")
    print("fermions, ce qui fait du choix de R128 une consequence, et elle demande la supersymetrie ou deux reglages exacts. Elle ne laisse aucun")
    print("element a caracteriser : le fil est fixe par sa tension seule, et sa taille apparente croit avec la resolution. Aucun nombre nouveau.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
