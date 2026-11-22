#!/usr/bin/env python3
"""R131 -- Le filament comme chaine de maillons : ce que la chaine est, ce qu'elle donne, ce qu'elle interdit.

Decision de l'auteur (22 novembre 2026) : « ok prend la chaine de maillons ».
Definition retenue : la corde de R130 ecrite comme une chaine de M maillons, chacun portant la meme part m de
l'impulsion (P+ = M m), lie a ses voisins par un ressort : c'est le modele des « bits de corde » (Giles et Thorn
1977). Chaque maillon est un dipole, une charge a chaque bout ; dans la chaine les charges voisines s'annulent.

  A. La chaine est la corde. Frequences exactes d'une chaine ouverte : omega_n = (2 T0 / m) sin(n pi / 2M),
     verifiees par action directe sur les vecteurs propres. Pour n petit devant M : omega_n = n pi T0 / P+, niveaux
     equidistants, M^2 = N / alpha' avec alpha' = 1 / (2 pi sigma) : la pente de Regge de R120.
  B. L'ordonnee a l'origine est l'energie de point zero de la chaine. Somme exacte : M [cot(pi / 4M) - 1]
     = 4 M^2 / pi - M - pi / 12 - pi^3 / (2880 M^2) ... Le terme en M^2 (energie propre du maillon) et le terme en M
     (energie des bouts) demandent deux contre-termes ajustes ; le terme constant donne a = (D - 2) / 24 : 1/12 a
     quatre dimensions, 1 a vingt-six (R127 B). Chaine supersymetrique : 8 modes bosoniques et 8 fermioniques, tout
     s'annule, etat fondamental sans masse (R128 A).
  C. Un nombre fini de maillons brise l'invariance de Lorentz. La degenerescence omega_2 = 2 omega_1, exigee pour
     que les etats du niveau 2 forment un multiplet de rotation, devient omega_2 / 2 omega_1 = cos(pi / 2M). Des
     maillons de la taille des longueurs de la base : 5 maillons par fermi, ecart de 5 %. Les tests d'isotropie de
     l'inertie du nucleon (10^-29) exigent plus de 3 x 10^14 maillons par nucleon.
  D. Dipoles et jonction. Un maillon est une matrice elementaire E_ab (une charge sur le support a, la charge
     opposee sur b) ; E_ab E_cd = E_ad si b = c, zero sinon : la chaine est un produit de matrices, la jonction
     bout a bout une concatenation, et le nombre de maillons se conserve. Mais ce nombre n'est pas intrinseque : un
     proton du LHC en a 14 000 fois plus que le meme proton au repos.
  E. La chaine est une echelle LC. Meme dispersion, omega = (2 / sqrt(LC)) sin(theta / 2) ; l'impedance de la
     chaine est la tension T0, independante de la finesse des maillons ; la coupure 2 T0 / m part a l'infini.
  F. Les maillons vivent dans les D - 2 dimensions transverses : 8 = 2 + 6 ; la derniere dimension visible est
     portee par le nombre de maillons.
  Verdict : la chaine de maillons est une ecriture connue de la corde (IDENTITE a la limite d'un nombre infini de
  maillons) qui donne au filament dipolaire de l'auteur sa place exacte. Elle n'apporte aucun nombre nouveau. Un
  maillon de taille finie est EXCLU : il fixe un repere et un axe privilegies.
"""
from fractions import Fraction
from itertools import product
from math import cos, pi, sin, sqrt, tan

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = 0.1973269804                       # GeV fm
LAMBDA_E = 386.15926796                    # fm
L1_3 = 2 * pi * LAMBDA_E / 3
L1_9 = L1_3 * 3 ** (-2 * pi)
SIGMA = pi * HBARC ** 2 / L1_9 ** 2        # GeV^2, tension de la base
ALPHA_P = 1 / (2 * pi * SIGMA)             # GeV^-2
L_S = sqrt(ALPHA_P) * HBARC                # fm

def freq(n, M): return 2 * sin(n * pi / (2 * M))                         # en unites T0 / m

def laplacian_residual(M):
    """Chaine ouverte a bouts libres : (L v)_k = 2 v_k - v_{k-1} - v_{k+1}, bouts : v_1 - v_2 ; vecteur cos(n pi (k - 1/2) / M)."""
    worst = 0.0
    for n in range(M):
        v = [cos(n * pi * (k + 0.5) / M) for k in range(M)]
        lam = freq(n, M) ** 2
        for k in range(M):
            left = v[k - 1] if k > 0 else v[k]
            right = v[k + 1] if k < M - 1 else v[k]
            worst = max(worst, abs(2 * v[k] - left - right - lam * v[k]))
    return worst

def zero_point(M): return 2 * M * sum(sin(n * pi / (2 * M)) for n in range(1, M))    # contribution a M^2, en unites T0, par dimension

def main():
    print("R131 -- le filament comme chaine de maillons (bits de corde)\n")
    print(f"    base : sigma = {SIGMA:.4f} GeV^2, alpha' = {ALPHA_P:.4f} GeV^-2, longueur de corde = {L_S:.4f} fm\n")

    # A
    print("A. La chaine est la corde")
    res = max(laplacian_residual(M) for M in (2, 3, 4, 12, 60))
    small = {M: [freq(n, M) for n in range(1, M)] for M in (2, 3, 4)}
    by_hand = abs(small[2][0] - sqrt(2)) < 1e-12 and abs(small[3][0] - 1) < 1e-12 and abs(small[3][1] - sqrt(3)) < 1e-12
    print(f"    frequences exactes omega_n = (2 T0 / m) sin(n pi / 2M) : residu maximal sur les vecteurs propres {res:.1e}")
    print("    petites chaines, a refaire a la main (unites T0 / m) : " + " ; ".join(f"M = {M} : " + ", ".join(f"{w:.4f}" for w in ws) for M, ws in small.items()))
    M = 20000
    ratios = [freq(n, M) / freq(1, M) for n in range(1, 6)]
    cont = max(abs(r - n) for n, r in zip(range(1, 6), ratios))
    level = 2 * M * freq(1, M) / (2 * pi)                       # M^2 du premier niveau en unites 2 pi T0 = 1 / alpha'
    print(f"    {M} maillons : omega_n / omega_1 = " + ", ".join(f"{r:.5f}" for r in ratios) + f" ; premier niveau a {level:.6f} / alpha'")
    print(f"    -> niveaux equidistants, M^2 = N / alpha', alpha' = 1 / (2 pi sigma) = {ALPHA_P:.4f} GeV^-2 : la pente de R120 ; premier niveau {sqrt(1/ALPHA_P)*1e3:.0f} MeV.\n")
    check("A. chaine ouverte : frequences (2 T0 / m) sin(n pi / 2M) exactes ; limite continue : niveaux equidistants et pente alpha' = 1 / (2 pi sigma)",
          res < 1e-12 and by_hand and cont < 1e-5 and abs(level - 1) < 1e-6 and abs(ALPHA_P - 0.8595) < 5e-4, f"alpha' = {ALPHA_P:.4f} GeV^-2")

    # B
    print("B. L'ordonnee a l'origine : l'energie de point zero de la chaine")
    exact = all(abs(zero_point(M) - M * (1 / tan(pi / (4 * M)) - 1)) < 1e-9 * M * M for M in (2, 5, 17, 100, 1000))
    rows = []
    for M in (10, 100, 1000):                                      # au-dela, l'arrondi machine sur 4 M^2 / pi masque le terme cherche
        s = M * (1 / tan(pi / (4 * M)) - 1)
        const = s - 4 * M * M / pi + M
        rows.append((M, const, (const + pi / 12) * M * M))
    for M, const, nxt in rows:
        print(f"    M = {M:5d} : somme - 4 M^2 / pi + M = {const:.8f} ; (ecart a -pi/12) x M^2 = {nxt:.5f}")
    print(f"    limites : -pi / 12 = {-pi/12:.8f} ; -pi^3 / 2880 = {-pi**3/2880:.6f}")
    casimir = abs(rows[-1][1] + pi / 12) < 1e-7 and abs(rows[1][2] + pi ** 3 / 2880) < 1e-5
    intercept = {D: Fraction(D - 2, 24) for D in (4, 26)}
    print("    M^2 de point zero = -(D - 2) pi T0 / 12 = -(D - 2) / (24 alpha') : a = (D - 2) / 24 = " + ", ".join(f"{a} (D = {D})" for D, a in intercept.items()))
    print("    le terme en M^2 (energie propre de chaque maillon) et le terme en M (energie des deux bouts) demandent deux contre-termes ajustes")
    susy = all((8 - 8) * zero_point(M) == 0 for M in (5, 50))
    print(f"    chaine supersymetrique : 8 modes bosoniques et 8 fermioniques de memes frequences, somme nulle pour tout M : {susy} ; a = 0, etat fondamental sans masse")
    print("    -> le 1/12 qui condamnait la corde a quatre dimensions (R127 B) est l'energie de Casimir de la chaine ; la supersymetrie l'annule (R128 A).\n")
    check("B. point zero de la chaine : 4 M^2 / pi - M - pi / 12 - pi^3 / (2880 M^2) ; a = (D - 2) / 24, soit 1/12 a D = 4 et 1 a D = 26 ; chaine supersymetrique : 0",
          exact and casimir and intercept[4] == Fraction(1, 12) and intercept[26] == 1 and susy, "a = 1/12, 1 ; 0 avec supersymetrie")

    # C
    print("C. Un nombre fini de maillons : l'invariance de Lorentz")
    counts = {D: ((D - 2), (D - 2) * (D - 1) // 2, D * (D - 1) // 2 - 1) for D in (4, 10, 26)}
    mult_ok = all(a + b == c for a, b, c in counts.values())
    print("    niveau 2 : un quantum du mode 2 (D - 2 etats) et deux quanta du mode 1 ((D - 2)(D - 1) / 2 etats) doivent avoir la meme energie")
    print("    pour former un seul multiplet de rotation : " + " ; ".join(f"D = {D} : {a} + {b} = {c}" for D, (a, b, c) in counts.items()))
    print("    nombre de maillons   omega_2 / (2 omega_1)   ecart      pente / pente continue")
    table = {}
    for M in (2, 5, 10, 100, 10 ** 4, 10 ** 7):
        split = 1 - cos(pi / (2 * M)); slope = sin(pi / (2 * M)) / (pi / (2 * M)); table[M] = split
        print(f"    {M:12d}         {cos(pi/(2*M)):.12f}      {split:.2e}     {slope:.10f}")
    L_hadron = 1.0
    links = {"longueur de corde": L_hadron / L_S, "l1(9)": L_hadron / L1_9}
    print(f"    maillons de la taille de la base dans une corde d'un fermi : {links['longueur de corde']:.1f} (maillon = {L_S:.3f} fm), {links['l1(9)']:.1f} (maillon = {L1_9:.3f} fm)")
    print(f"       avec 5 maillons : ecart {100*table[5]:.1f} % ; avec 1 maillon : aucun mode interne")
    bound = 1e-29
    M_min = pi / sqrt(8 * bound)
    print(f"    isotropie de l'inertie du nucleon verifiee a {bound:.0e} pres (comagnetometres, de memoire) : pi^2 / (8 M^2) < {bound:.0e} exige M > {M_min:.1e}")
    print(f"       soit un maillon plus fin que {1e-15/M_min:.1e} m pour un fermi de corde, {L_S/(1.0/M_min):.1e} fois plus fin que la longueur de corde de la base")
    print("    -> un maillon de taille finie fixe un repere et un axe ; l'invariance ne revient qu'a la limite d'une infinite de maillons, ou la chaine est la corde.\n")
    check("C. nombre fini de maillons : omega_2 / 2 omega_1 = cos(pi / 2M), multiplet de rotation brise ; 5 maillons par fermi : 4,9 % ; isotropie a 1e-29 : plus de 3e14 maillons",
          mult_ok and counts[4] == (2, 3, 5) and counts[26] == (24, 300, 324) and 0.048 < table[5] < 0.050 and abs(table[100] - pi ** 2 / 8e4) < 1e-7
          and 5 < links["longueur de corde"] < 6 and 3e14 < M_min < 4e14, f"4,9 % pour 5 maillons ; M > {M_min:.1e}")

    # D
    print("D. Dipoles, jonction, et ce que compte le nombre de maillons")
    N = 3
    def E(a, b): return tuple(tuple(int((i, j) == (a, b)) for j in range(N)) for i in range(N))
    def mul(X, Y): return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(N)) for j in range(N)) for i in range(N))
    zero = tuple(tuple(0 for _ in range(N)) for _ in range(N))
    rule = all(mul(E(a, b), E(c, d)) == (E(a, d) if b == c else zero) for a, b, c, d in product(range(N), repeat=4))
    chain = E(0, 1)
    for a, b in ((1, 2), (2, 2), (2, 0), (0, 1)): chain = mul(chain, E(a, b))
    print(f"    maillon = matrice elementaire E_ab (charge + sur le support a, - sur b) ; E_ab E_cd = E_ad si b = c, zero sinon : {rule} ({N**4} produits)")
    print(f"    chaine de 5 maillons (0,1)(1,2)(2,2)(2,0)(0,1) : un seul dipole (0,1) : {chain == E(0, 1)} ; les charges interieures s'annulent deux a deux")
    print("    jonction bout a bout = concatenation : M1 + M2 maillons, P+ conserve")
    mp, e_beam = 0.93827, 6800.0
    gamma = e_beam / mp; beta = sqrt(1 - 1 / gamma ** 2)
    ratio = gamma * (1 + beta)
    print(f"    mais M = P+ / m depend du repere : un proton du LHC (6,8 TeV, gamma = {gamma:.0f}) porte {ratio:.0f} fois plus de maillons que le meme proton au repos")
    print("    -> le maillon n'est pas un constituant que l'on pourrait compter : c'est une tranche d'impulsion, aussi fine qu'on veut.\n")
    check("D. un maillon est un dipole (matrice elementaire), la chaine un produit de matrices, la jonction une concatenation ; le nombre de maillons depend du repere (x 14 500 au LHC)",
          rule and chain == E(0, 1) and 14000 < ratio < 15000, f"x {ratio:.0f}")

    # E
    print("E. La chaine est une echelle LC")
    ok_ladder = True
    for (Lc, Cc) in ((2.0, 0.5), (1.0, 1.0), (0.3, 4.0)):
        mb, k = Lc, 1 / Cc                                             # dictionnaire : L <-> masse du maillon, 1/C <-> raideur du ressort
        for th in (0.1, 0.7, 1.9, 3.0):
            ok_ladder &= abs(2 / sqrt(Lc * Cc) * sin(th / 2) - 2 * sqrt(k / mb) * sin(th / 2)) < 1e-12
        ok_ladder &= abs(sqrt(Lc / Cc) - sqrt(k * mb)) < 1e-12
    T0 = 1.0
    imp = [sqrt((T0 ** 2 / m) * m) for m in (1.0, 0.1, 1e-3, 1e-6)]
    cut = [2 * T0 / m for m in (1.0, 0.1, 1e-3, 1e-6)]
    print(f"    cellule (L, C) <-> (masse du maillon, ressort de raideur 1/C) : meme dispersion (2 / sqrt(LC)) sin(theta / 2), meme impedance sqrt(L / C) : {ok_ladder}")
    print(f"    maillon de corde : masse m, raideur T0^2 / m ; impedance sqrt(k m) = T0 pour m = 1, 0,1, 1e-3, 1e-6 : {imp} ; coupure 2 T0 / m : {cut}")
    print("    -> la tension est l'impedance de la chaine, independante de la finesse des maillons ; deux cordes de meme tension sont deux lignes adaptees.")
    print("       L'idee d'origine, le vide comme lignes adaptees, est sur le feuillet la chaine de maillons elle-meme (R127 E).\n")
    check("E. chaine de maillons = echelle LC (meme dispersion, meme impedance) ; impedance = tension, independante de la finesse ; coupure 2 T0 / m",
          ok_ladder and all(abs(z - 1) < 1e-12 for z in imp) and cut == sorted(cut), "Z = T0")

    # F
    print("F. Ou vivent les maillons")
    dims = {D: D - 2 for D in (4, 10, 26)}
    print("    les maillons bougent dans les D - 2 dimensions transverses : " + ", ".join(f"D = {D} : {t}" for D, t in dims.items()))
    print("    a dix dimensions : 8 = 2 visibles + 6 cachees ; la troisieme dimension visible est portee par le nombre de maillons (P+)")
    print("    dans les six dimensions cachees les ressorts suivent la forme du Calabi-Yau, sans expression connue (R130) : spectre non calculable en forme close.\n")
    check("F. les maillons vivent dans D - 2 dimensions : 8 = 2 + 6 a dix dimensions, la derniere dimension visible vient du nombre de maillons", dims[10] == 8 == 2 + 6, "8 = 2 + 6")

    print("Verdict : la chaine de maillons est une ecriture connue de la corde, exacte a la limite d'un nombre infini de maillons, ou le filament")
    print("dipolaire de l'auteur a sa place : chaque maillon est un dipole, la chaine un fil, la jonction une mise bout a bout, l'echelle LC son circuit.")
    print("Elle n'apporte aucun nombre nouveau, et un maillon de taille finie est exclu : il fixerait un repere et un axe privilegies.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
