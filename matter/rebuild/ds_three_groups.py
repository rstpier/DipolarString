#!/usr/bin/env python3
"""R129 -- Neuf dimensions d'espace en trois groupes de trois, et « tores, non probablement pas ».

Precisions de l'auteur (22 novembre 2026) : « 9 dim : 3 groupes de 3 dimensions », soit le groupe visible (x, y, z),
le groupe « capacitif » (c1, c2, c3) et le groupe « inductif » (l1, l2, l3) ; puis, a la question « est-ce le meme
espace que les trois tores de R128 ? » : « tores, non probablement pas ». On calcule donc les deux lectures.

Premiere lecture : chaque groupe cache est ferme et plat.
  A. Meme espace. T^6 = (T^2)^3 = T^3_C x T^3_L : les paires (c_i, l_i) et les triplets (c1 c2 c3), (l1 l2 l3)
     partagent les memes six coordonnees.
  B. Objets qui remplissent des groupes entiers : p = 3k, soit D3, D6, D9. La dualite T sur un groupe entier echange
     les types IIA et IIB et, dans ce groupe, impulsion et enroulement.
  C. Generations. Deux branes alignees sur les groupes se coupent 0 ou 1 fois, jamais 3. Trois generations exigent un
     biais (enroulement d'au moins 2 melant une direction capacitive et une inductive) ou une symetrie d'ordre 3.
  D. Symetrie d'ordre 3 : 3 points fixes par tore hexagonal, donc 3, 9, 27 = 3^n.
  E. Meme groupe : I = 0, paires non chirales ; groupes differents : I non nul, fermion chiral.

Seconde lecture : les groupes ne sont pas des tores. Que peuvent-ils etre ?
  F. En trois dimensions la courbure est fixee entierement par le contenu : le tenseur de Riemann (6 composantes) est
     une fonction du tenseur de Ricci (6 composantes), identite verifiee ici exactement. Un groupe vide (ni flux ni
     gradient) est donc plat : un 3-tore ou l'un de ses neuf replis (dix varietes plates fermees, Bieberbach 1911 ;
     Hantzsche et Wendt 1935), les replis permis etant d'ordre 1, 2, 3, 4, 6. A partir de quatre dimensions ce n'est
     plus vrai (20 composantes contre 10 ; une onde gravitationnelle est courbe sans contenu) : un espace de six
     dimensions pris d'un bloc peut etre courbe pour rien, un produit de deux groupes de trois ne le peut pas.
  G. Un groupe traverse par un flux (champ de Kalb-Ramond) : courbure positive constante, donc une 3-sphere ou l'un
     de ses quotients, de rayon R^2 = k alpha' (k entier). Mais deux spheres cachees ne ferment pas le bilan : charge
     centrale 15 - 6 (1/k1 + 1/k2) < 15. Le manque se paie dans l'espace visible : soit un couplage qui varie d'un
     facteur e tous les 1/Q le long d'une direction d'espace, avec Q^2 = 1/R1^2 + 1/R2^2 (equations de la
     supergravite et bilan exact concordent ; un gradient dans le temps a le mauvais signe) ; soit une courbure
     negative de meme rayon, 1/L^2 = 1/R1^2 + 1/R2^2 : c'est la solution exacte connue a trois groupes de trois,
     AdS3 x S^3 x S^3 x S^1 (Elitzur, Feinerman, Giveon, Tsabar 1998), dont le groupe visible n'a que deux dimensions
     d'espace. Le compte exact de l'auteur avec trois spheres, temps + S^3 x S^3 x S^3, n'a pas de solution statique.
     Enonce general (secteur de Neveu-Schwarz, espace cache ferme, espace visible plat et constantes constantes) :
     l'equation du dilaton s'ecrit Laplacien(e^{-2 Phi}) = e^{-2 Phi} H^2 / 6, dont l'integrale est nulle : H = 0,
     donc espace cache sans courbure de Ricci, donc, pour un produit de groupes de trois, plat.
  H. Des groupes infinis : la gravite se dilue dans neuf dimensions, force en 1/r^8, sans orbite stable (sauf
     piegeage de la gravite sur notre groupe, mecanisme supplementaire non examine ici).
  I. L'issue sans cout : renoncer au produit. Une variete de Calabi-Yau vue comme une famille de 3-tores (un groupe)
     au-dessus d'une base de trois dimensions (l'autre groupe) : conjecture de Strominger, Yau et Zaslow (1996), ou
     la dualite T sur le groupe torique est la symetrie miroir (le B ci-dessus). Les generations y sont un nombre
     topologique : |chi|/2 dans la lecture heterotique. Quintique : chi = -200 ; variete de Tian et Yau : chi = -18,
     et son quotient libre d'ordre 3 : chi = -6, trois generations. S^3 x S^3 : b2 = 0, pas de Calabi-Yau.
  Verdict : « trois groupes de trois » en produit avec un espace visible plat impose des tores ou leurs replis
  (DERIVE, theoreme classique). « Pas des tores » : spheres tenues par un flux EXCLUES comme dimensions cachees,
  groupes infinis EXCLUS ; reste le Calabi-Yau, ou les deux groupes caches sont entrelaces et ou l'un reste torique :
  CONDITIONNEL, et le 3 des generations y est encore un nombre choisi avec l'espace. Parametrisation, comme R128.
  Hors champ : flux de Ramond et plans d'orientifold (compactifications a flux), non examines.
"""
from fractions import Fraction
from itertools import product
from math import comb, cos, gcd, pi, sqrt

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# ---------------------------------------------------------------- premiere lecture : branes et points fixes
def inter(a, b):
    I = 1
    for (n1, m1), (n2, m2) in zip(a, b): I *= n1 * m2 - m1 * n2
    return I

def z3_fixed_points():
    """Points fixes de la rotation d'ordre 3 sur le tore hexagonal : (theta - 1) z dans le reseau."""
    th = ((0, -1), (1, -1))                                   # theta dans la base (1, omega)
    a, b, c, d = th[0][0] - 1, th[0][1], th[1][0], th[1][1] - 1
    det = a * d - b * c
    pts = set()
    for k1 in range(-abs(det), abs(det) + 1):
        for k2 in range(-abs(det), abs(det) + 1):
            z1 = Fraction(d * k1 - b * k2, det); z2 = Fraction(-c * k1 + a * k2, det)
            assert a * z1 + b * z2 == k1 and c * z1 + d * z2 == k2
            pts.add((z1 % 1, z2 % 1))
    return len(pts), abs(det)

# ---------------------------------------------------------------- courbure exacte d'une metrique polynomiale
def P(terms):
    """Polynome : {exposants: coefficient}."""
    return {e: Fraction(c) for e, c in terms.items() if c != 0}

def pdiff(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e); f[i] -= 1
            out[tuple(f)] = out.get(tuple(f), 0) + c * e[i]
    return out

def peval(p, pt):
    s = Fraction(0)
    for e, c in p.items():
        t = c
        for x, k in zip(pt, e): t *= x ** k
        s += t
    return s

def inverse(m):
    n = len(m)
    a = [list(row) + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(m)]
    for col in range(n):
        piv = next(r for r in range(col, n) if a[r][col] != 0)
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]; a[col] = [v / pv for v in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                f = a[r][col]; a[r] = [v - f * w for v, w in zip(a[r], a[col])]
    return [row[n:] for row in a]

def curvature(metric, pt):
    """Riemann R_abcd, Ricci R_bd = R^a_bad, scalaire et metrique en un point rationnel : arithmetique exacte."""
    n = len(metric); rg = range(n)
    g = [[peval(metric[a][b], pt) for b in rg] for a in rg]
    dg = [[[peval(pdiff(metric[a][b], e), pt) for b in rg] for a in rg] for e in rg]
    ddg = [[[[peval(pdiff(pdiff(metric[a][b], e), f), pt) for b in rg] for a in rg] for f in rg] for e in rg]
    gi = inverse(g)
    c1 = [[[(dg[b][a][c] + dg[c][a][b] - dg[a][b][c]) / 2 for c in rg] for b in rg] for a in rg]          # [a; bc]
    dc1 = [[[[(ddg[e][b][a][c] + ddg[e][c][a][b] - ddg[e][a][b][c]) / 2 for c in rg] for b in rg] for a in rg] for e in rg]
    dgi = [[[-sum(gi[a][p] * dg[e][p][q] * gi[q][d] for p in rg for q in rg) for d in rg] for a in rg] for e in rg]
    G = [[[sum(gi[a][d] * c1[d][b][c] for d in rg) for c in rg] for b in rg] for a in rg]
    dG = [[[[sum(dgi[e][a][d] * c1[d][b][c] + gi[a][d] * dc1[e][d][b][c] for d in rg) for c in rg] for b in rg] for a in rg] for e in rg]
    Rup = [[[[dG[c][a][d][b] - dG[d][a][c][b] + sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in rg)
              for d in rg] for c in rg] for b in rg] for a in rg]
    R = [[[[sum(g[a][e] * Rup[e][b][c][d] for e in rg) for d in rg] for c in rg] for b in rg] for a in rg]
    Ric = [[sum(Rup[a][b][a][d] for a in rg) for d in rg] for b in rg]
    S = sum(gi[b][d] * Ric[b][d] for b in rg for d in rg)
    return g, R, Ric, S

def weyl(g, R, Ric, S):
    n = len(g); rg = range(n)
    return [[[[R[a][b][c][d] - (g[a][c] * Ric[b][d] - g[a][d] * Ric[b][c] - g[b][c] * Ric[a][d] + g[b][d] * Ric[a][c]) / (n - 2)
               + S * (g[a][c] * g[b][d] - g[a][d] * g[b][c]) / ((n - 1) * (n - 2)) for d in rg] for c in rg] for b in rg] for a in rg]

def flat4(t): return [v for a in t for b in a for c in b for v in c]

# ---------------------------------------------------------------- charges centrales (supercorde, total 15)
def c_flat(d): return Fraction(3 * d, 2)
def c_su2(k): return Fraction(9, 2) - Fraction(6) / k          # 3-sphere avec k quanta de flux
def c_sl2(k): return Fraction(9, 2) + Fraction(6) / k          # anti-de Sitter a trois dimensions

# ---------------------------------------------------------------- nombre d'Euler d'une intersection complete
def euler_cicy(dims, degrees):
    """chi = integrale de c3(X) sur X, X intersection complete dans un produit d'espaces projectifs."""
    m = len(dims)
    def mul(p, q):
        out = {}
        for e1, v1 in p.items():
            for e2, v2 in q.items():
                e = tuple(x + y for x, y in zip(e1, e2))
                if all(x <= n for x, n in zip(e, dims)): out[e] = out.get(e, 0) + v1 * v2
        return out
    one = {tuple([0] * m): 1}
    def lin(q): return {tuple(int(i == r) for i in range(m)): q[r] for r in range(m) if q[r]}
    c = dict(one)
    for r, n in enumerate(dims):
        h = dict(one); h[tuple(int(i == r) for i in range(m))] = 1
        for _ in range(n + 1): c = mul(c, h)
    top = dict(one)
    for q in degrees:
        L = lin(q); inv = dict(one); term = dict(one)
        for _ in range(sum(dims)):
            term = mul(term, {e: -v for e, v in L.items()})
            for e, v in term.items(): inv[e] = inv.get(e, 0) + v
        c = mul(c, inv); top = mul(top, L)
    c3 = {e: v for e, v in c.items() if sum(e) == 3}
    c1 = {e: v for e, v in c.items() if sum(e) == 1 and v}
    return mul(c3, top).get(tuple(dims), 0), not c1

def main():
    print("R129 -- neuf dimensions d'espace en trois groupes de trois ; « tores, non probablement pas »\n")
    print("=== Premiere lecture : chaque groupe cache est ferme et plat ===\n")

    # A
    print("A. Le meme espace, regroupe autrement")
    pairs = [("c1", "l1"), ("c2", "l2"), ("c3", "l3")]
    triples = [("c1", "c2", "c3"), ("l1", "l2", "l3")]
    s_pairs = sorted(x for p in pairs for x in p); s_triples = sorted(x for t in triples for x in t)
    print(f"    trois tores  : {pairs}")
    print(f"    deux 3-tores : {triples}")
    print("    memes six coordonnees ; avec le groupe visible (x, y, z) : 9 = 3 + 3 + 3 = 3 + 2 + 2 + 2\n")
    check("A. si les groupes sont fermes et plats : (T^2)^3 et T^3 x T^3 partagent les memes six dimensions",
          s_pairs == s_triples and 3 + 3 + 3 == 3 + 2 + 2 + 2 == 9, "9")

    # B
    print("B. Les branes qui remplissent des groupes entiers")
    group_branes = {k: 3 * k for k in (1, 2, 3)}
    iia = [p for p in range(0, 10) if p % 2 == 0]; iib = [p for p in range(0, 10) if p % 2 == 1]
    print(f"    k groupes remplis -> p = 3k : {group_branes} ; type IIA : p pair {iia} ; type IIB : p impair {iib}")
    print(f"    D6 (IIA) : visible + un groupe cache ; dualite T sur ce groupe : D{6-3} (IIB), sur l'autre groupe : D{6+3} (IIB)")
    print("    -> une dualite T sur un groupe entier echange les deux theories et, dans ce groupe, impulsion et enroulement.\n")
    check("B. D3, D6, D9 remplissent 1, 2, 3 groupes ; la dualite T sur un groupe (3 directions) echange IIA et IIB et envoie D6 sur D3 ou D9",
          group_branes == {1: 3, 2: 6, 3: 9} and 6 in iia and 3 in iib and 9 in iib and (6 - 3) % 2 != 6 % 2, "D3, D6, D9")

    # C
    print("C. Les generations pour des branes alignees sur les axes des groupes")
    axis = [tuple(choice) for choice in product([(1, 0), (0, 1)], repeat=3)]
    vals = sorted({abs(inter(a, b)) for a in axis for b in axis})
    n1 = sum(1 for a in axis for b in axis if abs(inter(a, b)) == 1)
    pure_c, pure_l = ((1, 0),) * 3, ((0, 1),) * 3
    print(f"    {len(axis)} branes alignees (direction c ou l sur chaque tore) : |I| prend les valeurs {vals} ; {n1} paires ordonnees sur {len(axis)**2} se coupent une fois")
    print(f"    groupe capacitif entier contre groupe inductif entier : I = {inter(pure_c, pure_l)} ; meme groupe : I = {inter(pure_c, pure_c)}")
    def frac3(rmax):
        vecs = [(n, m) for n in range(-rmax, rmax + 1) for m in range(-rmax, rmax + 1) if (n, m) != (0, 0) and gcd(abs(n), abs(m)) == 1]
        per = {}
        for a in vecs:
            for b in vecs:
                d = a[0] * b[1] - a[1] * b[0]; per[d] = per.get(d, 0) + 1
        hit = tot = 0
        for d1, c1 in per.items():
            for d2, c2 in per.items():
                for d3, c3 in per.items():
                    v = abs(d1 * d2 * d3)
                    if v: tot += c1 * c2 * c3
                    if v == 3: hit += c1 * c2 * c3
        return hit / tot
    tilted_a, tilted_b = ((1, 1), (1, 0), (1, 0)), ((1, -2), (0, 1), (0, 1))
    print(f"    enroulements limites a 1 en valeur absolue : part de |I| = 3 = {100*frac3(1):.0f} % ; il faut un enroulement d'au moins 2")
    print(f"    exemple de biais : {tilted_a} et {tilted_b} -> I = {inter(tilted_a, tilted_b)}")
    print("    -> alignees sur les groupes, les branes donnent une generation au plus ; trois exigent un biais entre une direction capacitive et une inductive.\n")
    check("C. branes alignees sur les groupes : |I| = 0 ou 1, jamais 3 ; avec des enroulements de 1 au plus, 3 est impossible ; un biais d'enroulement 2 donne 3",
          vals == [0, 1] and inter(pure_c, pure_l) in (1, -1) and inter(pure_c, pure_c) == 0 and frac3(1) == 0 and abs(inter(tilted_a, tilted_b)) == 3,
          f"|I| dans {vals} ; exemple en biais : {abs(inter(tilted_a, tilted_b))}")

    # D
    print("D. La symetrie d'ordre 3 et le « 3^n »")
    nfix, det = z3_fixed_points()
    print(f"    rotation de 120 degres sur un tore hexagonal : {nfix} points fixes (comptage direct ; |det(theta - 1)| = {det})")
    print(f"    un, deux, trois tores : {nfix}, {nfix**2}, {nfix**3} = 3^n points fixes")
    print("    -> le « 3^n » compte les points ou la matiere peut se loger dans l'orbifold Z_3, pas des cordes par particule.\n")
    check("D. orbifold Z_3 : 3 points fixes par tore, donc 3, 9, 27 = 3^n pour un, deux, trois tores", nfix == 3 == det and nfix ** 3 == 27, f"{nfix}, {nfix**2}, {nfix**3}")

    # E
    print("E. « Meme groupe » et annihilation")
    print("    branes paralleles (meme groupe, I = 0) : cordes en paires particule-antiparticule, annihilation et masse permises ;")
    print("    intersection (groupes differents, I non nul) : un fermion chiral, sans partenaire, sans masse permise.\n")
    check("E. meme groupe : I = 0 (matiere non chirale) ; groupes differents : I non nul (matiere chirale)", inter(pure_c, pure_c) == 0 and inter(pure_c, pure_l) != 0, "0 / non nul")

    print("=== Seconde lecture : « tores, non probablement pas » ===\n")

    # F
    print("F. Ce qu'un groupe de trois dimensions peut etre quand il est vide")
    riem = lambda n: n * n * (n * n - 1) // 12
    ricc = lambda n: n * (n + 1) // 2
    print("    composantes independantes (Riemann / Ricci) : " + " ; ".join(f"n = {n} : {riem(n)} / {ricc(n)}" for n in (3, 4, 6)))
    x = lambda i, j, k: (i, j, k)
    g3 = [[P({x(0,0,0): 2, x(2,0,0): 1, x(0,1,1): 1}), P({x(1,1,0): Fraction(1, 3)}), P({x(0,0,1): Fraction(1, 5)})],
          [P({x(1,1,0): Fraction(1, 3)}), P({x(0,0,0): 3, x(0,2,0): 1, x(1,0,1): 1}), P({x(1,0,0): Fraction(1, 7)})],
          [P({x(0,0,1): Fraction(1, 5)}), P({x(1,0,0): Fraction(1, 7)}), P({x(0,0,0): 4, x(0,0,2): 1, x(1,1,0): 1})]]
    pts = [(Fraction(1, 2), Fraction(1, 3), Fraction(-1, 4)), (Fraction(-2, 5), Fraction(3, 7), Fraction(1, 6)), (Fraction(1, 9), Fraction(-1, 2), Fraction(2, 3))]
    weyl_zero, curved = True, True
    for pt in pts:
        g, R, Ric, S = curvature(g3, pt)
        weyl_zero &= all(v == 0 for v in flat4(weyl(g, R, Ric, S)))
        curved &= any(v != 0 for v in flat4(R))
    print(f"    metrique quelconque a trois dimensions, {len(pts)} points rationnels : Riemann non nul : {curved} ; Riemann - (fonction de Ricci) = 0 exactement : {weyl_zero}")
    # 3-sphere unite au pole, coordonnees stereographiques tronquees a l'ordre 2 (la courbure au pole n'en demande pas plus)
    s3 = [[P({x(0,0,0): 4, x(2,0,0): -8, x(0,2,0): -8, x(0,0,2): -8}) if a == b else P({}) for b in range(3)] for a in range(3)]
    gs, Rs, Rics, Ss = curvature(s3, (Fraction(0),) * 3)
    sphere_ok = all(Rics[a][b] == 2 * gs[a][b] for a in range(3) for b in range(3)) and Ss == 6
    print(f"    3-sphere de rayon 1 : Ricci = 2 g : {sphere_ok} ; courbure scalaire = {Ss}")
    # quatre dimensions : onde gravitationnelle plane, ds^2 = -2 du dv + (x^2 - y^2) du^2 + dx^2 + dy^2
    y = lambda i, j, k, l: (i, j, k, l)
    z4 = P({})
    g4 = [[P({y(0,0,2,0): 1, y(0,0,0,2): -1}), P({y(0,0,0,0): -1}), z4, z4],
          [P({y(0,0,0,0): -1}), z4, z4, z4],
          [z4, z4, P({y(0,0,0,0): 1}), z4],
          [z4, z4, z4, P({y(0,0,0,0): 1})]]
    gw, Rw, Ricw, Sw = curvature(g4, (Fraction(1, 2), Fraction(1, 3), Fraction(1, 5), Fraction(1, 7)))
    wave_ricci_flat = all(v == 0 for row in Ricw for v in row)
    wave_curved = Rw[0][2][0][2] != 0
    print(f"    quatre dimensions, onde gravitationnelle plane : Ricci = 0 : {wave_ricci_flat} ; Riemann non nul : {wave_curved} (R_uxux = {Rw[0][2][0][2]})")
    orders = [n for n in range(1, 25) if abs(2 * cos(2 * pi / n) - round(2 * cos(2 * pi / n))) < 1e-12]
    print(f"    replis d'un reseau a trois dimensions (trace entiere, 1 + 2 cos(2 pi / n)) : ordres permis {orders}")
    print("    -> un groupe de trois dimensions vide est plat : un 3-tore ou l'un de ses neuf replis (dix varietes plates fermees,")
    print("       dont six orientables : holonomies 1, Z2, Z3, Z4, Z6, Z2 x Z2). Un espace de six dimensions d'un bloc echappe a la regle.\n")
    check("F. trois dimensions : Riemann est une fonction de Ricci (6 = 6, identite exacte), donc vide = plat = tore ou repli ; quatre dimensions et plus : faux (20 contre 10, onde plane)",
          riem(3) == ricc(3) == 6 and riem(4) == 20 and ricc(4) == 10 and riem(6) == 105 and weyl_zero and curved and sphere_ok and wave_ricci_flat and wave_curved and orders == [1, 2, 3, 4, 6],
          f"Riemann/Ricci : 6/6 a n = 3, 20/10 a n = 4, 105/21 a n = 6 ; replis {orders}")

    # G
    print("G. Un groupe traverse par un flux : la 3-sphere, et son prix")
    eps = {}
    for a, b, c in product(range(3), repeat=3):
        eps[(a, b, c)] = (a - b) * (b - c) * (c - a) // 2
    contr = [[sum(eps[(a, c, d)] * eps[(b, c, d)] for c in range(3) for d in range(3)) for b in range(3)] for a in range(3)]
    full = sum(v * v for v in eps.values())
    print(f"    flux H = h x (volume) : H_acd H_b^cd = {contr[0][0]} h^2 g_ab, H^2 = {full} h^2 ; Einstein : Ricci = (h^2 / 2) g, positif")
    print("    -> courbure positive constante : 3-sphere ou l'un de ses quotients, 2/R^2 = h^2/2, soit h^2 = 4/R^2 ; pas d'espace hyperbolique avec ce flux.")
    h2 = lambda R2: Fraction(4) / R2
    q2_sugra = lambda R2s: sum(Fraction(full) * h2(R2) for R2 in R2s) / 24            # equation du dilaton : Q^2 = H^2 / 24
    ns5 = all(c_flat(6) + Fraction(3, 2) + 6 * Fraction(1, k) + c_su2(k) == 15 for k in range(2, 60))
    deficit_ok = all(15 - (c_flat(4) + c_su2(k1) + c_su2(k2)) == 6 * (Fraction(1, k1) + Fraction(1, k2)) > 0 for k1 in range(2, 60) for k2 in range(2, 60))
    match = all(q2_sugra([Fraction(k1), Fraction(k2)]) == Fraction(1, k1) + Fraction(1, k2) for k1 in range(2, 60) for k2 in range(2, 60))
    ads = all(c_sl2(1 / (Fraction(1, k1) + Fraction(1, k2))) + c_su2(k1) + c_su2(k2) + c_flat(1) == 15 for k1 in range(2, 60) for k2 in range(2, 60))
    print(f"    etalonnage (gorge de la 5-brane, resultat d'ecole alpha' Q^2 = 1/k) : {ns5}")
    print(f"    espace visible plat + deux spheres : charge centrale 15 - 6 (1/k1 + 1/k2), manque strictement positif pour tout k : {deficit_ok}")
    print(f"    gradient de dilaton qui le comble : alpha' Q^2 = 1/k1 + 1/k2, soit Q^2 = 1/R1^2 + 1/R2^2 ; supergravite = bilan exact : {match}")
    print(f"    ou courbure negative : AdS3 x S^3 x S^3 x S^1 ferme le bilan avec 1/k = 1/k1 + 1/k2, soit 1/L^2 = 1/R1^2 + 1/R2^2 : {ads}")
    print("       compte 3 + 3 + 3 + 1 = 10, mais le temps est dans le premier groupe : deux dimensions d'espace visibles.")
    # le compte exact de l'auteur avec trois spheres : temps + S^3 x S^3 x S^3 ; gradient de dilaton V : c = D + 6 alpha' V.V
    need = lambda ks: (15 - (c_flat(1) + sum(c_su2(k) for k in ks))) / 6          # alpha' V.V requis
    static3 = all(need((k1, k2, k3)) == Fraction(1, k1) + Fraction(1, k2) + Fraction(1, k3) > 0
                  for k1 in range(2, 20) for k2 in range(2, 20) for k3 in range(2, 20))
    print(f"    temps + trois spheres (1 + 3 + 3 + 3) : il faudrait alpha' V.V = 1/k1 + 1/k2 + 1/k3 > 0, un gradient de genre espace : {static3} ;")
    print("       la seule direction ouverte est le temps (V.V < 0) : pas de solution statique.")
    R_hidden = 52e-6
    efold = R_hidden / sqrt(2)
    hubble = 1.37e26; dalpha = 1e-5; span = 1e26; omega_k = 0.005
    R_needed = sqrt(2) * span / dalpha
    L_curv = hubble / sqrt(omega_k)
    print(f"    spheres cachees sous la borne de la gravite a courte distance (R < 52 um) : couplage multiplie par e tous les {efold*1e6:.0f} um")
    print(f"    a l'inverse, des constantes stables a {dalpha:.0e} pres sur {span:.0e} m exigent R > {R_needed:.1e} m, soit {R_needed/hubble:.0e} rayons de Hubble")
    print(f"    et la variante courbe : rayon de courbure visible observe > {L_curv:.1e} m (|Omega_K| < 0,005), donc spheres plus grandes que l'univers observable")
    print("    -> des spheres cachees tenues par un flux : exclues. Le tore est plat, c'est pour cela qu'il ne coute rien.\n")
    check("G. flux -> 3-sphere (R^2 = k alpha') ; deux spheres + espace visible plat : manque 6 (1/k1 + 1/k2) ; Q^2 = 1/R1^2 + 1/R2^2 ou AdS de meme rayon : spheres cachees exclues",
          contr == [[2, 0, 0], [0, 2, 0], [0, 0, 2]] and full == 6 and ns5 and deficit_ok and match and ads and static3 and efold < 1e-4 and R_needed > 1e4 * hubble and L_curv > hubble,
          f"e-pli tous les {efold*1e6:.0f} um pour R = 52 um ; R > {R_needed:.0e} m sinon")

    # H
    print("H. Des groupes infinis")
    d = 9; p = d - 1
    print(f"    {d} dimensions d'espace ouvertes : les lignes de force se diluent sur une sphere de dimension {d-1}, force en 1/r^{p}")
    print(f"    doubler la distance divise la force par {2**p} (observe : 4) ; orbite circulaire stable si l'exposant est inferieur a 3 : 3 - {p} = {3-p}")
    print("    -> ni loi de Newton (verifiee de 52 um au systeme solaire) ni orbites : exclu, sauf piegeage de la gravite sur notre groupe,")
    print("       mecanisme supplementaire non examine ici.\n")
    check("H. trois groupes infinis : gravite en 1/r^8, pas d'orbite stable (exposant < 3 requis)", p == 8 and 2 ** p == 256 and 3 - p < 0, "1/r^8")

    # I
    print("I. L'issue sans cout : deux groupes entrelaces (Calabi-Yau)")
    chi_q, cy_q = euler_cicy([4], [[5]])
    chi_ty, cy_ty = euler_cicy([3, 3], [[3, 0], [0, 3], [1, 1]])
    chi_bc, cy_bc = euler_cicy([2, 2], [[3, 3]])
    poin_s3s3 = [1, 0, 0, 2, 0, 0, 1]
    poin_t6 = [comb(6, i) for i in range(7)]
    chi = lambda b: sum((-1) ** i * v for i, v in enumerate(b))
    print(f"    quintique : chi = {chi_q} (c1 = 0 : {cy_q}) -> {abs(chi_q)//2} generations ; bicubique : chi = {chi_bc} -> {abs(chi_bc)//2}")
    print(f"    variete de Tian et Yau : chi = {chi_ty} (c1 = 0 : {cy_ty}) ; quotient libre d'ordre 3 : chi = {chi_ty//3} -> {abs(chi_ty//3)//2} generations")
    print(f"    nombres de Betti : S^3 x S^3 {poin_s3s3} (b2 = 0 : pas de forme de Kahler, pas de Calabi-Yau, chi = {chi(poin_s3s3)}) ; T^6 {poin_t6} (chi = {chi(poin_t6)})")
    print("    -> le Calabi-Yau est courbe et ne coute rien, mais il n'est pas un produit : la forme du groupe torique varie le long de l'autre groupe")
    print("       (Strominger, Yau, Zaslow 1996), et la dualite T sur le groupe torique est la symetrie miroir. Le 3 des generations y est un nombre")
    print("       topologique de l'espace choisi (chi = +-6), pas un resultat.\n")
    check("I. Calabi-Yau : quintique chi = -200, Tian-Yau chi = -18 et quotient d'ordre 3 chi = -6 (trois generations) ; S^3 x S^3 a b2 = 0 et chi = 0",
          chi_q == -200 and cy_q and chi_ty == -18 and cy_ty and chi_ty // 3 == -6 and chi_bc == -162 and poin_s3s3[2] == 0 and chi(poin_s3s3) == 0 and chi(poin_t6) == 0,
          f"chi = {chi_q}, {chi_ty}, {chi_ty//3}")

    print("Verdict : en produit et avec un espace visible plat, trois groupes de trois imposent des tores ou leurs replis (theoreme classique).")
    print("« Pas des tores » : spheres tenues par un flux exclues comme dimensions cachees ; groupes infinis exclus ; reste le Calabi-Yau, ou les deux")
    print("groupes caches sont entrelaces, l'un restant torique. Le 3 des generations y est encore choisi avec l'espace : parametrisation, comme R128.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
