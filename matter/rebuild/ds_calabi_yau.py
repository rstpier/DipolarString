#!/usr/bin/env python3
"""R130 -- La CD sur une variete de Calabi-Yau : ce que le choix donne, ce qu'il impose, ce qu'il laisse libre.

Decision de l'auteur (22 novembre 2026), devant R129 : « ok prend Calabi-Yau ».
Definition retenue, la plus simple qui garde les deux groupes caches de l'auteur et le dipole comme corde ouverte :
supercorde de type IIA sur une variete de Calabi-Yau, vue comme une famille de 3-tores (le groupe torique) au-dessus
d'une base de trois dimensions (le groupe de base), avec des D6-branes enroulees sur des 3-cycles. Exemple de
travail : la quintique, la plus simple des varietes de Calabi-Yau.

  A. Ce que le choix donne. Six dimensions d'un bloc : 105 composantes de courbure pour 21 equations, donc courbe
     sans contenu, sans rien couter a l'espace visible (R129 F). L'holonomie SU(3) laisse invariant un spineur sur
     quatre : un quart des supersymetries (8 charges sur 32, 4 avec des branes), le minimum, celui qui permet une
     matiere chirale. Le tore les gardait toutes.
  B. Les entiers de l'espace. Nombre de formes h21 par comptage des deformations, nombre d'Euler par les classes de
     Chern, nombre de tailles h11 par difference : dix varietes classiques retrouvees, dont la quintique (1, 101),
     chi = -200, 204 cycles de dimension trois. Variete de Tian et Yau (14, 23), chi = -18 ; son quotient d'ordre 3
     est libre (verifie : aucun point fixe), d'ou chi = -6.
  C. Les generations avec des cordes ouvertes. Deux branes se coupent un nombre de fois egal a leur produit
     antisymetrique dans le reseau des cycles de dimension trois. Sur le 6-tore c'est la formule de R128 (un
     determinant 6 x 6, reseau de rang 20, unimodulaire) ; sur un Calabi-Yau, rang 2 + 2 h21. Tout entier est
     disponible : 3 se choisit. Bilan de charge (Gauss) sur chaque cycle : deux familles de branes seules ont un
     nombre d'intersection nul ; l'exemple a trois generations de R128 ne ferme donc pas seul son bilan.
  D. Le prix. Chaque parametre de taille ou de forme est un champ sans masse. 6-tore : 70 scalaires (controle : la
     supergravite N = 8), dont 22 de taille, de forme et de couplage ; quintique : 410, dont 204 ; quotient de Tian
     et Yau : 52, dont 25. Sans mecanisme supplementaire chacun porte une force de portee infinie ; la gravite
     verifiee jusqu'a 52 um exige une masse superieure a 3,8 meV.
  E. Les deux groupes. Un 3-tore et une 3-sphere ont un nombre d'Euler nul : une famille lisse de tores aussi. Tout
     le nombre d'Euler vient des endroits ou le tore se pince. Retourner le groupe torique (dualite T, l'echange
     L <-> C de R127 E) est la symetrie miroir : tailles et formes s'echangent, chi change de signe.
  F. Le Calabi-Yau le plus proche de R128 : les trois tores replies par leur symetrie d'ordre 3. (36, 0), chi = 72,
     courbure concentree aux 27 = 3^3 points fixes, aucun parametre de forme, deux cycles de dimension trois
     seulement. Mais la symetrie envoie la direction capacitive de chaque tore sur l'inductive : les deux groupes de
     l'auteur n'en font plus qu'un. Trois generations : 34 %, 16 %, 2 % des paires selon la plage, effet de plage.
  G. Le 3 d'Euler. La regle « generations = |chi| / 2 » est celle de la corde heterotique, qui n'a pas de corde
     ouverte, donc pas de dipole. Echantillon de 2362 configurations : chi entre -200 et 0, jamais -6.
  Hors calcul : la metrique d'un Calabi-Yau ferme existe (Yau 1977) mais n'a pas d'expression connue ; les masses,
  qui en dependent, ne se calculent pas en forme close. Nombre de varietes connues : 7890 intersections completes
  (Candelas, Dale, Lutken, Schimmrigk 1988), 473 800 776 polytopes reflexifs et 30 108 couples (h11, h21) (Kreuzer
  et Skarke 2000), de memoire.
  Verdict : choix coherent et sans cout pour l'espace visible (DERIVE). Aucun nombre mesure n'en sort : l'espace se
  choisit, le 3 des generations est un entier choisi, les masses ne sont pas calculables, et le choix ajoute des
  dizaines a des centaines de champs sans masse a fixer. PARAMETRISATION, plus large que R128.
"""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, gcd
import random

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# ---------------------------------------------------------------- intersections completes : Euler et deformations
def euler_cicy(dims, degrees):
    """chi = integrale de c3(X) ; X intersection complete dans un produit d'espaces projectifs. Series tronquees au degre 3."""
    m = len(dims)
    zero = tuple([0] * m)
    unit = lambda r: tuple(int(i == r) for i in range(m))
    def mul(p, q, cap):
        out = {}
        for e1, v1 in p.items():
            s1 = sum(e1)
            for e2, v2 in q.items():
                if cap is not None and s1 + sum(e2) > cap: continue
                e = tuple(x + y for x, y in zip(e1, e2))
                if any(x > n for x, n in zip(e, dims)): continue
                out[e] = out.get(e, 0) + v1 * v2
        return out
    c = {zero: 1}
    for r, n in enumerate(dims):
        h = {zero: 1, unit(r): 1}
        for _ in range(n + 1): c = mul(c, h, 3)
    top = {zero: 1}
    for q in degrees:
        L = {unit(r): q[r] for r in range(m) if q[r]}
        inv = {zero: 1}; term = {zero: 1}
        for _ in range(3):
            term = mul(term, {e: -v for e, v in L.items()}, 3)
            for e, v in term.items(): inv[e] = inv.get(e, 0) + v
        c = mul(c, inv, 3); top = mul(top, L, None)
    c1_zero = all(v == 0 for e, v in c.items() if sum(e) == 1)
    chi = sum(v * top.get(tuple(n - x for n, x in zip(dims, e)), 0) for e, v in c.items() if sum(e) == 3)
    return chi, c1_zero

def h21_count(dims, degrees):
    """Deformations de forme : monomes des equations, moins leurs redefinitions, moins les changements de coordonnees."""
    def nmono(q):
        if any(x < 0 for x in q): return 0
        r = 1
        for n, x in zip(dims, q): r *= comb(n + x, n)
        return r
    tot = sum(nmono(q) for q in degrees)
    red = sum(nmono(tuple(x - y for x, y in zip(qa, qb))) for qa in degrees for qb in degrees)
    aut = sum((n + 1) ** 2 - 1 for n in dims)
    return tot - red - aut

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest

def compositions(total, k):
    if k == 1: yield (total,); return
    for i in range(total + 1):
        for rest in compositions(total - i, k - 1): yield (i,) + rest

def sample(max_dim):
    """Configurations irreductibles (a permutation pres) dans des produits d'espaces projectifs de dimension totale <= max_dim."""
    seen = {}
    for total in range(4, max_dim + 1):
        K = total - 3
        for dims in partitions(total):
            groups = {}
            for n in dims: groups[n] = groups.get(n, 0) + 1
            ns = sorted(groups, reverse=True)
            dd = [n for n in ns for _ in range(groups[n])]
            choices = [list(combinations_with_replacement(list(compositions(n + 1, K)), groups[n])) for n in ns]
            for pick in product(*choices):
                rows = [row for grp in pick for row in grp]
                cols = list(zip(*rows))
                if any(sum(col) <= 1 for col in cols): continue            # ni equation vide ni equation lineaire dans un seul facteur
                ok = True
                for size in range(1, K + 1):                               # tout sous-ensemble d'equations laisse au moins une dimension
                    for E in combinations(range(K), size):
                        supp = [r for r in range(len(dd)) if any(rows[r][a] for a in E)]
                        if sum(dd[r] for r in supp) < size + 1: ok = False; break
                    if not ok: break
                if not ok: continue
                comp = {0}; frontier = [0]; eqs = set()                     # pas de produit de deux espaces
                while frontier:
                    r = frontier.pop()
                    for a in range(K):
                        if rows[r][a] and a not in eqs:
                            eqs.add(a)
                            for r2 in range(len(dd)):
                                if rows[r2][a] and r2 not in comp: comp.add(r2); frontier.append(r2)
                if len(comp) != len(dd) or len(eqs) != K: continue
                best = None
                for perm in permutations(range(K)):
                    key = []; i = 0
                    for n in ns:
                        key.append(tuple(sorted(tuple(rows[i + j][a] for a in perm) for j in range(groups[n])))); i += groups[n]
                    key = tuple(key)
                    if best is None or key < best: best = key
                best = (tuple(dd), best)
                if best not in seen: seen[best] = euler_cicy(dd, [list(col) for col in cols])[0]
    return seen

# ---------------------------------------------------------------- algebre lineaire exacte
def det_exact(M):
    M = [[Fraction(v) for v in row] for row in M]; n = len(M); d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fraction(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f: M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d

def rank_exact(M):
    M = [[Fraction(v) for v in row] for row in M]; r = 0
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r

# ---------------------------------------------------------------- 3-cycles du 6-tore
TRIPLES = list(combinations(range(6), 3))
def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s
PAIRING = [[(perm_sign(I + J) if not set(I) & set(J) else 0) for J in TRIPLES] for I in TRIPLES]

def one_cycles(w):
    out = []
    for i, (n, m) in enumerate(w):
        v = [0] * 6; v[2 * i] = n; v[2 * i + 1] = m; out.append(v)
    return out

def cycle_class(w):
    """Classe d'une brane factorisable dans le reseau des 20 cycles de dimension trois : les mineurs 3 x 3."""
    a = one_cycles(w)
    return [int(det_exact([[a[r][c] for c in I] for r in range(3)])) for I in TRIPLES]

def pair_classes(A, B): return sum(A[i] * PAIRING[i][j] * B[j] for i in range(20) for j in range(20))

def torus_formula(a, b):
    p = 1
    for (n1, m1), (n2, m2) in zip(a, b): p *= n1 * m2 - m1 * n2
    return p

# ---------------------------------------------------------------- orbifold d'ordre 3 des trois tores hexagonaux
def theta(u): n, m = u; return (-m, n - m)                    # rotation de 120 degres dans la base (1, omega)
def wedge(u, v): return u[0] * v[1] - u[1] * v[0]
def inter_orbifold(a, b):
    tot = 0; bb = list(b)
    for _ in range(3):
        p = 1
        for u, v in zip(a, bb): p *= wedge(u, v)
        tot += p; bb = [theta(v) for v in bb]
    return tot

def orbifold_distribution(rmax):
    pr = [(n, m) for n in range(-rmax, rmax + 1) for m in range(-rmax, rmax + 1) if (n, m) != (0, 0) and gcd(abs(n), abs(m)) == 1]
    per = {}
    for u in pr:
        for v in pr:
            k = (wedge(u, v), wedge(u, theta(v))); per[k] = per.get(k, 0) + 1
    two = {}
    for (x1, y1), c1 in per.items():
        for (x2, y2), c2 in per.items():
            k = (x1 * x2, y1 * y2, (x1 + y1) * (x2 + y2)); two[k] = two.get(k, 0) + c1 * c2
    cnt = {}
    for (X, Y, Z), c12 in two.items():
        for (x3, y3), c3 in per.items():
            v = abs(X * x3 + Y * y3 - Z * (x3 + y3)); cnt[v] = cnt.get(v, 0) + c12 * c3
    tot = sum(v for k, v in cnt.items() if k)
    order = sorted(((v, k) for k, v in cnt.items() if k), reverse=True)
    return cnt.get(3, 0) / tot, [k for v, k in order].index(3) + 1

def main():
    print("R130 -- la CD sur une variete de Calabi-Yau : ce que le choix donne, impose et laisse libre\n")

    # A
    print("A. Ce que le choix donne")
    riem = lambda n: n * n * (n * n - 1) // 12
    ricc = lambda n: n * (n + 1) // 2
    spinor = [w for w in product((1, -1), repeat=3) if w[0] * w[1] * w[2] == 1]            # poids du spineur de Spin(6), en demi-unites
    roots = [(1, -1, 0), (0, 1, -1)]                                                        # racines simples de SU(3)
    inv = [w for w in spinor if all(sum(a * b for a, b in zip(w, r)) == 0 for r in roots)]
    q_cy, q_torus = 32 * len(inv) // len(spinor), 32
    print(f"    six dimensions d'un bloc : {riem(6)} composantes de courbure pour {ricc(6)} equations ; {riem(6)-ricc(6)} restent libres : courbe sans contenu")
    print(f"    spineur a six dimensions : {len(spinor)} etats ; invariants sous SU(3) : {len(inv)} ; sous l'holonomie triviale du tore : {len(spinor)}")
    print(f"    supersymetries conservees : Calabi-Yau {q_cy} sur 32 (N = {q_cy//4} a quatre dimensions), {q_cy//2} avec des branes (N = {q_cy//8}) ;")
    print(f"                               6-tore {q_torus} sur 32 (N = {q_torus//4}), {q_torus//2} avec des branes (N = {q_torus//8})")
    print("    -> le Calabi-Yau garde le minimum, celui qui permet des fermions chiraux (N = 1 au plus) ; le tore gardait tout.\n")
    check("A. Calabi-Yau : courbe sans contenu (105 contre 21) ; un spineur invariant sur quatre, soit 8 supersymetries sur 32, 4 avec des branes",
          riem(6) == 105 and ricc(6) == 21 and len(spinor) == 4 and len(inv) == 1 and q_cy == 8 and q_cy // 8 == 1 and q_torus // 4 == 8, "1/4")

    # B
    print("B. Les entiers de l'espace")
    known = {"quintique": ([4], [[5]], (1, 101)), "P5[4,2]": ([5], [[4], [2]], (1, 89)), "P5[3,3]": ([5], [[3], [3]], (1, 73)),
             "P6[3,2,2]": ([6], [[3], [2], [2]], (1, 73)), "P7[2,2,2,2]": ([7], [[2]] * 4, (1, 65)), "bicubique": ([2, 2], [[3, 3]], (2, 83)),
             "P1xP3[2,4]": ([1, 3], [[2, 4]], (2, 86)), "tetraquadrique": ([1, 1, 1, 1], [[2, 2, 2, 2]], (4, 68)),
             "Tian-Yau": ([3, 3], [[3, 0], [0, 3], [1, 1]], (14, 23)), "Schoen": ([1, 2, 2], [[1, 3, 0], [1, 0, 3]], (19, 19))}
    hodge = {}; all_match = True
    print("    variete            chi (Chern)   h21 (deformations)   h11 = h21 + chi/2   b3 = 2 + 2 h21   publie")
    for name, (d, q, lit) in known.items():
        chi, cy = euler_cicy(d, q); h21 = h21_count(d, q); h11 = h21 + chi // 2
        hodge[name] = (h11, h21, chi); all_match &= cy and (h11, h21) == lit
        print(f"    {name:16s} {chi:8d} {h21:16d} {h11:19d} {2 + 2*h21:16d}      {lit}")
    # quotient d'ordre 3 de la variete de Tian et Yau : x -> (x0, a^2 x1, a x2, a x3), y -> (y0, a y1, a^2 y2, a^2 y3)
    px, py = (0, 2, 1, 1), (0, 1, 2, 2)
    invariant = all((3 * p) % 3 == 0 for p in px + py) and all((a + b) % 3 == 0 for a, b in zip(px, py))
    def zmul(u, v): return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0] - u[1] * v[1])        # entiers d'Eisenstein a + b omega
    cube_roots = [(1, 0), (0, 1), (-1, -1)]
    assert all(zmul(zmul(z, z), z) == (1, 0) for z in cube_roots)
    # points fixes : x et y dans un espace propre. Vecteur propre isole : x_i^3 = 1, pas sur la variete. Reste x = (0,0,a,b), y = (0,0,c,d)
    # avec b = -z a, d = -z' c (z, z' racines cubiques de 1) ; la troisieme equation vaut a c (1 + z z').
    free = all(tuple(p + q for p, q in zip((1, 0), zmul(z, zp))) != (0, 0) for z in cube_roots for zp in cube_roots)
    chi_ty = hodge["Tian-Yau"][2]
    print(f"    quotient d'ordre 3 de Tian-Yau : equations invariantes : {invariant} ; aucun point fixe sur la variete : {free} ; chi = {chi_ty} / 3 = {chi_ty//3}")
    print(f"    nombres publies du quotient : (6, 9), soit chi = {2*(6-9)}")
    print("    -> la quintique : 1 taille, 101 formes, 204 cycles de dimension trois, chi = -200. Ce sont des entiers ; ils se calculent a la main.\n")
    check("B. dix varietes classiques : h21 par comptage et chi par classes de Chern redonnent les couples publies ; quintique (1, 101), chi = -200 ; quotient de Tian-Yau libre, chi = -6",
          all_match and hodge["quintique"] == (1, 101, -200) and invariant and free and chi_ty // 3 == -6 == 2 * (6 - 9), "(1, 101), -200 ; -18 / 3 = -6")

    # C
    print("C. Les generations avec des cordes ouvertes")
    rng = random.Random(130)
    w = lambda r: tuple((rng.randint(-r, r), rng.randint(-r, r)) for _ in range(3))
    pairs = [(w(4), w(4)) for _ in range(200)]
    det_ok = all(det_exact(one_cycles(a) + one_cycles(b)) == -torus_formula(a, b) for a, b in pairs)
    lattice_ok = all(pair_classes(cycle_class(a), cycle_class(b)) == -torus_formula(a, b) for a, b in pairs[:25])
    antisym = all(PAIRING[i][j] == -PAIRING[j][i] for i in range(20) for j in range(20))
    unimod = det_exact(PAIRING)
    print(f"    6-tore : nombre d'intersection = determinant 6 x 6 des six cycles = (au signe pres) produit de R128 : {det_ok} sur {len(pairs)} paires")
    print(f"    reseau des cycles de dimension trois : rang {len(TRIPLES)}, accouplement antisymetrique : {antisym}, determinant {unimod} ; meme resultat dans le reseau : {lattice_ok}")
    b3 = {"6-tore": 20, "quintique": 2 + 2 * 101, "quotient de Tian-Yau": 2 + 2 * 9, "trois tores replies": 2 + 2 * 0}
    print("    sur un Calabi-Yau : meme accouplement, rang b3 = 2 + 2 h21 : " + ", ".join(f"{k} {v}" for k, v in b3.items()))
    every = all(pair_classes([1] + [0] * 19, [0] * 19 + [k]) == k for k in range(0, 7))
    print(f"    tout entier est disponible (un cycle contre k fois son cycle dual : k) : {every} -> le 3 des generations se choisit")
    A0, B0 = cycle_class(((1, 0), (1, 0), (1, 0))), cycle_class(((1, 3), (1, 1), (1, 1)))
    r128 = pair_classes(A0, B0)
    two_stack = all(pair_classes(cycle_class(a), [-k * x for x in cycle_class(a)]) == 0 for (a, _), k in zip(pairs[:25], range(1, 26)))
    print(f"    bilan de charge sur chaque cycle (Gauss, espace ferme) : somme des classes des branes = 0. Deux familles seules : classes opposees, intersection nulle : {two_stack}")
    print(f"    exemple de R128, (1,0)(1,0)(1,0) et (1,3)(1,1)(1,1) : |I| = {abs(r128)}, classes independantes (rang {rank_exact([A0, B0])}) : il ne ferme pas seul son bilan")
    stacks = [(rng.randint(1, 5), cycle_class(w(2))) for _ in range(5)]
    total = [sum(N * c[i] for N, c in stacks) for i in range(20)]
    stacks.append((1, [-t for t in total]))
    anomaly = all(sum(N * pair_classes(ca, cb) for N, cb in stacks) == 0 for _, ca in stacks)
    print(f"    bilan ferme -> pour chaque brane, somme des fermions chiraux comptes avec leur signe = 0 (pas d'anomalie) : {anomaly}\n")
    check("C. generations = produit antisymetrique dans le reseau des 3-cycles (rang 20 sur le tore, 204 sur la quintique), tout entier disponible ; deux familles de branes seules : zero ; l'exemple de R128 ne ferme pas seul",
          det_ok and lattice_ok and antisym and unimod == 1 and every and abs(r128) == 3 and two_stack and rank_exact([A0, B0]) == 2 and anomaly and b3["quintique"] == 204,
          "rang 20, 204, 20, 2 ; |I| = 3 choisi")

    # D
    print("D. Le prix : chaque parametre de taille ou de forme est un champ sans masse")
    t6_scalars = 21 + comb(6, 2) + 1 + 1 + 6 + comb(6, 3) + 6            # metrique, B, dilaton, dual de B, C1, C3, duaux de C3
    t6_vectors = 6 + 6 + 1 + comb(6, 2)
    def cy_scalars(h11, h21):
        by_fields = (h11 + 2 * h21) + h11 + 1 + 1 + (2 + 2 * h21)      # metrique, B, dilaton, dual de B, C3
        by_multiplets = 2 * h11 + 4 * (h21 + 1)
        assert by_fields == by_multiplets
        return by_fields, h11 + 2 * h21 + 1
    table = {"6-tore": (t6_scalars, 21 + 1), "quintique (1, 101)": cy_scalars(1, 101), "quotient de Tian-Yau (6, 9)": cy_scalars(6, 9),
             "trois tores replies (36, 0)": cy_scalars(36, 0)}
    for k, (s, geo) in table.items():
        print(f"    {k:30s} : {s:4d} scalaires sans masse, dont {geo:4d} de taille, de forme et de couplage")
    hbarc = 1.973269804e-7                                               # eV m
    m_min = hbarc / 52e-6
    print(f"    controle : 6-tore, {t6_scalars} scalaires et {t6_vectors} vecteurs, le contenu de la supergravite N = 8")
    print(f"    sans mecanisme supplementaire, chacun porte une force de portee infinie ; gravite verifiee jusqu'a 52 um : masse > {m_min*1e3:.1f} meV")
    print("    -> le choix de l'espace ne fixe ni ses tailles ni ses formes ; les fixer demande des flux, dont les choix entiers font le « paysage ».\n")
    check("D. champs sans masse : 6-tore 70 (dont 22 geometriques), quintique 410 (204), quotient de Tian-Yau 52 (25), trois tores replies 76 (37) ; masse requise > 3,8 meV",
          t6_scalars == 70 and t6_vectors == 28 and table["quintique (1, 101)"] == (410, 204) and table["quotient de Tian-Yau (6, 9)"] == (52, 25)
          and table["trois tores replies (36, 0)"] == (76, 37) and 3.7e-3 < m_min < 3.9e-3, "70, 410, 52, 76")

    # E
    print("E. Les deux groupes caches dans un Calabi-Yau")
    chi_b = lambda b: sum((-1) ** i * v for i, v in enumerate(b))
    t3, s3 = [1, 3, 3, 1], [1, 0, 0, 1]
    print(f"    3-tore : nombres de Betti {t3}, chi = {chi_b(t3)} ; 3-sphere : {s3}, chi = {chi_b(s3)} ; famille lisse de tores sur une base : chi = {chi_b(t3)} x chi(base) = 0")
    print("    quintique : chi = -200, donc le groupe torique se pince, et tout le nombre d'Euler loge aux pincements")
    mirror_ok = True
    for name in ("quintique", "bicubique", "Tian-Yau", "Schoen"):
        h11, h21, chi = hodge[name]
        even = 1 + h11 + h11 + 1
        mirror_ok &= 2 * (h21 - h11) == -chi and 2 + 2 * h11 == even
    print(f"    retourner le groupe torique (symetrie miroir) : (h11, h21) -> (h21, h11), chi -> -chi, cycles de dimension trois <-> cycles de dimension paire : {mirror_ok}")
    print("       quintique (1, 101) -> (101, 1) : 204 cycles de dimension trois deviennent 4 ; branes de dimension paire <-> impaire (R129 B)")
    print("    -> les deux groupes ne sont plus interchangeables : une base, et un tore dont la forme varie et se pince le long de la base.\n")
    check("E. chi(3-tore) = chi(3-sphere) = 0 : le nombre d'Euler vient des pincements ; symetrie miroir : (h11, h21) echanges, chi -> -chi",
          chi_b(t3) == 0 and chi_b(s3) == 0 and mirror_ok and hodge["quintique"][2] == -200, "0 ; -200 <-> +200")

    # F
    print("F. Le Calabi-Yau le plus proche de R128 : les trois tores replies par leur symetrie d'ordre 3")
    fixed_per_torus = len({(Fraction(a, 3), Fraction(b, 3)) for a in range(3) for b in range(3)
                           if ((-b - a) % 3 == 0 and (a - 2 * b) % 3 == 0)})          # (theta - 1) z dans le reseau, z = (a, b) / 3
    nfix = fixed_per_torus ** 3
    chi_orb = Fraction(0 + 8 * nfix, 3)                                                # (1/|G|) somme sur les paires (g, h) des chi des points fixes communs
    form = lambda p, q: comb(3, p) * comb(3, q) if (p - q) % 3 == 0 else 0             # formes dz^p dzbar^q invariantes
    h11_orb, h21_orb = form(1, 1) + nfix, form(2, 1)
    proper = form(1, 0) == 0 and form(2, 0) == 0 and form(3, 0) == 1
    print(f"    points fixes : {fixed_per_torus} par tore, {nfix} = 3^3 en tout ; nombre d'Euler (formule d'orbifold) : {chi_orb}")
    print(f"    formes invariantes : h11 = {form(1,1)} + {nfix} = {h11_orb}, h21 = {h21_orb} ; chi = 2 (h11 - h21) = {2*(h11_orb-h21_orb)} ; Calabi-Yau propre : {proper}")
    print("    -> la courbure est concentree aux 27 points ; aucun parametre de forme (la symetrie d'ordre 3 impose le reseau hexagonal).")
    pure_c, pure_l = ((1, 0),) * 3, ((0, 1),) * 3
    prim = [(n, m) for n in range(-2, 3) for m in range(-2, 3) if (n, m) != (0, 0) and gcd(abs(n), abs(m)) == 1]
    branes = [tuple(rng.choice(prim) for _ in range(3)) for _ in range(30)]
    orb_rank = rank_exact([[inter_orbifold(a, b) for b in branes] for a in branes])
    print(f"    cycles de dimension trois : rang de l'accouplement des orbites de branes = {orb_rank} = 2 + 2 h21")
    print(f"    la rotation envoie la direction capacitive (1,0) sur l'inductive {theta((1, 0))} : brane « tout capacitif » et brane « tout inductif » sont")
    print(f"       la meme brane, nombre d'intersection {inter_orbifold(pure_c, pure_l)} -> sur cet espace les deux groupes caches de l'auteur n'en font qu'un")
    shares = [orbifold_distribution(r) for r in (1, 2, 3)]
    print("    trois generations parmi les paires de branes qui se coupent, selon la plage d'enroulements (1, 2, 3) : "
          + ", ".join(f"{100*s:.1f} % (rang {k})" for s, k in shares))
    print("    -> 3 est frequent aux petits enroulements puis disparait : effet de plage, comme en R128 ; il se choisit.\n")
    check("F. trois tores replies : 27 points fixes, (36, 0), chi = 72 de deux facons, rang 2 ; capacitif et inductif identifies ; part de |I| = 3 : 34 %, 16 %, 2 % selon la plage",
          fixed_per_torus == 3 and nfix == 27 and chi_orb == 72 == 2 * (h11_orb - h21_orb) and (h11_orb, h21_orb) == (36, 0) and proper and orb_rank == 2
          and theta((1, 0)) == (0, 1) and inter_orbifold(pure_c, pure_l) == 0
          and 0.33 < shares[0][0] < 0.36 and 0.15 < shares[1][0] < 0.17 and shares[2][0] < 0.03 and shares[1][1] == 1 and shares[2][1] > 3,
          f"(36, 0), chi = 72 ; {100*shares[0][0]:.0f} %, {100*shares[1][0]:.0f} %, {100*shares[2][0]:.0f} %")

    # G
    print("G. Le 3 d'Euler")
    print("    regle « generations = |chi| / 2 » : corde heterotique, qui n'a pas de corde ouverte, donc pas de dipole ; lecture donnee pour memoire")
    print(f"    quintique {abs(hodge['quintique'][2])//2}, bicubique {abs(hodge['bicubique'][2])//2}, trois tores replies {int(chi_orb)//2}, quotient de Tian-Yau {abs(chi_ty//3)//2}")
    s = sample(7)
    vals = sorted(set(s.values()))
    near = sorted(v for v in vals if -12 <= v < 0)
    mult6 = sum(1 for v in s.values() if v and v % 6 == 0)
    print(f"    echantillon : {len(s)} configurations irreductibles (au plus 4 equations) ; {len(vals)} valeurs de chi, de {min(vals)} a {max(vals)} ; chi = -6 : {sum(1 for v in s.values() if v == -6)} fois")
    print(f"    valeurs non nulles les plus proches de zero : {near} ; multiples de 6 non nuls : {mult6}, qui demanderaient un quotient libre d'ordre |chi| / 6")
    print("    -> aucune ne donne trois generations directement ; il faut une symetrie libre du bon ordre, comme celle d'ordre 3 de Tian et Yau.\n")
    check("G. echantillon de 2362 configurations : chi entre -200 et 0, jamais -6 ; le 3 d'Euler demande un quotient choisi",
          len(s) == 2362 and min(vals) == -200 and max(vals) == 0 and -6 not in vals and len(vals) == 66 and abs(chi_ty // 3) // 2 == 3 and int(chi_orb) // 2 == 36,
          f"{len(s)} configurations, {len(vals)} valeurs, -6 absent")

    print("Verdict : le Calabi-Yau est un choix coherent et sans cout pour l'espace visible, qui garde deux groupes caches (une base, un tore).")
    print("Il ne donne aucun nombre mesure : l'espace se choisit, le 3 des generations est un entier choisi, les masses dependent d'une metrique")
    print("que l'on ne sait pas ecrire, et le choix ajoute des dizaines a des centaines de champs sans masse a fixer. Parametrisation, plus large que R128.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
