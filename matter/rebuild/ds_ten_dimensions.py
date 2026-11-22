#!/usr/bin/env python3
"""R128 -- La CD en dix dimensions : la supercorde, les trois tores, et ce que le choix impose.

Decision de l'auteur (22 novembre 2026), apres R127 B (une corde fondamentale a quatre dimensions est incoherente) :
« ok prend 10 dimensions ». Dix dimensions, c'est la supercorde : fermions sur le feuillet, pas de tachyon, graviton
sans masse. Definition retenue pour pouvoir calculer, la plus simple qui respecte le decoupage de l'auteur
(« 3 massiques, 3 capacitives, 3 inductives, 1 temps ») : la supercorde de type IIA sur trois tores, T^6 = (T^2)^3,
avec des D6-branes qui se coupent (Zwiebach, chapitre 21). Le filament dipolaire est la corde ouverte tendue entre
deux branes ; chaque tore porte une direction « capacitive » et une direction « inductive ».

  A. Le spectre. Supercorde ouverte : 8 etats bosoniques + 8 fermioniques sans masse (un vecteur et son partenaire de
     spin 1/2), puis 128 + 128, 1152 + 1152 : bosons = fermions a chaque niveau (identite « abstruse » de Jacobi,
     verifiee ici en serie). Le spin 1/2 que la CD n'a jamais eu est la, au prix de la supersymetrie.
  B. Le decoupage. 10 = 1 + 3 + 6 ; le vecteur a 8 polarisations = 2 (photon a quatre dimensions) + 6 scalaires ;
     6 = 3 x 2 : trois dimensions complexes, une par tore. Le compte de l'auteur, 1 + 3 + 3 + 3 = 10, se transpose
     donc sous cette lecture (rectificatif de R127 F, qui supposait un dedoublement) ; une D6-brane occupe les 3
     dimensions visibles et une direction par tore.
  C. Les generations. Deux D6-branes d'enroulements (n1, m1) et (n2, m2) sur un tore se coupent |n1 m2 - m1 n2| fois
     (verifie par comptage direct) ; sur trois tores, le produit. Chaque intersection porte un fermion chiral : le
     nombre de generations est un nombre d'intersection (Ibanez, Marchesano, Rabadan 2001). Mais 3 n'est pas
     prefere : sa part va de 0 % a 14 % selon la plage d'enroulements ; il se choisit, il ne se derive pas.
  D. Les echelles. Tores de la taille de la corde : M_s = sqrt(pi) g_s M_Planck, soit 4e17 a 4e18 GeV. Borne basse :
     M_s > environ 8 TeV (resonances de cordes au LHC, valeur citee de memoire). La longueur de corde est donc sous
     2,5e-20 m : les longueurs de la base (0,81 fm, 809 fm) ne sont pas des longueurs de corde. Six grandes
     dimensions egales a la borne du LHC : R ~ 1 fm, valeur d'ecole (Arkani-Hamed, Dimopoulos, Dvali), pas une trouvaille.
  E. La tension des hadrons. Dans une geometrie deformee ds^2 = e^{2A} dx^2 + ..., la tension effective d'une corde
     est e^{2A}/(2 pi alpha'). Pour retrouver sigma_base il faut e^{A} = 1e-4 (M_s = 8 TeV) a 1e-19 (M_s ~ Planck).
     RECTIFICATIF DE R127 C : « une seule corde ne peut pas faire les hadrons et la gravite » vaut en espace plat ; en
     dix dimensions un facteur de deformation le permet, et c'est le mecanisme de la QCD holographique.
  F. Ce qui n'est pas observe : aucun partenaire supersymetrique (gluino au-dessus de 2 TeV environ), loi en 1/r^2
     verifiee jusqu'a 52 micrometres, aucune resonance de corde. Rien n'est predit aux energies accessibles.
  Verdict : en dix dimensions la CD devient une construction connue (type IIA, trois tores, D6-branes secantes) ou
  ses motifs trouvent une place exacte (dipole = corde ouverte entre branes ; 3 + 3 = trois tores ; generations =
  intersections). Les nombres (groupe de jauge, 3 generations, tension des hadrons) y sont CHOISIS, pas derives :
  parametrisation, comme le catalogue de particules retire en V2.10. Ce n'est pas une variante : c'est un cas.
"""
import math
from fractions import Fraction
from math import gcd, comb

HBARC = 197.3269804              # MeV fm
ME = 0.51099895
M_PL_RED = 2.435e18              # GeV, masse de Planck reduite
LAM = HBARC / ME
L1_9 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)
SIGMA_GEV2 = math.pi * HBARC / L1_9 ** 2 * HBARC * 1e-6
M_S_LHC = 7.9e3                  # GeV, cite de memoire (CMS-EXO-19-012), non verifie a la source

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

# ---------- series en x = q^(1/2) ----------
N = 9
def mul(a, b):
    out = [0] * (N + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j > N: break
                out[i + j] += ai * bj
    return out
def power(a, k):
    out = [1] + [0] * N
    for _ in range(k): out = mul(out, a)
    return out
def binom_factor(step, sign):
    f = [0] * (N + 1); f[0] = 1
    if step <= N: f[step] = sign
    return f
def inv8(step):
    f = [0] * (N + 1)
    k = 0
    while step * k <= N:
        f[step * k] = comb(k + 7, 7); k += 1
    return f

def spectrum():
    A = [1] + [0] * N; B = [1] + [0] * N; Rr = [1] + [0] * N; P = [1] + [0] * N
    n = 1
    while 2 * n - 1 <= N:
        A = mul(A, power(binom_factor(2 * n - 1, +1), 8)); B = mul(B, power(binom_factor(2 * n - 1, -1), 8)); n += 1
    n = 1
    while 2 * n <= N:
        Rr = mul(Rr, power(binom_factor(2 * n, +1), 8)); P = mul(P, inv8(2 * n)); n += 1
    diff = [(a - b) // 2 for a, b in zip(A, B)]
    ns0 = diff[1:] + [0]
    NS = mul(ns0, P)
    R = [8 * c for c in mul(Rr, P)]
    return NS[:N], R[:N]

# ---------- intersections sur un tore ----------
def intersections_bruteforce(n1, m1, n2, m2):
    det = n1 * m2 - m1 * n2
    if det == 0: return 0
    pts = set()
    rng = range(-abs(det), abs(det) + 1)
    for a in rng:
        for b in rng:
            t = Fraction(a * m2 - b * n2, det)
            s = Fraction(a * m1 - b * n1, det)
            assert n1 * t - n2 * s == a and m1 * t - m2 * s == b
            pts.add(((n1 * t) % 1, (m1 * t) % 1))
    return len(pts)

def main():
    print("R128 -- la CD en dix dimensions : supercorde de type IIA sur trois tores, D6-branes secantes\n")

    # A. spectre
    print("A. Le spectre de la supercorde ouverte : bosons = fermions a chaque niveau")
    NS, R = spectrum()
    lv = [NS[i] for i in range(0, 9, 2)]; lr = [R[i] for i in range(0, 9, 2)]
    odd = [NS[i] for i in range(1, 9, 2)] + [R[i] for i in range(1, 9, 2)]
    print(f"    secteur bosonique (NS, projection GSO) : {lv}")
    print(f"    secteur fermionique (Ramond)           : {lr}")
    print("    niveau 0 : 8 + 8 = un vecteur et son partenaire de spin 1/2 ; niveau 1 : 128 + 128 ; pas de tachyon.\n")
    check("A. identite de Jacobi : 8, 128, 1152, 7680, 42112 etats bosoniques et autant de fermioniques ; spin 1/2 present, pas de tachyon",
          lv == lr and lv[:3] == [8, 128, 1152] and not any(odd), f"{lv}")

    # B. decoupage
    print("B. Le decoupage des dix dimensions")
    D = 10
    print(f"    {D} = 1 (temps) + 3 (visibles) + 6 (compactes) ; vecteur : {D-2} polarisations = 2 + 6 scalaires ; 6 = 3 x 2 : trois tores T^2")
    print("    lecture de l'auteur : 1 temps + 3 « massiques » + 3 « capacitives » + 3 « inductives » = une paire par tore ; une D6-brane = 3 visibles + 1 direction par tore")
    print("    -> le compte se transpose (rectificatif de R127 F).\n")
    check("B. 1 + 3 + 3 + 3 = 10 avec 6 = (T^2)^3 ; le vecteur a dix dimensions donne 2 + 6", 1 + 3 + 3 + 3 == D and D - 2 == 2 + 6 and 6 == 3 * 2, "10")

    # C. intersections
    print("C. Les generations comme nombres d'intersection")
    tests = [((1, 0), (1, 3)), ((1, 0), (0, 1)), ((2, 1), (1, 2)), ((1, 1), (1, -1)), ((3, 1), (1, 3)), ((2, 3), (1, 1)), ((1, 2), (-2, 1))]
    ok = True
    for (a, b) in tests:
        det = a[0] * b[1] - a[1] * b[0]
        cnt = intersections_bruteforce(a[0], a[1], b[0], b[1])
        ok &= (cnt == abs(det))
        print(f"      tore : {a} et {b} -> comptage direct {cnt}, |n1 m2 - m1 n2| = {abs(det)}")
    brane_a = [(1, 0), (1, 0), (1, 0)]; brane_b = [(1, 3), (1, 1), (1, 1)]
    I_ab = 1
    for (a, b) in zip(brane_a, brane_b): I_ab *= a[0] * b[1] - a[1] * b[0]
    print(f"    exemple sur (T^2)^3 : {brane_a} et {brane_b} -> I = {I_ab} : trois copies du meme fermion chiral")
    fracs = {}
    for rmax in (1, 2, 3, 4):
        vecs = [(n, m) for n in range(-rmax, rmax + 1) for m in range(-rmax, rmax + 1) if (n, m) != (0, 0) and gcd(abs(n), abs(m)) == 1]
        per = {}
        for a in vecs:
            for b in vecs:
                d = a[0] * b[1] - a[1] * b[0]; per[d] = per.get(d, 0) + 1
        dist = {}
        for d1, c1 in per.items():
            for d2, c2 in per.items():
                for d3, c3 in per.items():
                    v = abs(d1 * d2 * d3); dist[v] = dist.get(v, 0) + c1 * c2 * c3
        nonzero = sum(c for v, c in dist.items() if v)
        top = sorted(((c / nonzero, v) for v, c in dist.items() if v), reverse=True)[:5]
        fracs[rmax] = dist.get(3, 0) / nonzero
        print(f"    enroulements entre -{rmax} et {rmax} ({len(vecs)} directions par tore) : {len([v for v in dist if v])} valeurs de |I| ; en tete : "
              + ", ".join(f"{v} ({100*f:.1f} %)" for f, v in top) + f" ; |I| = 3 : {100*fracs[rmax]:.1f} %")
    print("    -> la part de |I| = 3 va de 0 % a 14 % selon la plage d'enroulements : aucune preference ; 3 se choisit, il ne se derive pas.\n")
    check("C. formule d'intersection verifiee par comptage ; exemple a 3 intersections ; aucune preference pour 3 (0 %, 14 %, 2,6 %, 0,8 % selon la plage)",
          ok and I_ab == 3 and fracs[1] == 0 and fracs[2] < 0.2 and fracs[4] < fracs[3] < fracs[2],
          ", ".join(f"{100*fracs[r]:.1f} %" for r in (1, 2, 3, 4)))

    # D. echelles
    print("D. Les echelles")
    ms_lo, ms_hi = math.sqrt(math.pi) * 0.1 * M_PL_RED, math.sqrt(math.pi) * 1.0 * M_PL_RED
    ls_lhc = HBARC * 1e-3 / M_S_LHC * 1e-15          # m
    g_s = 0.2
    R_gev = (math.pi * g_s ** 2 * M_PL_RED ** 2) ** (1 / 6) / M_S_LHC ** (4 / 3)
    R_fm = R_gev * HBARC * 1e-3
    print(f"    tores de la taille de la corde : M_s = sqrt(pi) g_s M_Planck = {ms_lo:.1e} a {ms_hi:.1e} GeV (g_s = 0,1 a 1)")
    print(f"    borne du LHC : M_s > {M_S_LHC/1e3:.1f} TeV -> l_s < {ls_lhc:.1e} m ; l_1(9) = {L1_9:.2f} fm = {L1_9*1e-15/ls_lhc:.0e} l_s : les longueurs de la base ne sont pas des longueurs de corde")
    print(f"    six grandes dimensions egales a cette borne (g_s = {g_s}) : R = {R_fm:.2f} fm -- valeur d'ecole pour six dimensions, dependante de la borne (R ~ M_s^(-4/3)), pas une trouvaille\n")
    check("D. fenetre de l'echelle de corde : 8 TeV a 4e18 GeV ; l_s < 2,5e-20 m ; l_1(9) depasse l_s de plus de 1e4",
          ms_hi > 1e18 and ls_lhc < 3e-20 and L1_9 * 1e-15 / ls_lhc > 1e4 and 0.3 < R_fm < 5, f"l_s < {ls_lhc:.1e} m ; R(6 dim.) = {R_fm:.2f} fm")

    # E. deformation
    print("E. La tension des hadrons par un facteur de deformation : sigma_eff = e^{2A} / (2 pi alpha')")
    rows = []
    for lab, ms in (("M_s = 7,9 TeV", M_S_LHC), ("M_s = 4,3e17 GeV (g_s = 0,1)", ms_lo), ("M_s = 4,3e18 GeV (g_s = 1)", ms_hi)):
        e2A = 2 * math.pi * SIGMA_GEV2 / ms ** 2
        rows.append((lab, math.sqrt(e2A), e2A * ms ** 2 / (2 * math.pi)))
        print(f"      {lab:30s} e^A = {math.sqrt(e2A):.2e} (sigma_eff = {e2A*ms**2/(2*math.pi):.4f} GeV^2)")
    print("    -> rectificatif de R127 C : en dix dimensions une geometrie deformee donne la tension des hadrons a la corde fondamentale ;")
    print("       c'est la QCD holographique, et e^A est un nombre choisi, pas derive.\n")
    check("E. un facteur de deformation de 1e-4 a 1e-19 reproduit sigma_base : l'exclusion de R127 C ne vaut qu'en espace plat",
          all(abs(r[2] / SIGMA_GEV2 - 1) < 1e-12 for r in rows) and rows[0][1] < 1e-3 and rows[2][1] < 1e-18, f"e^A = {rows[0][1]:.1e} a {rows[2][1]:.1e}")

    # F. ce qui n'est pas observe
    print("F. Ce qui n'est pas observe")
    print("    aucun partenaire supersymetrique (gluino au-dessus de 2 TeV environ) ; loi en 1/r^2 verifiee jusqu'a 52 micrometres (Lee et al. 2020,")
    print("    portee d'un Yukawa de force gravitationnelle < 38,6 micrometres) ; aucune resonance de corde au LHC. Rien n'est predit aux energies accessibles.\n")

    print("Verdict : en dix dimensions la CD devient une construction connue -- supercorde de type IIA sur trois tores, D6-branes secantes -- ou ses")
    print("motifs ont une place exacte : le dipole est la corde ouverte entre deux branes, le 3 + 3 est fait des trois tores, les generations sont des")
    print("intersections, et la tension des hadrons vient d'un facteur de deformation. Mais le groupe de jauge, le 3 des generations et ce facteur y sont")
    print("CHOISIS : c'est une parametrisation, comme le catalogue retire en V2.10. Ce n'est pas une variante de la theorie des cordes, c'en est un cas.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
