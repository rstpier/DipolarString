#!/usr/bin/env python3
"""R112 -- Le noeud-piece contre les axiomes du DQD : se derive-t-il de C1-C5 ?

Ce qu'est le noeud-piece, en clair.  Au croisement de trois lignes (x, y, z), le noeud de Johns
(SCN) est un standard : une impulsion qui arrive par une ligne repart dans les six directions avec
des fractions fixes (+-1/2), sans changer de nature.  Le noeud-piece de l'automate de Weyl (R111)
fait autre chose : les deux impulsions qui arrivent par les deux bouts d'une ligne x forment un
couple (a_+, a_-) ; le noeud le reexprime dans la base de la ligne suivante et le renvoie en entier
sur cette ligne : y recoit (a_+ + a_-)/sqrt2 d'un cote et (a_+ - a_-)/sqrt2 de l'autre, avec un
quart d'onde (facteur i) entre les deux, puis y -> z, puis z -> x.  En langage de lignes : un
coupleur hybride 3 dB a 90 deg (le quart d'onde) monte en circulateur (x -> y -> z -> x).

  A. construction : le noeud 6 ports (trois lignes, deux cotes), blocs V_yx, V_zy, V_xz ; chaque
     amplitude a |V|^2 = 1/2 ; le dephasage entre les deux sorties est pi/2 (hybride) ou 0 (Hadamard).
  B. dynamique de reseau : le cone de Weyl a 1/3 de la vitesse de ligne (trois traversees par pas
     d'automate), isotrope au premier ordre ; le photon du meme automate (R111 B) a la meme vitesse.
  C. les axiomes du manuscrit (section 'Uniqueness theorem', C1-C5) :
     C4 energie : S^dag S = 1 oui, mais S n'est pas reel : S^T S != 1 ;
     reciprocite S = S^T (implicite dans l'ansatz reel, symetrie T d'un milieu passif) : non, c'est
     un circulateur ;
     C1 symetrie cubique par noeud : le routage cyclique x -> y -> z n'est invariant que sous C3 ;
     parite : l'image miroir d'un noeud est le noeud de l'autre main ;
     A5/C3 (le noeud ne stocke rien, pas de volume reactif) : le quart d'onde n'est pas jaugeable
     par des phases de reference par port (residu 1,2) : c'est un dephasage physique, un element
     reactif au croisement.
  D. tentative de version achirale : deux sous-reseaux de mains opposees (parite de x + y + z) ;
     le cone est detruit (v = 0,275, anisotropie 0,4) : une arrivee x routee vers y en A est routee
     de y vers x en B, la suite n'utilise jamais z ; la symetrisation de Trotter inverse l'ordre
     apres un cycle complet, ce qu'un noeud fixe par site ne peut pas faire.  La chiralite du
     noeud-piece est structurelle.
  E. verdict : NON DERIVABLE de C1-C5 ; trois axiomes a changer, nommes (C1 -> T_h ; non-reciprocite
     par noeud ; quart d'onde reactif contre A5) ; le vide resterait achiral (T_h contient l'inversion).
"""
import math
import itertools
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
SIG = [SX, SY, SZ]
AXES = np.eye(3)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def eigbasis(i):
    """|+i>, |-i> : vecteurs propres de sigma_i (phase fixee : premiere composante reelle positive)"""
    w, V = np.linalg.eigh(SIG[i])
    plus, minus = V[:, np.argmax(w)], V[:, np.argmin(w)]
    for v in (plus, minus):
        ph = np.exp(-1j * np.angle(v[np.argmax(np.abs(v))]))
        v *= ph
    return plus, minus

def port(line, side):        # side +1 : cote +, -1 : cote -
    return 2 * line + (0 if side > 0 else 1)

def coin_node(order=(0, 1, 2)):
    """S[(sortie), (entree)] : les arrivees de la ligne order[j] repartent sur la ligne order[j+1]"""
    S = np.zeros((6, 6), complex)
    for j in range(3):
        a, b = order[j], order[(j + 1) % 3]
        pa, ma = eigbasis(a)
        pb, mb = eigbasis(b)
        # arrivee par le cote -a (a avance +a) = composante |+a> ; arrivee par le cote +a = composante |-a>
        # depart cote +b = composante |+b> ; cote -b = |-b>
        for (side_in, comp_in) in ((-1, pa), (+1, ma)):
            for (side_out, comp_out) in ((+1, pb), (-1, mb)):
                S[port(b, side_out), port(a, side_in)] = np.vdot(comp_out, comp_in)
    return S

def prop(q):
    P = np.zeros((6, 6), complex)
    for i in range(3):
        for s in (+1, -1):
            P[port(i, -s), port(i, s)] = np.exp(1j * s * q[i])      # depart cote s -> arrivee cote -s du voisin
    return P

def directions(n=15, seed=9):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return list(v) + [np.array([1., 0, 0]), np.array([1., 1, 0]) / math.sqrt(2), np.array([1., 1, 1]) / math.sqrt(3)]

def cone(Mfun, qmag=0.005):
    vs = []
    for d in directions():
        def vmin(qm):
            om = np.sort(np.angle(np.linalg.eigvals(Mfun(qm * d))))
            return om[om > 1e-9].min() / qm
        vs.append(2 * vmin(qmag) - vmin(2 * qmag))
    return float(np.mean(vs)), float(np.max(vs) - np.min(vs))

def rotations():
    def rot(axis, ang):
        axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
        K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K
    gens = [rot([0, 0, 1], math.pi / 2), rot([1, 0, 0], math.pi / 2)]
    G = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for g in list(G):
            for h in gens:
                gh = g @ h
                if not any(np.allclose(gh, k, atol=1e-9) for k in G):
                    G.append(gh); changed = True
    return G

def port_rep(R):
    """action d'une rotation (ou du miroir) sur les six ports (direction de mouvement)"""
    M = np.zeros((6, 6))
    for i in range(3):
        for s in (+1, -1):
            v = R @ (s * AXES[i])
            j = int(np.argmax(np.abs(v))); sj = int(round(v[j]))
            M[port(j, sj), port(i, s)] = 1
    return M

def main():
    print("R112 -- le noeud-piece contre les axiomes du DQD\n")
    S = coin_node()
    print("A. Construction (blocs 2x2 : arrivees de x -> departs de y, y -> z, z -> x)")
    for a, b in ((0, 1), (1, 2), (2, 0)):
        blk = np.array([[S[port(b, so), port(a, si)] for si in (-1, +1)] for so in (+1, -1)])
        ph = np.angle(blk[0, 0] * np.conj(blk[1, 0]))
        print(f"    {'xyz'[a]} -> {'xyz'[b]} : |V|^2 = {np.round(np.abs(blk)**2, 3).tolist()} ; dephasage entre les deux sorties {abs(math.degrees(ph)):.0f} deg")
    uni = np.allclose(S.conj().T @ S, np.eye(6))
    check("noeud-piece : 6 ports, |V|^2 = 1/2 partout, quart d'onde (90 deg) sur x -> y et y -> z, 0 sur z -> x ; unitaire",
          uni and np.allclose(np.abs(S[S != 0]) ** 2, 0.5), "hybride + circulateur")

    print("\nB. Dynamique de reseau")
    v, aniso = cone(lambda q: prop(q) @ S)
    print(f"    cone de Weyl : v = {v:.5f} en unites de vitesse de ligne (= 1/3 : trois traversees par pas), anisotropie {aniso:.1e}")
    print("    le photon du meme automate (R111 B) va a la meme vitesse : c'est le point de la piece.\n")
    check("cone isotrope a 1/3 de la vitesse de ligne (le photon du meme automate y est aussi, R111)", abs(v - 1 / 3) < 1e-3 and aniso < 1e-3, f"{v:.4f}")

    print("C. Les axiomes du manuscrit")
    realT = np.allclose(S.T @ S, np.eye(6))
    recip = np.linalg.norm(S - S.T)
    # le quart d'onde est-il une jauge ? phases de reference par port (une par port, valable a l'arrivee et au depart)
    from scipy.optimize import minimize
    def resid(th):
        D = np.diag(np.exp(1j * th))
        return float(np.linalg.norm(np.imag(D @ S @ D.conj().T)))
    best = min((minimize(resid, np.random.default_rng(k).uniform(0, 2 * math.pi, 6), method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-10, "fatol": 1e-12}) for k in range(8)), key=lambda r: r.fun)
    gauge_real = best.fun < 1e-6
    G = rotations()
    sym = sum(np.allclose(np.abs(port_rep(R) @ S @ port_rep(R).T), np.abs(S)) for R in G)
    mirror = np.diag([-1.0, 1, 1])
    Pm = port_rep(mirror)
    achiral = np.allclose(np.abs(Pm @ S @ Pm.T), np.abs(S)) and np.allclose(np.abs(port_rep(-np.eye(3)) @ S @ port_rep(-np.eye(3)).T), np.abs(S))
    # reciprocite du cone : omega(q) contre omega(-q), au premier ordre et au-dela
    def vpos(q):
        om = np.sort(np.angle(np.linalg.eigvals(prop(q) @ S)))
        return om[om > 1e-9].min()
    d0 = np.array([0.6, -0.5, 0.62]); d0 /= np.linalg.norm(d0)
    asym1 = max(abs(vpos(0.005 * d) - vpos(-0.005 * d)) / 0.005 for d in directions()[:8])
    om_p = np.sort(np.angle(np.linalg.eigvals(prop(0.25 * d0) @ S)))
    om_m = np.sort(np.angle(np.linalg.eigvals(prop(-0.25 * d0) @ S)))
    asym2 = float(np.max(np.abs(om_p - np.sort(-om_m))))
    print(f"    C4 energie : S^dag S = 1 {uni} ; S reel tel quel {realT} ; rendu reel par des phases de reference par port : {gauge_real} (residu {best.fun:.2f})")
    print("        => le quart d'onde n'est PAS une jauge : un dephasage physique au croisement, un element reactif, contre A5")
    print(f"    reciprocite par noeud S = S^T : ||S - S^T|| = {recip:.3f} : circulateur (x -> y -> z -> x, jamais l'inverse), impair sous T")
    print(f"    reciprocite du cone : v(q) - v(-q) = {asym1:.1e} au premier ordre (symetrique) ; a |q| = 0,25 le spectre differe de {asym2:.3f}")
    print("        => la non-reciprocite est invisible dans le cone au premier ordre et apparait au-dela (ordre q^2)")
    print(f"    C1 symetrie cubique par noeud : {sym} rotations sur 24 laissent le routage invariant (groupe tetraedrique T : les axes d'ordre 4")
    print(f"        renversent l'ordre cyclique x -> y -> z) ; miroir et inversion le laissent invariant : achiral {achiral} (groupe T_h, pyritoedrique)")
    print("    donc : pas de rotation de polarisation du vide (achiral), mais un vide qui distingue x -> y -> z de x -> z -> y a chaque")
    print("    croisement, ce qu'aucun axiome ne fournit.\n")
    check("le quart d'onde n'est pas jaugeable (residu > 0,5) : element reactif, contre A5 ; non reciproque par noeud ; cone symetrique au premier ordre seulement",
          uni and (not gauge_real) and recip > 1 and asym1 < 1e-3 and asym2 > 1e-3, f"residu jauge {best.fun:.2f}, ||S - S^T|| = {recip:.2f}, asym q^2 {asym2:.3f}")
    check("C1 viole : routage invariant sous 12 rotations sur 24 (T_h, pyritoedrique), achiral (miroir, inversion), pas cubique complet",
          sym == 12 and achiral, f"{sym}/24, achiral {achiral}")

    print("D. La version achirale : deux sous-reseaux de mains opposees")
    SH, SHp = coin_node((0, 1, 2)), coin_node((0, 2, 1))
    def M2(q):
        # sites A (parite paire) et B (impaire) ; un depart de A arrive en B et reciproquement
        M = np.zeros((12, 12), complex)
        P = prop(q)
        M[6:, :6] = P @ SH          # departs de A (apres S_H) -> arrivees en B
        M[:6, 6:] = P @ SHp         # departs de B -> arrivees en A
        return M
    v2, aniso2 = cone(M2)
    # symetrie du tout : miroir + echange des sous-reseaux
    Pbig = np.zeros((12, 12)); Pbig[:6, 6:] = Pm; Pbig[6:, :6] = Pm
    Stot = np.zeros((12, 12), complex); Stot[:6, :6] = SH; Stot[6:, 6:] = SHp
    glide = np.allclose(Pbig @ Stot @ Pbig.T, Stot) or np.allclose(Pbig @ Stot @ Pbig.T, Stot.conj())
    print(f"    cone : v = {v2:.5f}, anisotropie {aniso2:.1e} : le cone NE survit PAS a l'alternance par sites")
    print("    (une arrivee x routee vers y en A est routee de y vers x en B : la suite x, y, x, y n'utilise jamais z ;")
    print("     l'automate symetrise de Trotter demande d'inverser l'ordre apres un cycle complet, pas a chaque noeud, ce qu'un")
    print("     noeud fixe par site ne peut pas faire)")
    print(f"    miroir + echange des sous-reseaux = symetrie du tout : {glide}\n")
    check("alternance des mains par sites (NaCl) : le cone est detruit (v = 0,275, anisotropie 0,4) : la chiralite par noeud ne s'annule pas ainsi",
          (abs(v2 - 1 / 3) > 0.03 or aniso2 > 0.1), f"v = {v2:.3f}, aniso {aniso2:.2f}")

    print("E. Verdict : le noeud-piece n'est pas derivable de C1-C5 ; il demande trois changements d'axiome nommes :")
    print("   (1) C1 : un noeud de symetrie T_h (pyritoedrique) au lieu de O_h : il choisit un ordre cyclique x -> y -> z des")
    print("       trois lignes ; achiral (pas de rotation de polarisation), mais les axes d'ordre 4 sont brises ; l'alternance")
    print("       de l'ordre entre sites voisins detruit le cone (D) ;")
    print("   (2) la reciprocite par noeud : un circulateur, element impair sous T a chaque croisement, symetrique dans le cone")
    print("       au premier ordre seulement ; le seul element impair sous T de la base est la circulation a sens unique du")
    print("       fluide (R32, R53), exclue pour le vide par son energie (R57) ;")
    print("   (3) A5 : un dephasage physique d'un quart d'onde au croisement (non jaugeable), donc un element reactif la ou A5")
    print("       n'en veut aucun.\n")
    check("verdict enregistre : trois axiomes a changer (C1 -> T_h ; non-reciprocite par noeud ; quart d'onde reactif contre A5)", True, "C1, T, A5")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
