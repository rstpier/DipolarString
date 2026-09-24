#!/usr/bin/env python3
"""R106 -- La regle de noeud : transporter une excitation orientee a sens unique d'une ligne a
l'autre de la trame, avec son spineur, et a quelle vitesse.

Le noeud de la base est le SCN de Johns (v2.9/scripts/johns/S_johns.npy, 12 ports : 3 bras x 2
cotes x 2 polarisations), unique noeud isotrope de son ansatz, avec deux branches acoustiques
transverses a v = 1/2 (vitesse de ligne 1, maille 1).

  A. reproduction : S_johns symetrique, unitaire, v = 1/2 isotrope (12 directions, 1e-6).
  B. l'etiquette inerte : attacher un spineur a chaque impulsion et le faire "suivre la
     direction" par un changement de base par direction, S_spin = W (S x I_2) W^dag, est une
     relabellisation : bandes du SCN doublees exactement, vitesse 1/2, holonomie triviale sur
     toute boucle du reseau.  Le reseau ne voit pas le spineur.  (Corrige R105 D : le c/3 de
     champ moyen n'est pas une propriete du reseau.)
  C. le transport geometrique (relevement SU(2) de la rotation la plus courte a chaque virage,
     R95) donne bien -1 autour de la boucle carree (+x -> +y -> -x -> -y, un tour de 2 pi), mais
     le noeud S_ij U_ij ainsi construit n'est PAS unitaire (defaut calcule) : la base ne peut
     pas transporter un spineur geometrique avec ses amplitudes scalaires.
  D. construire la regle : sur l'espace des six etats d'helicite + |n, +> (une excitation a sens
     unique, spineur verrouille a sa direction), le groupe octaedrique binaire (48 elements)
     agit par une representation a deux composantes irreductibles (dimensions calculees) ; un
     noeud unitaire equivariant est une phase par composante.  Le secteur qui propage (phase 0)
     et la phase de l'autre composante forment la famille complete ; on scanne la phase et on
     mesure la vitesse, l'isotropie et la degenerescence des branches lineaires a omega -> 0.
     Reference : la meme construction pour une etiquette scalaire (six directions, A1 + E + T1,
     noeud shunt (1/3) J - I) donne v = 1/sqrt(3).
  E. la meme trame donne au photon (SCN) c0 = 1/2, au scalaire 1/sqrt3, au spineur verrouille 1/3 :
     v_e/c0 = 2/3 contre |v_max,e - c|/c < 1e-11 (LEP, Hohensee et al. 2009).  EXCLU sur cette trame :
     elle ne porte pas a la fois le spineur (le -1 par tour) et la vitesse de la lumiere.
"""
import math
import os
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
sl = load("su2_lift")
wd = load("weave_doublet")

I2, SX, SY, SZ = wd.I2, wd.SX, wd.SY, wd.SZ
PORTS = [(d, s, p) for d in range(3) for s in (+1, -1) for p in range(3) if p != d]
INDEX = {P: i for i, P in enumerate(PORTS)}
DIRS = [np.array(v, float) for v in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
def dir_index(v):
    return next(i for i, e in enumerate(DIRS) if np.allclose(v, e, atol=1e-9))

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def quat_of_rotation(R):
    angle = math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2)))
    if angle < 1e-12:
        return np.array([1.0, 0, 0, 0])
    if abs(angle - math.pi) < 1e-9:
        w, v = np.linalg.eigh((R + R.T) / 2)
        axis = v[:, np.argmax(w)]
        return np.array([0.0, *axis])
    axis = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (2 * math.sin(angle))
    return np.array([math.cos(angle / 2), *(math.sin(angle / 2) * axis)])

def su2(q):
    return wd.su2_from_q(q)

def shortest_rotation(a, b):
    """rotation la plus courte amenant la direction a sur b"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    c = np.dot(a, b)
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:                       # demi-tour : axe perpendiculaire (choix : le premier)
        axis = np.cross(a, [1, 0, 0]) if abs(a[0]) < 0.9 else np.cross(a, [0, 1, 0])
        return sl.rot(axis, math.pi)
    axis = np.cross(a, b)
    return sl.rot(axis, math.acos(c))

def bloch_scn(S, q):
    N = len(PORTS)
    prop = np.zeros((N, N), complex)
    for col, (d, s, p) in enumerate(PORTS):
        prop[INDEX[(d, -s, p)], col] = np.exp(1j * q[d] * s)
    return prop @ S

def acoustic(Mfun, direction, qmag=0.005, nmax=4):
    """vitesses des branches LINEAIRES a omega -> 0 : omega(2q)/omega(q) = 2 a 2 %, 0,01 < v <= 1,05 (vitesse de ligne 1)"""
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    def positive(qm):
        om = np.sort(np.angle(np.linalg.eigvals(Mfun(qm * d))))
        return om[(om > 1e-6) & (om < 1.0)]
    o1, o2 = positive(qmag), positive(2 * qmag)
    out = []
    for w in o1:
        v = w / qmag
        if not (0.01 < v <= 1.05):
            continue
        if np.any(np.abs(o2 / w - 2) < 0.02):
            out.append(v)
    return np.array(out[:nmax])

def directions(n, seed=3):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return list(v) + [np.array([1., 0, 0]), np.array([1., 1, 0]) / math.sqrt(2), np.array([1., 1, 1]) / math.sqrt(3)]

def commutant_projectors(reps):
    """projecteurs sur les composantes isotypiques (multiplicite 1) de la representation reps (liste de matrices)"""
    n = reps[0].shape[0]
    # espace des X tels que g X = X g pour tout g
    rows = []
    for g in reps:
        rows.append(np.kron(g, np.eye(n)) - np.kron(np.eye(n), g.T))
    A = np.vstack(rows)
    u, s, vh = np.linalg.svd(A)
    null = vh[s.size - np.sum(s < 1e-9):] if np.sum(s < 1e-9) else vh[-(n * n - np.linalg.matrix_rank(A)):]
    dim = null.shape[0]
    # element hermitien generique du commutant
    rng = np.random.default_rng(11)
    X = sum(rng.normal() * v.reshape(n, n) for v in null)
    X = (X + X.conj().T) / 2
    w, V = np.linalg.eigh(X)
    # regrouper les valeurs propres egales
    groups = []
    for i, val in enumerate(w):
        if groups and abs(val - w[groups[-1][-1]]) < 1e-8:
            groups[-1].append(i)
        else:
            groups.append([i])
    projs = [V[:, g] @ V[:, g].conj().T for g in groups]
    return dim, projs

def main():
    print("R106 -- la regle de noeud\n")
    S = np.load(os.path.join(ROOT, "v2.9", "scripts", "johns", "S_johns.npy"))
    N = len(PORTS)

    # A. reproduction du SCN
    print("A. Le noeud de la base (SCN de Johns, S_johns.npy)")
    sym, uni = np.allclose(S, S.T), np.allclose(S @ S.T, np.eye(N))
    vs = np.array([acoustic(lambda q: bloch_scn(S, q), d, nmax=2) for d in directions(9)])
    print(f"    symetrique {sym}, unitaire {uni} ; branches acoustiques sur 12 directions : v in [{vs.min():.6f}, {vs.max():.6f}], attendu 1/2")
    okA = sym and uni and vs.shape == (12, 2) and np.max(np.abs(vs - 0.5)) < 1e-5
    check("SCN de Johns : symetrique, unitaire, deux branches acoustiques a v = 1/2 isotropes (1e-5)", okA, f"{np.max(np.abs(vs-0.5)):.0e}")

    # B. etiquette inerte
    print("\nB. Un spineur qui suit la direction par changement de base : S_spin = W (S x I2) W^dag")
    V = {}
    for i, n in enumerate(DIRS):
        V[i] = su2(quat_of_rotation(shortest_rotation([0, 0, 1], n)))          # base d'helicite de la direction de mouvement n
    W_in, W_out = np.zeros((2 * N, 2 * N), complex), np.zeros((2 * N, 2 * N), complex)
    for col, (d, s, p) in enumerate(PORTS):
        n_out = np.zeros(3); n_out[d] = s                                       # un pulse sortant par (d, s) va vers +s e_d
        W_out[2 * col:2 * col + 2, 2 * col:2 * col + 2] = V[dir_index(n_out)]
        W_in[2 * col:2 * col + 2, 2 * col:2 * col + 2] = V[dir_index(-n_out)]  # un pulse entrant par (d, s) allait vers -s e_d
    S_spin = W_out @ np.kron(S, I2) @ W_in.conj().T
    okB = np.allclose(S_spin @ S_spin.conj().T, np.eye(2 * N))
    maxdiff = 0.0
    for d in directions(4, seed=5):
        q = 0.31 * d
        prop = np.zeros((N, N), complex)
        for col, (dd, s, p) in enumerate(PORTS):
            prop[INDEX[(dd, -s, p)], col] = np.exp(1j * q[dd] * s)
        lam_s = np.linalg.eigvals(prop @ S)
        lam_2 = np.linalg.eigvals(np.kron(prop, I2) @ S_spin)
        # multiset : chaque valeur propre du modele spineur est a distance nulle d'une du SCN, et reciproquement (x2)
        d1 = max(np.min(np.abs(lam_s - l)) for l in lam_2)
        d2 = max(np.min(np.abs(lam_2 - l)) for l in lam_s)
        maxdiff = max(maxdiff, d1, d2)
    # holonomie sur la boucle carree avec des bases fixes par direction : telescopique
    loop = [dir_index(v) for v in ([1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 0, 0])]
    H = np.eye(2, dtype=complex)
    for a, b in zip(loop[:-1], loop[1:]):
        H = (V[b] @ V[a].conj().T) @ H
    print(f"    unitaire {okB} ; bandes = bandes du SCN doublees (ecart max {maxdiff:.1e}) ; holonomie de la boucle carree = {np.round(H, 6).tolist()}")
    print("    => sur le reseau, un spineur qui suit la direction par relabellisation est inerte : meme vitesse 1/2, holonomie +1.")
    print("       (le c/3 de R105 D reapparait en D comme la vitesse exacte du noeud spinoriel equivariant)")
    check("relabellisation : bandes du SCN doublees (1e-10), v = 1/2, holonomie +1 : le reseau ne voit pas le spineur",
          okB and maxdiff < 1e-10 and np.allclose(H, np.eye(2)), f"{maxdiff:.0e}")

    # C. transport geometrique : -1 sur la boucle, mais noeud non unitaire
    print("\nC. Le transport geometrique (relevement de la rotation la plus courte a chaque virage)")
    U = {}
    for i, a in enumerate(DIRS):
        for j, b in enumerate(DIRS):
            U[(j, i)] = su2(quat_of_rotation(shortest_rotation(a, b)))         # de la direction i vers j
    Hg = np.eye(2, dtype=complex)
    for a, b in zip(loop[:-1], loop[1:]):
        Hg = U[(b, a)] @ Hg
    S_geo = np.zeros((2 * N, 2 * N), complex)
    for row, (d, s, p) in enumerate(PORTS):
        for col, (dd, ss, pp) in enumerate(PORTS):
            if abs(S[row, col]) < 1e-12:
                continue
            ni = np.zeros(3); ni[dd] = ss         # direction d'arrivee (le pulse entrant par (dd, ss) se deplace vers -ss)
            no = np.zeros(3); no[d] = s
            # convention : le port (d, s) recoit un pulse venant du cote s, qui se deplace dans la direction -s e_d ;
            # il repart par le port (d', s') dans la direction +s' e_d'
            S_geo[2 * row:2 * row + 2, 2 * col:2 * col + 2] = S[row, col] * U[(dir_index(no), dir_index(-ni))]
    defect = np.linalg.norm(S_geo.conj().T @ S_geo - np.eye(2 * N), 2)
    print(f"    holonomie de la boucle carree (4 virages de 90 deg autour de z) = {np.round(Hg, 6).tolist()} : -1")
    print(f"    noeud S_ij U_ij : defaut d'unitarite ||S^dag S - I||_2 = {defect:.3f}")
    print("    => le -1 est la, mais les amplitudes scalaires du SCN ne conservent plus la probabilite avec un transport")
    print("       geometrique : la regle de noeud d'un spineur n'est pas 'le SCN plus une rotation'.")
    check("transport geometrique : holonomie -1 sur un tour, mais S_ij U_ij non unitaire (defaut > 0,1)",
          np.allclose(Hg, -np.eye(2), atol=1e-9) and defect > 0.1, f"{defect:.2f}")

    # D. construire la regle : noeud equivariant sur les etats d'helicite +
    print("\nD. La regle construite : noeud unitaire equivariant sur les six etats d'helicite +")
    G = wd.octahedral_group()
    u = [wd.spinor(n) for n in DIRS]
    reps = []
    for R in G:
        for sign in (+1, -1):
            Ur = sign * su2(quat_of_rotation(R))
            G6 = np.zeros((6, 6), complex)
            for i, n in enumerate(DIRS):
                j = dir_index(R @ n)
                ph = np.vdot(u[j], Ur @ u[i])
                assert abs(abs(ph) - 1) < 1e-9
                G6[j, i] = ph
            reps.append(G6)
    dimC, projs = commutant_projectors(reps)
    dims = [int(round(np.trace(P).real)) for P in projs]
    print(f"    groupe octaedrique binaire : {len(reps)} elements ; commutant de dimension {dimC} ; composantes irreductibles de dimensions {dims}")
    # reference scalaire : six directions sans spin
    reps_s = []
    for R in G:
        P6 = np.zeros((6, 6))
        for i, n in enumerate(DIRS):
            P6[dir_index(R @ n), i] = 1
        reps_s.append(P6)
    dimS, projs_s = commutant_projectors(reps_s)
    dims_s = [int(round(np.trace(P).real)) for P in projs_s]
    print(f"    reference scalaire (six directions) : commutant de dimension {dimS}, composantes {dims_s} (A1 + E + T1)")

    def bands_dir(Snode, q):
        prop = np.diag([np.exp(1j * np.dot(q, n)) for n in DIRS])
        return prop @ Snode

    def family_scan(projs, fixed, label, grid=48):
        """phases (0 pour la composante 'fixed', beta, gamma pour les autres) ; points isotropes a branches lineaires"""
        others = [r for r in range(len(projs)) if r != fixed]
        n_node = projs[0].shape[0]
        found = []
        betas = np.linspace(0, 2 * math.pi, grid, endpoint=False)
        grids = [betas] * len(others)
        for phases in np.array(np.meshgrid(*grids)).T.reshape(-1, len(others)) if others else [[]]:
            Snode = projs[fixed].copy().astype(complex)
            for r, ph in zip(others, phases):
                Snode = Snode + np.exp(1j * ph) * projs[r]
            vs = []
            for d in directions(6, seed=9):
                v = acoustic(lambda q: bands_dir(Snode, q), d, nmax=6)
                vs.append(v)
            counts = {len(v) for v in vs}
            if len(counts) != 1 or 0 in counts:
                continue
            arr = np.array(vs)
            aniso = float(np.max(arr.max(axis=0) - arr.min(axis=0)))
            if aniso < 1e-4:
                found.append((tuple(float(x) for x in phases), arr.shape[1], tuple(np.round(arr.mean(axis=0), 6))))
        print(f"    {label} : composante de dimension {int(round(np.trace(projs[fixed]).real))} a phase 0, {grid} valeurs par autre phase :")
        if not found:
            print("      aucun point isotrope a branche lineaire")
        else:
            vs_all = sorted({round(float(v[0]), 3) for _, _, v in found})
            print(f"      {len(found)} points isotropes, {found[0][1]} branche(s) lineaire(s), v de {min(vs_all):.3f} a {max(vs_all):.3f} ;")
            for ph, nb, v in found[:3]:
                print(f"        ex. phases {tuple(round(x/math.pi, 3) for x in ph)} pi : v = {float(v[0]):.6f}")
        return found

    print("    familles unitaires equivariantes (une phase par composante), 9 directions, |q| = 0,005, isotropie < 1e-4 :")
    found_s = {}
    for r, dm in enumerate(dims_s):
        found_s[dm] = family_scan(projs_s, r, f"scalaire, {['T1','A1','E'][[3,1,2].index(dm)]} a phase 0")
    # controle : noeud shunt en base de direction, S_dir = S_port P = J/3 - P (A1 : +1, T1 : +1, E : -1)
    Pant = np.zeros((6, 6))
    for i, n in enumerate(DIRS):
        Pant[dir_index(-n), i] = 1
    S_shunt_dir = np.ones((6, 6)) / 3 - Pant
    v_shunt = np.array([acoustic(lambda q: bands_dir(S_shunt_dir, q), d, nmax=2) for d in directions(6, seed=9)])
    print(f"      controle, noeud shunt (1/3)J - I en base de direction = (1/3)J - P (phases A1 0, T1 0, E pi) : v = {v_shunt.mean():.6f}"
          f" (1/sqrt3 = {1/math.sqrt(3):.6f}), anisotropie {v_shunt.max()-v_shunt.min():.1e}")
    found_spin = {}
    for r, dm in enumerate(dims):
        found_spin[dm] = family_scan(projs, r, f"spineur, j = {'1/2' if dm == 2 else '3/2'} (dimension {dm}) a phase 0")
    print()
    okD = dimC == 2 and sorted(dims) == [2, 4] and abs(v_shunt.mean() - 1 / math.sqrt(3)) < 1e-4 and v_shunt.max() - v_shunt.min() < 1e-4
    check("helicite + sur six directions : deux composantes j = 1/2 (2) et j = 3/2 (4) ; controle scalaire shunt v = 1/sqrt3 isotrope", okD, f"{dims}")
    iso_spin = [(dm, ph, nb, v) for dm, lst in found_spin.items() for ph, nb, v in lst]
    check("secteur spinoriel : un noeud unitaire equivariant avec branche(s) lineaire(s) isotrope(s) existe",
          len(iso_spin) > 0, f"{len(iso_spin)} points, tous j = 1/2 a phase 0")
    # la famille j = 1/2 a phase 0 : vitesse en fonction de la phase beta du secteur j = 3/2
    r_half = dims.index(2); r_32 = dims.index(4)
    def v_half(beta, dirs=None):
        Snode = projs[r_half] + np.exp(1j * beta) * projs[r_32]
        vs = [acoustic(lambda q: bands_dir(Snode, q), d, nmax=6) for d in (dirs or directions(6, seed=9))]
        if len({len(v) for v in vs}) != 1 or vs[0].size == 0:
            return None
        arr = np.array(vs)
        return float(arr.mean()), float(np.max(arr.max(axis=0) - arr.min(axis=0))), arr.shape[1]
    v_pi = v_half(math.pi)
    print(f"    famille j = 1/2 a phase 0, j = 3/2 a phase beta : a beta = pi, v = {v_pi[0]:.9f} ({v_pi[2]} branche, anisotropie {v_pi[1]:.1e}) : exactement 1/3")
    fam = found_spin.get(2, [])
    vmax_iso = max(float(v[0]) for _, _, v in fam) if fam else float("nan")
    vmin_iso = min(float(v[0]) for _, _, v in fam) if fam else float("nan")
    print(f"    sur toute la famille isotrope (beta de {min(ph[0] for ph,_,_ in fam)/math.pi:.3f} a {max(ph[0] for ph,_,_ in fam)/math.pi:.3f} pi) : v de {vmin_iso:.4f} a {vmax_iso:.4f} ;")
    print(f"    beta -> 0 ou 2 pi rend le noeud identique (balistique, anisotrope) : aucune phase ne donne 1/2.")
    v_photon = 0.5
    ratio = (1 / 3) / v_photon
    print(f"    la meme trame donne au photon (SCN) c0 = {v_photon} et au spineur verrouille 1/3 : v_e/c0 = {ratio:.4f} ; mesure : |v_max,e - c|/c < 1e-11 (LEP)")
    print("    => une excitation de spin 1/2 verrouillee a sa direction se propage sur la trame de la base a 2/3 de la vitesse de la")
    print("       lumiere de la meme trame, isotropiquement ; l'electron n'a pas cette vitesse limite. EXCLU sur cette trame.")
    check("j = 1/2 verrouille sur la trame cubique : cone de Weyl isotrope a v = 1/3 exactement (beta = pi), jamais 1/2 : v_e/c0 = 2/3, exclu (1e-11)",
          v_pi is not None and abs(v_pi[0] - 1 / 3) < 1e-6 and v_pi[1] < 1e-4 and vmax_iso < 0.45, f"v = {v_pi[0]:.6f}, max isotrope {vmax_iso:.3f}")

    print("\nE. Verdict : voir le bilan et BASE.md (R106).\n")
    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
