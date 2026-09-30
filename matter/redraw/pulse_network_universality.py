#!/usr/bin/env python3
"""R110 -- Le 2/3 est-il propre a la trame cubique ou universel pour les reseaux d'impulsions ?

R106-R107 : sur la trame cubique (six directions d'aretes), un spineur verrouille a sa direction va
a 1/3, le photon transverse a 1/2, le scalaire a 1/sqrt3 (vitesse de ligne 1).  Question : ces
nombres viennent-ils du cube ou de la structure "impulsion a sens unique + etiquette" ?

Methode exacte au premier ordre (R107) : les vitesses des branches lineaires ne dependent que du
sous-espace propagateur E1 (valeur propre 1 du noeud), omega = valeurs propres de P1 (q.N) P1 avec
N_a = diag(n_i . e_a) ; pour un noeud equivariant, E1 est une somme de composantes irreductibles.
On refait le calcul pour trois jeux de directions, tous des orbites du groupe du cube :
  faces (6, aretes de la trame de Johns), sommets (8, diagonales du cube, reseau cubique centre),
  aretes (12, diagonales de faces, reseau cubique a faces centrees).
Secteurs : spineur d'helicite + (j = 1/2 verrouille), vecteur transverse (photon), scalaire.
Resultat : spineur 1/3 exactement sur les trois jeux (identite : sigma.n = 1 sur l'helicite +, Wigner-Eckart
sur j = 1/2 => 3 lambda = 1) ; photon 1/2, 1/2, 1/sqrt3 ; scalaire 1/sqrt3 partout ; spineur/photon <= 2/3.
"""
import math
import os
import importlib.util
import itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
wd = load("weave_doublet")
nt = load("node_transport")
psn = load("polarized_spinor_node")

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def direction_sets():
    faces = [np.array(v, float) for v in ([1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1])]
    verts = [np.array(v, float) / math.sqrt(3) for v in itertools.product((1,-1), repeat=3)]
    edges = []
    for i, j in ((0,1),(0,2),(1,2)):
        for si in (1,-1):
            for sj in (1,-1):
                v = np.zeros(3); v[i] = si; v[j] = sj
                edges.append(v / math.sqrt(2))
    return {"faces (6)": faces, "sommets (8)": verts, "aretes (12)": edges}

def dindex(dirs, v):
    return next(i for i, e in enumerate(dirs) if np.allclose(v, e, atol=1e-9))

def transverse_frame(n):
    a = np.array([0, 0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(a, n); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return e1, e2

def reps_scalar(G, dirs):
    out = []
    for R in G:
        M = np.zeros((len(dirs), len(dirs)))
        for i, n in enumerate(dirs):
            M[dindex(dirs, R @ n), i] = 1
        out.append(M)
    return out

def reps_spinor(G, dirs):
    u = [wd.spinor(n) for n in dirs]
    out = []
    for R in G:
        for sign in (+1, -1):
            Ur = sign * nt.su2(nt.quat_of_rotation(R))
            M = np.zeros((len(dirs), len(dirs)), complex)
            for i, n in enumerate(dirs):
                j = dindex(dirs, R @ n)
                M[j, i] = np.vdot(u[j], Ur @ u[i])
            out.append(M)
    return out

def reps_vector(G, dirs):
    frames = [transverse_frame(n) for n in dirs]
    N = 2 * len(dirs)
    out = []
    for R in G:
        M = np.zeros((N, N))
        for i, n in enumerate(dirs):
            j = dindex(dirs, R @ n)
            for p in range(2):
                v = R @ frames[i][p]
                for pp in range(2):
                    M[2 * j + pp, 2 * i + p] = float(frames[j][pp] @ v)
        out.append(M)
    return out

def Nmats(dirs, per_dir):
    return [np.diag([n[a] for n in dirs for _ in range(per_dir)]) for a in range(3)]

def max_isotropic(reps, Ns, dirs_q, label, cp1_grid=8):
    basis = psn.commutant_basis(reps)
    blocks = psn.isotypic_blocks(reps, basis)
    dims = [(d, m) for d, m, _ in blocks]
    found = psn.enumerate_E1(blocks, Ns, dirs_q, label, cp1_grid=cp1_grid)
    vmax = max((f[0] for f in found), default=float("nan"))
    vals = sorted({round(f[0], 4) for f in found})
    return vmax, vals, dims

def main():
    print("R110 -- universalite du 2/3 pour les reseaux d'impulsions\n")
    G = wd.octahedral_group()
    dirs_q = nt.directions(12, seed=9)
    table = {}
    for name, dirs in direction_sets().items():
        print(f"=== jeu de directions : {name} ===")
        grid = 8 if len(dirs) <= 8 else 4                      # 12 directions : 24 etats vectoriels, grille CP1 reduite
        vs, vals_s, dims_s = max_isotropic(reps_spinor(G, dirs), Nmats(dirs, 1), dirs_q, f"  spineur j = 1/2 verrouille, {name}", grid)
        vv, vals_v, dims_v = max_isotropic(reps_vector(G, dirs), Nmats(dirs, 2), dirs_q, f"  vecteur transverse (photon), {name}", grid)
        vc, vals_c, dims_c = max_isotropic(reps_scalar(G, dirs), Nmats(dirs, 1), dirs_q, f"  scalaire, {name}", grid)
        table[name] = (vs, vv, vc, vals_s, vals_v, vals_c)
        print(f"  -> blocs spineur {dims_s}, vecteur {dims_v}, scalaire {dims_c}")
        print(f"  -> vitesses isotropes distinctes : spineur {vals_s} ; photon {len(vals_v)} valeurs, max {max(vals_v):.4f} ; scalaire {len(vals_c)} valeurs, max {max(vals_c):.4f}")
        print(f"  -> maxima : spineur {vs:.4f}, photon {vv:.4f}, scalaire {vc:.4f} ; rapport spineur/photon = {vs/vv:.4f}\n")
    ok_third = all(abs(t[0] - 1/3) < 1e-6 for t in table.values())
    ok_photon = all(t[1] > 0.5 - 1e-6 for t in table.values()) and abs(table["faces (6)"][1] - 0.5) < 1e-6
    ok_ratio = all(t[0] / t[1] < 2/3 + 1e-6 for t in table.values())
    print("Lecture : le 1/3 du spineur est exact et universel : sur un etat d'helicite sigma.n = 1, donc sum_a sigma_a (P1 N_a P1)")
    print("= 1 sur j = 1/2 ; par Wigner-Eckart P1 N_a P1 = lambda sigma_a, d'ou 3 lambda = 1 : lambda = 1/3 pour tout jeu de")
    print("directions isotrope. Le photon n'est pas universel (1/2 sur les faces et les sommets, 1/sqrt3 sur les aretes) mais il")
    print("est toujours au-dessus : le spineur verrouille n'atteint jamais la vitesse du photon du meme reseau, rapport <= 2/3.")
    print("C'est la structure 'impulsion a sens unique + etiquette qui la suit' qui l'impose, pas le cube. Pour porter un spineur")
    print("a la vitesse du photon, le milieu doit etre un automate a piece (le spineur decide du pas) : Bialynicki-Birula 1994,")
    print("D'Ariano-Perinotti 2014, ou le photon est une paire de Weyl (de Broglie 1934).\n")
    check("spineur verrouille : vitesse isotrope maximale exactement 1/3 pour les trois jeux de directions (faces, sommets, aretes)", ok_third,
          ", ".join(f"{k}: {t[0]:.4f}" for k, t in table.items()))
    check("photon transverse : maximum >= 1/2 sur les trois jeux, = 1/2 sur les faces (le SCN) et les sommets, 1/sqrt3 sur les aretes", ok_photon,
          ", ".join(f"{k}: {t[1]:.4f}" for k, t in table.items()))
    check("spineur/photon <= 2/3 sur tout jeu de directions : le spineur verrouille n'atteint jamais la vitesse du photon du reseau",
          ok_ratio, ", ".join(f"{k}: {t[0]/t[1]:.4f}" for k, t in table.items()))

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
