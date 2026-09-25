#!/usr/bin/env python3
"""R107 -- Le dernier balayage sur la trame de la base : le spineur avec la polarisation du SCN.

R106 : un spineur verrouille a sa direction, sur les six directions de la trame cubique, va a 1/3
(photon du SCN : 1/2).  Le photon n'atteint 1/2 que grace a sa polarisation (deux etats par
direction).  Question : une excitation qui porte a la fois la polarisation du SCN et le spineur
(douze etats d'helicite +) peut-elle atteindre 1/2 sur la meme trame ?

  A. representations : le photon (douze etats |n, e_p>) se decompose sous O en 2 T1 + 2 T2
     (commutant de dimension 8, centre de dimension 2) ; le spineur polarise (|n, e_p, +>) sous le
     groupe octaedrique binaire en 2 G + E1/2 + E'1/2 (commutant 6, centre 3).  Le SCN en base
     de direction est unitaire, equivariant, et va a 1/2 (reproduit).
  B. la methode exacte au premier ordre : les branches lineaires a omega -> 0 viennent du sous-espace
     E1 de valeur propre 1 du noeud (omega = 0 a q = 0) ; au premier ordre en q, omega = valeurs
     propres de P1 (q.N) P1 avec N_a = diag(n_i . e_a) (perturbation degeneree exacte), ce qui ne
     depend que de E1, pas des phases des autres composantes.  Pour un noeud equivariant, E1 est
     une somme de composantes irreductibles (et, dans un bloc de multiplicite 2, une copie
     quelconque, un point de CP^1, ou les deux).  On enumere donc TOUS les E1 possibles, on calcule
     les vitesses et l'isotropie de chacun : c'est la famille complete, sans balayage de phases.
     Calibration : le E1 du SCN donne deux branches isotropes a v = 1/2.
  C. spineur polarise : la plus grande vitesse isotrope de tous les E1 possibles, contre 1/2 ; et
     la reference du spineur sans polarisation (R106) retrouvee : E1 = j = 1/2 -> 1/3.
"""
import math
import os
import importlib.util
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
wd = load("weave_doublet")
nt = load("node_transport")
DIRS = nt.DIRS

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

STATES = []
for i, n in enumerate(DIRS):
    d = int(np.argmax(np.abs(n)))
    for p in range(3):
        if p != d:
            STATES.append((i, p))
IDX = {s: k for k, s in enumerate(STATES)}
NS = len(STATES)

def rotate_state(R, i, p):
    j = nt.dir_index(R @ DIRS[i])
    e2 = R @ np.eye(3)[p]
    p2 = int(np.argmax(np.abs(e2)))
    return j, p2, int(round(e2[p2]))

def reps_photon(G):
    out = []
    for R in G:
        M = np.zeros((NS, NS))
        for k, (i, p) in enumerate(STATES):
            j, p2, sg = rotate_state(R, i, p)
            M[IDX[(j, p2)], k] = sg
        out.append(M)
    return out

def reps_spinor(G):
    u = [wd.spinor(n) for n in DIRS]
    out = []
    for R in G:
        for sign in (+1, -1):
            Ur = sign * nt.su2(nt.quat_of_rotation(R))
            M = np.zeros((NS, NS), complex)
            for k, (i, p) in enumerate(STATES):
                j, p2, sg = rotate_state(R, i, p)
                M[IDX[(j, p2)], k] = sg * np.vdot(u[j], Ur @ u[i])
            out.append(M)
    return out


def bands_dir(Snode, q):
    prop = np.diag([np.exp(1j * np.dot(q, DIRS[i])) for (i, p) in STATES])
    return prop @ Snode

def evaluate(Snode, dirs):
    vs = [nt.acoustic(lambda q: bands_dir(Snode, q), d, nmax=6) for d in dirs]
    if len({len(v) for v in vs}) != 1 or vs[0].size == 0:
        return None
    arr = np.array(vs)
    return float(arr.mean()), float(np.max(arr.max(axis=0) - arr.min(axis=0))), arr.shape[1], arr.mean(axis=0)

def commutant_basis(reps):
    n = reps[0].shape[0]
    A = np.vstack([np.kron(g, np.eye(n)) - np.kron(np.eye(n), g.T) for g in reps])
    u, s, vh = np.linalg.svd(A)
    k = int(np.sum(s < 1e-9)) + (n * n - s.size)
    null = vh[-k:].conj() if k else np.zeros((0, n * n))   # noyau droit : conjugue des lignes de vh (cas complexe)
    mats = [v.reshape(n, n) for v in null]
    herm = []
    for X in mats:
        for H in (X + X.conj().T, 1j * (X - X.conj().T)):
            if np.linalg.norm(H) > 1e-9:
                herm.append(H / np.linalg.norm(H, 2))
    basis = []
    for H in herm:
        for B in basis:
            H = H - np.trace(B.conj().T @ H).real * B
        if np.linalg.norm(H) > 1e-9:
            basis.append(H / np.linalg.norm(H))
    return basis

def center_dim(basis):
    Mat = np.array([[(B @ C - C @ B).reshape(-1) for C in basis] for B in basis])
    Lin = np.transpose(Mat, (0, 2, 1)).reshape(-1, len(basis))
    return len(basis) - np.linalg.matrix_rank(Lin, tol=1e-8)

def isotypic_blocks(reps, basis):
    """composantes isotypiques (par un element central generique) et, pour chacune, une base alignee des copies"""
    n = reps[0].shape[0]
    # element central generique : combinaison de la base qui commute avec toute la base
    Mat = np.array([[(B @ C - C @ B).reshape(-1) for C in basis] for B in basis])
    Lin = np.transpose(Mat, (0, 2, 1)).reshape(-1, len(basis))
    u, s, vh = np.linalg.svd(Lin)
    k = int(np.sum(s < 1e-8)) + (len(basis) - s.size)
    cvecs = vh[-k:]
    rng = np.random.default_rng(5)
    Z = sum(rng.normal() * sum(c[j] * basis[j] for j in range(len(basis))) for c in cvecs)   # element central generique
    Z = (Z + Z.conj().T) / 2
    assert all(np.linalg.norm(Z @ B - B @ Z) < 1e-8 for B in basis), "element non central"
    w, Vz = np.linalg.eigh(Z)
    blocks = []
    i = 0
    while i < n:
        j = i
        while j < n and abs(w[j] - w[i]) < 1e-7:
            j += 1
        blocks.append(Vz[:, i:j])
        i = j
    out = []
    Y = sum(rng.normal() * B for B in basis)          # element generique du commutant (pour aligner les copies)
    for Wb in blocks:
        dimb = Wb.shape[1]
        # multiplicite : dimension du commutant restreint au bloc
        Yb = Wb.conj().T @ Y @ Wb
        # copies : espaces propres d'un element hermitien generique du commutant restreint
        Hb = Wb.conj().T @ sum(rng.normal() * B for B in basis) @ Wb
        Hb = (Hb + Hb.conj().T) / 2
        wb, Vb = np.linalg.eigh(Hb)
        copies = []
        ii = 0
        while ii < dimb:
            jj = ii
            while jj < dimb and abs(wb[jj] - wb[ii]) < 1e-7:
                jj += 1
            copies.append(Wb @ Vb[:, ii:jj])
            ii = jj
        d = copies[0].shape[1]
        m = len(copies)
        aligned = [copies[0]]
        for c in copies[1:]:
            T = c.conj().T @ Y @ copies[0]                  # entrelaceur copie 0 -> copie c (Schur : proportionnel a un unitaire)
            uT, sT, vT = np.linalg.svd(T)
            aligned.append(c @ (uT @ vT))                   # base de la copie alignee sur la copie 0
        out.append((d, m, aligned))
    return out

def velocity_spectrum(P1, Ns, dirs):
    """valeurs propres non nulles de q.V, V_a = P1 N_a P1, sur les directions : (isotropie, vitesses moyennes)"""
    rows = []
    for d in dirs:
        Vq = sum(d[a] * (P1 @ Ns[a] @ P1) for a in range(3))
        ev = np.sort(np.abs(np.linalg.eigvalsh((Vq + Vq.conj().T) / 2)))
        rows.append(ev[ev > 1e-9])
    if len({len(r) for r in rows}) != 1:
        return None, None
    arr = np.array(rows)
    if arr.size == 0:
        return 0.0, np.array([])
    return float(np.max(arr.max(axis=0) - arr.min(axis=0))), arr.mean(axis=0)

def enumerate_E1(blocks, Ns, dirs, label, cp1_grid=24):
    """toutes les sommes de composantes ; dans un bloc de multiplicite 2 : rien, une copie (CP1), les deux"""
    import itertools
    n = Ns[0].shape[0]
    options = []
    for (d, m, aligned) in blocks:
        opts = [("0", None)]
        if m == 1:
            opts.append(("1", aligned[0]))
        else:
            for t in np.linspace(0, math.pi, cp1_grid // 2 + 1):
                for ph in np.linspace(0, 2 * math.pi, cp1_grid, endpoint=False):
                    if (t == 0 or abs(t - math.pi) < 1e-12) and ph > 0:
                        continue
                    W = math.cos(t / 2) * aligned[0] + np.exp(1j * ph) * math.sin(t / 2) * aligned[1]
                    opts.append((f"copie({t/math.pi:.2f}pi,{ph/math.pi:.2f}pi)", W))
            opts.append(("les deux", np.hstack(aligned)))
    # produit cartesien
        options.append(opts)
    found = []
    for combo in itertools.product(*options):
        Ws = [W for _, W in combo if W is not None]
        if not Ws:
            continue
        W = np.hstack(Ws)
        P1 = W @ W.conj().T
        aniso, vel = velocity_spectrum(P1, Ns, dirs)
        if aniso is None or vel is None or vel.size == 0:
            continue
        if aniso < 1e-6:
            found.append((float(vel.max()), tuple(np.round(vel, 6)), " + ".join(f"[{d}]{name}" for (d, m, _), (name, _) in zip(blocks, combo) if name != "0")))
    found.sort(key=lambda f: -f[0])
    print(f"    {label} : {len(found)} sous-espaces E1 isotropes a branche(s) lineaire(s)")
    seen = set()
    for vmax, vel, desc in found:
        key = vel
        if key in seen:
            continue
        seen.add(key)
        print(f"      v = {vel} : E1 = {desc}")
        if len(seen) >= 8:
            break
    return found

def main():
    print("R107 -- le spineur avec la polarisation du SCN sur la trame de la base\n")
    G = wd.octahedral_group()
    dirs = nt.directions(12, seed=9)
    print("A. Representations et decomposition")
    rp = reps_photon(G)
    rs = reps_spinor(G)
    bp, bs = commutant_basis(rp), commutant_basis(rs)
    cp, cs = center_dim(bp), center_dim(bs)
    blP, blS = isotypic_blocks(rp, bp), isotypic_blocks(rs, bs)
    print(f"    photon (12 etats |n, e_p>) : commutant {len(bp)}, centre {cp}, blocs (dim, mult) {[(d, m) for d, m, _ in blP]} -> 2 T1 + 2 T2")
    print(f"    spineur polarise (12 etats |n, e_p, +>) : commutant {len(bs)}, centre {cs}, blocs {[(d, m) for d, m, _ in blS]} -> 2 G + E1/2 + E'1/2")
    # alignement des copies : le projecteur sur une copie melangee commute avec le groupe
    okalign = True
    for reps, blocks in ((rp, blP), (rs, blS)):
        for d, m, aligned in blocks:
            if m == 2:
                W = (aligned[0] + 1j * aligned[1]) / math.sqrt(2)
                P = W @ W.conj().T
                okalign &= all(np.allclose(g @ P, P @ g, atol=1e-8) for g in reps)
    print(f"    copies alignees (le projecteur sur une copie melangee commute avec le groupe) : {okalign}")
    S = np.load(os.path.join(ROOT, "v2.9", "scripts", "johns", "S_johns.npy"))
    Sdir = np.zeros((NS, NS))
    for row, (d, s, p) in enumerate(nt.PORTS):
        for col, (dd, ss, pp) in enumerate(nt.PORTS):
            n_out = np.zeros(3); n_out[d] = s
            n_in = np.zeros(3); n_in[dd] = -ss
            Sdir[IDX[(nt.dir_index(n_out), p)], IDX[(nt.dir_index(n_in), pp)]] = S[row, col]
    eq = all(np.allclose(g @ Sdir @ g.T, Sdir) for g in rp)
    r_scn = evaluate(Sdir, nt.directions(6, seed=9))
    print(f"    SCN en base de direction : unitaire {np.allclose(Sdir @ Sdir.T, np.eye(NS))}, equivariant {eq}, v = {r_scn[0]:.6f} ({r_scn[2]} branches, anisotropie {r_scn[1]:.0e})\n")
    check("photon 2 T1 + 2 T2 (commutant 8, centre 2) ; spineur polarise 2 G + E1/2 + E'1/2 (commutant 6, centre 3) ; copies alignees ; SCN equivariant a v = 1/2",
          len(bp) == 8 and cp == 2 and len(bs) == 6 and cs == 3 and okalign and eq and abs(r_scn[0] - 0.5) < 1e-4, f"{len(bp)}, {cp}; {len(bs)}, {cs}")

    Ns = [np.diag([DIRS[i][a] for (i, p) in STATES]) for a in range(3)]
    print("B. Premier ordre exact : E1 du SCN, puis tous les E1 possibles de la famille du photon (calibration)")
    w, V = np.linalg.eig(Sdir)
    E1 = V[:, np.abs(w - 1) < 1e-9]
    Q, _ = np.linalg.qr(E1)
    P1 = Q @ Q.conj().T
    aniso, vel = velocity_spectrum(P1, Ns, dirs)
    print(f"    E1(SCN) de dimension {E1.shape[1]} : vitesses au premier ordre {np.round(vel, 6)}, anisotropie {aniso:.0e}")
    foundP = enumerate_E1(blP, Ns, dirs, "photon")
    vmaxP = max((f[0] for f in foundP), default=float("nan"))
    check("calibration : le premier ordre redonne v = 1/2 pour le E1 du SCN, et l'enumeration des E1 du photon contient 1/2",
          aniso < 1e-6 and vel is not None and abs(vel.max() - 0.5) < 1e-6 and any(abs(f[0] - 0.5) < 1e-6 for f in foundP), f"max {vmaxP:.4f}")

    print("\nC. Tous les E1 possibles de la famille du spineur polarise")
    foundS = enumerate_E1(blS, Ns, dirs, "spineur polarise")
    vmaxS = max((f[0] for f in foundS), default=float("nan"))
    # reference : spineur sans polarisation (R106), E1 = j = 1/2
    rs6 = []
    u6 = [wd.spinor(n) for n in DIRS]
    for R in G:
        for sign in (+1, -1):
            Ur = sign * nt.su2(nt.quat_of_rotation(R))
            M = np.zeros((6, 6), complex)
            for i, n in enumerate(DIRS):
                j = nt.dir_index(R @ n)
                M[j, i] = np.vdot(u6[j], Ur @ u6[i])
            rs6.append(M)
    b6 = commutant_basis(rs6)
    bl6 = isotypic_blocks(rs6, b6)
    N6 = [np.diag([n[a] for n in DIRS]) for a in range(3)]
    found6 = enumerate_E1(bl6, N6, dirs, "reference R106, spineur sans polarisation")
    print(f"    plus grande vitesse isotrope : photon {vmaxP:.4f} ; spineur polarise {vmaxS:.4f} ; spineur seul {max((f[0] for f in found6), default=float('nan')):.4f}")
    check("spineur polarise : la plus grande vitesse isotrope de toute la famille reste sous 1/2 (ecart > 5 %)",
          not math.isnan(vmaxS) and vmaxS < 0.475, f"{vmaxS:.4f}")
    check("reference R106 : E1 = j = 1/2 donne 1/3 au premier ordre",
          any(abs(f[0] - 1 / 3) < 1e-6 for f in found6), "1/3")

    print("\nVerdict : voir BASE.md (R107).\n")
    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
