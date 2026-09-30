#!/usr/bin/env python3
"""R111 -- L'automate de Weyl minimal sur les DQD, et son photon-paire a la meme vitesse.

R110 : aucun reseau d'impulsions (etiquette qui suit la direction) ne porte un spineur a la
vitesse du photon.  Les modeles qui le font sont des automates a piece (Bialynicki-Birula 1994) :
un pas de temps = trois deplacements successifs le long de x, y, z, chacun commande par la
composante du spineur sur l'axe correspondant :
    W(k) = [cos k_z + i sigma_z sin k_z][cos k_y + i sigma_y sin k_y][cos k_x + i sigma_x sin k_x]
(pas a = 1, temps tau = 1).  Le facteur d'un axe, P+ S+ + P- S-, est exactement la paire du
telegraphiste sur une ligne de cet axe (R104 A) : la composante + avance, la composante - recule.

  A. Weyl : W unitaire ; dispersion omega(k) = phases propres ; cone isotrope a petit k avec
     v = 1 (un pas par axe par tau) ; anisotropie d'ordre k^3.
  B. Maxwell (Riemann-Silberstein, F = E + iB, 3 composantes) : meme automate avec les matrices
     de spin 1, W_M(k) = prod [P+ e^{ik} + P0 + P- e^{-ik}] ; deux modes transverses a +-omega,
     un mode longitudinal a omega = 0 (a retirer par la contrainte div F = 0) ; MEME vitesse v = 1.
  C. Le photon comme paire de Weyl (de Broglie 1934, D'Ariano-Perinotti 2014) : le generateur de
     deux spineurs, (sigma x 1 + 1 x sigma)/2, restreint au triplet symetrique, EST la matrice de
     spin 1 ; deux Weyl paralleles de moment K/2 ont omega_1 + omega_2 = |K| : la paire va a la
     vitesse du spineur, exactement ce que la paire de Maxwell demande.
  D. Ce que cela veut dire pour la base : la ligne DQD est deja le facteur d'un axe ; le noeud
     doit etre une PIECE, pas le SCN : une impulsion arrivant sur une ligne z est reexprimee dans
     la base des mouvants de la ligne x (etats propres de sigma_x, amplitudes 1/sqrt2) puis de la
     ligne y (etats propres de sigma_y, amplitudes (1, +-i)/sqrt2 : un dephasage de quart d'onde
     entre les deux mouvants).  Le i de R105 B est ce quart d'onde.  Trois familles de lignes,
     trois bases de mouvants liees par l'algebre de Pauli.
  E. Les prix : doublement (modes lineaires supplementaires aux coins de la zone, omega = pi),
     et l'isotropie seulement au premier ordre (anisotropie ~1 % a k = 0,5).
"""
import math
import itertools
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
SIG = [SX, SY, SZ]
# spin 1 (base cartesienne, (S_i)_jk = -i eps_ijk)
S1 = [np.zeros((3, 3), complex) for _ in range(3)]
for i, j, k in itertools.permutations(range(3)):
    eps = np.linalg.det(np.eye(3)[[i, j, k]])
    S1[i][j, k] = -1j * eps

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def W_weyl(k):
    W = I2.copy()
    for i in range(3):
        W = (math.cos(k[i]) * I2 + 1j * math.sin(k[i]) * SIG[i]) @ W
    return W

def projectors_spin1(i):
    w, V = np.linalg.eigh(S1[i])
    P = {}
    for val in (-1, 0, 1):
        cols = V[:, np.abs(w - val) < 1e-9]
        P[val] = cols @ cols.conj().T
    return P

PROJ1 = [projectors_spin1(i) for i in range(3)]

def W_maxwell(k):
    W = np.eye(3, dtype=complex)
    for i in range(3):
        F = PROJ1[i][1] * np.exp(1j * k[i]) + PROJ1[i][0] + PROJ1[i][-1] * np.exp(-1j * k[i])
        W = F @ W
    return W

def phases(W):
    return np.sort(np.angle(np.linalg.eigvals(W)))

def directions(n=15, seed=9):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return list(v) + [np.array([1., 0, 0]), np.array([1., 1, 0]) / math.sqrt(2), np.array([1., 1, 1]) / math.sqrt(3)]

def cone_velocity(Wfun, qmag=0.005):
    """vitesse du cone extrapolee (Richardson, 2 v(q) - v(2q)) : le produit ordonne des facteurs a une correction O(k)"""
    vs = []
    for d in directions():
        v1 = phases(Wfun(qmag * d)); v1 = v1[v1 > 1e-9].min() / qmag
        v2 = phases(Wfun(2 * qmag * d)); v2 = v2[v2 > 1e-9].min() / (2 * qmag)
        vs.append(2 * v1 - v2)
    return float(np.mean(vs)), float(np.max(vs) - np.min(vs))

def main():
    print("R111 -- l'automate de Weyl minimal et son photon-paire\n")

    print("A. Weyl")
    k = np.array([0.3, -0.2, 0.7])
    W = W_weyl(k)
    uni = np.allclose(W @ W.conj().T, I2)
    v, aniso = cone_velocity(W_weyl)
    # anisotropie a k fini
    big = [phases(W_weyl(0.5 * d))[-1] / 0.5 for d in directions()]
    print(f"    unitaire {uni} ; cone a petit k (extrapole) : v = {v:.6f} (pas/tau), anisotropie {aniso:.1e} sur 18 directions ;")
    print(f"    a |k| = 0,5 : v de {min(big):.4f} a {max(big):.4f} (anisotropie {100*(max(big)-min(big)):.1f} %, ordre k^3)")
    # le facteur d'un axe = paire du telegraphiste
    kz = 0.4
    Fz = math.cos(kz) * I2 + 1j * math.sin(kz) * SZ
    print(f"    facteur z seul : phases propres {np.round(np.angle(np.linalg.eigvals(Fz)), 4)} = +-k_z : la composante + avance, la - recule (R104 A)\n")
    check("Weyl : W unitaire ; cone isotrope a v = 1 pas/tau (1e-4) ; facteur d'un axe = paire du telegraphiste",
          uni and abs(v - 1) < 1e-4 and aniso < 1e-4 and np.allclose(np.sort(np.angle(np.linalg.eigvals(Fz))), [-kz, kz]), f"v = {v:.5f}")

    print("B. Maxwell (Riemann-Silberstein) avec le meme automate et les matrices de spin 1")
    WM = W_maxwell(k)
    uniM = np.allclose(WM @ WM.conj().T, np.eye(3))
    zero_modes = all(np.min(np.abs(phases(W_maxwell(0.005 * d)))) < 1e-9 for d in directions())
    vM, anisoM = cone_velocity(W_maxwell)
    print(f"    unitaire {uniM} ; un mode longitudinal a omega = 0 exactement (a retirer par div F = 0) : {zero_modes} ;")
    print(f"    modes transverses : v = {vM:.6f} pas/tau, anisotropie {anisoM:.1e} : la MEME vitesse que le spineur.\n")
    check("Maxwell : W_M unitaire ; deux modes transverses a v = 1, un mode longitudinal a 0 : meme vitesse que Weyl",
          uniM and abs(vM - 1) < 1e-4 and anisoM < 1e-4 and zero_modes, f"v = {vM:.5f}")

    print("C. Le photon comme paire de Weyl")
    Ssum = [(np.kron(s, I2) + np.kron(I2, s)) / 2 for s in SIG]
    # triplet symetrique
    sym = np.array([[1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]], float).T
    sym[:, 1] /= math.sqrt(2)
    S_trip = [sym.T @ S @ sym for S in Ssum]
    # algebre de spin 1 : [S_x, S_y] = i S_z, S^2 = 2
    comm = np.allclose(S_trip[0] @ S_trip[1] - S_trip[1] @ S_trip[0], 1j * S_trip[2])
    casimir = np.allclose(sum(S @ S for S in S_trip), 2 * np.eye(3))
    K = np.array([0.016, 0.004, -0.012]); Kn = np.linalg.norm(K)
    om1 = phases(W_weyl(K / 2)).max()
    print(f"    (sigma x 1 + 1 x sigma)/2 sur le triplet : [S_x, S_y] = i S_z {comm}, S^2 = 2 {casimir} : c'est le spin 1 de B.")
    print(f"    deux Weyl paralleles de moment K/2 : omega_1 + omega_2 = {2*om1:.5f} contre |K| = {Kn:.5f} : la paire va a la vitesse du spineur.\n")
    check("paire : le generateur de deux spineurs sur le triplet est le spin 1 ; deux Weyl paralleles a K/2 donnent omega = |K| (1 %)",
          comm and casimir and abs(2 * om1 / Kn - 1) < 1e-2, f"{2*om1:.5f} vs {Kn:.5f}")

    print("D. Ce que cela demande a la base : un noeud-piece")
    H = np.array([[1, 1], [1, -1]], complex) / math.sqrt(2)          # sigma_z -> sigma_x
    Q = np.array([[1, 1], [1j, -1j]], complex) / math.sqrt(2)         # sigma_z -> sigma_y (quart d'onde sur le second mouvant)
    okD = np.allclose(H.conj().T @ SX @ H, SZ) and np.allclose(Q.conj().T @ SY @ Q, SZ)
    print(f"    ligne z -> ligne x : base des mouvants tournee par Hadamard (amplitudes 1/sqrt2) : {np.allclose(H.conj().T @ SX @ H, SZ)}")
    print(f"    ligne z -> ligne y : Hadamard + quart d'onde (i) sur un mouvant : {np.allclose(Q.conj().T @ SY @ Q, SZ)} ; c'est le i de R105 B.")
    print("    Le SCN transmet une impulsion avec son etiquette ; la piece la reexprime dans la base de la ligne suivante.")
    print("    A ecrire dans la base : un noeud a trois familles de lignes, diviseur 1/sqrt2, dephasage pi/2 sur la famille y.\n")
    check("noeud-piece : changement de base z -> x (Hadamard) et z -> y (Hadamard + quart d'onde) realisent sigma_x, sigma_y", okD, "H, Q")

    print("E. Les prix")
    corners = [np.array(c, float) * math.pi for c in itertools.product((0, 1), repeat=3)]
    extra = 0
    for c in corners[1:]:
        om = phases(W_weyl(c + 0.005 * np.array([1., 1, 1]) / math.sqrt(3)))
        # cone lineaire autour de omega = pi (mod 2 pi) ?
        dist = np.min(np.abs(np.abs(om) - math.pi))
        extra += int(dist < 0.01)
    print(f"    doublement : {extra} coins de zone (sur 7) portent un cone lineaire a omega = pi ; l'isotropie n'est qu'au premier ordre")
    print(f"    (anisotropie {100*(max(big)-min(big)):.1f} % a |k| = 0,5). Connu (Bialynicki-Birula 1994) ; non resolu ici.\n")
    check("doublement present aux coins de la zone (enregistre, non resolu)", extra >= 1, f"{extra}/7")

    print("Verdict : l'automate a piece porte Weyl et Maxwell a la meme vitesse, et Maxwell y est la paire symetrique")
    print("de deux Weyl. La ligne DQD est deja le facteur d'un axe ; ce qui manque a la base est le noeud-piece :")
    print("trois familles de lignes dont les bases de mouvants sont liees par Hadamard et un quart d'onde. CONDITIONNEL.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
