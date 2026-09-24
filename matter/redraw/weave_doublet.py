#!/usr/bin/env python3
"""R105 -- Le doublet du milieu sur les trois directions de la trame : l'algebre de Pauli.

R104 F : pour que le quantum du milieu soit l'electron de Dirac, son doublet doit tourner comme un
spineur avec la direction de propagation, avec l'algebre de Pauli sigma_x sigma_y = i sigma_z.
Question : les trois directions de la trame cubique (noeud de Johns, O_h) peuvent-elles porter
cette algebre ?  Reponse calculee : non, et voici ce qui la porte.

  A. la trame reelle ne peut pas : (i) une trame reciproque (renversement du temps T^2 = +1,
     champ (V, I) reel) a des hamiltoniens effectifs reels symetriques, dont la composante
     sigma_y est identiquement nulle : au plus deux matrices de Pauli, jamais trois ;
     (ii) les six demi-lignes d'un noeud O_h se decomposent en A1 + E + T1 : un scalaire, un
     doublet reel E (representation de O dans SO(3), 2 pi -> +1, bosonique) et un vecteur
     (spin 1).  Pas de doublet spinoriel.  C'est le "carrier neutre" de la phase A, redit en
     representations.
  B. ce qui fournit le i : une circulation a sens unique.  Les deux quadratures d'une onde
     progressive tournent sous la translation du temps avec un generateur J, J^2 = -1 ; T le
     renverse (T J T^-1 = -J), donc une onde stationnaire (T-symetrique) n'a pas de structure
     complexe.  L'anneau a sens unique (R53) l'a ; la trame reciproque ne l'a pas.
  C. ce qui porte l'algebre de Pauli : l'orientation d'une unite qui circule (R95).  Les
     coordonnees spinorielles psi(n) = (cos(theta/2), e^{i phi} sin(theta/2)) avec la phase de
     circulation forment un C^2 (fibration de Hopf, psi^dag sigma psi = n) ; une rotation de
     l'orientation agit par l'element SU(2) du relevement de R95, U psi(n) = psi(R n) a une
     phase pres ; le commutateur de deux petites rotations est la rotation de leur commutateur
     et son relevement vaut exp(-i delta^2 sigma_z/2) : [sigma_x/2, sigma_y/2] = i sigma_z/2
     realise par la composition des rotations.  L'algebre est celle de l'orientation, pas des
     directions du reseau.
  D. l'operateur moyen d'un milieu de telles unites : chaque unite porte la paire 1+1 D le long de
     son axe n, h = c (sigma.n)(n.p) dans son repere ; la moyenne sur les orientations
     (isotrope, ou les trois axes cubiques) donne <(sigma.n)(n.p)> = (1/3) sigma.p : la FORME de
     Weyl sort, avec une vitesse c/3.  L'electron a c dans sa relation de dispersion : la regle
     de noeud doit restaurer c pour le secteur spinoriel comme le theoreme 3 du manuscrit le fait
     pour le secteur scalaire ; non calcule ici.
"""
import math
import os
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
sl = load("su2_lift")

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
SIG = [SX, SY, SZ]

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def su2_from_q(q):
    w, x, y, z = q
    return w * I2 - 1j * (x * SX + y * SY + z * SZ)

def spinor(n):
    n = np.asarray(n, float) / np.linalg.norm(n)
    theta = math.acos(max(-1, min(1, n[2])))
    phi = math.atan2(n[1], n[0])
    return np.array([math.cos(theta / 2), np.exp(1j * phi) * math.sin(theta / 2)])

def octahedral_group():
    gens = [sl.rot([0, 0, 1], math.pi / 2), sl.rot([1, 0, 0], math.pi / 2)]
    group = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for g in list(group):
            for h in gens:
                gh = g @ h
                if not any(np.allclose(gh, k, atol=1e-9) for k in group):
                    group.append(gh)
                    changed = True
    return group

def class_of(R):
    angle = math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2)))
    if angle < 1e-9:
        return "E"
    if abs(angle - 2 * math.pi / 3) < 1e-6:
        return "C3"
    if abs(angle - math.pi / 2) < 1e-6:
        return "C4"
    # angle pi : axe de coordonnee (C2) ou diagonale de face (C2')
    w, v = np.linalg.eig(R)
    axis = np.real(v[:, np.argmin(np.abs(w - 1))])
    axis = np.abs(axis / np.linalg.norm(axis))
    return "C2" if np.sum(axis > 0.9) == 1 else "C2p"

CHAR_TABLE = {   # O : E, 8C3, 3C2, 6C4, 6C2'
    "A1": {"E": 1, "C3": 1, "C2": 1, "C4": 1, "C2p": 1},
    "A2": {"E": 1, "C3": 1, "C2": 1, "C4": -1, "C2p": -1},
    "E":  {"E": 2, "C3": -1, "C2": 2, "C4": 0, "C2p": 0},
    "T1": {"E": 3, "C3": 0, "C2": -1, "C4": 1, "C2p": -1},
    "T2": {"E": 3, "C3": 0, "C2": -1, "C4": -1, "C2p": 1},
}

def main():
    print("R105 -- le doublet du milieu sur les trois directions de la trame\n")

    # A(i) : reel symetrique => pas de sigma_y
    print("A. La trame reelle")
    rng = np.random.default_rng(1)
    max_sy = 0.0
    for _ in range(1000):
        a, b, c = rng.normal(size=3)
        H = np.array([[a, b], [b, c]], complex)
        max_sy = max(max_sy, abs(np.trace(H @ SY) / 2))
    print(f"    (i) 1000 hamiltoniens 2x2 reels symetriques : composante sigma_y maximale = {max_sy:.1e} ;")
    print("        un champ (V, I) reel avec T^2 = +1 n'a que sigma_x et sigma_z : deux matrices de Pauli, jamais trois.")
    # A(ii) : les six demi-lignes du noeud O_h
    G = octahedral_group()
    dirs = [np.array(v, float) for v in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    chars = {}
    counts = {}
    for R in G:
        P = np.zeros((6, 6))
        for j, d in enumerate(dirs):
            img = R @ d
            i = next(k for k, e in enumerate(dirs) if np.allclose(img, e, atol=1e-9))
            P[i, j] = 1
        cl = class_of(R)
        chars[cl] = int(round(np.trace(P)))
        counts[cl] = counts.get(cl, 0) + 1
    mult = {}
    for irr, tab in CHAR_TABLE.items():
        mult[irr] = sum(counts[cl] * chars[cl] * tab[cl] for cl in counts) / len(G)
    print(f"    (ii) noeud O_h, {len(G)} rotations, classes {dict(sorted(counts.items()))} ; caracteres de la representation des six")
    print(f"        demi-lignes : {dict(sorted(chars.items()))} ; decomposition : " + " + ".join(f"{int(round(m))} {k}" for k, m in mult.items() if round(m)))
    print("        A1 (scalaire) + E (doublet reel de O dans SO(3) : 2 pi -> +1, bosonique) + T1 (vecteur, spin 1).")
    print("        Aucun doublet spinoriel : le carrier de la phase A (neutre, reel) redit en representations.\n")
    okA = max_sy < 1e-12 and len(G) == 24 and all(abs(mult[k] - v) < 1e-9 for k, v in {"A1": 1, "A2": 0, "E": 1, "T1": 1, "T2": 0}.items())
    check("trame reelle : sigma_y = 0 (T^2 = +1) ; six demi-lignes O_h = A1 + E + T1, pas de spineur", okA, "A1+E+T1")

    # B : le i de la circulation a sens unique
    print("B. Ce qui fournit le i")
    J = np.array([[0, -1], [1, 0]], float)             # generateur de rotation des quadratures (cos, sin)
    T = np.array([[1, 0], [0, -1]], float)             # t -> -t : cos garde, sin change de signe
    print(f"    quadratures (cos wt, sin wt) d'une onde progressive : generateur J, J^2 = {np.diag(J @ J)} ; T J T^-1 = {'-J' if np.allclose(T @ J @ T, -J) else '?'}")
    print("    une onde progressive (sens unique) porte une structure complexe ; une onde stationnaire, T-symetrique, n'en a")
    print("    pas de definie. L'anneau a sens unique (R53) a le i ; la trame reciproque ne l'a pas.\n")
    check("structure complexe = circulation a sens unique : J^2 = -1, T J T^-1 = -J", np.allclose(J @ J, -np.eye(2)) and np.allclose(T @ J @ T, -J), "J")

    # C : l'algebre de Pauli portee par l'orientation (R95)
    print("C. Ce qui porte l'algebre de Pauli : l'orientation d'une unite qui circule (R95)")
    okC = True
    for n in ([0, 0, 1], [1, 0, 0], [0.3, -0.4, 0.866], [-0.6, 0.5, -0.62]):
        psi = spinor(n)
        nb = np.array([np.vdot(psi, S @ psi).real for S in SIG])
        okC &= np.allclose(nb, np.asarray(n) / np.linalg.norm(n), atol=1e-12)
    print(f"    psi^dag sigma psi = n pour quatre orientations : {okC} (fibration de Hopf, R96)")
    for axis, ang in (([0, 1, 0], 0.7), ([1, 1, 0], 2.1), ([1, 0, 0], math.pi)):
        R = sl.rot(axis, ang)
        U = su2_from_q(sl.track(lambda t, R=R, axis=axis, ang=ang: sl.rot(axis, ang * t)))
        for n in ([0, 0, 1], [0.3, -0.4, 0.866]):
            ov = abs(np.vdot(spinor(R @ np.asarray(n, float)), U @ spinor(n)))
            okC &= abs(ov - 1) < 1e-6
    print(f"    U psi(n) = psi(R n) a une phase pres (|recouvrement| = 1 a 1e-6) pour trois rotations : {okC}")
    d = 0.05
    Ux, Uy = su2_from_q(sl.small_step_quaternion(sl.rot([1, 0, 0], d))), su2_from_q(sl.small_step_quaternion(sl.rot([0, 1, 0], d)))
    Uc = np.linalg.inv(Uy) @ np.linalg.inv(Ux) @ Uy @ Ux
    Rc = sl.rot([0, 1, 0], -d) @ sl.rot([1, 0, 0], -d) @ sl.rot([0, 1, 0], d) @ sl.rot([1, 0, 0], d)
    q_net = sl.small_step_quaternion(Rc)
    coef = {k: (np.trace(Uc @ S) / 2) for k, S in zip("xyz", SIG)}
    print(f"    commutateur de rotations de {d} rad autour de x puis y : rotation nette de {2*math.acos(min(1,q_net[0])):.5f} rad (delta^2 = {d*d:.5f})")
    print(f"    relevement SU(2) : coefficients sigma_x {coef['x'].imag:+.2e}, sigma_y {coef['y'].imag:+.2e}, sigma_z {coef['z'].imag:+.5f} (attendu -+ delta^2/2 = {d*d/2:.5f})")
    okC &= abs(abs(coef['z'].imag) - d * d / 2) < 1e-4 and abs(coef['x']) < 1e-4 and abs(coef['y']) < 1e-4 and np.allclose(Uc, su2_from_q(q_net), atol=1e-6)
    print("    => [sigma_x/2, sigma_y/2] = i sigma_z/2 est la loi de composition des rotations de l'orientation, relevee par R95.")
    print("       L'algebre de Pauli est portee par l'orientation de l'unite qui circule, pas par les directions du reseau.\n")
    check("orientation + circulation = C^2 (Hopf) ; U psi(n) = psi(R n) ; commutateur -> -i delta^2 sigma_z/2 (algebre de Pauli)", okC, f"{abs(coef['z'].imag):.5f}")

    # D : operateur moyen d'un milieu d'unites
    print("D. L'operateur moyen d'un milieu de telles unites : h = c (sigma.n)(n.p) par unite")
    rng = np.random.default_rng(7)
    v = rng.normal(size=(200000, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    iso = (v[:, :, None] * v[:, None, :]).mean(axis=0)
    cubic = sum(np.outer(e, e) for e in np.eye(3)) / 3
    err_iso = np.max(np.abs(iso - np.eye(3) / 3))
    print(f"    <n n^T> isotrope (2e5 tirages) = I/3 a {err_iso:.1e} ; trois axes cubiques a poids egaux : exactement I/3.")
    print("    => <(sigma.n)(n.p)> = (1/3) sigma.p : la forme de Weyl sort de la moyenne, avec une vitesse c/3.")
    print("       L'electron a c dans E^2 = p^2 c^2 + m^2 c^4 : la regle de noeud doit restaurer c pour le secteur spinoriel,")
    print("       comme le theoreme 3 du manuscrit le fait pour le secteur scalaire. Non calcule : c'est R106.\n")
    check("moyenne d'orientation : <(sigma.n)(n.p)> = (1/3) sigma.p (forme de Weyl, vitesse c/3 en champ moyen)",
          err_iso < 5e-3 and np.allclose(cubic, np.eye(3) / 3), f"{err_iso:.0e}")

    print("Verdict : les trois directions de la trame donnent A1 + E + T1, du spin 0 et 1, jamais un spineur, et la")
    print("trame reelle n'a pas de sigma_y. Le i vient d'une circulation a sens unique, et l'algebre de Pauli de")
    print("l'orientation d'une unite qui circule (R95). Le champ electronique universel est donc le champ des")
    print("excitations orientees a sens unique du milieu (le secteur antiperiodique de R104), pas un mode de la trame ;")
    print("sa forme moyenne est bien sigma.p, a vitesse c/3 tant que la regle de noeud n'est pas ecrite. CONDITIONNEL.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
