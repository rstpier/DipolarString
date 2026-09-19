#!/usr/bin/env python3
"""R85 -- Le nombre de Chern du mode transverse transporte sur toutes les
orientations de l'anneau (test propose par l'auteur apres R83/R84) :
A = i <psi|d psi>, F = dA, c_1 = (1/2 pi) int_{S^2} F.

Si le mode deplace de R84, transporte sur la sphere des orientations n de
l'anneau, donne |c_1| = 1, c'est le fibre de spin 1/2 et le -1 sous 2 pi est
une propriete globale, pas une regle de torsion choisie.

Methode : etats sur S^2 construits par rotation rigide de l'anneau,
psi(n) = D(R_n) psi_0, ou psi_0 a une charge axiale k definie (comportement
e^{i k alpha} sous la rotation de l'anneau autour de son propre axe) ; courbure
de Berry par la methode des plaquettes (Fukui-Hatsugai-Suzuki), invariante de
jauge, sommee sur la sphere.

  A. verification : k = 0, 1/2, 1 via les matrices de Wigner -> c_1 = 0, 1, 2.
  B. le coeur deplace de R84, un vecteur du plan perpendiculaire a n tournant
     avec la circulation a sens unique (R53) : charge axiale k = 1 (vecteur de
     polarisation circulaire) -> c_1 = 2, calcule directement sur le champ.
     Ce n'est pas un spineur : c'est un objet d'helicite 1.
  C. avec le demi-tour de section par circuit (t = 1/2), le mode dans le repere
     de section a charge axiale 1/2 -> c_1 = 1 : le fibre de spin 1/2.
  D. donc c_1 mesure t, il ne le derive pas : c_1 = 2 t' avec t' la charge
     axiale du mode ; le -1 exige t = 1/2.  Le candidat de la base pour
     t = 1/2 : la brisure (R4) ou chaque fille herite la moitie de l'enroulement
     de la mere (phase A, spin_from_breaking.py, re-execute).
"""
import math
import os
import subprocess
import sys
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
BREAKING = os.path.join(HERE, "..", "scripts", "spin_from_breaking.py")

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def spin_matrices(j):
    dim = int(round(2 * j + 1))
    m = np.array([j - i for i in range(dim)])
    Jz = np.diag(m).astype(complex)
    Jp = np.zeros((dim, dim), dtype=complex)
    for i in range(1, dim):
        Jp[i - 1, i] = math.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    Jy = (Jp - Jp.conj().T) / (2j)
    return Jy, Jz, m

def wigner_states(j, k, n_theta, n_phi):
    Jy, Jz, m = spin_matrices(j)
    idx = int(np.argmin(abs(m - k)))
    psi0 = np.zeros(len(m), dtype=complex)
    psi0[idx] = 1.0
    thetas = np.linspace(0, math.pi, n_theta + 1)
    phis = np.linspace(0, 2 * math.pi, n_phi, endpoint=False)
    states = np.zeros((n_theta + 1, n_phi, len(m)), dtype=complex)
    for a, th in enumerate(thetas):
        Ry = expm(-1j * th * Jy)
        for b, ph in enumerate(phis):
            states[a, b] = expm(-1j * ph * Jz) @ Ry @ psi0
    return states

def chern_plaquettes(states):
    """somme des phases de plaquette / 2 pi (methode FHS), phi periodique"""
    nt, nph, _ = states.shape
    total = 0.0
    for a in range(nt - 1):
        for b in range(nph):
            b2 = (b + 1) % nph
            p1, p2, p3, p4 = states[a, b], states[a, b2], states[a + 1, b2], states[a + 1, b]
            u = np.vdot(p1, p2) * np.vdot(p2, p3) * np.vdot(p3, p4) * np.vdot(p4, p1)
            total += -np.angle(u)
    return total / (2 * math.pi)

def polarization_states(n_theta, n_phi, helicity=+1):
    """vecteur de polarisation circulaire (e_theta + i h e_phi)/sqrt2 du plan perpendiculaire a n"""
    thetas = np.linspace(0, math.pi, n_theta + 1)
    phis = np.linspace(0, 2 * math.pi, n_phi, endpoint=False)
    states = np.zeros((n_theta + 1, n_phi, 3), dtype=complex)
    for a, th in enumerate(thetas):
        for b, ph in enumerate(phis):
            e_th = np.array([math.cos(th) * math.cos(ph), math.cos(th) * math.sin(ph), -math.sin(th)])
            e_ph = np.array([-math.sin(ph), math.cos(ph), 0.0])
            states[a, b] = (e_th + 1j * helicity * e_ph) / math.sqrt(2)
    return states

def main():
    print("R85 -- nombre de Chern du mode transverse sur la sphere des orientations\n")
    N_T, N_P = 60, 120

    # A. verification Wigner
    print("A. Etats de charge axiale k transportes rigidement : c_1 par plaquettes")
    c1 = {}
    for j, k in ((1, 0), (0.5, 0.5), (1, 1)):
        c = chern_plaquettes(wigner_states(j, k, N_T, N_P))
        c1[k] = c
        print(f"  j = {j}, charge axiale k = {k:4.1f} : c_1 = {c:+.4f}")
    print("  -> c_1 = 2k : un scalaire donne 0, un spineur 1, un vecteur 2.\n")
    check("c_1 = 0, 1, 2 pour k = 0, 1/2, 1 (a 1e-3)",
          all(abs(c1[k] - 2 * k) < 1e-3 for k in c1), f"{[round(float(c1[k]), 4) for k in c1]}")

    # B. le coeur deplace de R84 comme champ vectoriel
    print("B. Le coeur deplace de R84 : un vecteur du plan perpendiculaire a l'axe, tournant avec la")
    print("   circulation a sens unique (R53) -> polarisation circulaire, calculee sur le champ lui-meme")
    c_vec = chern_plaquettes(polarization_states(N_T, N_P, +1))
    c_vec_m = chern_plaquettes(polarization_states(N_T, N_P, -1))
    print(f"  helicite +1 : c_1 = {c_vec:+.4f} ; helicite -1 : c_1 = {c_vec_m:+.4f}")
    print("  -> |c_1| = 2 : le coeur deplace seul est un objet d'helicite 1, pas un spineur.\n")
    check("coeur deplace seul : |c_1| = 2 (helicite 1)", abs(abs(c_vec) - 2) < 1e-3 and abs(abs(c_vec_m) - 2) < 1e-3,
          f"{c_vec:+.4f}, {c_vec_m:+.4f}")

    # C. avec le demi-tour de section
    print("C. Avec le demi-tour de section par circuit (t = 1/2, R78/R82)")
    print("  dans le repere de section le mode a la charge axiale 1/2 : c_1 = 2 x 1/2 = 1 (calcule en A).")
    print("  -> c'est le fibre de spin 1/2 ; le -1 sous 2 pi devient global. Mais seulement avec t = 1/2.\n")
    check("avec t = 1/2 : c_1 = 1 (spin 1/2)", abs(c1[0.5] - 1) < 1e-3, f"{c1[0.5]:+.4f}")

    # D. c_1 mesure t ; candidat de la base pour t = 1/2
    print("D. c_1 = 2 t' : le nombre de Chern mesure la charge axiale, il ne la derive pas.")
    print("  Candidat de la base pour t = 1/2 : la brisure R4 (spin_from_breaking.py, phase A) ou chaque")
    print("  fille herite la moitie de l'enroulement de la mere (enroulement 1 sur six brins -> 1/2 sur trois).")
    try:
        r = subprocess.run([sys.executable, BREAKING], capture_output=True, text=True, timeout=300)
        tail = [l for l in r.stdout.strip().splitlines() if l.strip()][-2:]
        for l in tail:
            print("   ", l[:160])
        ok = r.returncode == 0
    except Exception as exc:                        # noqa: BLE001
        ok, tail = False, [str(exc)]
    print("  -> a tester ensuite : que l'enroulement 1 de la mere se partage en deux demi-torsions de")
    print("     SECTION (une par fille), et non en deux demi-enroulements du triple.\n")
    check("spin_from_breaking.py (phase A) passe encore", ok, f"code {r.returncode if ok else 'erreur'}")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
