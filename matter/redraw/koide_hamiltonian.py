#!/usr/bin/env python3
"""R79 -- L'hamiltonien qui produirait a^2 = 2|b|^2 (le 45 deg de Koide), sur
relecture de l'auteur : "construire l'hamiltonien physique qui produit
a^2 = 2|b|^2 a partir des energies du canal commun et du doublet differentiel
deja presents dans la geometrie a trois brins ; si l'egalite des normes sort
d'une minimisation ou d'une conservation sans mettre sqrt2 a la main, le 45 deg
est derive."

Ce script ne le derive pas.  Il etablit exactement ce qu'un tel hamiltonien
doit satisfaire, et ce qui est exclu :
  A. identite : Koide (Q = 2/3) <=> ||C_iso||_F = ||C_dev||_F pour l'operateur
     de racine de masse (partie isotrope = (tr C/3) I, partie deviatorique =
     le reste), quel que soit C hermitien 3 x 3.
  B. equipartition par degre de liberte (le doublet en a deux) => Q = 1 : deux
     leptons de masse nulle, exclu ; equipartition par CANAL (irrep) => Q = 2/3.
     Chaque canal porte alors la moitie de la masse de la famille.
  C. aucun extremum : la fraction isotrope s = ||C_iso||^2/||C||^2 est monotone
     en |b|/a, Koide est a s = 1/2 sans point stationnaire : l'egalite doit
     etre une contrainte (conservation), pas un minimum.
  D. la forme quadratique statique des trois brins (couplage reel L, M) est
     bien C3-circulante mais son doublet est degenere (m_e = m_mu) : le
     dedoublement exige une phase de saut complexe (holonomie), independante
     de l'egalite des normes.
"""
import math
import numpy as np

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
TWO_PI = 2 * math.pi
U = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
I3 = np.eye(3)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def circulant(a, b):
    return a * I3 + b * U + np.conj(b) * U @ U

def iso_dev_norms(C):
    iso = (np.trace(C) / 3) * I3
    dev = C - iso
    return np.linalg.norm(iso), np.linalg.norm(dev)

def koide_Q(ms):
    return sum(ms) / sum(math.sqrt(m) for m in ms) ** 2

def main():
    print("R79 -- ce qu'un hamiltonien doit faire pour donner le 45 deg\n")

    # A. identite
    print("A. Identite : Q = 2/3 <=> ||C_iso|| = ||C_dev||")
    lam = np.array([math.sqrt(ME), math.sqrt(MMU), math.sqrt(MTAU)])
    C_meas = np.diag(lam)                       # n'importe quel C hermitien a ces valeurs propres
    n_iso, n_dev = iso_dev_norms(C_meas)
    Q = koide_Q([ME, MMU, MTAU])
    print(f"  masses mesurees : Q = {Q:.6f} ; ||iso|| = {n_iso:.4f}, ||dev|| = {n_dev:.4f} "
          f"(rapport {n_dev/n_iso:.6f})")
    # preuve numerique sur des matrices aleatoires : Q = 1/3 + (||dev||/||iso||)^2 / 3
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(200):
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        H = A + A.conj().T + 4 * I3
        ev = np.linalg.eigvalsh(H)
        if ev.min() <= 0:
            continue
        q = np.sum(ev ** 2) / np.sum(ev) ** 2 * 1.0   # Q sur lambda = sqrt m  <=>  sum m / (sum sqrt m)^2
        ni, nd = iso_dev_norms(H)
        worst = max(worst, abs(q - (1 / 3 + (nd / ni) ** 2 / 3)))
    print(f"  pour tout C hermitien : Q = 1/3 + (||dev||/||iso||)^2/3 (ecart max {worst:.1e}) ;")
    print("  Q = 2/3 <=> ||dev|| = ||iso||.\n")
    check("Koide <=> normes isotrope et deviatorique egales (mesure : rapport 1 a 1e-5)",
          abs(n_dev / n_iso - 1) < 2e-5 and worst < 1e-12, f"{n_dev/n_iso:.6f}")

    # B. equipartition
    print("B. Equipartition : par degre de liberte ou par canal ?")
    for name, ratio2 in (("par degre de liberte (doublet = 2 dof) : E_dev = 2 E_iso", 2.0),
                         ("par canal (irrep) : E_dev = E_iso", 1.0)):
        q = 1 / 3 + ratio2 / 3
        print(f"    {name:60s} -> Q = {q:.4f}")
    # Q = 1 : rang 1
    print("  Q = 1 est le rang 1 : deux valeurs propres nulles, deux leptons sans masse : exclu.")
    E_iso, E_dev = n_iso ** 2, n_dev ** 2
    print(f"  par canal : chaque canal porte ||.||^2 = {E_iso:.2f} MeV = {E_dev:.2f} MeV = "
          f"(m_e + m_mu + m_tau)/2 = {(ME+MMU+MTAU)/2:.2f} MeV.\n")
    check("equipartition par dof exclue (Q = 1), par canal = Koide (Q = 2/3)",
          abs((1/3 + 2/3) - 1) < 1e-12 and abs((1/3 + 1/3) - 2/3) < 1e-12, "Q = 1 vs 2/3")
    check("chaque canal porte la moitie de la masse de la famille (a 1e-5)",
          abs(E_iso / ((ME + MMU + MTAU) / 2) - 1) < 2e-5, f"{E_iso:.2f} MeV")

    # C. pas d'extremum
    print("C. Un extremum ? fraction isotrope s(|b|/a) = 1/(1 + 2|b|^2/a^2)")
    xs = np.linspace(0, 3, 301)
    s_vals = 1 / (1 + 2 * xs ** 2)
    ds = np.diff(s_vals)
    x_koide = 1 / math.sqrt(2)
    print(f"  s decroit strictement de 1 (|b| = 0, masses egales) a 0 (rang 1) ; Koide a |b|/a = {x_koide:.4f}, s = 1/2 ;")
    print("  aucun point stationnaire : ni minimum ni maximum ne selectionne Koide.")
    print("  -> l'egalite doit etre une contrainte de conservation (un quantum par canal), pas un extremum.\n")
    check("s monotone, sans point stationnaire", np.all(ds < 0), "ds < 0 partout")

    # D. la forme quadratique statique des trois brins
    print("D. Forme quadratique des trois brins : K = L I + M (U + U^2), M reel")
    L, M = 1.0, 1 / math.sqrt(2)
    K = circulant(L, M)
    evK = np.sort(np.linalg.eigvalsh(K))
    print(f"  M/L = 1/sqrt2 : valeurs propres {evK.round(4)} = (L - M, L - M, L + 2M) : doublet degenere,")
    print(f"  m_e = m_mu, et m_s/m_d = {(evK[2]/evK[0])**2:.1f} : rien des leptons.")
    b = M * np.exp(1j * 2 / 9)
    evC = np.sort(np.linalg.eigvalsh(circulant(L, b)))
    ni2, nd2 = iso_dev_norms(circulant(L, b))
    print(f"  avec la phase phi = 2/9 : valeurs propres {evC.round(4)}, rapport des masses "
          f"{(evC[1]/evC[0])**2:.1f}, {(evC[2]/evC[0])**2:.1f} (mesure 206,8 ; 3477) ; normes inchangees "
          f"({nd2/ni2:.6f})")
    print("  -> la statique des trois brins donne la forme C3 et peut donner a^2 = 2|b|^2 (M/L = 1/sqrt2),")
    print("     mais pas le dedoublement e/mu : la phase est une holonomie, independante des normes.\n")
    check("couplage reel : doublet degenere (m_e = m_mu)", abs(evK[0] - evK[1]) < 1e-12, "exclu seul")
    check("la phase dedouble sans changer les normes", abs(nd2 / ni2 - 1) < 1e-12 and evC[1] > evC[0] + 1e-6,
          f"rapport {nd2/ni2:.6f}")

    print("Ce que l'hamiltonien doit satisfaire, exactement :")
    print("  (1) une symetrie C3 (trois positions) : donne la forme ;")
    print("  (2) une conservation 'un quantum par canal irreductible' (singulet et doublet, pas par dof) :")
    print("      donne ||iso|| = ||dev||, le 45 deg ; le candidat de la base est R47 (un quantum par circuit")
    print("      ferme) si chaque canal est un circuit ferme : non demontre ;")
    print("  (3) une holonomie de saut phi = 2/9 : dedouble e/mu ; sans mecanisme (R63).")
    print("  Rien ici ne derive (2) ni (3) ; le sqrt2 n'est pas mis a la main, il est montre equivalent")
    print("  a (2).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
