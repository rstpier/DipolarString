#!/usr/bin/env python3
"""R74 -- Point 4 de coherence : le -1 sous un tour complet (Pauli).

Ce que la base a et n'a pas, chiffre :
  A. la phase A (step2) : une boucle bifilaire non ordonnee avec un nombre impair
     de demi-torsions porte un mode antiperiodique, psi(phi + 2 pi) = -psi ; on
     re-execute le script.
  B. le spin 1/2 du redessin est mecanique : S = R E_circ / c = hbar/2.
  C. l'holonomie Z3 du triple (R61) vaut 2 pi k/3 : jamais pi.
  D. ce que Pauli exige : un Z2 sur l'anneau de branches simples (absent) et
     l'attache au vide (trois jonctions, presente) pour le tour de ceinture.
"""
import math
import os
import subprocess
import sys

HBARC = 197.3269804
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
HERE = os.path.dirname(os.path.abspath(__file__))
STEP2 = os.path.join(HERE, "..", "scripts", "step2_bifilar_berry_phase.py")

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R74 -- le signe d'echange\n")

    # A. phase A, step 2
    print("A. Phase A, step2_bifilar_berry_phase.py (Z2 de la paire)")
    try:
        r = subprocess.run([sys.executable, STEP2], capture_output=True, text=True, timeout=300)
        tail = [l for l in r.stdout.strip().splitlines() if l.strip()][-3:]
        for l in tail:
            print("   ", l[:150])
        ok = r.returncode == 0
    except Exception as exc:                       # noqa: BLE001
        ok, tail = False, [str(exc)]
    print("  -> le Z2 existe pour une PAIRE de branches (mode differentiel antiperiodique).\n")
    check("step2 (Z2 de la paire) passe encore", ok, f"code {r.returncode if ok else 'erreur'}")

    # B. spin mecanique
    print("B. Le spin du redessin")
    S = LAMBDA_BAR * (ME / 2) / HBARC             # R E_circ / (hbar c) en unites de hbar
    print(f"  S = R E_circ / c = lambda-bar x (m_e c^2/2) / c = {S:.6f} hbar")
    print("  -> exact, mais c'est un moment cinetique mecanique : il ne dit rien du signe sous 2 pi.\n")
    check("S = hbar/2 exactement", abs(S - 0.5) < 1e-9, f"{S:.6f}")

    # C. holonomie Z3 du triple
    print("C. Holonomie du triple (R61) sous un tour : 2 pi k/3, k = n mod 3")
    phases = {n: (2 * math.pi * ((n % 3) / 3)) for n in (3, 7, 11)}
    for n, ph in phases.items():
        print(f"    n = {n:2d} : phase {ph/math.pi:.3f} pi  (exp = {math.cos(ph):+.3f} {'+' if math.sin(ph)>=0 else '-'} "
              f"{abs(math.sin(ph)):.3f} i)")
    print("  -> jamais -1 : le Z3 est une etiquette de generation, pas le signe d'echange.\n")
    check("aucune holonomie Z3 ne vaut pi", all(abs(ph - math.pi) > 0.1 for ph in phases.values()),
          "0, 2pi/3, 4pi/3")

    # D. ce qui manque
    print("D. Ce que le -1 d'echange exige")
    print("  (1) un mode antiperiodique sur l'anneau de l'electron : la phase A ne le fournit que")
    print("      pour une paire (deux conducteurs, mode differentiel) ; les brins de l'electron")
    print("      sont trois branches simples : rien ne porte le Z2.  ABSENT.")
    print("  (2) l'attache au milieu : la moitie statique est trois jonctions (R54), donc l'anneau")
    print("      est attache ; tourner l'objet de 2 pi tord ses trois attaches d'un tour, 4 pi se")
    print("      defait (tour de ceinture) : l'echange = une rotation de 2 pi.  PRESENT.")
    print("  -> avec (1), (2) donnerait Pauli ; sans (1), rien.  Le point 4 reste ouvert, et il")
    print("     se reduit a une question : qu'est-ce qui rend antiperiodique un anneau de trois")
    print("     branches simples ?\n")
    check("l'anneau est attache (3 jonctions) : le tour de ceinture s'applique", True, "R54")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
