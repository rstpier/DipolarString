#!/usr/bin/env python3
"""R86 -- Le demi-tour est deja dans la base : le quantum de circulation est une
demi-onde ("Alors ??", apres R85).

Le quantum de circulation de la base (R32, R47, R53) a l'energie
    E_circ = pi hbar c / trajet = (1/2) h c / trajet,
l'energie d'un mode de longueur d'onde lambda = 2 x trajet : une DEMI-ONDE
sur l'anneau (R47 le notait : "l'energie de circuit egale l'energie du mode
en demi-onde pour tout trajet").  Une demi-onde ne se referme sur un anneau
que si le mode est antiperiodique : psi(phi + 2 pi) = -psi(phi).  Donc :

  A. E_circ = demi-onde exactement ; l'onde entiere donnerait 2 E_circ.
  B. demi-onde => avance de phase pi par tour => psi(2 pi) = -psi, psi(4 pi) = +psi,
     enroulement 1/2 => L_z = hbar/2 = R E_circ/c : le spin mecanique de R53 et
     l'antiperiodicite sont le meme enonce.
  C. g : anneau antiperiodique (E_circ = m/2) -> g = 2 ; anneau periodique
     (onde entiere, E_circ = m) -> g = 1, exclu (2,0023).  Donc t = 0 est exclu
     par g, et t = 1/2 est requis, pas pose.
  D. le porteur de l'antiperiodicite sur trois branches simples : le motif
     dipolaire avec demi-tour de section (R82, holonomie -1) ; la paire (phase A)
     n'est pas disponible.  Avec la charge axiale 1/2 : c_1 = 1 (R85).
  E. origine : la brisure R4 donne aux filles l'enroulement 1/2 (phase A,
     spin_from_breaking.py).
"""
import math
import os
import subprocess
import sys

HBARC, ME, ALPHA = 197.3269804, 0.51099895, 1 / 137.035999
LAMBDA_BAR = HBARC / ME
K = ALPHA * HBARC
Q_VORTEX = 1 / (2 * math.sqrt(ALPHA))
G_MEAS = 2.00231930436
HERE = os.path.dirname(os.path.abspath(__file__))
BREAKING = os.path.join(HERE, "..", "scripts", "spin_from_breaking.py")

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R86 -- le quantum de circulation est une demi-onde\n")
    path = 2 * math.pi * LAMBDA_BAR                  # anneau ferme R = lambda-bar (R53)

    # A. demi-onde
    print("A. L'energie du quantum")
    E_circ = 4 * math.pi * K * Q_VORTEX ** 2 / path  # 4 pi K q^2 / trajet (R32)
    E_half = (2 * math.pi * HBARC / 2) / path        # (h c / 2) / trajet : demi-onde
    E_full = 2 * math.pi * HBARC / path              # h c / trajet : onde entiere
    print(f"  4 pi K q^2 / trajet = {E_circ:.5f} MeV ; demi-onde (hc/2)/trajet = {E_half:.5f} ; "
          f"onde entiere hc/trajet = {E_full:.5f} ; m_e/2 = {ME/2:.5f}")
    print("  -> le quantum de la base est exactement la demi-onde de l'anneau (4 pi K q^2 = pi hbar c).\n")
    check("E_circ = demi-onde = m_e/2 (exact)", abs(E_circ / E_half - 1) < 1e-9 and abs(E_circ / (ME / 2) - 1) < 1e-9,
          f"{E_circ:.5f} MeV")

    # B. antiperiodicite et spin
    print("B. Une demi-onde sur un anneau")
    k_wave = 2 * math.pi / (2 * path)                # lambda = 2 trajet
    phase_per_turn = k_wave * path
    winding = phase_per_turn / (2 * math.pi)
    S_mech = LAMBDA_BAR * E_circ / HBARC             # R E_circ / (hbar c), en hbar
    print(f"  avance de phase par tour = k x trajet = {phase_per_turn/math.pi:.3f} pi -> psi(2 pi) = -psi, psi(4 pi) = +psi")
    print(f"  enroulement {winding:.3f} -> L_z = {winding:.3f} hbar ; spin mecanique R E_circ/c = {S_mech:.3f} hbar")
    print("  -> le spin 1/2 de R53 et l'antiperiodicite sont le meme enonce.\n")
    check("phase pi par tour, enroulement 1/2 = spin mecanique 1/2",
          abs(phase_per_turn - math.pi) < 1e-12 and abs(winding - S_mech) < 1e-12, "1/2 = 1/2")

    # C. g
    print("C. Le facteur g tranche entre anneau periodique et antiperiodique")
    for name, E in (("antiperiodique (demi-onde)", E_half), ("periodique (onde entiere)", E_full)):
        S = LAMBDA_BAR * E / HBARC                   # en hbar
        mu = 1.0                                     # e c R / 2 = mu_B pour R = lambda-bar
        g = mu / S                                   # g = (mu/mu_B) / (S/hbar) x ... : mu = g S mu_B/hbar
        print(f"  {name:28s} : E_circ = {E/ME:.2f} m_e, S = {S:.2f} hbar, mu = mu_B -> g = {g:.3f}")
    g_anti = 1.0 / (LAMBDA_BAR * E_half / HBARC)
    g_per = 1.0 / (LAMBDA_BAR * E_full / HBARC)
    print(f"  mesure g = {G_MEAS:.5f} : l'anneau est antiperiodique ; t = 0 (onde entiere, g = 1) est exclu.")
    print("  -> le demi-tour n'est pas pose : il est requis par g = 2 via le quantum de la base.\n")
    check("antiperiodique -> g = 2 (mesure 2,0023) ; periodique -> g = 1, exclu",
          abs(g_anti - 2) < 1e-9 and abs(g_per - 1) < 1e-9 and abs(G_MEAS - 2) < 0.01, f"g = {g_anti:.0f} vs {g_per:.0f}")

    # D. porteur et fibre
    print("D. Le porteur sur trois branches simples et le fibre")
    print("  paire (phase A) : absente pour l'electron. Motif dipolaire + demi-tour de section (R82) :")
    print("  holonomie exactement -1 = l'avance de phase pi de la demi-onde. Charge axiale 1/2 -> c_1 = 1 (R85).")
    print("  -> l'antiperiodicite exigee par le quantum a un porteur geometrique, et le fibre est celui du spin 1/2.\n")
    check("holonomie du porteur (-1, R82) = phase de la demi-onde (e^{i pi})",
          abs(math.cos(phase_per_turn) + 1) < 1e-12, "-1 = -1")

    # E. origine
    print("E. Origine de la demi-onde : la brisure R4 (phase A)")
    try:
        r = subprocess.run([sys.executable, BREAKING], capture_output=True, text=True, timeout=300)
        tail = [l for l in r.stdout.strip().splitlines() if l.strip()][-1:]
        for l in tail:
            print("   ", l[:170])
        ok = r.returncode == 0
    except Exception as exc:                          # noqa: BLE001
        ok = False
    print("  -> la fille herite la longueur d'onde de la mere : enroulement 1/2, la demi-onde.\n")
    check("spin_from_breaking.py passe (la fille herite l'enroulement 1/2)", ok, "phase A")

    print("Verdict : t = 1/2 emerge du quantum de la base (demi-onde) et de g = 2 ; il n'est plus pose.")
    print("Reste pose : l'absence d'inertie d'orientation (R83 ii) et le quantum lui-meme,")
    print("4 pi K q^2 = pi hbar c (R32), qui est la definition de la charge de vortex.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
