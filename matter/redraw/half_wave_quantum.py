#!/usr/bin/env python3
"""R86 -- Le quantum de circulation de la base correspond a une action de
demi-tour, pi hbar (corrige apres relecture de l'auteur ; verdict CONDITIONNEL).

Version initiale : "le quantum est une demi-onde, donc l'anneau est
antiperiodique, et g = 2 l'exige".  Deux fautes, relevees par l'auteur :
  (1) R32/R47 ont ecarte la lecture onde ("le fluide n'est pas une onde, c'est
      un vortex") ; "meme energie qu'une demi-onde" n'implique pas "est une
      demi-onde" ;
  (2) circularite : q = e/(2 sqrt alpha) a ete fixe (vortex.py, "spin condition
      restated") pour que E_circ = m/2, d'ou S = hbar/2 et g = 2 ; g ne peut pas
      servir une seconde fois pour prouver l'antiperiodicite.
Formulation retenue (action, sans onde) :
  A. E_circ = 4 pi K q^2 / L = pi hbar c / L exactement (identite R32).
  B. avec p = E_circ/c, l'action du quantum autour du circuit vaut
     oint p dl = E_circ L / c = pi hbar pour tout L ; SI la phase physique du
     quantum est e^{iS/hbar}, l'holonomie vaut e^{i pi} = -1.  Conditionnel.
  C. la chaine q -> E = m/2 -> S = hbar/2 -> g = 2 est exacte a chaque pas :
     g = 2 est la meme condition que le quantum, pas une preuve independante.
  D. porteur geometrique (R82, holonomie -1) et fibre (R85, c_1 = 1) : coherents
     avec l'action pi hbar, ne la derivent pas.
  E. spin_from_breaking.py (phase A) : arithmetique coherente, mais suppose
     l'enroulement 1/2 dans sa lecture onde ; pas une derivation independante.
Etabli : le quantum R32 correspond exactement a une action de demi-tour.
Non etabli : que R32 derive l'antiperiodicite.  Cible : deriver 4 pi K q^2 =
pi hbar c de la structure du DQD, independamment de g (R87).
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
    print("B. L'action du quantum autour du circuit (formulation sans onde)")
    k_wave = 2 * math.pi / (2 * path)                # lambda = 2 trajet
    phase_per_turn = k_wave * path
    winding = phase_per_turn / (2 * math.pi)
    S_mech = LAMBDA_BAR * E_circ / HBARC             # R E_circ / (hbar c), en hbar
    action = E_circ * path / HBARC                   # oint p dl / hbar, p = E/c
    print(f"  oint p dl = E_circ L / c = {action/math.pi:.6f} pi hbar (pour tout L : pi hbar c/L x L/c)")
    print(f"  SI la phase physique est e^(iS/hbar) : holonomie e^(i pi) = {math.cos(action):+.0f} ; L_z = S/(2 pi) = {action/(2*math.pi):.3f} hbar")
    print(f"  spin mecanique R E_circ/c = {S_mech:.3f} hbar : le meme nombre, par la meme identite.")
    print("  -> etabli : action de demi-tour ; conditionnel : que cette action soit la phase du quantum.\n")
    check("action du quantum = pi hbar exactement (identite R32)", abs(action - math.pi) < 1e-12, f"{action/math.pi:.6f} pi")

    # C. g
    print("C. Circularite : g = 2 est la meme condition que le quantum")
    for name, E in (("antiperiodique (demi-onde)", E_half), ("periodique (onde entiere)", E_full)):
        S = LAMBDA_BAR * E / HBARC                   # en hbar
        mu = 1.0                                     # e c R / 2 = mu_B pour R = lambda-bar
        g = mu / S                                   # g = (mu/mu_B) / (S/hbar) x ... : mu = g S mu_B/hbar
        print(f"  {name:28s} : E_circ = {E/ME:.2f} m_e, S = {S:.2f} hbar, mu = mu_B -> g = {g:.3f}")
    g_anti = 1.0 / (LAMBDA_BAR * E_half / HBARC)
    g_per = 1.0 / (LAMBDA_BAR * E_full / HBARC)
    print(f"  chaine : q = e/(2 sqrt alpha) -> E_circ = m/2 -> S = hbar/2 -> g = {g_anti:.0f} (mesure {G_MEAS:.5f}) ;")
    print(f"  l'alternative E_circ = m donnerait g = {g_per:.0f}. Mais q a ete choisi pour E = m/2 (vortex.py) :")
    print("  -> g = 2 ne prouve pas l'action pi hbar, il la restate. La version initiale de R86 etait circulaire.\n")
    check("chaine q -> E = m/2 -> S = 1/2 -> g = 2 exacte (donc g n'est pas une preuve independante)",
          abs(g_anti - 2) < 1e-9 and abs(E_circ / (ME / 2) - 1) < 1e-9, "circulaire")

    # D. porteur et fibre
    print("D. Le porteur sur trois branches simples et le fibre")
    print("  paire (phase A) : absente pour l'electron. Motif dipolaire + demi-tour de section (R82) :")
    print("  holonomie exactement -1 = e^(i pi) de l'action. Charge axiale 1/2 -> c_1 = 1 (R85).")
    print("  -> coherents avec l'action pi hbar ; ils ne la derivent pas.\n")
    check("holonomie du porteur (-1, R82) = e^(i x action) de R32", abs(math.cos(action) + 1) < 1e-12, "-1 = -1")

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
    print("  -> arithmetique coherente, mais elle suppose l'enroulement 1/2 (lecture onde) : pas independante.\n")
    check("spin_from_breaking.py passe (coherent, non independant)", ok, "phase A")

    print("Verdict (corrige) : CONDITIONNEL. Etabli : le quantum R32 est exactement une action de")
    print("demi-tour pi hbar. Non etabli : que R32 derive l'antiperiodicite. Cible : deriver")
    print("4 pi K q^2 = pi hbar c de la structure du DQD, independamment de g (R87).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
