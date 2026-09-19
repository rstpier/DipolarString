#!/usr/bin/env python3
"""R72 -- Point 3 de coherence : un nucleon avec un seul compte (spin, moments,
masse du neutron).

Regle testee : les trois circuits de quark sont statiques (763 MeV, R42), un
seul anneau circule et porte a lui seul le spin 1/2, le moment et le reste de
la masse.  Pour les moments, la circulation est partagee comme la charge :
la paire de quarks identiques en porte 2/3 chacun, le quark impair -1/3
(ce que le modele des quarks appelle SU(6)).  Pour le neutron, la forme du
pole est celle de plus basse energie dans la section (capacite maximale).
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
L1_9 = L1_3 * (3 / 9) ** (2 * math.pi)
LAMBDA_P = HBARC / MP                      # 0.2103 fm
MU_P, MU_N = 2.79284734, -1.91304273       # en mu_N
DM_OBS = MN - MP

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R72 -- un seul compte pour le nucleon\n")

    # A. spin une seule fois : l'anneau
    print("A. Trois circuits statiques + un anneau qui porte S = 1/2")
    E_static = math.pi * HBARC / L1_9                  # 763 MeV, R42
    E_ring = MP - E_static
    R = HBARC / (2 * E_ring)                           # S = R E/c = hbar/2
    print(f"  statique pi hbar c / l1(9) = {E_static:.1f} MeV ; anneau = m_p - statique = {E_ring:.1f} MeV "
          f"({100*E_ring/MP:.1f} % de m_p)")
    print(f"  S = R E_anneau / c = hbar/2  ->  R = {R:.3f} fm (R25 : 0,587 fm par mu_p avec Q = e ; "
          f"R49 : 0,671)")
    print("  -> le spin est compte une fois, par l'anneau ; les circuits n'en portent pas.\n")
    check("R de l'anneau par le spin seul a 10 % du R25 (0,587 fm)", abs(R / 0.587 - 1) < 0.10,
          f"R = {R:.3f} fm")
    check("part circulante 15-25 % de m_p (R25 : 18 %)", 0.15 < E_ring / MP < 0.25,
          f"{100*E_ring/MP:.1f} %")

    # B. moments : quelle charge circule ?
    print("B. Moments magnetiques, mu/mu_N = Q x R/lambda_p")
    q_p = MU_P / (R / LAMBDA_P)
    q_n = MU_N / (R / LAMBDA_P)
    print(f"  Q requis : proton {q_p:+.3f} e, neutron {q_n:+.3f} e ; rapport mu_p/mu_n = {MU_P/MU_N:.4f}")
    rules = {
        "toutes les charges circulent (p : +1, n : 0)": (1.0, 0.0),
        "paire en avant, impair en arriere a egalite (p : 5/3, n : -4/3)": (5/3, -4/3),
        "partage comme la charge : paire 2/3 chacun, impair -1/3 (p : 1, n : -2/3)": (1.0, -2/3),
    }
    for name, (qp, qn) in rules.items():
        ratio = "inf" if qn == 0 else f"{qp/qn:+.3f}"
        print(f"    {name:70s} rapport {ratio}")
    print(f"  -> le partage comme la charge donne -3/2 ({100*((1/(-2/3))/(MU_P/MU_N)-1):+.1f} %) et, "
          f"avec R du spin, Q_p = 1 ({100*(q_p-1):+.0f} %), Q_n = -2/3 ({100*(q_n/(-2/3)-1):+.0f} %).\n")
    check("partage comme la charge : Q_p = 1 et Q_n = -2/3 a 10 %",
          abs(q_p - 1) < 0.10 and abs(q_n / (-2/3) - 1) < 0.10, f"{q_p:.3f}, {q_n:.3f}")
    check("rapport mu_p/mu_n = -3/2 a 3 %", abs((MU_P / MU_N) / (-1.5) - 1) < 0.03,
          f"{MU_P/MU_N:.4f}")

    # C. un seul nombre pour le neutron : le pole de plus basse energie
    print("C. Le neutron : forme du pole = plus basse energie dans la section w x w")
    shapes = {   # c = C/(4 pi eps0 w) (pole_shape.py, 5/5) ; m_n - m_p avec la regle R44 (meme script)
        "cube / bouchon": (0.6607, 1.270),
        "sphere d = w":   (0.5000, 1.324),
        "plaque carree":  (0.3667, 1.380),
        "disque":         (0.3183, 1.410),
    }
    best = max(shapes.items(), key=lambda kv: kv[1][0])
    for name, (c, dm) in shapes.items():
        print(f"    {name:16s} c = {c:.4f}  F = 1/(2c) = {1/(2*c):.3f}  m_n - m_p = {dm:.3f} MeV "
              f"({100*(dm/DM_OBS-1):+.1f} %)")
    print(f"  -> capacite maximale = energie minimale : {best[0]}, m_n - m_p = {best[1][1]:.3f} MeV "
          f"({100*(best[1][1]/DM_OBS-1):+.1f} %), la meme forme que R5 (bouchon).\n")
    check("la forme de plus basse energie est le bouchon (R5) : un seul nombre, -1,8 %",
          best[0] == "cube / bouchon" and abs(best[1][1] / DM_OBS - 1) < 0.03,
          f"{best[1][1]:.3f} MeV")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
