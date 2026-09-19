#!/usr/bin/env python3
"""R73 -- Point 2 de coherence : une seule echelle de masse.

L'echelle (n/3)^(2 pi) donne mu et tau a 1 % ; Koide les donne a 1e-5.  Peut-on
corriger l'echelle avec un seul nombre ?  Sinon : les masses sont Koide, n
n'est qu'une etiquette, et 2 pi l'approximation a 1 %.
"""
import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
TWO_PI = 2 * math.pi

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def ladder(n, p=TWO_PI, delta=0.0):
    return ME * ((n + delta) / (3 + delta)) ** p

def koide_masses(theta):
    """Brannen : sqrt(m_k) = A [1 + sqrt2 cos(theta + 2 pi k/3)], A fixe par m_e"""
    f = [1 + math.sqrt(2) * math.cos(theta + TWO_PI * k / 3) for k in range(3)]
    f_sorted = sorted(f)
    A = math.sqrt(ME) / f_sorted[0]
    return [(A * x) ** 2 for x in f_sorted]

def main():
    print("R73 -- l'echelle contre Koide\n")

    # (i) l'echelle nue
    mu0, tau0 = ladder(7), ladder(11)
    print(f"(i)   echelle 2 pi : mu {mu0:.2f} ({100*(mu0/MMU-1):+.2f} %), tau {tau0:.1f} "
          f"({100*(tau0/MTAU-1):+.2f} %)")

    # (ii) un decalage delta ajuste sur le muon
    target = (MMU / ME) ** (1 / TWO_PI)          # (7+d)/(3+d)
    delta = (7 - 3 * target) / (target - 1)
    tau_d = ladder(11, delta=delta)
    print(f"(ii)  decalage n -> n + delta, delta = {delta:+.4f} (muon exact) : tau {tau_d:.1f} "
          f"({100*(tau_d/MTAU-1):+.2f} %)")

    # (iii) un exposant p ajuste sur le muon
    p = math.log(MMU / ME) / math.log(7 / 3)
    tau_p = ladder(11, p=p)
    print(f"(iii) exposant p = {p:.4f} (muon exact) : tau {tau_p:.1f} ({100*(tau_p/MTAU-1):+.2f} %)")
    print("  -> un seul nombre de plus rend le tau pire (+2 %) que 2 pi nu (+1 %) : l'echelle")
    print("     n'est pas une loi a corriger, c'est une approximation.\n")
    check("aucune correction a un parametre ne fait mieux que 2 pi sur le tau",
          abs(tau_d / MTAU - 1) > abs(tau0 / MTAU - 1) and abs(tau_p / MTAU - 1) > abs(tau0 / MTAU - 1),
          f"delta : {100*(tau_d/MTAU-1):+.1f} %, p : {100*(tau_p/MTAU-1):+.1f} %, 2 pi : {100*(tau0/MTAU-1):+.1f} %")

    # (iv) Koide + 2/9 + m_e
    m_e_k, m_mu_k, m_tau_k = koide_masses(2 / 9)
    print(f"(iv)  Koide + theta = 2/9 + m_e : mu {m_mu_k:.4f} ({100*(m_mu_k/MMU-1):+.4f} %), "
          f"tau {m_tau_k:.3f} ({100*(m_tau_k/MTAU-1):+.4f} %)")
    check("Koide + 2/9 predit mu et tau a 0,01 %",
          abs(m_mu_k / MMU - 1) < 1e-4 and abs(m_tau_k / MTAU - 1) < 1e-4,
          f"mu {100*(m_mu_k/MMU-1):+.4f} %, tau {100*(m_tau_k/MTAU-1):+.4f} %")

    # (v) ce que l'echelle vaut comme Koide
    ms = [ME, mu0, tau0]
    Q = sum(ms) / sum(math.sqrt(m) for m in ms) ** 2
    Q_true = (ME + MMU + MTAU) / (math.sqrt(ME) + math.sqrt(MMU) + math.sqrt(MTAU)) ** 2
    print(f"(v)   Q = sum m / (sum sqrt m)^2 : echelle {Q:.4f}, mesure {Q_true:.6f}, Koide 2/3")
    check("l'echelle est Koide a 1 %", abs(Q / (2 / 3) - 1) < 0.01, f"Q = {Q:.4f}")

    # (vi) consequence sur ce qui utilise l'echelle
    l1_ratio = (3 / 9) ** TWO_PI
    print(f"(vi)  ce qui garde l'echelle : le barreau des quarks, l1(9)/l1(3) = 3^(-2 pi) = {l1_ratio:.3e} ;")
    print("      son incertitude de 1 % est sous celles du reseau (sigma : 2 %, n-p fort : 12 %).")
    print("  -> decision de coherence : masses des leptons = Koide ; n = etiquette ; 2 pi = 1 %.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
