#!/usr/bin/env python3
"""R117 -- Les deux meilleurs nombres du depot contre les determinations de 2024-2026.

Le depot compare sa partie forte de m_n - m_p (2,474 MeV) et sa part QED implicite (-1,20 MeV) au
reseau BMW 2015 (2,52 +- 0,29 ; -1,00 +- 0,16), et sa tension de corde (430,3 MeV) a la bande
420-440 MeV (Necco-Sommer 2002). Nouvelles references :
  - Bulava, Knechtli, Koch, Morningstar, Peardon, PLB 854 (2024) 138754 : sqrt(sigma) = 445(3)(6) MeV,
    QCD a 2+1 saveurs extrapolee aux masses physiques (hamiltonien modele sur le spectre statique) ;
  - regle de somme dispersive pour les polarisabilites isovectorielles, arXiv:2609.15154 (sept. 2026) :
    (m_p - m_n)_QED = 0,71 (+0,03/-0,06) MeV ;
  - Gasser, Leutwyler, Rusetsky, PLB 814 (2021) 136087 (Cottingham) : 0,58 +- 0,16 MeV ;
  - Walker-Loud, Carlson, Miller, PRL 108 (2012) 232301 (Cottingham) : 1,30 +- 0,47 MeV.

  A. sqrt(sigma) : 430,3 contre 445 +- 6,7 : 2,2 sigma ; la fenetre de l'exposant (427-435) reste
     a 1,5 sigma du nouveau nombre : tension, pas exclusion (< 3 sigma).
  B. part QED : la base (1,20 bouchon, 1,15 sphere) est dans l'etalement des determinations
     (0,58 a 1,30), a 1,3 sigma du reseau, mais a 3,9 sigma de GLR et ~8 sigma de la regle de somme
     2026 ; le reseau et le dispersif se contredisent entre eux (1,00 +- 0,16 contre 0,71 +0,03/-0,06).
  C. partie QCD impliquee par chaque determination (1,2933 + QED) : 2,29 (BMW), 1,87 (GLR), 2,59
     (WCM), 2,00 (2026) ; BMW direct 2,52. La formule 2,474 est a 0,2 sigma de BMW direct, 1,1 sigma
     de BMW derive, 3,7 sigma de GLR, 7,9 sigma de la regle de somme.
  D. le total m_n - m_p de la base (1,27-1,32) ne depend pas de la coupe ; la masse du tau de
     Belle II (1777,09 +- 0,14) deplace p(e, tau) de moins de 1e-4.
  Verdict : EN TENSION. La coupe de la base suit le reseau et contredit le dispersif ; le critere
  de HADRONS 8 est precise : une valeur reseau au point physique de la partie QCD sous 2,2 MeV a
  +- 0,1, ou une confirmation reseau du 0,71 dispersif, tue la formule ; 445 +- 7 est a surveiller.
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MMU, MTAU_PDG, MTAU_BELLE2 = 0.51099895, 105.6583755, 1776.86, 1777.09
DM_EXP = 1.29333
F_PLUG = 0.757
LAM = HBARC / ME
L1 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)
W9 = (4 * LAM / math.pi ** 2) * 3 ** (-2 * math.pi)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def base_numbers():
    E_dqd = 4 * math.pi * K / (9 * L1)
    M_mut = K / (math.sqrt(3) * L1)
    def coulomb(F):
        a = F * W9 / 4
        S_u = K / (2 * L1) * (math.log(4 * L1 / a) - 1)
        S_d = K / L1 * (math.log(2 * L1 / a) - 1)
        return (4 / 9) * S_u - (1 / 9) * S_d + M_mut / 3
    sigma = math.pi * HBARC / L1 ** 2
    return dict(strong=E_dqd, qed_plug=coulomb(F_PLUG), qed_sphere=coulomb(1.0), sqrt_sigma=math.sqrt(sigma * HBARC))

def pull(x, ref, err):
    return (x - ref) / err

def main():
    b = base_numbers()
    print("R117 -- les deux meilleurs nombres du depot contre 2024-2026\n")
    print(f"  base : partie forte {b['strong']:.3f} MeV ; part QED {b['qed_plug']:.3f} (bouchon) / {b['qed_sphere']:.3f} (sphere) ; total {b['strong']-b['qed_plug']:.3f} / {b['strong']-b['qed_sphere']:.3f} ; sqrt(sigma) = {b['sqrt_sigma']:.1f} MeV\n")

    # A. sqrt(sigma)
    print("A. La tension de corde")
    SIG_NEW, SIG_ERR = 445.0, math.hypot(3.0, 6.0)
    pw = pull(b["sqrt_sigma"], SIG_NEW, SIG_ERR)
    p_lo, p_hi = 6.2758, 6.2925
    def chain_sigma(p):
        l1 = (2 * math.pi * LAM / 3) * (1 / 3) ** p
        return math.sqrt(math.pi * HBARC / l1 ** 2 * HBARC)
    win = sorted([chain_sigma(p_lo), chain_sigma(p_hi)])
    edge = pull(win[1], SIG_NEW, SIG_ERR)
    print(f"    Bulava et al. 2024 (2+1 saveurs, masses physiques) : {SIG_NEW:.0f} +- {SIG_ERR:.1f} MeV ; ancienne bande 420-440 (Necco-Sommer, trempe et non trempe)")
    print(f"    base {b['sqrt_sigma']:.1f} : {pw:+.2f} sigma ; fenetre de l'exposant [{win[0]:.1f} ; {win[1]:.1f}] : bord a {edge:+.2f} sigma")
    print("  -> tension (2,2 sigma), pas exclusion ; la definition (hamiltonien modele sur le spectre statique avec brisure de corde)")
    print("     n'est pas celle de la bande trempee : a comparer avec une seconde determination au point physique.\n")
    check("A. sqrt(sigma) = 430 contre 445 +- 7 (2024) : 2,2 sigma, le bord de la fenetre a 1,5 sigma : tension, pas exclusion",
          2.0 < abs(pw) < 2.5 and abs(edge) < 2.0, f"{pw:+.2f} sigma, bord {edge:+.2f} sigma")

    # B. part QED
    print("B. La part electromagnetique (m_p - m_n)_QED")
    dets = {"BMW 2015 (reseau)": (1.00, 0.16), "Gasser-Leutwyler-Rusetsky 2021 (Cottingham)": (0.58, 0.16),
            "Walker-Loud-Carlson-Miller 2012 (Cottingham)": (1.30, 0.47), "regle de somme 2026 (arXiv:2609.15154)": (0.71, 0.06)}
    pulls = {}
    for k, (v, e) in dets.items():
        pulls[k] = (pull(b["qed_plug"], v, e), pull(b["qed_sphere"], v, e))
        print(f"    {k:46s} {v:.2f} +- {e:.2f} MeV : base bouchon {pulls[k][0]:+.1f} sigma, sphere {pulls[k][1]:+.1f} sigma")
    lo, hi = min(v for v, _ in dets.values()), max(v for v, _ in dets.values())
    print(f"    etalement des determinations : {lo:.2f} a {hi:.2f} MeV ; reseau contre dispersif : {pull(1.00, 0.71, math.hypot(0.16, 0.06)):+.1f} sigma entre eux")
    print("  -> la base est dans l'etalement, du cote du reseau ; les deux valeurs dispersives sont au-dessous d'elle.\n")
    check("B. part QED 1,20 : dans l'etalement 0,58-1,30, a 1,3 sigma du reseau, > 3 sigma des deux determinations dispersives",
          lo < b["qed_plug"] < hi and abs(pulls["BMW 2015 (reseau)"][0]) < 1.5
          and abs(pulls["Gasser-Leutwyler-Rusetsky 2021 (Cottingham)"][0]) > 3 and abs(pulls["regle de somme 2026 (arXiv:2609.15154)"][0]) > 3,
          f"BMW {pulls['BMW 2015 (reseau)'][0]:+.1f}, GLR {pulls['Gasser-Leutwyler-Rusetsky 2021 (Cottingham)'][0]:+.1f}, 2026 {pulls['regle de somme 2026 (arXiv:2609.15154)'][0]:+.1f} sigma")

    # C. partie QCD
    print("C. La partie QCD de m_n - m_p impliquee (1,2933 + QED), et BMW direct")
    qcd = {k: (DM_EXP + v, e) for k, (v, e) in dets.items()}
    qcd["BMW 2015 direct"] = (2.52, math.hypot(0.17, 0.24))
    worst_disp = 0.0
    for k, (v, e) in qcd.items():
        pq = pull(b["strong"], v, e)
        print(f"    {k:46s} {v:.2f} +- {e:.2f} MeV : formule 2,474 a {pq:+.1f} sigma")
        if "Cottingham" in k or "2026" in k:
            worst_disp = max(worst_disp, abs(pq))
    p_bmw = pull(b["strong"], *qcd["BMW 2015 direct"])
    p_sr = pull(b["strong"], *qcd["regle de somme 2026 (arXiv:2609.15154)"])
    print("  -> la formule (2 alpha/3) 3^(2 pi) m_e c^2 vit avec le reseau et meurt avec le dispersif ; le desaccord est entre eux.\n")
    check("C. partie QCD 2,474 : 0,2 sigma de BMW direct, 3,7 sigma de GLR, ~8 sigma de la regle de somme 2026",
          abs(p_bmw) < 0.5 and abs(pull(b["strong"], *qcd["Gasser-Leutwyler-Rusetsky 2021 (Cottingham)"])) > 3 and abs(p_sr) > 6,
          f"BMW {p_bmw:+.1f}, 2026 {p_sr:+.1f} sigma")

    # D. ce qui ne bouge pas
    print("D. Ce qui ne bouge pas")
    tot = (b["strong"] - b["qed_plug"], b["strong"] - b["qed_sphere"])
    p_tau_pdg = math.log(MTAU_PDG / ME) / math.log(11 / 3)
    p_tau_b2 = math.log(MTAU_BELLE2 / ME) / math.log(11 / 3)
    print(f"    total m_n - m_p de la base : {tot[0]:.3f} (bouchon) / {tot[1]:.3f} (sphere) MeV contre {DM_EXP:.4f} : {100*(tot[0]/DM_EXP-1):+.1f} % / {100*(tot[1]/DM_EXP-1):+.1f} %, independant de la coupe")
    print(f"    exposant depuis le tau : PDG {MTAU_PDG} -> p = {p_tau_pdg:.5f} ; Belle II 2023 {MTAU_BELLE2} -> p = {p_tau_b2:.5f} ; ecart {abs(p_tau_b2-p_tau_pdg):.1e}")
    print("    rayon du proton : 0,8409(4) fm inchange (PRad-II en analyse) ; electron : g-2 2023, alpha Rb/Cs a 5,5 sigma non resolu.\n")
    check("D. le total de la base (1,27-1,32) et l'exposant (Belle II : dp ~ 1e-4) ne bougent pas",
          abs(tot[0] / DM_EXP - 1) < 0.03 and abs(p_tau_b2 - p_tau_pdg) < 2e-4, f"dp = {abs(p_tau_b2-p_tau_pdg):.1e}")

    print("Verdict : EN TENSION. sqrt(sigma) a 2,2 sigma d'une determination 2+1 saveurs au point physique ; la coupe")
    print("QCD/QED de la base suit le reseau (0,2 sigma) et contredit les determinations dispersives (3,7 et ~8 sigma), qui")
    print("contredisent aussi le reseau. Critere precise : une valeur reseau au point physique de la partie QCD sous")
    print("2,2 MeV a +- 0,1, ou une confirmation reseau du 0,71 dispersif, tue la formule sans nombre libre ; une seconde")
    print("determination de sqrt(sigma) au point physique hors de 427-435 de plus de 3 sigma tue l'autre.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
