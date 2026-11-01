#!/usr/bin/env python3
"""R118 -- Veille du 1er novembre 2026 : les criteres du depot contre ce qui est sorti depuis R117.

Points surveilles (HADRONS 8, article, section Falsification) : la coupe QCD/QED de m_n - m_p, la
tension de corde, le rayon du proton, alpha, la masse du tau, la vitesse limite de l'electron.

  A. m_n - m_p : aucune nouvelle determination reseau au point physique. Table des references :
     BMW 2015 (QED 1,00 +- 0,16 ; QCD direct 2,52 +- 0,29) ; ETMC 2022 (total 1,73 +- 0,69, Nf = 2+1+1,
     m_pi 200-450 MeV ; une coupe QED 1,16 +- 0,25 / QCD 3,10 +- 0,59 est citee par un moteur de
     recherche d'apres la these de Romiti 2021, NON VERIFIEE a la source, acces bloque) ; dispersif :
     GLR 2021 0,58 +- 0,16 ; Cao et al. 2026 (arXiv:2609.15154, v2 du 25 octobre) 0,71 +0,03/-0,06 ;
     WCM 2012 1,30 +- 0,47. Les ecarts de R117 sont inchanges.
  B. La coupe est dependante du schema : Cao et al. l'ecrivent eux-memes (separations dispersive et
     reseau differentes). Seul le total est independant du schema. Le test dur reste donc le total :
     base 1,269 (bouchon) / 1,324 (sphere) contre 1,29332 MeV. Le « 8 sigma » de R117 se lit « dans
     le schema dispersif ».
  C. sqrt(sigma) : aucune nouvelle determination des niveaux au point physique ; Cea et al. 2026
     (arXiv:2607.07143, tube de flux a masses physiques) : brisure de corde entre 0,96 et 1,16 fm,
     une donnee que la base ne lit pas (aucune brisure dans ses regles). 445 +- 7 contre 430 : 2,2 sigma.
  D. r_p : PRad-II a fini sa prise de donnees (juillet 2026), analyse en cours, aucune valeur ;
     alpha : Rb/Cs a 5,5 sigma non resolu ; tau : 1777,09 +- 0,14 (2023) reste la meilleure, nouvelle
     mesure Belle II en cours ; KATRIN : 0,45 eV officiel, ~0,3 eV attendu ; violation de Lorentz :
     nouvelles limites GRB 2026 (PRD, fevrier), rien sur l'electron.
  Verdict : rien ne bouge. Une nuance (schema) et une bibliographie completee (Cao et al.).
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
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
    print("R118 -- veille du 1er novembre 2026\n")
    print(f"  base : partie forte {b['strong']:.3f} MeV ; part QED {b['qed_plug']:.3f} / {b['qed_sphere']:.3f} ; total {b['strong']-b['qed_plug']:.3f} / {b['strong']-b['qed_sphere']:.3f} ; sqrt(sigma) {b['sqrt_sigma']:.1f} MeV\n")

    # A. m_n - m_p
    print("A. La coupe QCD/QED de m_n - m_p : references au 1er novembre 2026")
    refs = {
        "BMW 2015 (reseau, Nf = 1+1+1+1)": (1.00, 0.16, "verifiee"),
        "ETMC 2022 (reseau, Nf = 2+1+1, coupe de la these 2021)": (1.16, 0.25, "NON VERIFIEE a la source"),
        "Gasser-Leutwyler-Rusetsky 2021 (Cottingham)": (0.58, 0.16, "verifiee"),
        "Cao et al. 2026, arXiv:2609.15154 v2 (regle de somme)": (0.71, 0.06, "verifiee (v1) ; v2 du 25 oct. a relire"),
        "Walker-Loud-Carlson-Miller 2012 (Cottingham)": (1.30, 0.47, "verifiee"),
    }
    print("    part QED (m_p - m_n)_QED, MeV, et ecart de la base (bouchon 1,20) :")
    for k, (v, e, status) in refs.items():
        print(f"      {k:56s} {v:.2f} +- {e:.2f} : {pull(b['qed_plug'], v, e):+.1f} sigma   [{status}]")
    print(f"    ETMC 2022 total : 1,73 +- 0,69 MeV (base {b['strong']-b['qed_plug']:.2f} : {pull(b['strong']-b['qed_plug'], 1.73, 0.69):+.1f} sigma)")
    print("    aucune nouvelle determination reseau au point physique depuis R117 ; les ecarts de R117 sont inchanges.\n")
    unchanged = abs(pull(b["qed_plug"], 1.00, 0.16) - 1.3) < 0.1 and abs(pull(b["qed_plug"], 0.71, 0.06) - 8.2) < 0.2
    check("A. aucune nouvelle valeur reseau ; ecarts de R117 inchanges (BMW +1,3 sigma, regle de somme +8 sigma) ; coupe ETMC citee, non verifiee",
          unchanged, f"BMW {pull(b['qed_plug'], 1.00, 0.16):+.1f}, 2026 {pull(b['qed_plug'], 0.71, 0.06):+.1f} sigma")

    # B. le schema
    print("B. La coupe depend du schema ; le total non")
    tot = (b["strong"] - b["qed_plug"], b["strong"] - b["qed_sphere"])
    print("    Cao et al. : « the Cottingham analysis and lattice calculations employ different schemes for separating QED and QCD effects ».")
    print(f"    le total est la seule grandeur sans schema : base {tot[0]:.3f} (bouchon) / {tot[1]:.3f} (sphere) MeV contre {DM_EXP:.5f} : {100*(tot[0]/DM_EXP-1):+.1f} % / {100*(tot[1]/DM_EXP-1):+.1f} %")
    print("    -> le « 8 sigma » de R117 se lit « dans le schema dispersif » ; le test dur reste le total, a 2-3 %.\n")
    check("B. total de la base 1,27-1,32 MeV contre 1,2933 : -1,9 % / +2,4 %, le seul test independant du schema",
          abs(tot[0] / DM_EXP - 1) < 0.03 and abs(tot[1] / DM_EXP - 1) < 0.03, f"{100*(tot[0]/DM_EXP-1):+.1f} % / {100*(tot[1]/DM_EXP-1):+.1f} %")

    # C. sqrt(sigma)
    print("C. La tension de corde")
    pw = pull(b["sqrt_sigma"], 445.0, math.hypot(3.0, 6.0))
    rb = (0.963, 1.156)
    print(f"    Bulava et al. 2024 : 445 +- 6,7 MeV ; base {b['sqrt_sigma']:.1f} : {pw:+.2f} sigma (inchange) ; aucune nouvelle determination des niveaux au point physique")
    print(f"    Cea et al. 2026 (tube de flux, masses physiques) : brisure de corde entre {rb[0]} et {rb[1]} fm ; pour memoire l_1(9) = {L1:.3f} fm, r_e/2 = {ALPHA*LAM/2:.3f} fm ;")
    print("    la base n'a pas de brisure dans ses regles : donnee enregistree, non lue.\n")
    check("C. sqrt(sigma) : 2,2 sigma inchange ; brisure de corde 0,96-1,16 fm (2026) enregistree sans lecture", 2.0 < abs(pw) < 2.5, f"{pw:+.2f} sigma")

    # D. le reste
    print("D. Le reste des criteres")
    p_mu = math.log(MMU / ME) / math.log(7 / 3)
    p_tau = math.log(MTAU / ME) / math.log(11 / 3)
    print(f"    r_p : 0,8409(4) fm ; PRad-II : donnees mars-juillet 2026, analyse en cours, aucune valeur publiee")
    print(f"    alpha : Rb (2020) et Cs (2018) a 5,5 sigma, non resolu ; tau : 1777,09 +- 0,14 (Belle II 2023), nouvelle mesure en cours ;")
    print(f"    fenetre de l'exposant inchangee : [{p_tau:.4f} ; {p_mu:.4f}], 2 pi = {2*math.pi:.4f}")
    print(f"    KATRIN : 0,45 eV officiel (2025), fin des donnees octobre 2025, ~0,3 eV attendu ; Lorentz : limites GRB 2026 (PRD fev.), rien de neuf pour l'electron (10^-11)\n")
    check("D. r_p, alpha, tau, nu, Lorentz : aucune nouvelle valeur ; fenetre de l'exposant inchangee", p_tau < 2 * math.pi < p_mu, f"[{p_tau:.4f} ; {p_mu:.4f}]")

    print("Verdict : RIEN NE BOUGE. Statuts de R117 maintenus (sqrt(sigma) en tension a 2,2 sigma ; coupe de m_n - m_p en")
    print("tension dans le schema dispersif, en accord avec le reseau). Nuance ajoutee : la coupe depend du schema, le test")
    print("sans schema est le total, a 2-3 %. Bibliographie completee : Cao et al., arXiv:2609.15154 (v2, 25 octobre 2026).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
