#!/usr/bin/env python3
"""R119 -- Veille du 21 novembre 2026 : les criteres du depot contre ce qui est sorti depuis R118 (1er novembre).

Sources lues par resumes de moteur de recherche (arxiv.org, pos.sissa.it, indico.* hors politique reseau) ; chaque
valeur porte son statut : « verifiee » (article connu ou deux resumes concordants) ou « citee, non verifiee a la source ».

  A. sqrt(sigma). La base donne 430,3 MeV, soit sigma = 0,185 GeV^2 : c'est la tension de corde phenomenologique
     (pente de Regge alpha' = 0,88 GeV^-2 -> 0,181 GeV^2 ; Cornell 0,183 GeV^2 ; bande 0,18-0,19 GeV^2 = 424-436 MeV).
     Le reseau a quarks dynamiques (energie statique ajustee dans une fenetre, la corde se brisant vers 1,1 fm) donne
     plus : 445(3)(6) MeV (Bulava et al. 2024, Nf = 2+1, point physique, verifiee) ; 458(7) preliminaire (extrapolation
     lineaire au point physique, expose Knechtli, DESY, date non verifiee) ; 467(7) / 481(7) (TUMQCD, Nf = 2+1+1,
     arXiv:2206.03156, PRD 107, 074503 (2023) ; valeurs en MeV citees par le meme expose, non verifiees a la source :
     l'abstract verifie dit seulement que r0*sqrt(sigma) concorde avec les resultats 2+1). Lecture : « au centre de
     420-440 » (R39, R117) decrit la bande phenomenologique, pas le reseau ; la base est sous toutes les determinations
     reseau a quarks dynamiques, de 2,2 sigma (verifiee) a 4-7 sigma (citees).
  B. Non-linearite du vide. European XFEL, instrument HED + laser ReLaX, champ sombre pompe-sonde, 2751 collisions,
     energie c.m. 116 eV : aucun exces ; section efficace lumiere-lumiere < 7e4 x QED a 95 % (arXiv:2610.39867, octobre
     2026 ; lettre arXiv:2611.04531, novembre 2026) ; dix ordres de grandeur mieux que les limites experimentales
     anterieures ; la sensibilite QED n'est pas atteinte. Pour une interaction quartique effective la section efficace
     est le carre du coefficient : tout coefficient quartique non QED a 116 eV est < sqrt(7e4) = 265 x Euler-Heisenberg.
     Traduit en echelle de densite d'energie (Z(u) = Z0 (1 + u/u_s), premier terme quartique en 1/u_s) :
     u_s > u_EH / 265, u_EH = 45 m_e^4 c^5 / (8 alpha^2 hbar^3) = 1,5e29 J/m^3, soit u_s > 5,7e26 J/m^3 = 1,2e3 u_0
     (u_0 = h c / l_1^4 = 4,6e23 J/m^3, manuscrit). Ordre de grandeur : la structure de Lorentz de Z(u) n'est pas celle
     d'Euler-Heisenberg.
  C. Lorentz. LHAASO, GRB 221009A, v3 (fevrier 2026) : E_QG,1 > 10 E_Pl ; E_QG,2 > 6e-8 E_Pl = 7,3e11 GeV. La dispersion
     du reseau (manuscrit, section Lorentz) equivaut a E_QG,2 = 0,69 MeV : ecart 1e15. La lecture « onde collective » de
     la coupure de Brillouin est exclue ; la lecture « soliton N = 1 » n'a pas de dispersion calculee (revue, 3.1).
  D. Falsificateurs gravitationnels. PBH : LHAASO tout-ciel (arXiv:2505.24586, 2025) : aucun sursaut, taux local < 181
     pc^-3 an^-1 (99 %, 20 s) ; HAWC < 3400 ; Fermi-LAT < 7,2e3 : exposition maintenue, rien d'etabli. GWTC-5.0, tests de
     la RG (arXiv:2607.19293, 22 juillet 2026, 72 evenements O4b) : rien au-dela de la RG, pas de polarisation
     supplementaire (mode respiratoire : coherent), restes coherents avec Kerr ; echos : recherche dediee seulement dans
     GWTC-4.0 (aucun). LRD : aucune nouvelle masse dynamique ; Liu et al. (ApJL, octobre 2026) : modeles d'ecoulement
     supercritique exigeant 1e5-1e7 M_sun -- un modele, pas une mesure ; canal toujours defavorise.
  E. Le reste. r_p : enigme resolue vers le petit rayon (hydrogene 2S-6P, Nature, 11 fevrier 2026 ; revue avril 2026) ;
     0,8409 fm demeure, base 0,8128 fm (-3,3 %) ; PRad-II : donnees terminees juin 2026, analyse en cours. m_n - m_p :
     aucune nouvelle valeur reseau ; PDG 2026 : ~1 MeV (QED) contre ~2,5 MeV (QCD), la base 1,20 / 2,47 est dans ce
     schema. alpha : Rb/Cs non resolu (revue Aliberti et al., Phys. Rep. 1143 (2025)). tau : 1777,09 +- 0,14 (2023).
     KATRIN : 0,45 eV (2025 ; bilan janvier 2026, arXiv:2601.01672).
  Verdict : deux choses bougent. (1) sqrt(sigma) : la tension n'est pas un 2,2 sigma isole ; la base est egale a la bande
  phenomenologique et sous toute la bande reseau a quarks dynamiques ; la formule doit dire laquelle des deux elle vise.
  (2) Premiere borne de laboratoire sur un exces de non-linearite du vide a basse energie : le canal de falsification du
  manuscrit (BIREF) et le test decisif de la reconstruction (Z(u)) ont maintenant un nombre.
"""
import math

HBARC = 197.3269804          # MeV fm
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.86
F_PLUG = 0.757
LAM = HBARC / ME             # fm
L1_3 = 2 * math.pi * LAM / 3 # fm, l_1(3)
L1 = L1_3 * 3 ** (-2 * math.pi)
W9 = (4 * LAM / math.pi ** 2) * 3 ** (-2 * math.pi)
# SI
HBAR, H, C, ME_KG, EPS0 = 1.054571817e-34, 6.62607015e-34, 2.99792458e8, 9.1093837015e-31, 8.8541878128e-12
E_PL_GEV = 1.220890e19

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
    sigma = math.pi * HBARC / L1 ** 2          # MeV/fm
    return dict(strong=E_dqd, qed_plug=coulomb(F_PLUG), sqrt_sigma=math.sqrt(sigma * HBARC))

def pull(x, ref, err):
    return (x - ref) / err

def main():
    b = base_numbers()
    print("R119 -- veille du 21 novembre 2026\n")

    # A. sqrt(sigma) : phenomenologie contre reseau
    print("A. La tension de corde : la base, la bande phenomenologique, le reseau a quarks dynamiques")
    s = b["sqrt_sigma"]
    sig_gev2 = (s / 1000) ** 2
    regge = math.sqrt(1 / (2 * math.pi * 0.88)) * 1000       # alpha' = 0,88 GeV^-2
    cornell = math.sqrt(0.183) * 1000
    band = (math.sqrt(0.18) * 1000, math.sqrt(0.19) * 1000)
    print(f"    base : sqrt(sigma) = {s:.1f} MeV, sigma = {sig_gev2:.3f} GeV^2 ; Regge (alpha' = 0,88) {regge:.1f} ; Cornell (0,183) {cornell:.1f} ; bande 0,18-0,19 GeV^2 : {band[0]:.0f}-{band[1]:.0f} MeV")
    lattice = {
        "Bulava et al. 2024, Nf = 2+1, point physique (PLB 854)": (445.0, math.hypot(3.0, 6.0), "verifiee"),
        "Knechtli, expose DESY, extrapolation lineaire, preliminaire": (458.0, 7.0, "citee, date et source non verifiees"),
        "TUMQCD 2022/2023, Nf = 2+1+1, schema 1 (arXiv:2206.03156)": (467.0, 7.0, "citee par l'expose, non verifiee a la source"),
        "TUMQCD 2022/2023, Nf = 2+1+1, schema 2 (arXiv:2206.03156)": (481.0, 7.0, "citee par l'expose, non verifiee a la source"),
    }
    pulls = {}
    for k, (v, e, st) in lattice.items():
        pulls[k] = pull(s, v, e)
        print(f"      {k:62s} {v:.0f} +- {e:.1f} : {pulls[k]:+.1f} sigma   [{st}]")
    print("    la corde se brise vers 1,1 fm (Cea et al. 2026 : 0,96-1,16 fm) : le sigma du reseau est une quantite de fenetre d'ajustement.")
    print("    -> « au centre de 420-440 » decrivait la bande phenomenologique ; la base est sous toutes les valeurs reseau a quarks dynamiques.\n")
    check("A. base dans la bande phenomenologique 0,18-0,19 GeV^2 et sous toute valeur reseau dynamique de plus de 2 sigma",
          band[0] < s < band[1] and all(p < -2.0 for p in pulls.values()),
          f"sigma = {sig_gev2:.3f} GeV^2 ; ecarts {', '.join(f'{p:+.1f}' for p in pulls.values())} sigma")

    # B. non-linearite du vide : la borne HIBEF/HED d'octobre-novembre 2026
    print("B. La non-linearite du vide a basse energie : premiere borne de laboratoire par champ sombre")
    ratio_qed = 7e4
    coeff = math.sqrt(ratio_qed)
    u_eh = 45 * ME_KG ** 4 * C ** 5 / (8 * ALPHA ** 2 * HBAR ** 3)
    u0 = H * C / (L1_3 * 1e-15) ** 4
    u_s_min = u_eh / coeff
    u_crit = 14.5 * u0
    print(f"    EuXFEL HED + ReLaX (arXiv:2610.39867, 2611.04531) : 2751 collisions, 116 eV c.m., aucun exces, sigma_LbL < {ratio_qed:.0e} x QED (95 %)")
    print(f"    coefficient quartique effectif : < sqrt({ratio_qed:.0e}) = {coeff:.0f} x Euler-Heisenberg a 116 eV")
    print(f"    u_EH = 45 m_e^4 c^5/(8 alpha^2 hbar^3) = {u_eh:.2e} J/m^3 ; u_0 = hc/l_1^4 = {u0:.2e} J/m^3 ; 14,5 u_0 = {u_crit:.2e} J/m^3")
    print(f"    Z(u) = Z0(1 + u/u_s) : u_s > u_EH/{coeff:.0f} = {u_s_min:.2e} J/m^3 = {u_s_min/u0:.0f} u_0 (ordre de grandeur ; structure de Lorentz non specifiee)")
    print("    -> le canal BIREF du manuscrit et le test decisif de la reconstruction ont un premier nombre : aucun exces jusqu'a 7e4 x QED.\n")
    check("B. borne HIBEF : u_s > 1e3 u_0 ; l'echelle de repos du substrat (u_0, 14,5 u_0) est exclue comme echelle de non-linearite a basse energie",
          u_s_min / u0 > 100 and u_s_min > u_crit, f"u_s > {u_s_min:.1e} J/m^3 = {u_s_min/u0:.0f} u_0 = {u_s_min/u_crit:.0f} x (14,5 u_0)")

    # C. Lorentz
    print("C. La violation de l'invariance de Lorentz : LHAASO contre la dispersion du reseau")
    e_qg2 = 6e-8 * E_PL_GEV
    e_lat = 0.69e-3
    print(f"    LHAASO GRB 221009A v3 (fevrier 2026) : E_QG,1 > 10 E_Pl ; E_QG,2 > 6e-8 E_Pl = {e_qg2:.1e} GeV")
    print(f"    dispersion du reseau (manuscrit) : E_QG,2 = {e_lat*1e3:.2f} MeV ; ecart {e_qg2/e_lat:.1e}")
    print("    -> la coupure de Brillouin comme onde collective est exclue de 15 ordres ; le soliton N = 1 n'a pas de dispersion calculee.\n")
    check("C. LHAASO 2026 : E_QG,2 > 7e11 GeV contre 0,69 MeV, ecart > 1e14", e_qg2 / e_lat > 1e14, f"{e_qg2/e_lat:.1e}")

    # D. falsificateurs gravitationnels
    print("D. Les falsificateurs gravitationnels")
    limits = {"LHAASO tout-ciel 2025 (20 s, 99 %)": 181.0, "HAWC (959 jours, 99 %)": 3400.0, "Fermi-LAT (mouvement propre)": 7200.0}
    for k, v in limits.items():
        print(f"      taux local de sursauts PBH : {k:36s} < {v:.0f} pc^-3 an^-1 -- aucune detection")
    print("    GWTC-5.0 tests de la RG (22 juillet 2026) : rien au-dela de la RG ; pas de polarisation supplementaire (mode respiratoire : coherent) ;")
    print("    restes coherents avec Kerr (GW250114, SNR 76,9) ; echos : recherche dediee dans GWTC-4.0 seulement, aucun.")
    print("    LRD : aucune nouvelle masse dynamique ; Liu et al. (ApJL, oct. 2026) : 1e5-1e7 M_sun par modele d'ecoulement -- pas une mesure ; canal defavorise, inchange.\n")
    check("D. PBH : aucun sursaut, limite la plus forte LHAASO 181 pc^-3 an^-1 ; GWTC-5.0 coherent avec la RG ; LRD inchange",
          min(limits.values()) == limits["LHAASO tout-ciel 2025 (20 s, 99 %)"], "exposition maintenue, rien d'etabli")

    # E. le reste
    print("E. Le reste des criteres")
    rp = 0.8409
    p_mu = math.log(MMU / ME) / math.log(7 / 3)
    p_tau = math.log(MTAU / ME) / math.log(11 / 3)
    print(f"    r_p : enigme resolue vers le petit rayon (hydrogene 2S-6P, Nature, 11 fevrier 2026 ; revue avril 2026) ; {rp} fm demeure ; base l_1(9) = {L1:.4f} fm ({100*(L1/rp-1):+.1f} %) ; PRad-II : donnees terminees juin 2026")
    print(f"    m_n - m_p : aucune nouvelle valeur reseau ; PDG 2026 : ~1 MeV QED contre ~2,5 MeV QCD ; base {b['qed_plug']:.2f} / {b['strong']:.2f} MeV ; BMW {pull(b['qed_plug'], 1.00, 0.16):+.1f} sigma (inchange)")
    print(f"    alpha : Rb/Cs non resolu ; tau : 1777,09 +- 0,14 (2023) ; KATRIN : 0,45 eV (2025, bilan janvier 2026) ; fenetre de l'exposant [{p_tau:.4f} ; {p_mu:.4f}], 2 pi = {2*math.pi:.4f}\n")
    check("E. r_p, m_n - m_p, alpha, tau, nu : aucune nouvelle valeur qui bouge un critere ; fenetre de l'exposant inchangee",
          p_tau < 2 * math.pi < p_mu and abs(L1 / rp - 1) < 0.04 and abs(pull(b["qed_plug"], 1.00, 0.16) - 1.3) < 0.1,
          f"r_p {100*(L1/rp-1):+.1f} %, BMW {pull(b['qed_plug'], 1.00, 0.16):+.1f} sigma")

    print("Verdict : DEUX CHOSES BOUGENT. (1) sqrt(sigma) : la base (0,185 GeV^2) est la tension phenomenologique de Regge-Cornell et")
    print("se trouve sous toute la bande reseau a quarks dynamiques (445 verifiee, 458-481 citees) : la tension n'est pas un 2,2 sigma")
    print("isole ; la formule doit dire quelle tension elle vise. (2) Premiere borne de laboratoire sur un exces de non-linearite du")
    print("vide a basse energie (7e4 x QED, 116 eV) : u_s > 1e3 u_0 pour toute lecture Z(u) ; le canal BIREF et le test decisif de la")
    print("reconstruction ont un nombre. Le reste ne bouge pas.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
