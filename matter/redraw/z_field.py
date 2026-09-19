#!/usr/bin/env python3
"""R68 -- "Le champ Z ?"

Deux sens : le boson Z0 du modele standard (91 GeV, quantum du champ faible
neutre) et le champ d'impedance Z de la base (Z0, Z(r), Z(I) de R66).  On
verifie si l'echelle manquante de R64-R67 (M = 1,62 TeV) est celle du champ
faible, et si la base a un objet a cette echelle.
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ALPHA_MZ = 1 / 127.95
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
K = ALPHA * HBARC                    # 1.44 MeV fm
D0 = 2 * math.cosh(math.pi) * LAMBDA_BAR / 37.1   # 241 fm, ecart du DQD (manuscrit)

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

def m_ladder(n):
    return ME * (n / 3) ** (2 * math.pi)

MMU, TAU_MU = 105.6583755, 2.1969811e-6
HBAR_MEV_S = 6.582119569e-22
GF = 1.1663788e-5 * 1e-6              # MeV^-2
MW, MZ = 80369.2, 91187.6             # MeV
SIN2_W = 0.23122
M_R64 = 1.62e6                        # MeV, echelle ajustee en R64

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R68 -- le champ Z ?\n")

    # A. le tour de l'echelle est la periode de Compton
    print("A. Le 'tour' du redessin")
    T_turn = 3 * l1(7) * 1e-15 / 2.99792458e8
    T_compton_ladder = 2 * math.pi * HBAR_MEV_S / m_ladder(7)
    T_compton = 2 * math.pi * HBAR_MEV_S / MMU
    print(f"  3 l1(n) = 2 pi lambda-bar_e m_e/m(n) : T_tour(7) = {T_turn:.4e} s, "
          f"2 pi hbar/(m(7) c^2) = {T_compton_ladder:.4e} s")
    print(f"  (avec la vraie masse du muon : {T_compton:.4e} s, l'ecart de 0,8 % de l'echelle "
          f"a Koide, R59)")
    print("  -> le tour est exactement la periode de Compton de la particule.\n")
    check("le tour est la periode de Compton (exact sur l'echelle)",
          abs(T_turn / T_compton_ladder - 1) < 1e-6, f"rapport {T_turn/T_compton_ladder:.7f}")

    # B. l'echelle manquante est la constante de Fermi
    print("B. L'echelle manquante contre la constante de Fermi")
    # Gamma_mu = G_F^2 m^5 / (192 pi^3) ; fuite/tour = Gamma T = G_F^2 m^4 / (96 pi^2)
    # = (m/M)^4  =>  M = (96 pi^2)^(1/4) / sqrt(G_F)
    M_fermi = (96 * math.pi ** 2) ** 0.25 / math.sqrt(GF)
    gamma_sm = GF ** 2 * MMU ** 5 / (192 * math.pi ** 3) / HBAR_MEV_S
    print(f"  fuite/tour = Gamma T = G_F^2 m^4/(96 pi^2) = (m/M)^4 avec "
          f"M = (96 pi^2)^(1/4)/sqrt(G_F) = {M_fermi/1e6:.3f} TeV")
    print(f"  R64 avait ajuste M = {M_R64/1e6:.2f} TeV ({100*(M_R64/M_fermi-1):+.1f} %)")
    print(f"  (controle : 1/Gamma_mu par cette formule = {1/gamma_sm*1e6:.3f} us, "
          f"mesure {TAU_MU*1e6:.4f})")
    print("  -> l'objet manquant est le champ faible, sans autre nombre que G_F.\n")
    check("M de R64 = (96 pi^2)^(1/4)/sqrt(G_F) a 1 %",
          abs(M_R64 / M_fermi - 1) < 0.01, f"{100*(M_R64/M_fermi-1):+.2f} %")

    # C. G_F en fonction du boson W : deux nombres de plus
    print("C. Ce que le champ faible demande en plus de alpha")
    gf_tree = math.pi * ALPHA / (math.sqrt(2) * MW ** 2 * SIN2_W)
    gf_run = math.pi * ALPHA_MZ / (math.sqrt(2) * MW ** 2 * SIN2_W)
    print(f"  G_F = pi alpha / (sqrt2 M_W^2 sin^2 theta_W) : {gf_tree*1e6:.3e} "
          f"({100*(gf_tree/GF-1):+.1f} %) avec alpha, {gf_run*1e6:.3e} "
          f"({100*(gf_run/GF-1):+.1f} %) avec alpha(M_Z)")
    print(f"  M_W/M_Z = cos theta_W : {MW/MZ:.4f} contre {math.sqrt(1-SIN2_W):.4f}")
    print(f"  -> deux nombres nouveaux : une masse (~80-90 GeV) et un angle (sin^2 = {SIN2_W}).\n")
    check("G_F = pi alpha/(sqrt2 M_W^2 s^2) avec alpha(M_Z) a 1 %",
          abs(gf_run / GF - 1) < 0.01, f"{100*(gf_run/GF-1):+.2f} %")

    # D. un barreau de l'echelle a M_W ou M_Z ?
    print("D. L'echelle l1(n) a-t-elle un barreau au W ou au Z ?")
    best = []
    for n in range(15, 30):
        m = m_ladder(n)
        tag = "paire (3+4k)" if (n - 3) % 4 == 0 else ""
        dw, dz = 100 * (m / MW - 1), 100 * (m / MZ - 1)
        if abs(dw) < 30 or abs(dz) < 30:
            print(f"  n = {n:2d} : {m/1e3:7.1f} GeV  (W {dw:+.1f} %, Z {dz:+.1f} %) {tag}")
        best.append(min(abs(dw), abs(dz)))
    n_w = 3 * (MW / ME) ** (1 / (2 * math.pi))
    n_z = 3 * (MZ / ME) ** (1 / (2 * math.pi))
    print(f"  n exact : W -> {n_w:.2f}, Z -> {n_z:.2f} ; les barreaux appareilles sont 19 et 23")
    print("  -> aucun barreau a moins de 4 % ; le W et le Z ne sont pas des leptons lourds.\n")
    check("aucun barreau entier a moins de 4 % de M_W ou M_Z", min(best) > 4,
          f"le plus proche : {min(best):.1f} %")

    # E. l'objet neutre de spin 1 de la base : le DQD ; a quelle echelle vaut-il 91 GeV ?
    print("E. Le seul objet neutre de spin 1 de la base est le DQD (R3, R9)")
    e_dqd_vac = 4 * math.pi * K / (9 * D0)
    e_dqd_9 = 4 * math.pi * K / (9 * l1(9))
    d_z = 4 * math.pi * K / (9 * MZ)
    print(f"  energie bifilaire (e/3)^2/(eps0 d) : {e_dqd_vac*1e3:.1f} keV a l'ecart du vide "
          f"D0 = {D0:.0f} fm, {e_dqd_9:.2f} MeV a l1(9) (R42)")
    print(f"  91 GeV demanderait un ecart d = {d_z:.1e} fm, soit {l1(9)/d_z:.0f} fois sous le quark")
    print("  -> la base n'a pas d'excitation neutre de spin 1 a 10^2 GeV ; le DQD du vide est a Z0,")
    print("     sans echelle de masse. Le 'champ Z' de la base (impedance) et le Z faible ne")
    print("     coincident que si l'impedance sature a I_c (R66) : meme mur.\n")
    check("aucune energie de la base a moins d'un facteur 10 de M_Z",
          e_dqd_9 * 10 < MZ and HBARC / (0.157) * 10 < MZ,
          f"max base ~ {HBARC/0.157:.0f} MeV (hbar c / w_9)")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
