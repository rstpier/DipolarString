#!/usr/bin/env python3
"""R65 -- "C'est le cas le plus probable [un defaut fixe, R64]. Sauf si resonance
a Z != Z0."

Alternative testee : la particule instable n'est pas un circuit adapte avec un
petit defaut, mais une resonance dont l'impedance propre differe de Z0.  Quatre
lectures d'un Z != Z0, chacune confrontee au nombre de tours que vit le muon
(5,6e16) et l'electron (> 2,6e56).
"""
import math

HBARC = 197.3269804
C_FM = 2.99792458e23
HBAR_MEV_S = 6.582119569e-22
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
W_E = 4 * LAMBDA_BAR / math.pi ** 2          # 156.5 fm, R17

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

MMU, TAU_MU = 105.6583755, 2.1969811e-6
MDELTA, GAMMA_DELTA = 1232.0, 117.0
TAU_E_BOUND = 6.6e28 * 3.156e7               # s (Borexino, > 6.6e28 ans)
ZE_V29 = 0.73                                # Z_e / Z0 calibre dans v2.9 (retire)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R65 -- resonance a Z != Z0 ?\n")
    T_mu = 3 * l1(7) / C_FM
    N_mu = TAU_MU / T_mu
    leak_mu = 1 / N_mu
    T_e = 3 * L1_3 / C_FM
    N_e = TAU_E_BOUND / T_e
    print(f"  muon : {N_mu:.2e} tours, fuite/tour {leak_mu:.2e}")
    print(f"  electron : > {N_e:.1e} tours, fuite/tour < {1/N_e:.1e}\n")

    # 1. lecture transmission : resonateur ferme sur le milieu par sa desadaptation
    #    fuite/tour = 4 Z Z0 / (Z + Z0)^2  (petite seulement si Z << Z0 ou Z >> Z0)
    print("1. Resonance confinee par la desadaptation (fuite = 4ZZ0/(Z+Z0)^2)")
    z_ratio = leak_mu / 4                    # Z/Z0 << 1
    print(f"  muon : il faut Z/Z0 = {z_ratio:.1e} ou {1/z_ratio:.1e}")
    for zr in (0.5, 0.73, 0.1, 0.01):
        f = 4 * zr / (1 + zr) ** 2
        print(f"  Z/Z0 = {zr:5.2f} : fuite/tour {f:.3f}, vie {1/f:.0f} tours")
    print("  -> aucune impedance n'est a 10^17 de Z0 : exclue.\n")
    check("lecture transmission : Z/Z0 requis < 1e-15 (aucune impedance)",
          z_ratio < 1e-15, f"Z/Z0 = {z_ratio:.1e}")

    # 2. lecture reflexion : petite desadaptation, fuite = (dZ/2Z0)^2
    print("2. Petite desadaptation (fuite = (dZ/2Z0)^2)")
    dz_mu = 2 * math.sqrt(leak_mu)
    dz_e = 2 * math.sqrt(1 / N_e)
    print(f"  muon : dZ/Z0 = {dz_mu:.1e} ; electron : dZ/Z0 < {dz_e:.1e}")
    print("  -> a 1e-8 pres c'est Z0 : c'est le cas du defaut fixe de R64, pas une "
          "autre impedance.\n")
    check("lecture reflexion : dZ/Z0 du muon < 1e-7 (= Z0, cas R64)",
          dz_mu < 1e-7, f"{dz_mu:.1e}")
    check("l'electron est adapte a mieux que 1e-27",
          dz_e < 1e-27, f"dZ/Z0 < {dz_e:.1e}")

    # 3. la calibration v2.9 : Z_e = 0.73 Z0
    print("3. La calibration de v2.9, Z_e = 0,73 Z0, lue comme impedance de l'anneau")
    f_v29 = ((1 - ZE_V29) / (1 + ZE_V29)) ** 2
    tau_v29 = T_e / f_v29
    print(f"  fuite/tour {f_v29:.4f}, vie {1/f_v29:.0f} tours = {tau_v29:.1e} s "
          f"(borne mesuree > {TAU_E_BOUND:.0e} s)")
    print("  -> Z_e = 0,73 Z0 ne peut pas etre l'impedance du circuit : "
          "l'adaptation exacte de R10 est requise.\n")
    check("Z_e = 0.73 Z0 tuerait l'electron en < 1e-15 s",
          tau_v29 < 1e-15, f"{tau_v29:.1e} s")

    # 4. resonance ouverte (onde stationnaire) : rayonne
    print("4. Resonance ouverte (onde stationnaire de charge e, R47) : temps de Larmor")
    tau_rad_mu = 3 * HBAR_MEV_S / (ALPHA * MMU)          # 3 hbar / (alpha m c^2)
    tau_rad_delta = 3 * HBAR_MEV_S / (ALPHA * MDELTA)
    tau_delta = HBAR_MEV_S / GAMMA_DELTA
    print(f"  muon : rayonne en {tau_rad_mu:.1e} s, vit {TAU_MU:.1e} s "
          f"(rapport {TAU_MU/tau_rad_mu:.0e})")
    print(f"  Delta : rayonnerait en {tau_rad_delta:.1e} s, vit {tau_delta:.1e} s "
          f"(se defait {tau_rad_delta/tau_delta:.0f} x plus vite : rupture, pas rayonnement)")
    print("  -> une resonance ouverte ne vit pas plus de 1e-21 s : exclue pour mu, tau, n, pi.\n")
    check("resonance ouverte : le muon vit > 1e10 fois le temps de Larmor",
          TAU_MU / tau_rad_mu > 1e10, f"rapport {TAU_MU/tau_rad_mu:.0e}")

    # 5. lecture dispersion : Z != Z0 => v != c, le paquet glisse et cesse de se fermer
    print("5. Dispersion (Z != Z0 => v != c) : le paquet glisse d'un ruban en N tours")
    w7 = W_E * l1(7) / L1_3
    eps = w7 / (3 * l1(7) * N_mu)
    print(f"  muon : 1 - v/c = {eps:.1e}")
    print("  -> encore Z0 a 1e-18 pres.\n")
    check("lecture dispersion : 1 - v/c requis < 1e-10", eps < 1e-10, f"{eps:.1e}")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
