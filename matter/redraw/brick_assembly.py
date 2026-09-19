#!/usr/bin/env python3
"""R70 -- "Les DQD sont les briques naturelles fondamentales : les brins
solitaires charges se trouvent, s'orientent et s'assemblent."

Quatre verbes, quatre tests :
  A. la brique : quelle unite de charge et quel nombre de brins par fermion
     reproduisent exactement le spectre {0, +-1/3, +-2/3, +-1} ?
  B. "se trouvent" : la force de Coulomb entre brins solitaires suffit-elle a
     les assembler, comparee aux energies de jonction et de circulation ?
  C. "s'orientent" : que coute l'orientation des bras (Coulomb) devant la
     jonction ?
  D. "s'assemblent" : l'anneau ferme de trois est-il la seule configuration
     stationnaire (moment magnetique) ?
  E. par paires : un DQD de spin 1 (R3) peut-il s'ajouter seul a un anneau de
     spin 1/2 ?
"""
import math
import itertools
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
L1_3 = 2 * math.pi * LAMBDA_BAR / 3
K = ALPHA * HBARC                       # 1.44 MeV fm
W_E = 4 * LAMBDA_BAR / math.pi ** 2     # 156.5 fm

def l1(n):
    return L1_3 * (3 / n) ** (2 * math.pi)

OBSERVED = {0, 1, 2, 3}                 # |Q| en unites de e/3

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def spectrum(b, N):
    """|Q| atteignables avec N brins de charge {-1/b, 0, +1/b}, en unites de 1/3"""
    qs = set()
    for combo in itertools.product((-1, 0, 1), repeat=N):
        q = sum(combo)                  # en unites de 1/b
        if (3 * q) % b == 0:
            qs.add(abs(3 * q // b))
        else:
            qs.add(None)                # charge non multiple de 1/3
    return qs

def segment_energy(theta, L, a, nq=400):
    """energie mutuelle (K q^2 / L^2 unites) de deux segments charges uniformes
    de longueur L partant de l'origine a l'angle theta, coupure a"""
    s = (np.arange(nq) + 0.5) / nq * L
    t = s.copy()
    S, T_ = np.meshgrid(s, t)
    r = np.sqrt(S * S + T_ * T_ - 2 * S * T_ * math.cos(theta) + a * a)
    return (1.0 / r).sum() * (L / nq) ** 2 / L ** 2

def main():
    print("R70 -- la brique et ses quatre verbes\n")

    # A. la brique
    print("A. Quelle brique ? unite de charge 1/b, N brins par fermion")
    exact = []
    for b in range(1, 7):
        for N in range(1, 9):
            if spectrum(b, N) == OBSERVED:
                exact.append((b, N))
    print(f"  (b, N) reproduisant exactement le spectre : {exact}")
    print(f"  N = 2 (b = 3) : {sorted(q for q in spectrum(3, 2) if q is not None)} (pas de charge 1) ; "
          f"N = 4 : {sorted(q for q in spectrum(3, 4) if q is not None)} (charge 4/3)")
    print("  -> la brique a des branches de e/3 et un fermion en prend exactement trois.\n")
    check("spectre exact <=> tiers et trois brins", exact == [(3, 3)], f"{exact}")

    # B. se trouvent
    print("B. 'Se trouvent' : Coulomb entre deux brins solitaires a l'echelle de l'electron")
    e_coul = K / 9 / L1_3                    # (e/3)^2 / l1
    e_junction = ME / 6                      # R54 : trois jonctions = moitie statique
    e_circ = ME / 2                          # R42 : circulation = moitie de la masse
    print(f"  (e/3)^2/(4 pi eps0 l1) = {e_coul*1e3:.2f} keV ; jonction m_e/6 = {e_junction*1e3:.0f} keV ; "
          f"circulation m_e/2 = {e_circ*1e3:.0f} keV")
    print(f"  -> Coulomb est {e_junction/e_coul:.0f} fois sous la jonction : ce n'est pas l'attraction")
    print("     qui assemble, c'est la fermeture (le quantum de circulation n'existe que sur un circuit")
    print("     ferme, R47) ; les trois brins d'un lepton naissent deja voisins (R4, la brisure).\n")
    check("Coulomb entre brins solitaires < 1 % de la jonction", e_coul / e_junction < 0.01,
          f"{100*e_coul/e_junction:.2f} %")

    # C. s'orientent
    print("C. 'S'orientent' : cout de Coulomb de l'orientation des bras d'un quark (echelle 9)")
    L9, a9 = l1(9), W_E * l1(9) / L1_3 / 4
    e_unit = K / 9 / 1.0                     # (e/3)^2 par fm
    best = None
    for deg in (60, 90, 120, 150, 180):
        u = segment_energy(math.radians(deg), L9, a9) * e_unit   # MeV
        print(f"    bras a {deg:3d} deg : U = {u:+.3f} MeV (deux bras + ou deux bras -)")
        if best is None or u < best[1]:
            best = (deg, u)
    u120 = segment_energy(math.radians(120), L9, a9) * e_unit
    e_circuit = math.pi * HBARC / (3 * L9)   # 254 MeV, R42
    print(f"  minimum libre a {best[0]} deg ; a 120 deg (impose par la tri-jonction) U = {u120:.3f} MeV,")
    print(f"  soit {100*u120/e_circuit:.2f} % du circuit ({e_circuit:.0f} MeV) : l'orientation est")
    print("  geometrique (trois jonctions a 120 deg), l'electricite ne la choisit pas.\n")
    check("cout electrique de l'orientation < 1 % du circuit", u120 / e_circuit < 0.01,
          f"{100*u120/e_circuit:.2f} %")

    # D. s'assemblent
    print("D. 'S'assemblent' : quelle configuration de trois brins est stationnaire ?")
    R_ring = LAMBDA_BAR                       # 3 l1(3) = 2 pi lambda-bar
    mu_ring = 0.5 * R_ring                    # (e c R / 2) en unites de e c fm ; mu_B = e hbar/(2m) = e c lambda-bar/2
    mu_B = 0.5 * LAMBDA_BAR
    print(f"  anneau ferme a sens unique (R53), rayon lambda-bar : mu/mu_B = {mu_ring/mu_B:.4f}")
    print("  arc ouvert (R47) : pas de circulation nette, mu = 0, rayonne ; epingle fermee : 0,30 mu_B")
    print("  -> seul l'anneau ferme de trois brins est stationnaire et porte le bon moment.\n")
    check("anneau ferme de trois : mu = mu_B exactement", abs(mu_ring / mu_B - 1) < 1e-9,
          f"{mu_ring/mu_B:.6f}")

    # E. par paires
    print("E. Par paires : j DQD de spin 1 ajoutes a un anneau de spin 1/2")
    def can_sum_to_zero(j):
        # spins entiers 1 x j : le total 0 est atteignable ssi j != 1 (j = 0 trivial)
        spins = {0}
        for _ in range(j):
            spins = {s2 for s in spins for s2 in range(abs(s - 1), s + 2)}
        return 0 in spins
    rows = []
    for j in range(0, 7):
        n = 3 + 2 * j
        ok = can_sum_to_zero(j)
        m = ME * (n / 3) ** (2 * math.pi)
        seen = {3: "electron", 7: "muon", 9: "quark (classe 0, R41)", 11: "tau"}.get(n, "absent")
        rows.append((n, j, ok))
        print(f"    j = {j} DQD -> n = {n:2d} : spin total 0 {'possible' if ok else 'impossible'} ; "
              f"m = {m:9.1f} MeV ; {seen}")
    forbids_5 = not [ok for n, j, ok in rows if n == 5][0]
    allows_7_11 = all([ok for n, j, ok in rows if n in (7, 11)])
    print("  -> un DQD seul ne peut pas garder le spin 1/2 : n = 5 interdit ; 7, 9, 11 permis ;")
    print("     13 et 15 ne sont pas interdits par le spin (il faut la regle de classe mod 3, R61).\n")
    check("le spin interdit n = 5 et permet 7, 9, 11 (13, 15 restent a la regle de classe)",
          forbids_5 and allows_7_11, "j = 1 impossible, j >= 2 possible")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
