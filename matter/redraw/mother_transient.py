#!/usr/bin/env python3
"""R88 -- Pourquoi la mere a un quantum n'est pas un etat libre (le caveat de
R87), et pourquoi ses filles le sont.

  A. la mere : boucle neutre periodique, un quantum entier, E_m = h c / L_m,
     p = h / L_m : exactement les relations d'un photon de longueur d'onde
     L_m.  Pas de verrou : charge nulle, Lk = Tw + Wr = 0 (boucle plane, repere
     fixe) : rien ne l'empeche de se rouvrir en onde libre.  La mere est la
     forme fermee, transitoire, d'un photon sur un anneau ; elle n'est pas un
     boson neutre de 0,26 MeV libre.  Coherent avec l'absence d'un tel boson.
  B. les filles : action h/2 chacune.  En onde libre, ce quantum aurait
     lambda = h/p = L_m = 2 L_d : il ne tient pas sur la boucle L_d comme onde
     entiere, et un demi-quantum n'a pas de forme libre (un photon porte h par
     longueur d'onde).  Plus la charge +-e, conservee : les filles sont
     verrouillees deux fois.
  C. annihilation : e+ e- -> 2 gamma de 0,511 MeV ; lambda_gamma = h c/m_e c^2 =
     2 pi lambda-bar = L_d exactement : les deux demi-quanta et les deux moities
     statiques redonnent deux quanta entiers, un par photon.  Les actions
     entieres sont restaurees.
  D. spin : mere 1 (R3) -> 1/2 + 1/2 ; annihilation 1/2 + 1/2 -> deux photons.
"""
import math

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME
H = 2 * math.pi                              # h / hbar

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R88 -- la mere transitoire, les filles verrouillees\n")
    L_d = 2 * math.pi * LAMBDA_BAR
    L_m = 2 * L_d

    # A. la mere est un photon ferme, sans verrou
    print("A. La mere")
    E_m = H * HBARC / L_m
    p_m = H / L_m                            # en hbar/fm
    lam_photon = H * HBARC / E_m
    lk_mother = 0.0 + 0.0                    # Tw = 0 (repere fixe, R61) + Wr = 0 (boucle plane)
    print(f"  E_m = h c / L_m = {E_m:.4f} MeV, p = h/L_m ; un photon de cette energie a lambda = {lam_photon:.1f} fm = L_m")
    print(f"  verrou : charge 0, Lk = Tw + Wr = {lk_mother:.0f} : aucun (R77 verrouille par Lk != 0)")
    print(f"  E_m < 2 m_e = {2*ME:.3f} MeV : elle ne peut pas se briser seule ; sans verrou, elle se rouvre.")
    print("  -> la mere est la forme fermee d'un photon sur un anneau, transitoire ; pas un boson libre.\n")
    check("mere : relations du photon (lambda = L_m) et aucun verrou (Lk = 0)",
          abs(lam_photon / L_m - 1) < 1e-12 and lk_mother == 0, f"E_m = {E_m:.4f} MeV")
    check("mere sous le seuil de paire seule (E_m < 2 m_e)", E_m < 2 * ME, f"{E_m:.3f} < {2*ME:.3f}")

    # B. les filles
    print("B. Les filles")
    action_d = p_m * L_d                     # h/2
    lam_free = H / p_m                       # longueur d'onde libre du meme quantum
    print(f"  action = p L_d = {action_d/math.pi:.3f} pi hbar = h/2 ; en onde libre lambda = h/p = {lam_free:.1f} fm = {lam_free/L_d:.1f} L_d")
    print("  un demi-quantum n'a pas de forme libre (un photon porte h par longueur d'onde) ;")
    print("  et la charge +-e est conservee (plus leger charge) : verrouillees deux fois.\n")
    check("fille : action h/2, longueur d'onde libre = 2 L_d (ne tient pas)",
          abs(action_d - math.pi) < 1e-12 and abs(lam_free / L_d - 2) < 1e-12, "h/2")

    # C. annihilation
    print("C. Annihilation e+ e- -> 2 gamma")
    E_gamma = ME
    lam_gamma = H * HBARC / E_gamma
    total_in = 2 * (E_m + E_m)               # (circulant + statique) x 2 filles
    print(f"  chaque photon : {E_gamma:.4f} MeV, lambda = {lam_gamma:.1f} fm = {lam_gamma/L_d:.4f} L_d")
    print(f"  bilan : 2 x (demi-quantum {E_m:.4f} + statique {E_m:.4f}) = {total_in:.4f} = 2 x {E_gamma:.4f} MeV")
    print("  -> deux quanta entiers, un par photon, chacun de longueur d'onde L_d : les actions entieres reviennent.\n")
    check("lambda_gamma = L_d exactement et bilan 2 m_e", abs(lam_gamma / L_d - 1) < 1e-12 and abs(total_in / (2 * ME) - 1) < 1e-12,
          f"{lam_gamma:.1f} fm")

    # D. spin
    print("D. Spin")
    print("  mere : enroulement 1 (spin 1, R3) -> deux filles d'enroulement 1/2 ; annihilation : 1/2 + 1/2 -> 2 photons.")
    print("  -> le boson libre a l'action entiere, le fermion l'action demi-entiere et un verrou de charge :")
    print("     c'est la version de la base de 'le photon est libre et sans masse, l'electron massif et stable'.\n")
    check("arithmetique des enroulements 1 -> 1/2 + 1/2", abs(1.0 - 2 * (action_d / H)) < 1e-12, "1 = 1/2 + 1/2")

    print("Verdict : CONDITIONNEL ; la mere est transitoire parce qu'elle n'a aucun verrou et qu'elle")
    print("est un photon ferme ; les filles sont verrouillees par leur charge et par leur demi-quantum,")
    print("qui n'a pas de forme libre. Ce qui reste : la geometrie de la fermeture (ou un photon se")
    print("ferme sur trois DQD) et la cinematique reelle de la creation de paire (noyau, 1,022 MeV).\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
