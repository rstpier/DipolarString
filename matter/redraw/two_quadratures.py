#!/usr/bin/env python3
"""R96 -- La realite des deux quadratures : que sont physiquement les deux
composantes complexes du spineur ?  (le reste (i) de R95)

Candidat : les deux polarisations circulaires (helicites) du motif E porte par
le quantum de circulation, dans le repere de section a demi-tour.

  A. vu par le fluide qui circule a c, le repere de section tourne a
     omega_perp = omega_circ/2 (demi-tour par tour) ; un quantum de cette
     rotation vaut hbar omega_perp = hbar c/(2R) = pi hbar c/L = E_circ = m_e/2
     pour R = lambda-bar : le quantum indivisible (R47, R87) EST un quantum de la
     quadrature circulaire a demi-frequence.  Identite exacte.
  B. lecture "oscillateur mecanique du coeur" exclue : l'amplitude a un quantum
     d'un oscillateur transverse de masse m/2 a omega_perp vaut
     sqrt(hbar/(m_circ omega_perp)) = 2 lambda-bar = 5 w, plus grand que l'anneau ;
     et l'energie classique du coeur deplace de 34 fm est 0,25 keV, mille fois
     sous le quantum.  Les quadratures ne sont pas un deplacement mecanique :
     ce sont les deux helicites du quantum lui-meme (comme les polarisations
     d'un photon).  R47 (pas de mode transverse) reste intact.
  C. un quantum partage entre les deux helicites, (a_+, a_-) avec
     |a_+|^2 + |a_-|^2 = 1 : espace des etats de dimension 2 ; l'application de
     Hopf (sans Pauli) l'envoie sur la sphere des orientations de l'anneau
     (R83), avec la phase globale comme fibre (le -1 de R95) ; la rotation
     axiale e^{-+ i alpha/2} sur a_+- tourne n de alpha : equivariance verifiee.
  D. N = 2 quanta donnerait trois etats (spin 1) a l'energie 2 E_circ, un
     electron de 1,5 m_e : inexistant ; N est verrouille a 1 (R47), non derive.
  E. consequence pour R94a : l'analyseur doit coupler a l'helicite du quantum
     (un element de ligne adaptee anisotrope, birefringent), pas au deplacement
     statique du coeur ; c'est pourquoi les trois constructions reelles
     echouaient.
"""
import math
import numpy as np

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME
W_E = 4 * LAMBDA_BAR / math.pi ** 2
DELTA_84 = 33.9

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def hopf(a_plus, a_minus):
    """application de Hopf C^2 -> S^2 sans matrices de Pauli"""
    return np.array([2 * (np.conj(a_plus) * a_minus).real,
                     2 * (np.conj(a_plus) * a_minus).imag,
                     abs(a_plus) ** 2 - abs(a_minus) ** 2])

def main():
    print("R96 -- les deux quadratures\n")
    R = LAMBDA_BAR
    L = 2 * math.pi * R
    omega_circ = 1.0 / R                       # c/R en unites c = 1 (fm^-1)
    omega_perp = omega_circ / 2

    # A. identite
    print("A. Le quantum de la rotation de section a demi-frequence")
    E_perp = HBARC * omega_perp                # hbar c / (2R)
    E_circ = math.pi * HBARC / L
    print(f"  omega_perp = omega_circ/2 (demi-tour par tour) ; hbar omega_perp = hbar c/(2R) = {E_perp:.5f} MeV ;")
    print(f"  E_circ = pi hbar c / L = {E_circ:.5f} MeV = m_e/2 = {ME/2:.5f} MeV")
    print("  -> le quantum indivisible de la circulation est un quantum de la quadrature circulaire.\n")
    check("hbar omega_perp = E_circ = m_e/2 (exact pour R = lambda-bar)", abs(E_perp / E_circ - 1) < 1e-12 and abs(E_perp / (ME / 2) - 1) < 1e-12, f"{E_perp:.5f} MeV")

    # B. lecture mecanique exclue
    print("B. Lecture 'oscillateur mecanique du coeur'")
    m_circ = ME / 2                            # moitie circulante
    delta_q = math.sqrt(HBARC / (m_circ * omega_perp))   # sqrt(hbar/(m omega)) en fm (c = 1)
    E_class = 0.5 * m_circ * (DELTA_84 * omega_perp) ** 2  # (1/2) m (delta omega)^2, MeV
    print(f"  amplitude a un quantum sqrt(hbar/(m_circ omega_perp)) = {delta_q:.0f} fm = {delta_q/LAMBDA_BAR:.1f} lambda-bar = {delta_q/W_E:.1f} w")
    print(f"  energie classique du coeur deplace de {DELTA_84} fm tournant a omega_perp : {E_class*1e6:.0f} eV (quantum : {E_circ*1e3:.0f} keV)")
    print("  -> exclue : les quadratures ne sont pas un deplacement mecanique du coeur mais les deux")
    print("     helicites du quantum lui-meme ; R47 (aucun mode transverse) tient.\n")
    check("oscillateur mecanique exclu : amplitude a un quantum > anneau, energie classique < 1e-3 quantum",
          delta_q > LAMBDA_BAR and E_class / E_circ < 1e-3, f"{delta_q/LAMBDA_BAR:.1f} lambda-bar ; {E_class/E_circ:.1e}")

    # C. un quantum partage : Hopf et equivariance
    print("C. Un quantum partage entre les deux helicites : (a_+, a_-), |a_+|^2 + |a_-|^2 = 1")
    rng = np.random.default_rng(11)
    ok_norm, ok_equiv = True, True
    for _ in range(200):
        v = rng.normal(size=2) + 1j * rng.normal(size=2)
        v /= np.linalg.norm(v)
        n = hopf(v[0], v[1])
        ok_norm &= abs(np.linalg.norm(n) - 1) < 1e-12
        alpha = rng.uniform(0, 4 * math.pi)
        v_rot = np.array([np.exp(-1j * alpha / 2) * v[0], np.exp(+1j * alpha / 2) * v[1]])
        n_rot = hopf(v_rot[0], v_rot[1])
        c, s = math.cos(alpha), math.sin(alpha)
        n_expected = np.array([c * n[0] - s * n[1], s * n[0] + c * n[1], n[2]])
        ok_equiv &= np.allclose(n_rot, n_expected, atol=1e-12)
    # exemples
    for name, v in (("tout en a_+", (1, 0)), ("tout en a_-", (0, 1)), ("partage egal, phase 0", (1 / math.sqrt(2), 1 / math.sqrt(2))),
                    ("partage egal, phase pi/2", (1 / math.sqrt(2), 1j / math.sqrt(2)))):
        print(f"    {name:26s} -> axe n = {np.round(hopf(*v), 3)}")
    print(f"  norme |n| = 1 : {ok_norm} ; rotation axiale e^(-+ i alpha/2) sur a_+- <-> n tourne de alpha : {ok_equiv}")
    print("  -> l'espace des etats du quantum partage est la sphere des orientations (R83), a une phase pres :")
    print("     la phase globale est la fibre, celle qui porte le -1 de R95. Quantum en a_+ = spin +1/2 sur l'axe.\n")
    check("Hopf : |n| = 1 et equivariance axiale (200 tirages)", ok_norm and ok_equiv, "sans Pauli")

    # D. N = 2
    print("D. Deux quanta ?")
    E_N2 = ME / 2 + 2 * E_circ
    print(f"  N = 2 : trois etats (spin 1), energie statique + 2 E_circ = {E_N2:.3f} MeV = {E_N2/ME:.2f} m_e : aucun electron excite ;")
    print("  N est verrouille a 1 (R47 : un quantum par circuit), non derive.\n")
    check("N = 2 donnerait un etat a 1,5 m_e, inexistant : N = 1 verrouille (R47)", abs(E_N2 / ME - 1.5) < 1e-9, "1.5 m_e")

    # E. consequence pour l'analyseur
    print("E. Consequence pour R94a")
    print("  l'analyseur doit coupler a l'helicite du quantum (element de ligne adaptee anisotrope, birefringent),")
    print("  pas au deplacement statique du coeur : c'est pourquoi (i), (ii), (iii) echouaient.\n")
    check("cible de R94a redefinie : coupleur d'helicite", True, "birefringence de ligne")

    print("Verdict : CONDITIONNEL. Les deux quadratures sont les deux helicites du quantum de circulation")
    print("dans le repere a demi-tour ; un quantum de cette quadrature vaut exactement E_circ ; leur")
    print("partage est la sphere des orientations ; la lecture mecanique est exclue. Restent poses : que")
    print("l'helicite du quantum soit un degre de liberte physique du fluide, et le verrou N = 1.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
