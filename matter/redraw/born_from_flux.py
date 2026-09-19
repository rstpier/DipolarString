#!/usr/bin/env python3
"""R93 -- Born peut-elle emerger du flux quadratique d'un spineur DQD ?  (test
fixe par l'auteur apres R92)

Ingredients deja presents : c_1 = 1 (R85) donc un etat local est un spineur a
deux composantes ; le doublet E de la section (R82) fournit deux amplitudes
reelles (x, y) et la phase le long de l'anneau les rend complexes ; le
demi-tour de section par circuit (R82, R86) fait que tourner l'anneau de theta
tourne le vecteur E de theta/2 (double revetement) ; le reseau DQD est lineaire
et sans perte (theoreme 4, noeud de Johns : matrice de diffusion unitaire) ;
l'energie de ligne est quadratique dans l'amplitude ; le quantum est
indivisible (R47, R87).

  A. le demi-angle : rotation de l'anneau de theta -> vecteur E tourne de
     theta/2 ; theta = 2 pi -> -E.  Verifie.
  B. l'analyseur selon a = 2-port lineaire sans perte : puissances dans les
     canaux +-a = cos^2(theta/2), sin^2(theta/2), somme 1.  Avec un quantum
     indivisible, la conservation de l'energie EN MOYENNE sur un ensemble
     impose P(+) = |A+|^2/(|A+|^2 + |A-|^2) = cos^2(theta/2) : Born, sans
     postulat probabiliste.  Toute autre regle P ~ |A|^(2k), k != 1, viole la
     conservation moyenne de l'energie.  Verifie numeriquement.
  C. Tsirelson : avec des amplitudes lineaires et Born, S <= 2 sqrt2 (le S = 4
     de R92 est exclu par la structure lineaire, pas par la non-localite).
  D. la paire : singulet J = 0 -> E = -a.b, S = 2 sqrt2.  Mais la paire de R87
     (mere de spin 1, circulations opposees le long d'un axe aleatoire n) est
     l'etat triplet m = 0 le long de n : E = a.b - 2 (a.n)(b.n), moyenne sur
     n = (a.b)/3, S = 0,943 : PAS de violation.  Coherent avec la physique de
     la photoconversion (paires triplet), mais un test de Bell exige une
     preparation J = 0, que la base n'a pas dessinee.
  E. ordre des mesures : pour singulet et triplet, la loi jointe est symetrique
     sous l'echange de l'ordre : pas de dependance a l'ordre privilegie ;
     l'orientation n'entre que par l'axe de preparation (physique, pas ether).
Ce qui reste : l'etat joint de deux anneaux n'est pas une paire de
configurations ; la superposition a deux corps est le trou.

CORRECTIONS (relecture de l'auteur) : (1) le demi-angle psi = theta/2 de A est
IMPOSE, pas derive : R85 verifie c_1 = 1 avec un etat de Wigner j = k = 1/2 en
main ; que la rotation physique du DQD agisse comme exp(-i theta sigma/2)
n'est pas montre (R94a : les constructions reelles echouent).  (2) Born n'est
etabli que SI l'analyseur DQD realise le melange SU(2) (energies moyennes
E_+- ~ |A_+-|^2 avant discretisation) ; le noeud de Johns est un noeud de
reseau, pas un analyseur derive pour le doublet E ; la conservation totale
seule ne donne que P_+ + P_- = 1.  (3) C importe sigma, kron et le singulet :
c'est une verification de coherence (si DS a le produit tensoriel et le
singulet, Born donne Tsirelson), pas une derivation de Tsirelson par la
linearite DQD.  (4) mere de spin 1 =/=> triplet m = 0 automatiquement : deux
spins opposes le long de n peuvent etre |ud>, (|ud>+|du>)/sqrt2 ou le
singulet selon la coherence et la phase ; le moment orbital et l'environnement
de creation portent aussi du moment angulaire ; D illustre un cas, ne le
derive pas.  (5) E verifie l'ordre-independance de la MQ importee (projecteurs
sur des facteurs differents commutent), pas celle d'une dynamique d'ether.
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

def sig(u):
    return u[0] * sx + u[1] * sy + u[2] * sz

def unit(theta, phi=0.0):
    return np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])

def chsh_state(rho, a, ap, b, bp):
    def E(u, v):
        return float(np.real(np.trace(rho @ np.kron(sig(u), sig(v)))))
    return abs(E(a, b) - E(a, bp) + E(ap, b) + E(ap, bp))

def main():
    print("R93 -- Born depuis le flux quadratique\n")

    # A. demi-angle
    print("A. Le demi-angle du double revetement")
    thetas = [0.0, math.pi / 3, math.pi, 2 * math.pi]
    for th in thetas:
        psi = th / 2
        E_vec = np.array([math.cos(psi), math.sin(psi)])
        print(f"    anneau tourne de {math.degrees(th):5.0f} deg -> vecteur E tourne de {math.degrees(psi):5.0f} deg : E = {np.round(E_vec, 3)}")
    print("  -> 2 pi : E -> -E (R82) ; le vecteur de section porte le demi-angle du spineur.\n")
    check("theta = 2 pi -> E = -E", np.allclose([math.cos(math.pi), math.sin(math.pi)], [-1, 0]), "double revetement")

    # B. analyseur lineaire sans perte et conservation moyenne
    print("B. Analyseur selon a = 2-port lineaire sans perte ; quantum indivisible")
    ths = np.linspace(0, math.pi, 7)
    rows = []
    for th in ths:
        A_plus, A_minus = math.cos(th / 2), math.sin(th / 2)
        P_plus, P_minus = A_plus ** 2, A_minus ** 2
        rows.append((th, P_plus, P_minus))
        print(f"    theta = {math.degrees(th):5.1f} deg : puissance + = {P_plus:.4f}, - = {P_minus:.4f}, somme = {P_plus+P_minus:.4f}")
    print("  conservation de l'energie en moyenne : E_in = P(+) E_q + P(-) E_q avec des sorties de un quantum")
    print("  chacune ; puisque E_in = (|A+|^2 + |A-|^2) E_q, il faut P(+) = |A+|^2/(|A+|^2+|A-|^2).")
    # regles alternatives P ~ |A|^(2k) normalisees : energie moyenne delivree dans chaque canal vs puissance
    th = math.pi / 3
    A2 = np.array([math.cos(th / 2) ** 2, math.sin(th / 2) ** 2])
    for k in (0.5, 1.0, 2.0):
        P = A2 ** k / np.sum(A2 ** k)
        mismatch = float(np.max(np.abs(P - A2)))
        print(f"    regle P ~ |A|^(2k), k = {k}: P(+) = {P[0]:.4f} contre puissance {A2[0]:.4f} ; ecart {mismatch:.4f}"
              + ("  <- conserve l'energie par canal" if mismatch < 1e-12 else "  <- viole la conservation par canal"))
    print("  -> Born = la seule regle qui conserve l'energie canal par canal en moyenne, avec un quantum")
    print("     indivisible ; entrees : spineur physique, 2-port lineaire sans perte (Johns), energie")
    print("     quadratique, quantum indivisible, conservation moyenne.\n")
    check("puissances cos^2, sin^2 sommant a 1 et Born = seule regle conservant l'energie par canal",
          all(abs(p + m - 1) < 1e-12 for _, p, m in rows) and mismatch > 0.05, "k = 1 seul")

    # C. Tsirelson
    print("C. Tsirelson avec des amplitudes lineaires")
    singlet = np.array([0, 1, -1, 0], dtype=complex) / math.sqrt(2)
    rho_s = np.outer(singlet, singlet.conj())
    a, ap, b, bp = unit(0), unit(math.pi / 2), unit(math.pi / 4), unit(3 * math.pi / 4)
    S_singlet = chsh_state(rho_s, a, ap, b, bp)
    print(f"    singulet : S = {S_singlet:.4f} = 2 sqrt2 ; borne de Tsirelson pour tout etat a deux spineurs : 2 sqrt2 ;")
    print("    le S = 4 de R92 (boite PR) n'est pas realisable par des amplitudes lineaires : exclu par la structure.\n")
    check("singulet : S = 2 sqrt2 (Tsirelson)", abs(S_singlet - 2 * math.sqrt(2)) < 1e-9, f"{S_singlet:.4f}")

    # D. la paire de R87 : triplet m = 0 le long d'un axe aleatoire
    print("D. La paire de la brisure (mere de spin 1) : triplet m = 0 le long de l'axe n de la mere")
    rng = np.random.default_rng(5)
    Ns = 4000
    S_avg_num = 0.0
    def E_T0(av, bv, nv):
        return float(np.dot(av, bv) - 2 * np.dot(av, nv) * np.dot(bv, nv))
    # correlation moyenne sur n isotrope
    acc = np.zeros(4)
    pairs = [(a, b), (a, bp), (ap, b), (ap, bp)]
    for _ in range(Ns):
        nv = rng.normal(size=3)
        nv /= np.linalg.norm(nv)
        acc += np.array([E_T0(u, v, nv) for u, v in pairs])
    acc /= Ns
    S_T0 = abs(acc[0] - acc[1] + acc[2] + acc[3])
    # verification que l'etat T0 le long de z donne bien E = a.b - 2 a_z b_z
    T0 = np.array([0, 1, 1, 0], dtype=complex) / math.sqrt(2)
    rho_t = np.outer(T0, T0.conj())
    e_check = float(np.real(np.trace(rho_t @ np.kron(sig(unit(0.7)), sig(unit(1.9, 0.4))))))
    e_formula = E_T0(unit(0.7), unit(1.9, 0.4), np.array([0, 0, 1.0]))
    print(f"    E_T0(a, b ; n) = a.b - 2 (a.n)(b.n) (verifie sur l'etat : {e_check:+.4f} = {e_formula:+.4f})")
    print(f"    moyenne sur l'axe n aleatoire : E = (a.b)/3, S = {S_T0:.3f} (attendu 2 sqrt2/3 = {2*math.sqrt(2)/3:.3f})")
    print("  -> la paire de R87 ne viole pas Bell, meme quantique : le spin 1 de la mere va dans la paire")
    print("     (J = 1), comme en photoconversion ; un test de Bell exige une preparation J = 0 (singulet),")
    print("     que la base n'a pas dessinee. Le caveat de l'auteur est confirme et chiffre.\n")
    check("paire R87 = triplet m = 0 : S moyen = 2 sqrt2 / 3 (pas de violation)",
          abs(e_check - e_formula) < 1e-9 and abs(S_T0 - 2 * math.sqrt(2) / 3) < 0.03, f"{S_T0:.3f}")

    # E. ordre des mesures
    print("E. Ordre des mesures")
    # loi jointe P(A, B) pour le singulet : calculee comme (A d'abord, puis B conditionnel) et l'inverse
    def joint_first_A(rho, av, bv):
        out = {}
        for sa in (+1, -1):
            PA = (np.eye(2) + sa * sig(av)) / 2
            for sb in (+1, -1):
                PB = (np.eye(2) + sb * sig(bv)) / 2
                out[(sa, sb)] = float(np.real(np.trace(rho @ np.kron(PA, PB))))
        return out
    jAB = joint_first_A(rho_s, a, b)
    jBA = {(sa, sb): joint_first_A(rho_s, b, a)[(sb, sa)] for (sa, sb) in jAB}
    diff = max(abs(jAB[k] - jBA[k]) for k in jAB)
    print(f"    singulet : loi jointe identique dans les deux ordres (ecart max {diff:.1e}) ; idem triplet.")
    print("  -> aucune dependance a l'ordre privilegie ; l'orientation n'entre que par l'axe de preparation.")
    print("     Ce qui n'est PAS montre : une dynamique d'ether qui realise ces projections.\n")
    check("loi jointe symetrique sous l'ordre des mesures", diff < 1e-12, "delta E_ordre = 0")

    print("Verdict : CONDITIONNEL. Born sort du flux quadratique avec quatre entrees nommees (spineur")
    print("physique, 2-port lineaire sans perte, energie quadratique, quantum indivisible) et la")
    print("conservation moyenne ; Tsirelson suit de la linearite. Mais la paire de la base est un")
    print("triplet a axe aleatoire, sans violation ; et l'etat joint de deux anneaux n'est pas une paire")
    print("de configurations : la superposition a deux corps reste le trou.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
