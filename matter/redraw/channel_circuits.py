#!/usr/bin/env python3
"""R80 -- "Chaque canal C3 (singulet, doublet) est-il un circuit ferme portant un
quantum ?"  (le point a demontrer laisse par R79)

  A. topologie dessinee (R42, R54) : un anneau ferme de trois brins en serie,
     trois jonctions -> un triangle ; son espace des cycles est de dimension 1
     et porte le seul singulet : le doublet n'est pas un circuit, R47 ("un
     quantum par circuit ferme") ne peut pas lui donner de quantum.
  B. topologie alternative compatible avec R53 et R54 : trois boucles fermees
     coaxiales (une par brin), chacune refermee sur sa propre jonction pole a
     pole ; espace des cycles de dimension 3 = singulet + doublet ; mu_B, S =
     hbar/2 et les trois jonctions de R54 sont conserves.
  C. dans cette topologie, la matrice de couplage des trois boucles est le
     circulant L I + M (U + U^2) ; M/L avec les ecarts de la base (brins qui se
     touchent, ecart de jonction D = 6 lambda/pi^2) contre 1/sqrt2 (Koide).
  D. energies par canal avec un quantum egal par canal : E_s/E_d = (L+2M)/(L-M),
     jamais 1 pour M > 0 : "un quantum par circuit" ne donne pas l'egalite des
     canaux dans la lecture couplage.
"""
import math
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
ME = 0.51099895
LAMBDA_BAR = HBARC / ME
W_E = 4 * LAMBDA_BAR / math.pi ** 2          # 156.5 fm
D_J = 6 * LAMBDA_BAR / math.pi ** 2          # 235 fm, ecart de jonction (R54)
GMD_SQUARE = 0.44705                         # distance moyenne geometrique d'un carre / cote

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def cycle_space(n_vertices, edges):
    """dimension de l'espace des cycles et composantes connexes"""
    parent = list(range(n_vertices))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    comps = len({find(v) for v in range(n_vertices)})
    return len(edges) - n_vertices + comps, comps

def z3_character(perm_matrix):
    """caractere (identite, rotation, rotation^2) et decomposition singulet/doublet"""
    P = perm_matrix
    chars = [int(round(float(np.trace(np.linalg.matrix_power(P, k))))) for k in range(3)]
    n_singlet = round(sum(chars) / 3)
    n_doublet = round((chars[0] - n_singlet) / 2)
    return chars, n_singlet, n_doublet

def ring_L(R, r):
    return math.log(8 * R / r) - 2          # en unites mu0 R

def main():
    print("R80 -- les canaux C3 sont-ils des circuits ?\n")

    # A. topologie dessinee : triangle
    print("A. Anneau de trois brins en serie, trois jonctions (R42, R54) : un triangle")
    dim, comps = cycle_space(3, [(0, 1), (1, 2), (2, 0)])
    U = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    # Z3 agit sur l'unique cycle par l'identite
    print(f"  E = 3, V = 3 : dimension de l'espace des cycles = {dim} ; le cycle est invariant sous Z3")
    print("  -> un seul circuit ferme, le singulet ; le doublet n'a pas de circuit : R47 ne lui donne")
    print("     aucun quantum. Dans cette topologie, 'un quantum par canal' est impossible.\n")
    check("triangle : un seul cycle (singulet), pas de circuit pour le doublet", dim == 1, f"dim = {dim}")

    # B. topologie alternative : trois boucles coaxiales
    print("B. Trois boucles fermees coaxiales, une jonction chacune")
    dim3, comps3 = cycle_space(3, [(0, 0), (1, 1), (2, 2)])
    chars, ns, nd = z3_character(U)
    print(f"  E = 3 boucles, V = 3 jonctions, {comps3} composantes : dimension = {dim3}")
    print(f"  caractere de Z3 sur les cycles : {chars} = {ns} singulet + {nd} doublet")
    mu = 3 * (1 / 3) * LAMBDA_BAR / 2 / (LAMBDA_BAR / 2)        # (sum q_i) c R / 2 en mu_B
    S = LAMBDA_BAR * (ME / 2) / HBARC                          # en hbar
    print(f"  trois boucles de -e/3 a c sur R = lambda-bar : mu = {mu:.3f} mu_B, S = {S:.3f} hbar ;")
    print("  trois jonctions pole a pole (une par boucle) : la moitie statique de R54 inchangee.")
    print("  -> ici les deux canaux sont des circuits ; c'est la topologie ou l'enonce a un sens.\n")
    check("trois boucles : cycles = singulet + doublet, mu_B et S conserves",
          dim3 == 3 and ns == 1 and nd == 1 and abs(mu - 1) < 1e-12 and abs(S - 0.5) < 1e-12,
          f"dim {dim3}, {ns} + {nd}x2")

    # C. couplage des trois boucles
    print("C. Couplage des trois boucles : L = ln(8R/r_g) - 2, M = ln(8R/d) - 2 (unites mu0 R)")
    R = LAMBDA_BAR
    r_g = GMD_SQUARE * W_E
    L = ring_L(R, r_g)
    target = 1 / math.sqrt(2)
    for name, d in (("brins qui se touchent, d = w", W_E), ("ecart de jonction, d = 6 lambda/pi^2", D_J),
                    ("d = 2w", 2 * W_E)):
        M = ring_L(R, d)
        print(f"    {name:38s} M/L = {M/L:.3f}  (||dev||/||iso|| = sqrt2 M/L = {math.sqrt(2)*M/L:.3f})")
    # d qui donnerait 1/sqrt2
    d_star = 8 * R / math.exp(target * L + 2)
    print(f"  M/L = 1/sqrt2 demanderait d = {d_star:.0f} fm = {d_star/W_E:.2f} w : les brins se recouvriraient.")
    print("  -> le couplage magnetique statique des trois boucles ne donne pas le 45 deg.\n")
    M_touch = ring_L(R, W_E)
    check("M/L aux ecarts de la base loin de 1/sqrt2 (> 15 %)",
          abs(M_touch / L / target - 1) > 0.15 and abs(ring_L(R, D_J) / L / target - 1) > 0.15,
          f"{M_touch/L:.3f}, {ring_L(R, D_J)/L:.3f} vs {target:.3f}")

    # D. un quantum egal par canal
    print("D. Un quantum egal par canal dans la lecture couplage")
    for name, d in (("d = w", W_E), ("d = D", D_J)):
        M = ring_L(R, d)
        ratio = (L + 2 * M) / (L - M)
        print(f"    {name:8s} E_s/E_d = (L+2M)/(L-M) = {ratio:.2f}")
    print("  -> jamais 1 tant que M > 0 : 'un quantum par circuit' (R47) n'egalise pas les canaux ;")
    print("     l'egalite des normes de R79 n'est pas une egalite d'energies de circuits.\n")
    check("E_s/E_d != 1 pour tout M > 0", (L + 2 * M_touch) / (L - M_touch) > 1.5, "M > 0")

    print("Verdict : le canal doublet n'est un circuit que si l'electron est trois boucles coaxiales")
    print("(compatible avec R53, R54) ; meme alors, ni le couplage statique (M/L = 0,3-0,6) ni un")
    print("quantum egal par circuit ne donnent le 45 deg. La conservation 'un quantum par canal' de")
    print("R79 n'est pas celle de R47 : elle reste a inventer, ou le 45 deg vient d'ailleurs.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
