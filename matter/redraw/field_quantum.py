#!/usr/bin/env python3
"""R104 -- L'electron comme quantum d'un champ universel : la paire (V, I) de la ligne, l'anneau
ferme, et ce qu'il manque pour que son quantum soit l'electron de Dirac.

Lecture testee (proposee par l'auteur apres R103) : l'electron n'est pas un objet localise mais
le quantum du champ du milieu ; il est etale comme le mode tant qu'il n'interagit pas, ponctuel
quand il est absorbe ou diffuse.  Le manuscrit le dit deja : "l'electron est le quantum de la
ligne, pas son soliton".

  A. les equations du telegraphiste sur une ligne adaptee sont un champ a deux composantes du
     premier ordre : les deux vecteurs propres du symbole sont les deux sens de propagation,
     V +- Z0 I, a la vitesse +-c : la paire de Weyl 1+1 D.  Identite, verifiee.
  B. sur un anneau ferme de longueur L, le spectre du champ est hbar c |k| avec k = 2 pi n/L
     (secteur periodique) ou 2 pi (n + 1/2)/L (secteur antiperiodique) : le quantum entier
     2 pi hbar c/L de la mere (R87) et le demi-quantum pi hbar c/L de la fille (R47, R86) sont
     les etats fondamentaux des deux secteurs de conditions aux limites.  Le -1 sous un tour est
     la condition antiperiodique du champ sur l'anneau.  Identite, verifiee.
  C. le mode fondamental antiperiodique porte L_z = hbar/2 a tout rayon ; son energie
     hbar c/(2R) fixe le rayon.  Avec toute l'energie de repos dans le mode, R = lambda-bar/2,
     tour a 2 m c^2/hbar, phase a m c^2/hbar : les trois nombres du Zitterbewegung de Dirac
     (R103 B), sans moitie statique.  La lecture R53 (R = lambda-bar, m/2 dans le mode) est le
     meme mode a un rayon double, qui doit loger m/2 ailleurs (R54, coincidence).
  D. le vertex de charge d'un quantum de champ : <k+q| rho(q) |k> = e pour tout q, quelle que soit
     l'etendue du mode ; F = 1 par construction, contre j_0(qR) pour une charge classique etalee
     (R100).  C'est ce que la lecture champ achete.
  E. le moment sans boucle : pour un champ a deux composantes couple au premier ordre,
     (sigma.pi)^2 = pi^2 - e hbar sigma.B (verifie dans la representation de Landau), donc
     H = (sigma.pi)^2/2m contient -mu_B sigma.B : mu = mu_B, g = 2, sans courant de charge
     spatial.  Levy-Leblond (1967) : la linearisation de Schrodinger suffit, la relativite n'est
     pas necessaire.
  F. ce qui manque a la base : son doublet (V_+, V_-) est attache a UNE direction de ligne
     (sigma_z d/dx) ; pour sigma.grad il faut que le doublet tourne comme un spineur avec la
     direction de propagation (le relevement SU(2) de R95 applique au mode du milieu, pas a un
     anneau), avec l'algebre de Pauli entre les trois directions de la trame.  V2.10 n'a que
     Phi et theta : c'est l'addition structurelle, une seule, nommee.
  G. consequence sur l'echelle des hadrons : l'ancre l_1(3) = 2 pi lambda-bar/3 etait "le circuit
     de l'electron divise par trois" ; dans la lecture champ le circuit de l'electron mesure
     pi lambda-bar, et un circuit de quark "electron a l'echelle 9" pesent 509 MeV et non 254.
     Les nombres hadroniques (763 MeV, sigma, n - p) ne changent que si l'on garde
     l'ancre 2 pi lambda-bar comme longueur posee ; ils ne sont plus derives de la geometrie de
     l'electron.
"""
import math
import numpy as np

HBARC, ME = 197.3269804, 0.51099895
LAMBDA_BAR = HBARC / ME

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R104 -- l'electron comme quantum d'un champ\n")

    # A. telegraphiste = paire de Weyl
    print("A. Les equations du telegraphiste : d_x V = -L' d_t I, d_x I = -C' d_t V")
    Lp, Cp = 2.0, 0.5                      # unites arbitraires : c = 1/sqrt(L'C') = 1, Z0 = sqrt(L'/C') = 2
    c, Z0 = 1 / math.sqrt(Lp * Cp), math.sqrt(Lp / Cp)
    okA = True
    for k in (0.3, 1.0, 7.0):
        M = np.array([[0, -1j * k / Cp], [-1j * k / Lp, 0]])       # d_t (V, I) = M (V, I) pour e^{ikx}
        lam, vec = np.linalg.eig(M)
        for l, v in zip(lam, vec.T):
            ratio = Z0 * v[1] / v[0]
            okA &= abs(abs(l) - c * k) < 1e-12 and abs(abs(ratio) - 1) < 1e-12 and abs(l.real) < 1e-12
        print(f"    k = {k:4.1f} : valeurs propres {lam[0]:+.3f}, {lam[1]:+.3f} = -+ i c k ; Z0 I/V = {Z0*vec[1,0]/vec[0,0]:+.3f}, {Z0*vec[1,1]/vec[0,1]:+.3f}")
    print("  => deux modes propres, V +- Z0 I, a +-c, sans dispersion : un champ a deux composantes du premier")
    print("     ordre, i d_t psi = -i c sigma_z d_x psi avec psi = (V + Z0 I, V - Z0 I) : la paire de Weyl 1+1 D.\n")
    check("telegraphiste = paire de Weyl 1+1 D : modes propres V +- Z0 I a +-c (identite)", okA, "sigma_z d_x")

    # B. spectre sur l'anneau : periodique / antiperiodique
    print("B. Le champ sur un anneau ferme de longueur L (hbar = c = 1)")
    L = 2 * math.pi
    N = 64
    def spectrum(theta):
        # operateur -i sigma_z d_x dans la base de Fourier tordue k_n = 2 pi (n + theta)/L
        ks = 2 * math.pi * (np.arange(-N, N + 1) + theta) / L
        return np.sort(np.abs(np.concatenate([ks, -ks])))
    per, anti = spectrum(0.0), spectrum(0.5)
    e_per = per[per > 1e-12][0]
    e_anti = anti[0]
    print(f"    periodique   : plus basse energie non nulle = {e_per:.6f} = 2 pi/L ({2*math.pi/L:.6f}) [plus un mode zero]")
    print(f"    antiperiodique : plus basse energie = {e_anti:.6f} = pi/L ({math.pi/L:.6f})")
    print("  => la mere (R87, quantum entier h c/L) est le fondamental periodique ; la fille (R47, R86, demi-quantum")
    print("     pi hbar c/L) est le fondamental antiperiodique. Le -1 sous un tour est la condition aux limites du champ.\n")
    check("anneau : periodique -> 2 pi hbar c/L (mere), antiperiodique -> pi hbar c/L (fille) (identite)",
          abs(e_per - 2 * math.pi / L) < 1e-12 and abs(e_anti - math.pi / L) < 1e-12, "R87 / R47")

    # C. le mode antiperiodique fondamental : L_z, rayon, frequences
    print("C. Le mode fondamental antiperiodique comme objet 3D")
    phi = np.linspace(0, 2 * math.pi, 2001)[:-1]
    u = np.exp(1j * phi / 2)                                   # e^{i phi/2}
    Lz = (np.vdot(u, -1j * np.gradient(u, phi)) / np.vdot(u, u)).real   # -i d/dphi
    R_field = LAMBDA_BAR / 2                                   # E = hbar c/(2R) = m c^2
    E_field = HBARC / (2 * R_field)
    omega_turn = 2.99792458e23 / R_field                       # rad/s, c/R
    omega_dB = ME / (HBARC / 2.99792458e23)                    # m c^2/hbar, rad/s
    print(f"    L_z = -i d/dphi sur e^(i phi/2) : {Lz:.6f} hbar (a tout rayon)")
    print(f"    E = hbar c/(2R) = m c^2  =>  R = lambda-bar/2 = {R_field:.2f} fm ; E = {E_field:.4f} MeV")
    print(f"    frequence de tour c/R = {omega_turn:.4e} rad/s = 2 m c^2/hbar ({2*omega_dB:.4e}) ; phase E/hbar = m c^2/hbar")
    print(f"    R103 B (Zitterbewegung de Dirac) : amplitude lambda-bar/2, omega = 2 m c^2/hbar, action pi hbar : les trois.")
    print(f"    lecture R53 : le meme mode a R = lambda-bar, E = m/2 dans le mode, phase a m c^2/(2 hbar), et m/2 a loger ailleurs (R54).\n")
    check("mode antiperiodique : L_z = hbar/2, E = m c^2 => R = lambda-bar/2, tour a 2 m c^2/hbar : Dirac, sans moitie statique",
          abs(Lz - 0.5) < 1e-6 and abs(E_field - ME) < 1e-9 and abs(omega_turn / (2 * omega_dB) - 1) < 1e-9, f"R = {R_field:.1f} fm")

    # D. vertex de charge d'un quantum de champ
    print("D. Le vertex de charge d'un quantum de champ sur l'anneau (R = lambda-bar, modes antiperiodiques)")
    R = LAMBDA_BAR
    Lr = 2 * math.pi * R
    x = np.linspace(0, Lr, 4001)[:-1]
    dx = x[1] - x[0]
    def mode(n):
        k = 2 * math.pi * (n + 0.5) / Lr
        return np.exp(1j * k * x) / math.sqrt(Lr)
    okD = True
    rows = []
    for nq in (1, 2, 5, 20, 200):
        q = 2 * math.pi * nq / Lr                              # transfert permis, fm^-1
        amp = np.sum(np.conj(mode(3 + nq)) * mode(3) * np.exp(1j * q * x)) * dx
        classical = math.sin(q * R) / (q * R)                  # anneau classique, moyenne d'orientation
        okD &= abs(abs(amp) - 1) < 1e-9
        rows.append((q * HBARC, abs(amp), classical))
        print(f"    q = {q*HBARC:8.3f} MeV/c : |<k+q| rho(q) |k>|/e = {abs(amp):.9f} ; charge classique sur l'anneau : j_0(qR) = {classical:+.4f}")
    print("  => le quantum porte sa charge en entier a chaque transfert, quelle que soit l'etendue du mode : F = 1 par")
    print("     construction. La 'taille' d'un quantum est la longueur d'onde de son paquet, pas une distribution de charge.\n")
    check("quantum de champ : |<k+q|rho(q)|k>| = e pour tout q (F = 1), contre j_0(qR) classique", okD, "1e-9")

    # E. (sigma.pi)^2 = pi^2 - e hbar sigma.B : g = 2 sans boucle
    print("E. Le moment sans boucle : (sigma.pi)^2 dans la representation de Landau (e hbar B = 1)")
    nmax = 80
    a = np.diag(np.sqrt(np.arange(1, nmax)), 1).astype(complex)
    ad = a.conj().T
    px = (a + ad) / math.sqrt(2)
    py = 1j * (ad - a) / math.sqrt(2)
    pz = 0.37 * np.eye(nmax)
    comm = px @ py - py @ px
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
    I2, In = np.eye(2), np.eye(nmax)
    sp = np.kron(sx, px) + np.kron(sy, py) + np.kron(sz, pz)
    lhs = sp @ sp
    rhs = np.kron(I2, px @ px + py @ py + pz @ pz) - np.kron(sz, In)
    cut = nmax - 3
    idx = np.r_[0:cut, nmax:nmax + cut]                        # etats loin de la troncature
    err = np.max(np.abs((lhs - rhs)[np.ix_(idx, idx)]))
    print(f"    [pi_x, pi_y] = i e hbar B : {comm[0,0]:.6f} ; max |(sigma.pi)^2 - (pi^2 - e hbar B sigma_z)| = {err:.1e} (loin de la troncature)")
    print("  => H = (sigma.pi)^2/2m = pi^2/2m - (e hbar/2m) sigma.B = ... - 2 mu_B S.B/hbar : mu = mu_B, g = 2.")
    print("     Levy-Leblond (1967) : la linearisation de Schrodinger (equation du premier ordre a deux composantes)")
    print("     donne exactement cela sans relativite. Le moment est dans la structure du champ, pas dans une boucle.\n")
    check("(sigma.pi)^2 = pi^2 - e hbar sigma.B (1e-10) : mu_B et g = 2 par la structure du premier ordre, sans courant",
          err < 1e-10 and abs(comm[0, 0] - 1j) < 1e-12, f"{err:.0e}")

    # F. ce qui manque
    print("F. Ce qui manque a la base")
    print("    le doublet (V_+, V_-) de la ligne est sigma_z d_x : une direction. Pour sigma.grad il faut que le doublet")
    print("    tourne comme un spineur avec la direction de propagation (R95 applique au mode du milieu, pas a un anneau),")
    print("    et que les trois directions de la trame portent l'algebre de Pauli (sigma_x sigma_y = i sigma_z).")
    print("    V2.10 n'a que Phi et theta : une addition structurelle, une seule, a ecrire (R105 candidat).\n")

    # G. consequence sur l'ancre des hadrons
    print("G. Consequence sur l'echelle des hadrons")
    l1_9_base = (2 * math.pi * LAMBDA_BAR / 3) * 3 ** (-2 * math.pi)
    l1_9_field = (math.pi * LAMBDA_BAR / 3) * 3 ** (-2 * math.pi)
    Eq_base, Eq_field = math.pi * HBARC / (3 * l1_9_base), math.pi * HBARC / (3 * l1_9_field)
    print(f"    ancre de la base : circuit de l'electron 2 pi lambda-bar (R = lambda-bar, m/2 dans le mode) -> quark {Eq_base:.1f} MeV, nucleon statique {3*Eq_base:.0f} MeV")
    print(f"    lecture champ    : circuit de l'electron pi lambda-bar (R = lambda-bar/2, m dans le mode)  -> quark {Eq_field:.1f} MeV, nucleon statique {3*Eq_field:.0f} MeV (> m_p)")
    print("  => les nombres hadroniques (R41-R45) tiennent seulement si l'ancre 2 pi lambda-bar est une longueur posee ;")
    print("     'le quark est l'electron a l'echelle 9' n'est plus une geometrie, c'est E_quark = (m_e/2) 3^(2 pi).\n")
    check("ancre des hadrons : la lecture champ (R = lambda-bar/2) doublerait le quark (509 MeV) ; 2 pi lambda-bar devient une longueur posee",
          abs(Eq_field / Eq_base - 2) < 1e-12 and 3 * Eq_field > 938.3, f"{Eq_base:.0f} -> {Eq_field:.0f} MeV")

    print("Verdict : la lecture 'quantum d'un champ universel' est coherente avec R99-R103 et les force : F = 1 par")
    print("construction, lambda-bar cinematique, le demi-quantum = secteur antiperiodique, la mere = secteur periodique,")
    print("Dirac (lambda-bar/2, 2 m c^2/hbar, L_z = hbar/2) sans moitie statique, mu_B et g = 2 par la structure du premier")
    print("ordre sans boucle. Son prix : le doublet du milieu doit etre un spineur sous les rotations (une addition, F),")
    print("et l'ancre des hadrons 2 pi lambda-bar devient une longueur posee (G). CONDITIONNEL.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
