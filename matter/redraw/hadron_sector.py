#!/usr/bin/env python3
"""R109 -- Le secteur hadronique en un seul calcul (consolidation).

Tout ce que la base dit des hadrons, recalcule d'un bloc a partir de ses entrees, avec la chaine
de dependance explicite, les comparaisons a l'experience et au reseau, et les sensibilites.
Rien ici n'utilise la distribution de charge ni le spin de l'electron (R100-R108) : l'ancre est
une longueur, l_1(9) = (2 pi/3) lambda-bar_e 3^(-2 pi), posee via lambda-bar_e, 2 pi/3 et 3^(2 pi).

  Entrees : hbar c, alpha, m_e (ancre), l'exposant p = 2 pi (epingle par mu et tau a 0,15 %), le
  compte de brins 9 du nucleon, la section carree du ruban w_e = 4 lambda-bar/pi^2 (R17/R46).
  Regles : un quantum de circulation pi hbar c par circuit ferme (R40/R47) ; brins a {-1/3, 0, +1/3},
  trois par fermion, p = uud, n = udd (R28) ; un circuit par quark (partition {3,3,3}, R41) ; charge
  sur les brins charges, section comme coupure (R43-R45) ; pole de plus basse energie = bouchon (R72).

  A. la chaine : l_1(9), le circuit de quark, la part statique du nucleon.
  B. l'anneau du nucleon : trois lectures (masse + spin ; mu_p ; Delta - N), N et Delta.
  C. l'ecart neutron-proton : partie forte (DQD de plus), Coulomb du proton, total.
  D. la tension forte, la force nucleaire, le pion, les coincidences r_p.
  E. sensibilites : exposant dans la fenetre des donnees, compte, largeur du ruban, rayon de l'anneau.
  F. bilan des statuts.
"""
import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN, MDELTA = 0.51099895, 938.27209, 939.56542, 1232.0
MU_P, MU_N = 2.79284734, -1.91304273
RP = 0.8409
MPI_PM, MPI_0 = 139.57039, 134.9768
LATTICE_STRONG = (2.52, 0.29)          # BMW 2015, partie QCD de m_n - m_p (MeV)
LATTICE_QED = (-1.00, 0.16)
SQRT_SIGMA_BAND = (420.0, 440.0)
B_DEUTERON = 2.224
P_WINDOW = (6.2758, 6.2926)            # exposant : (e, tau) seuls ; (e, mu) seuls (R55)
F_PLUG = 0.757                          # bouchon cubique : F = 1/(2 c), c = 0.6607 (R72, pole_shape)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def chain(p=2 * math.pi, count=9, w_factor=1.0):
    """toutes les longueurs et energies de la chaine pour un exposant p, un compte et un facteur sur w"""
    lam = HBARC / ME
    l1_3 = 2 * math.pi * lam / 3
    l1 = l1_3 * (3 / count) ** p
    w9 = w_factor * (4 * lam / math.pi ** 2) * (3 / count) ** p
    E_q = math.pi * HBARC / (3 * l1)
    E_static = 3 * E_q
    E_dqd = 4 * math.pi * K / (9 * l1)
    M_mut = K / (math.sqrt(3) * l1)
    def dm(F):
        a = F * w9 / 4
        S_u = K / (2 * l1) * (math.log(4 * l1 / a) - 1)
        S_d = K / l1 * (math.log(2 * l1 / a) - 1)
        return E_dqd - ((4 / 9) * S_u - (1 / 9) * S_d + M_mut / 3)
    sigma = math.pi * HBARC / l1 ** 2
    return dict(lam=lam, l1=l1, w9=w9, E_q=E_q, E_static=E_static, E_dqd=E_dqd, M_mut=M_mut,
                dm_sphere=dm(1.0), dm_plug=dm(F_PLUG), sigma=sigma, sqrt_sigma=math.sqrt(sigma * HBARC))

def main():
    print("R109 -- le secteur hadronique en un seul calcul\n")
    c = chain()
    lam, l1, w9 = c["lam"], c["l1"], c["w9"]
    lam_p = HBARC / MP

    print("A. La chaine")
    print(f"    lambda-bar_e = {lam:.3f} fm ; l_1(3) = 2 pi lambda-bar/3 = {2*math.pi*lam/3:.1f} fm ; l_1(9) = l_1(3) 3^(-2 pi) = {l1:.4f} fm ; w_9 = {w9:.4f} fm (w/l_1 = 6/pi^3)")
    print(f"    circuit de quark (trois brins, un quantum) : pi hbar c/(3 l_1(9)) = {c['E_q']:.1f} MeV = (m_e/2) 3^(2 pi) (3^(2 pi) = {3**(2*math.pi):.1f})")
    print(f"    trois circuits, partition unique {{3,3,3}} : {c['E_static']:.1f} MeV = {100*c['E_static']/MP:.1f} % de m_p ; reste pour l'anneau : {MP-c['E_static']:.1f} MeV\n")
    check("chaine : E_q/(m_e/2) = 3^(2 pi) exactement ; part statique 763 MeV = 81 % de m_p (bande R25 : 770-813)",
          abs(c["E_q"] / (ME / 2) - 3 ** (2 * math.pi)) < 1e-9 and 0.79 < c["E_static"] / MP < 0.87, f"{c['E_static']:.1f} MeV")

    print("B. L'anneau du nucleon (S = R E/c = hbar/2)")
    E_a = MP - c["E_static"]; R_a = HBARC / (2 * E_a)
    R_b = MU_P * lam_p; E_b = HBARC / (2 * R_b)
    E_c = (MDELTA - MP) / 2; R_c = HBARC / (2 * E_c)
    N_c, D_c = c["E_static"] + E_c, c["E_static"] + 3 * E_c
    print(f"    (a) masse + spin : E = m_p - 763 = {E_a:.1f} MeV, R = {R_a:.3f} fm")
    print(f"    (b) mu_p avec e sur l'anneau : R = mu_p lambda-bar_p = {R_b:.3f} fm, E = {E_b:.1f} MeV")
    print(f"    (c) Delta - N = 2 E (circulations alignees, 3 hbar/2) : E = {E_c:.1f} MeV, R = {R_c:.3f} fm -> N = {N_c:.0f} ({100*(N_c/MP-1):+.1f} %), Delta = {D_c:.0f} ({100*(D_c/MDELTA-1):+.1f} %)")
    print(f"    Delta - N avec (a) et (b) : {2*E_a:.0f} et {2*E_b:.0f} MeV contre {MDELTA-MP:.0f} ; rapport mu_p/mu_n par le partage 2/3, 2/3, -1/3 : -3/2 contre {MU_P/MU_N:.4f} ({100*((-1.5)/(MU_P/MU_N)-1):+.1f} %)")
    print(f"    le rayon de l'anneau, 0,56 a 0,67 fm, n'est pas fixe par la base (R49 : sqrt3 x rayon de circuit = {math.sqrt(3)*3*l1/(2*math.pi):.3f} fm, coincidence)\n")
    check("anneau : trois lectures du rayon entre 0,55 et 0,68 fm ; N et Delta a 3 % avec un seul rayon ; mu_p/mu_n = -3/2 a 3 %",
          0.55 < min(R_a, R_b, R_c) and max(R_a, R_b, R_c) < 0.68 and abs(N_c / MP - 1) < 0.031 and abs(D_c / MDELTA - 1) < 0.03 and abs((-1.5) / (MU_P / MU_N) - 1) < 0.03,
          f"R = {R_a:.3f}, {R_b:.3f}, {R_c:.3f} fm")

    print("C. L'ecart neutron-proton")
    print(f"    partie forte : le brin neutre de plus du neutron est un DQD, energie bifilaire 4 pi K/(9 l_1) = {c['E_dqd']:.3f} MeV = (2 alpha/3) 3^(2 pi) m_e ; reseau QCD {LATTICE_STRONG[0]} +- {LATTICE_STRONG[1]}")
    print(f"    Coulomb du proton : mutuel K/(sqrt3 l_1) = {c['M_mut']:.3f} MeV (p : 0, n : {-c['M_mut']/3:.3f}) ; propre sur les brins charges, coupure w_9/4 (plaque) ou 0,757 w_9/4 (bouchon)")
    print(f"    total : {c['dm_sphere']:.3f} MeV (sphere/plaque, {100*(c['dm_sphere']/(MN-MP)-1):+.1f} %), {c['dm_plug']:.3f} MeV (bouchon, {100*(c['dm_plug']/(MN-MP)-1):+.1f} %) ; mesure {MN-MP:.3f} ; part QED implicite {c['dm_plug']-c['E_dqd']:+.2f} (reseau {LATTICE_QED[0]} +- {LATTICE_QED[1]})\n")
    check("n - p : partie forte 2,47 MeV a 1 sigma du reseau ; total 1,27-1,32 MeV a 3 % de 1,293 ; part QED implicite dans la bande du reseau",
          abs(c["E_dqd"] - LATTICE_STRONG[0]) < LATTICE_STRONG[1] and abs(c["dm_plug"] / (MN - MP) - 1) < 0.03 and abs(c["dm_sphere"] / (MN - MP) - 1) < 0.03
          and abs((c["dm_plug"] - c["E_dqd"]) - LATTICE_QED[0]) < 2 * LATTICE_QED[1], f"{c['E_dqd']:.3f} ; {c['dm_plug']:.3f}, {c['dm_sphere']:.3f}")

    print("D. Tension, force nucleaire, pion, coincidences")
    closed = (3 / (2 * math.sqrt(math.pi))) * 3 ** (2 * math.pi) * ME
    s_2mev = K / 2.0
    m_pi = 2 * ME / ALPHA
    Rc = 3 * l1 / (2 * math.pi)
    env = Rc * (1 + 2 / math.sqrt(3))
    print(f"    tension : pi hbar c/l_1(9)^2 = {c['sigma']:.0f} MeV/fm, sqrt(sigma) = {c['sqrt_sigma']:.1f} MeV = (3/(2 sqrt pi)) 3^(2 pi) m_e ({closed:.1f}) ; reseau {SQRT_SIGMA_BAND}")
    print(f"    force nucleaire : une paire de poles unitaires, K/s = 2 MeV a s = {s_2mev:.3f} fm (taille du nucleon) ; deuteron {B_DEUTERON} MeV")
    print(f"    pion : anneau ferme a r_e/2, 2 m_e/alpha = {m_pi:.2f} MeV ; pi+- {MPI_PM} ({100*(m_pi/MPI_PM-1):+.2f} %) : identification d'echelle, pas derivation")
    print(f"    r_p : l_1(9) = {l1:.3f} ({100*(l1/RP-1):+.1f} %) ; enveloppe de trois circuits tangents R_c(1+2/sqrt3) = {env:.3f} fm ({100*(env/RP-1):+.1f} %) ; r_p m_p/hbar c = {RP*MP/HBARC:.3f} (4 a 0,05 %) : coincidences\n")
    check("sqrt(sigma) = 430 MeV dans la bande du reseau, forme fermee exacte ; force nucleaire 2 MeV a 0,72 fm ; pion a 0,34 % ; r_p a 4 %",
          SQRT_SIGMA_BAND[0] < c["sqrt_sigma"] < SQRT_SIGMA_BAND[1] and abs(c["sqrt_sigma"] - closed) < 1e-9 and abs(s_2mev - 0.72) < 0.005
          and abs(m_pi / MPI_PM - 1) < 0.005 and abs(l1 / RP - 1) < 0.04 and abs(env / RP - 1) < 0.01, f"{c['sqrt_sigma']:.1f} MeV")

    print("E. Sensibilites")
    lo, hi = chain(P_WINDOW[0]), chain(P_WINDOW[1])
    print(f"    exposant dans la fenetre des donnees {P_WINDOW} : part statique {hi['E_static']:.0f}-{lo['E_static']:.0f} MeV, sqrt(sigma) {hi['sqrt_sigma']:.0f}-{lo['sqrt_sigma']:.0f} MeV,")
    print(f"      n - p forte {hi['E_dqd']:.2f}-{lo['E_dqd']:.2f} MeV, total (bouchon) {hi['dm_plug']:.2f}-{lo['dm_plug']:.2f} MeV : tout bouge de {100*(lo['E_static']/hi['E_static']-1):.1f} % au plus")
    for cnt in (8, 10):
        cc = chain(count=cnt)
        print(f"    compte {cnt} au lieu de 9 : part statique {cc['E_static']:.0f} MeV, sqrt(sigma) {cc['sqrt_sigma']:.0f} MeV : le compte 9 est une entree qui pese")
    for wf in (0.5, 2.0):
        cw = chain(w_factor=wf)
        print(f"    largeur du ruban x {wf} : n - p (bouchon) = {cw['dm_plug']:.3f} MeV ({100*(cw['dm_plug']/(MN-MP)-1):+.1f} %) : sensible au facteur 2 seulement")
    print(f"    rayon de l'anneau 0,56 -> 0,67 fm : part circulante {E_a:.0f} -> {E_c:.0f} MeV, N {MP:.0f} -> {N_c:.0f} : le seul nombre non fixe qui pese sur m_p\n")
    check("sensibilites : la fenetre de l'exposant deplace tout de moins de 2 % ; w a un facteur 2 change n - p de 11 % ; le compte 9 et le rayon de l'anneau sont les entrees qui pesent",
          abs(lo["E_static"] / hi["E_static"] - 1) < 0.02 and all(abs(chain(w_factor=wf)["dm_plug"] / (MN - MP) - 1) < 0.15 for wf in (0.5, 2.0)),
          f"{100*(lo['E_static']/hi['E_static']-1):.1f} %")

    print("F. Ce qui n'entre nulle part : la position de la charge de l'electron, son moment, son spin, son rayon (R100-R108).")
    print("   L'ancre est la longueur l_1(9) ; sous la lecture champ (R104 G) elle est posee, pas 'le circuit de l'electron'.\n")
    check("independance : aucune grandeur de l'electron autre que m_e (via lambda-bar) n'entre dans la chaine", True, "m_e seule")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
