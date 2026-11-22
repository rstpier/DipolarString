#!/usr/bin/env python3
"""R127 -- La CD comme variante de la theorie des cordes : definition, dictionnaire, premieres consequences.

Decision de l'auteur (22 novembre 2026) : « On va transformer dipolar string comme une variante de la TdC ».
Definition retenue pour pouvoir calculer : le filament dipolaire recoit la loi elastique d'une corde relativiste,
l'action de Nambu-Goto (tension sigma, plus un terme de rigidite optionnel), avec des charges +q et -q a ses
extremites (c'est une corde ouverte a facteurs de Chan-Paton opposes : un dipole), et il est quantifie (l'axiome
d'amplitudes est alors d'office). R127 etait reserve a « la loi elastique du filament » : c'est elle.

  A. Identites. Avec la tension de la base sigma = pi hbar c / l_1(9)^2 : alpha' = 1/(2 pi sigma) = 0,859 GeV^-2,
     longueur de corde l_s = sqrt(hbar c / (2 pi sigma)) = 0,183 fm. La part statique de la base est hbar c/(sqrt2 l_s)
     = 763 MeV ; le rayon de la boucle fermee a un quantum (R120) est l_s ; l'anneau rigide (kappa_b = sigma l_s^2) a
     le rayon l_s/sqrt2, celui du demi-quantum. Toutes les longueurs de la CD « corde » sont l_s : aucun nombre nouveau.
  B. Coherence quantique. Corde ouverte bosonique : M^2 = (N - a)/alpha', a = (D-2)/24. A D = 4 : a = 1/12 ; fondamental
     tachyonique (|M| = 311 MeV) ; premier niveau massif (1033 MeV) avec D-2 = 2 etats, alors qu'un vecteur massif en
     exige D-1 = 3 : l'invariance de Lorentz est perdue. Un vecteur sans masse exige a = 1, soit D = 26 ; avec fermions
     sur le feuillet a = (D-2)/16 = 1/2, soit D = 10. Une corde FONDAMENTALE en D = 4 est incoherente. Une corde
     EFFECTIVE (longue devant l_s) est coherente en toute dimension (Polchinski-Strominger 1991) ; sa limite de validite
     a D = 4 est le rayon d'Arvis r_c = (pi/sqrt3) l_s = 0,33 fm.
  C. Gravitation. Le graviton est le niveau sans masse de la corde fermee (D = 26 ou 10) et son couplage est
     G ~ g_s^2 alpha' : avec la tension de la base il faut g_s = l_Planck / l_s = 9e-20 ; a couplage d'ordre 1 la
     gravite serait 1e38 fois trop forte. En espace plat, la meme corde ne peut pas faire les hadrons et la
     gravitation (rectificatif R128 : en dix dimensions une geometrie deformee reduit la tension effective et le permet).
  D. Le dipole. Corde ouverte a charges +q, -q : champ critique E_c = sigma/q au-dela duquel le champ arrache les
     extremites. q = e : 9,4e23 V/m ; q = e/3 : 2,8e24 V/m ; champ de Schwinger : 1,3e18 V/m. Le vide reel cree des
     paires e+e- 7e5 fois avant.
  E. Ligne de transmission et feuille d'univers. Chaque coordonnee de la corde est un champ sans masse a 1+1
     dimensions, c'est-a-dire une ligne de transmission ; bout ouvert = Neumann, court-circuit = Dirichlet. Pour une
     ligne d'impedance Z et de quantum de charge q, la phase compacte (manuscrit, eq. phase_action) est un boson
     compact de rayon R^2/alpha' = h/(q^2 Z). Vide adapte, q = e : R^2/alpha' = R_K/Z_0 = 1/(2 alpha) = 68,5.
     Dualite T, R -> alpha'/R : Z -> (h/q^2)^2 / Z, l'echange L <-> C ; point auto-dual Z = h/q^2 = 25,8 kOhm.
     Identite alpha = Z_0/(2 R_K) : la constante de structure fine lue comme un rayon, pas derivee.
  Remarque (R121). « 10 = 3 + 3 + 3 + 1 » : en cordes, capacitif/inductif = enroulement/impulsion, deux descriptions
  duales des MEMES dimensions compactes ; 3 + 3 dedouble trois dimensions, il n'en fait pas six (1 + 3 + 3 = 7, pas 10).
  (Rectificatif R128 : sous dix dimensions le compte se transpose autrement, 3 + 3 = trois tores.)
  Verdict : le dictionnaire est exact (IDENTITES) ; la CD comme corde fondamentale a quatre dimensions est EXCLUE ;
  une seule corde pour les hadrons et la gravite est EXCLUE en espace plat (1e38 ; voir R128) ; il reste une corde effective hadronique dont le
  seul nombre propre a la CD est la tension (0,185 GeV^2, deja 2 a 7 sigma sous le reseau, R119). « Variante de la
  TdC » designe donc aujourd'hui un dictionnaire et une corde effective : les chemins A et B de la feuille de route.
"""
import math

HBARC = 197.3269804            # MeV fm
ME = 0.51099895
ALPHA = 1 / 137.035999084
H, E_CH, C = 6.62607015e-34, 1.602176634e-19, 2.99792458e8
HBAR = H / (2 * math.pi)
ME_KG = 9.1093837015e-31
MU0 = 4 * math.pi * 1e-7 * (1 + 5.5e-10)      # CODATA 2018, a 1e-9 pres
Z0 = MU0 * C
L_PLANCK = 1.616255e-35        # m

LAM = HBARC / ME
L1 = (2 * math.pi * LAM / 3) * 3 ** (-2 * math.pi)       # l_1(9), fm
SIGMA = math.pi * HBARC / L1 ** 2                        # MeV/fm
SIGMA_GEV2 = SIGMA * HBARC * 1e-6
ALPHA_P = 1 / (2 * math.pi * SIGMA_GEV2)                 # GeV^-2
L_S = math.sqrt(HBARC / (2 * math.pi * SIGMA))           # fm

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def main():
    print("R127 -- la CD comme variante de la theorie des cordes : dictionnaire et premieres consequences\n")
    print(f"  tension de la base : sigma = {SIGMA:.1f} MeV/fm = {SIGMA_GEV2:.4f} GeV^2 ; alpha' = {ALPHA_P:.4f} GeV^-2 ; l_s = {L_S:.4f} fm\n")

    # A. identites
    print("A. Les identites : toutes les longueurs de la CD « corde » sont l_s")
    static = math.pi * HBARC / L1
    static_ls = HBARC / (math.sqrt(2) * L_S)
    r_loop_1 = math.sqrt(1.0 * HBARC / (2 * math.pi * SIGMA))
    r_loop_half = math.sqrt(0.5 * HBARC / (2 * math.pi * SIGMA))
    kappa_b = SIGMA * L_S ** 2                               # rigidite de courbure, MeV fm
    r_rigid = math.sqrt(kappa_b / (2 * SIGMA))               # minimum de 2 pi sigma R + pi kappa_b / R
    print(f"    part statique : pi hbar c / l_1(9) = {static:.2f} MeV ; hbar c/(sqrt2 l_s) = {static_ls:.2f} MeV")
    print(f"    boucle fermee, un quantum (R120) : R* = {r_loop_1:.4f} fm = l_s ; demi-quantum : {r_loop_half:.4f} fm = l_s/sqrt2")
    print(f"    anneau rigide, kappa_b = sigma l_s^2 = {kappa_b:.1f} MeV fm : R* = sqrt(kappa_b/2 sigma) = {r_rigid:.4f} fm = l_s/sqrt2\n")
    check("A. identites : part statique = hbar c/(sqrt2 l_s) ; rayons de boucle l_s et l_s/sqrt2 ; l'anneau rigide retombe sur l_s/sqrt2",
          abs(static / static_ls - 1) < 1e-12 and abs(r_loop_1 / L_S - 1) < 1e-12 and abs(r_rigid / r_loop_half - 1) < 1e-12,
          f"{static:.2f} MeV ; {r_loop_1:.4f}, {r_rigid:.4f} fm")

    # B. coherence quantique
    print("B. La coherence quantique : l'ordonnee a = (D-2)/24 et la dimension critique")
    def a_bos(D): return (D - 2) / 24
    D = 4
    m2_0 = (0 - a_bos(D)) / ALPHA_P; m2_1 = (1 - a_bos(D)) / ALPHA_P
    D_bos = 2 + 24 * 1                 # a = 1
    D_sup = 2 + 16 * 0.5               # a = (D-2)/16 = 1/2
    r_c = math.sqrt(math.pi * (D - 2) / (12 * SIGMA / HBARC))     # Arvis, fm
    print(f"    D = 4 : a = {a_bos(D):.4f} ; N = 0 : M^2 = {m2_0:.4f} GeV^2 (tachyon, |M| = {1000*math.sqrt(-m2_0):.0f} MeV) ; N = 1 : M = {1000*math.sqrt(m2_1):.0f} MeV, {D-2} etats pour {D-1} requis")
    print(f"    vecteur sans masse : a = 1 -> D = {D_bos} (bosonique) ; avec fermions de feuillet a = (D-2)/16 = 1/2 -> D = {D_sup:.0f}")
    print(f"    corde effective : valable pour r >> l_s ; rayon d'Arvis a D = 4 : r_c = {r_c:.4f} fm = {r_c/L_S:.4f} l_s (pi/sqrt3 = {math.pi/math.sqrt(3):.4f})\n")
    check("B. corde fondamentale a D = 4 incoherente (a = 1/12, vecteur massif a 2 etats) ; D critique 26 et 10 ; corde effective valable au-dela de 0,33 fm",
          abs(a_bos(4) - 1 / 12) < 1e-15 and D_bos == 26 and D_sup == 10 and (D - 2) < (D - 1) and abs(r_c / L_S - math.pi / math.sqrt(3)) < 1e-12,
          f"a(4) = {a_bos(4):.4f} ; r_c = {r_c:.3f} fm")

    # C. gravitation
    print("C. La gravitation : G ~ g_s^2 alpha'")
    g_s = L_PLANCK / (L_S * 1e-15)
    ratio = (L_S * 1e-15 / L_PLANCK) ** 2
    print(f"    l_Planck / l_s = {g_s:.2e} : couplage de corde requis ; a g_s ~ 1 la gravite serait {ratio:.1e} fois trop forte")
    print("    -> en espace plat, la meme corde ne peut pas porter la tension des hadrons et le graviton (rectificatif R128 : une geometrie deformee a dix dimensions le permet).\n")
    check("C. tension hadronique + graviton en espace plat : exclu de 1e38 (ou g_s = 9e-20) ; rectificatif R128 pour dix dimensions deformees", g_s < 1e-18 and ratio > 1e37, f"g_s = {g_s:.1e} ; {ratio:.1e}")

    # D. le dipole
    print("D. Le dipole : corde ouverte a charges +q et -q, champ critique E_c = sigma/q")
    e_c = SIGMA * 1e6 / 1e-15                   # V/m pour q = e
    e_c3 = 3 * e_c
    e_s = ME_KG ** 2 * C ** 3 / (E_CH * HBAR)
    print(f"    q = e : E_c = {e_c:.2e} V/m ; q = e/3 : {e_c3:.2e} V/m ; Schwinger : {e_s:.2e} V/m ; rapport {e_c/e_s:.1e}")
    print("    -> le vide reel cree des paires bien avant que la corde hadronique ne casse ; en cordes, E_c est le champ maximal de Born-Infeld.\n")
    check("D. champ critique du dipole-corde 7e5 fois au-dessus de Schwinger (q = e)", 1e5 < e_c / e_s < 1e7, f"{e_c/e_s:.1e}")

    # E. ligne de transmission <-> boson compact
    print("E. La ligne de transmission comme boson compact : R^2/alpha' = h/(q^2 Z)")
    r_k = H / E_CH ** 2
    r2_e = r_k / Z0
    r2_e3 = 9 * r_k / Z0
    z_dual = r_k ** 2 / Z0
    print(f"    R_K = h/e^2 = {r_k:.1f} ohm ; Z_0 = {Z0:.3f} ohm ; alpha = Z_0/(2 R_K) = {Z0/(2*r_k):.9f} (1/{2*r_k/Z0:.6f})")
    print(f"    vide adapte, q = e : R^2/alpha' = {r2_e:.3f} = 1/(2 alpha) = {1/(2*ALPHA):.3f} ; q = e/3 : {r2_e3:.1f}")
    print(f"    dualite T : Z -> R_K^2/Z ; point auto-dual Z = R_K = {r_k/Z0:.2f} Z_0 ; dual du vide : {z_dual:.3e} ohm = {z_dual/Z0:.0f} Z_0 = Z_0/(2 alpha)^2\n")
    check("E. identite : la ligne du vide (q = e, Z_0) est un boson compact de rayon R^2/alpha' = 1/(2 alpha) = 68,5 ; auto-dual a Z = R_K",
          abs(r2_e * 2 * ALPHA - 1) < 1e-8 and abs(z_dual / Z0 - (1 / (2 * ALPHA)) ** 2) / (1 / (2 * ALPHA)) ** 2 < 1e-8,
          f"R^2/alpha' = {r2_e:.3f}")

    # F. le compte des dimensions
    print("F. « 10 = 3 + 3 + 3 + 1 » relu en cordes")
    phys = 1 + 3 + 3
    print(f"    capacitif / inductif = enroulement / impulsion : deux descriptions duales des memes dimensions compactes ; 3 + 3 dedouble 3 dimensions")
    print(f"    dimensions physiques : 1 + 3 + 3 = {phys}, pas 10 ; la supercorde en exige 10 dont 6 compactes : sous cette lecture le compte ne se transpose pas (R121) ; sous dix dimensions il devient trois tores (R128).\n")
    check("F. le dedoublement 3 + 3 donne 7 dimensions physiques, pas 10 : sous cette lecture le compte ne se transpose pas (R128 : trois tores)", phys == 7 and phys != 10, f"{phys}")

    print("Verdict : le dictionnaire CD <-> cordes est exact (identites A, E) et aucun nombre nouveau n'en sort (tout est l_s). La CD")
    print("comme corde FONDAMENTALE a quatre dimensions est EXCLUE (B) ; une seule corde pour les hadrons et la gravite est EXCLUE en espace plat (C ; rectificatif R128).")
    print("Reste une corde EFFECTIVE hadronique, dont le seul nombre propre a la CD est la tension, deja en tension avec le reseau (R119).")
    print("« Variante de la theorie des cordes » designe aujourd'hui un dictionnaire et une corde effective : chemins A et B de la feuille de route.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
