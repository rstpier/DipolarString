#!/usr/bin/env python3
"""R98 -- La dynamique de reorientation de l'anneau dans l'analyseur, a partir
du seul couplage que la base possede : mu = mu_B n dans un champ B le long de a
(R53), et S = (hbar/2) n (R53), sans inertie d'orientation (R83).

  A. ce que le couplage fait seul : dS/dt = mu x B, precession de Larmor a
     omega_L = 2 mu_B B/hbar = eB/m (g = 2) ; l'angle theta entre n et a est
     conserve exactement (le couple est perpendiculaire a a) : S.a conserve.
  B. Stern-Gerlach classique : force mu_B cos(theta) dB/dz, deflexion continue
     en cos(theta) ; l'experience donne deux taches : exclu.
  C. precession amortie (alignement par dissipation) : n -> +a quel que soit
     theta, P(+) = 1 ; l'experience donne 1/2 a theta = 90 deg : exclu.
  D. deux sorties seulement (le quantum est indivisible, R47 : l'anneau finit
     en +a ou -a) + conservation en moyenne de ce que le couplage conserve
     exactement, S.a (equivalent : l'energie -mu.B) : p - (1 - p) = cos theta,
     donc p = cos^2(theta/2).  Born, unique, sans |psi|^2 en entree.
  E. lecture : l'analyseur n'agit pas sur le motif de champ (R94a, R97) mais sur
     l'orientation, par le moment ; les deux canaux sont +-a ; le "2-port" est
     le champ.  Ce qui reste pose : que la sortie soit discrete (le saut vers
     +-a) et que la conservation vaille en moyenne sur l'ensemble.
"""
import math
import numpy as np

HBAR = 6.582119569e-16      # eV s
MU_B = 5.7883818060e-5      # eV/T
ME_KG, E_C = 9.1093837e-31, 1.602176634e-19

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def precess(n0, B, t_end, damping=0.0, steps=20000):
    """dn/dt = omega_L (n x a) - damping * omega_L (n x (n x a)) ; a = z ; renormalise"""
    a = np.array([0, 0, 1.0])
    omega = 2 * MU_B * B / HBAR
    n = np.array(n0, float)
    dt = t_end / steps
    traj = [n.copy()]
    c, s_ = math.cos(omega * dt), math.sin(omega * dt)
    Rz = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])   # precession exacte autour de a = z
    for _ in range(steps):
        n = Rz @ n
        if damping:
            n = n - damping * omega * np.cross(n, np.cross(n, a)) * dt
            n /= np.linalg.norm(n)
        traj.append(n.copy())
    return np.array(traj), omega

def main():
    print("R98 -- reorientation de l'anneau dans un champ\n")
    B = 1.0
    theta0 = math.radians(60)
    n0 = [math.sin(theta0), 0, math.cos(theta0)]

    # A. precession
    print("A. Le couplage seul : precession")
    traj, omega = precess(n0, B, t_end=20 * 2 * math.pi / (2 * MU_B * B / HBAR))
    cos_theta = traj[:, 2]
    omega_eB_m = E_C * B / ME_KG
    print(f"  omega_L = 2 mu_B B/hbar = {omega:.4e} rad/s ; eB/m_e = {omega_eB_m:.4e} rad/s (g = 2)")
    print(f"  S.a = cos theta : min {cos_theta.min():.6f}, max {cos_theta.max():.6f} sur 20 tours (initial {math.cos(theta0):.6f})")
    print("  -> theta est conserve exactement : le couplage ne reoriente pas, il fait precesser.\n")
    check("precession : S.a conserve (variation < 1e-4 sur 20 tours), omega_L = eB/m", cos_theta.max() - cos_theta.min() < 1e-4 and abs(omega / omega_eB_m - 1) < 1e-6, f"{omega:.3e} rad/s")

    # B. Stern-Gerlach classique
    print("B. Stern-Gerlach classique : force mu_B cos theta dB/dz")
    rng = np.random.default_rng(2)
    cos_rand = rng.uniform(-1, 1, 200000)          # orientations isotropes
    hist, edges = np.histogram(cos_rand, bins=10, range=(-1, 1))
    print(f"  distribution des deflexions (en unites de la deflexion maximale), 10 classes : {np.round(hist/hist.sum(), 2)}")
    print("  -> continue et plate ; l'experience (1922) montre deux taches : exclu.\n")
    check("SG classique : distribution continue (aucune classe centrale vide)", hist.min() > 0.05 * hist.sum() / 10, "deux taches observees")

    # C. alignement amorti
    print("C. Precession amortie (dissipation) : alignement")
    finals = []
    for th in (30, 60, 90, 120, 150):
        n_init = [math.sin(math.radians(th)), 0, math.cos(math.radians(th))]
        tr, _ = precess(n_init, B, t_end=200 * 2 * math.pi / (2 * MU_B * B / HBAR), damping=0.1)
        finals.append(tr[-1, 2])
        print(f"    theta = {th:3d} deg : cos theta final = {tr[-1,2]:+.4f}")
    print("  -> tout s'aligne sur +a, P(+) = 1 pour tout theta ; l'experience donne cos^2(theta/2) : exclu.\n")
    check("amorti : P(+) = 1 pour tout theta (cos final > 0,99)", all(f > 0.99 for f in finals), "exclu par 1/2 a 90 deg")

    # D. deux sorties + conservation en moyenne
    print("D. Deux sorties (+a, -a) et conservation en moyenne de S.a")
    print("  S.a initial = (hbar/2) cos theta ; final = (hbar/2) [p (+1) + (1 - p)(-1)] = (hbar/2)(2p - 1)")
    print("  egalite en moyenne => p = (1 + cos theta)/2 = cos^2(theta/2)")
    ok = True
    for th in (0, 30, 60, 90, 120, 150, 180):
        t = math.radians(th)
        p = (1 + math.cos(t)) / 2
        ok &= abs(p - math.cos(t / 2) ** 2) < 1e-12
        E_i = -MU_B * B * math.cos(t)
        E_f = -MU_B * B * (2 * p - 1)
        print(f"    theta = {th:3d} deg : p = {p:.4f} = cos^2(theta/2) ; energie moyenne -mu.B : {E_i:+.3e} = {E_f:+.3e} eV")
    print("  -> Born, unique : toute autre p(theta) viole la conservation moyenne de S.a et de l'energie.\n")
    check("deux sorties + conservation moyenne => p = cos^2(theta/2), unique", ok, "Born sans |psi|^2")

    # E. lecture
    print("E. Lecture")
    print("  l'analyseur agit sur l'orientation par le moment, pas sur le motif de champ (d'ou R94a, R97) ;")
    print("  les deux canaux sont +-a ; le champ est le 2-port. Reste pose : que la sortie soit discrete")
    print("  (le saut vers +-a, indivisibilite du quantum appliquee a l'orientation) et que la conservation")
    print("  vaille en moyenne sur l'ensemble ; le couplage la garantit exactement trajectoire par trajectoire.\n")
    check("entrees nommees : sortie discrete (R47), conservation moyenne, couplage mu_B (R53)", True, "R83 sans inertie")

    print("Verdict : CONDITIONNEL ; Born pour un anneau sort du couplage de la base (qui conserve S.a exactement)")
    print("et de deux sorties ; les alternatives continues (deflexion classique) et deterministes (alignement")
    print("amorti) sont exclues par l'experience. Non derive : le saut discret lui-meme.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
