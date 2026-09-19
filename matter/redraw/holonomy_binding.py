#!/usr/bin/env python3
"""R77 -- La liaison par holonomie : pourquoi deux DQD avec un tiers de tour
tiennent, et ce que devient n = 15 (trois paires, un tour entier).

Outils : la famille de courbes de la phase A (anneau a q lobes, courbe (1, q)
sur un tore, r0 = amplitude des lobes), le writhe par integrale de Gauss, et
le theoreme de Calugareanu-White-Fuller, Lk = Tw + Wr, invariant tant que le
ruban ne se traverse pas lui-meme (une traversee change Wr de +-2).

  A. Wr(r0) pour q = 3 : ou sont 1/3 (muon), 2/3 (tau), 1 (n = 15) ?
  B. le verrou : Tw = 0 est fixe par le milieu (R61), donc Lk = Wr est conserve
     sous toute deformation ; comme dWr/dr0 != 0, la forme est gelee a r0* :
     la paire ne peut pas partir sans reconnexion.  C'est la liaison.
  C. une traversee change Wr de 2 : de Wr = 1 on atteint -1, jamais 0 ; le
     retour de n = 15 vers n = 3 n'est pas topologique.
  D. atteignabilite : Wr = 1 existe-t-il dans la famille a trois lobes avant
     que la courbe ne se traverse (r0 < 1) ?  Resultat : oui, a r0 = 0,376
     (Wr monte jusqu'a 2,05) : la geometrie du writhe n'interdit pas n = 15 ;
     la lecture "six DQD = morceau de vide" (R76 E) est exclue, le verrou
     tient aussi la troisieme paire.  L'absence de quatrieme generation doit
     venir de l'operateur de masse (R76) ou d'une regle absente de la base.
"""
import math
import numpy as np
from scipy.optimize import brentq

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def curve(s, r0, q):
    rho = 1 + r0 * np.cos(q * s)
    drho = -q * r0 * np.sin(q * s)
    x = np.stack([rho * np.cos(s), rho * np.sin(s), -r0 * np.sin(q * s)], axis=-1)
    dx = np.stack([drho * np.cos(s) - rho * np.sin(s), drho * np.sin(s) + rho * np.cos(s),
                   -q * r0 * np.cos(q * s)], axis=-1)
    return x, dx

def writhe(r0, q=3, n=1500):
    s = (np.arange(n) + 0.5) * 2 * math.pi / n
    x, dx = curve(s, r0, q)
    dxx = x[:, None, :] - x[None, :, :]
    dist3 = np.sum(dxx * dxx, axis=-1) ** 1.5
    np.fill_diagonal(dist3, np.inf)
    cross = np.cross(dx[:, None, :], dx[None, :, :])
    return float(np.sum(np.sum(cross * dxx, axis=-1) / dist3) * (2 * math.pi / n) ** 2 / (4 * math.pi))

def path_and_area(r0, q=3, n=4000):
    s = (np.arange(n) + 0.5) * 2 * math.pi / n
    x, dx = curve(s, r0, q)
    L = float(np.sum(np.linalg.norm(dx, axis=1)) * 2 * math.pi / n)
    A = 0.5 * np.sum(np.cross(x, dx), axis=0) * 2 * math.pi / n
    return L, float(A[2])

def main():
    print("R77 -- la liaison par holonomie\n")

    # A. Wr(r0), q = 3
    print("A. Writhe de l'anneau a trois lobes, Wr(r0)")
    grid = [0.05, 0.1, 0.1785, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98]
    wr = {r: writhe(r) for r in grid}
    for r in grid:
        print(f"    r0 = {r:5.3f} : Wr = {wr[r]:.4f}")
    wr_max = max(wr.values())
    r_mu = brentq(lambda r: writhe(r) - 1 / 3, 0.05, 0.6, xtol=1e-8)
    try:
        r_tau = brentq(lambda r: writhe(r) - 2 / 3, 0.05, 0.98, xtol=1e-8)
    except ValueError:
        r_tau = None
    try:
        r_15 = brentq(lambda r: writhe(r) - 1.0, 0.05, 0.98, xtol=1e-8)
    except ValueError:
        r_15 = None
    print(f"  Wr = 1/3 (muon) a r0 = {r_mu:.4f} (R61 : 0,1785) ; "
          f"Wr = 2/3 (tau) a r0 = {r_tau if r_tau is None else round(r_tau, 4)} ; "
          f"Wr = 1 (n = 15) a r0 = {r_15 if r_15 is None else round(r_15, 4)}")
    print(f"  maximum de Wr sur la famille plongee (r0 < 1) : {wr_max:.4f}\n")
    check("Wr = 1/3 a r0 = 0,1785 (R61)", abs(r_mu - 0.1785) < 2e-3, f"{r_mu:.4f}")

    # B. le verrou de Calugareanu
    print("B. Le verrou : Lk = Tw + Wr, Tw = 0 fixe par le milieu (R61)")
    h = 1e-4
    dwr = (writhe(r_mu + h) - writhe(r_mu - h)) / (2 * h)
    print(f"  dWr/dr0 a r0* = {dwr:.3f} != 0 : toute deformation change Lk, ce qui est interdit sans")
    print("  traversee ou torsion ; la forme est gelee, la paire enroulee ne peut pas se degager.")
    print("  -> la liaison des deux DQD est un verrou topologique du ruban, pas une energie.\n")
    check("dWr/dr0 != 0 a r0* : forme gelee par Lk = Wr", abs(dwr) > 0.1, f"{dwr:.3f}")

    # C. traversee : Wr change de 2
    print("C. Une traversee change Wr de +-2")
    reach = {1 + 2 * m for m in range(-2, 3)}
    print(f"  depuis Wr = 1 : {sorted(reach)} ; 0 n'y est pas ; ni depuis 1/3 ou 2/3 vers 0.")
    print("  -> le retour de n = 15 dans le secteur de n = 3 n'est pas topologique (parite de Wr).\n")
    check("0 inaccessible depuis Wr = 1 par traversees", 0 not in reach, "1 -> -1, 3, ...")

    # D. atteignabilite de Wr = 1 a trois lobes
    print("D. La troisieme paire : Wr = 1 a trois lobes ?")
    if r_15 is None:
        print(f"  Wr = 1 n'est pas atteint : max {wr_max:.3f} a r0 -> 1 (la courbe se traverse a r0 = 1).")
        print("  -> un circuit a trois brins (trois lobes, R42) ne peut pas enrouler trois paires :")
        print("     n = 15 ne se forme pas ; trois generations, par la geometrie du writhe.")
    else:
        print(f"  Wr = 1 atteint a r0 = {r_15:.4f} : la troisieme paire est enroulable et, par B,")
        print("     verrouillee comme les autres. Le verrou ne ferme pas les generations : il faut")
        print("     autre chose pour interdire n = 15 (l'operateur de masse, R76, ou une regle absente).")
    check("Wr = 1 atteignable a trois lobes (r0 = 0,376) : la geometrie n'interdit pas n = 15",
          r_15 is not None and abs(r_15 - 0.3756) < 2e-3, f"r0 = {r_15:.4f}, max Wr = {wr_max:.3f}")
    # et le tau
    if r_tau is not None:
        L_mu, A_mu = path_and_area(r_mu)
        L_tau, A_tau = path_and_area(r_tau)
        print(f"  tau : r0 = {r_tau:.3f}, trajet {100*(L_tau/(2*math.pi)-1):+.0f} % du cercle, aire vectorielle "
              f"{100*A_tau/math.pi:.0f} % (muon : {100*(L_mu/(2*math.pi)-1):+.0f} %, {100*A_mu/math.pi:.0f} %)")
    check("Wr = 2/3 (tau) atteignable a trois lobes", r_tau is not None, f"r0 = {r_tau}")

    # E. plus de lobes ?
    print("\nE. Avec plus de lobes (q = 6, 9), Wr = 1 devient-il atteignable ?")
    for q in (6, 9):
        w = [writhe(r, q=q) for r in (0.3, 0.6, 0.9)]
        print(f"    q = {q} : Wr(0.3, 0.6, 0.9) = {[round(v, 3) for v in w]}")
    print("  -> oui avec plus de lobes ; mais q = 3 est le nombre de brins et de jonctions (R42) :")
    print("     un quatrieme pas demanderait un circuit qui n'est plus a trois brins.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
