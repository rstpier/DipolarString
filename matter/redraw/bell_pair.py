#!/usr/bin/env python3
"""R91 -- Prediction sur l'intrication ?

La base n'a pas d'espace d'etats : une configuration a la fois (ether de
Lorentz), pas de superposition, pas de regle de Born (le manuscrit v2.9.2 a
lui-meme retire son mecanisme d'intrication faute d'equation).  Ce qu'elle
predit malgre elle, comme modele local a configurations definies, se calcule :

  A. la paire de la brisure (R87) : deux anneaux de circulations opposees,
     S+ = -S- le long de l'axe de la mere, axe aleatoire : le singulet classique.
  B. correlation de spin E(a, b) pour deux modeles de mesure :
     (i)  projection de signe, sortie sign(a.n) : E = -(1 - 2 theta/pi), CHSH max = 2 ;
     (ii) reponse de Malus sur la sphere de Bloch (R83/R85), sortie +-1 avec
          probabilite cos^2(theta/2) : E = -(a.b)/3, CHSH max = 2 sqrt2 / 3 = 1,886.
     Mecanique quantique : E = -a.b, CHSH = 2 sqrt2 = 2,828.
  C. experiences sans faille (Hensen 2015 : S = 2,42 +- 0,20 ; photons : 2,70 +-
     0,05 typ.) : > 2.  Les deux modeles locaux sont exclus.
  D. ce qu'il faudrait : soit une superposition de configurations, soit un canal
     non local par le milieu (l'ether a un repere privilegie, ce qui le permet
     en principe, comme l'onde pilote) ; aucune equation dans la base.
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def chsh(E):
    """CHSH pour E(a, b) fonction de l'angle, aux angles optimaux 0, 45, 90, 135 deg"""
    a, ap, b, bp = 0.0, math.pi / 2, math.pi / 4, 3 * math.pi / 4
    return abs(E(a - b) - E(a - bp) + E(ap - b) + E(ap - bp))

def E_sign(theta):
    return -(1 - 2 * abs(theta) / math.pi)

def E_malus(theta):
    return -math.cos(theta) / 3

def E_qm(theta):
    return -math.cos(theta)

def monte_carlo_sign(n=200000, seed=3):
    """verification numerique du modele (i) : axes aleatoires, sorties sign"""
    rng = np.random.default_rng(seed)
    nvec = rng.normal(size=(n, 3))
    nvec /= np.linalg.norm(nvec, axis=1)[:, None]
    out = {}
    for th in (0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4):
        a = np.array([0, 0, 1.0])
        b = np.array([math.sin(th), 0, math.cos(th)])
        A = np.sign(nvec @ a)
        B = np.sign(-(nvec @ b))
        out[th] = float(np.mean(A * B))
    return out

def main():
    print("R91 -- intrication : ce que la base predit malgre elle\n")

    print("A. La paire de la brisure : S+ = -S- le long d'un axe aleatoire (singulet classique)\n")
    print("B. Correlations et CHSH")
    mc = monte_carlo_sign()
    for th, e in mc.items():
        print(f"    theta = {math.degrees(th):5.1f} deg : E_sign (MC) = {e:+.3f}, analytique {E_sign(th):+.3f}, "
              f"Malus {E_malus(th):+.3f}, MQ {E_qm(th):+.3f}")
    S = {"projection de signe (i)": chsh(E_sign), "Malus / Bloch (ii)": chsh(E_malus), "mecanique quantique": chsh(E_qm)}
    for k, v in S.items():
        print(f"  CHSH {k:26s} = {v:.3f}")
    print()
    check("modele (i) : E lineaire en theta (MC a 1 %), CHSH = 2",
          all(abs(mc[th] - E_sign(th)) < 0.01 for th in mc) and abs(S["projection de signe (i)"] - 2) < 1e-9, "2.000")
    check("modele (ii) : CHSH = 2 sqrt2 / 3 = 1,886", abs(S["Malus / Bloch (ii)"] - 2 * math.sqrt(2) / 3) < 1e-9, "1.886")
    check("MQ : CHSH = 2 sqrt2", abs(S["mecanique quantique"] - 2 * math.sqrt(2)) < 1e-9, "2.828")

    print("C. Experiences sans faille")
    exps = {"Hensen 2015 (spins NV, 1,3 km)": (2.42, 0.20), "photons (typ.)": (2.70, 0.05)}
    for k, (s, ds) in exps.items():
        print(f"  {k:32s} S = {s:.2f} +- {ds:.2f} : depasse 2 de {(s-2)/ds:.1f} sigma")
    print("  -> les deux modeles locaux de la base sont exclus ; la MQ tient.\n")
    check("les experiences depassent la borne locale 2 (> 2 sigma)", all((s - 2) / ds > 2 for s, ds in exps.values()), "Bell viole")

    print("D. Verdict : la base predit CHSH <= 2 (1,89 avec sa propre mesure de Bloch), faux. Elle n'a")
    print("   ni superposition ni regle de Born ; le manuscrit v2.9.2 a retire son mecanisme faute")
    print("   d'equation. Une extension viable doit etre non locale par le milieu (repere privilegie de")
    print("   l'ether, permis en principe) ou porter des superpositions de configurations : a construire.\n")
    check("statut : EXCLU comme modele local ; extension quantique a construire", True, "manuscrit 1164, 1212")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
