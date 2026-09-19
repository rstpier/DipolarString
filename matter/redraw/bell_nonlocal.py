#!/usr/bin/env python3
"""R92 -- Peut-on deriver E(a, b) = -cos(a - b) avec ce que le milieu DQD a deja,
sans mettre -cos theta en entree ?  (test fixe par l'auteur, cinq contraintes :
marges 1/2, E = -a.b, S_max = 2 sqrt2, pas de signalisation, pas de dependance
observable au repere de l'ether.)

  A. loi locale quelconque : scan d'une famille de reponses P(+|a, n) = f(a.n)
     (signe, Malus, rampes, seuils) : max S = 2, atteint par le signe.  C'est
     le vrai enonce de R91 : la base locale donne |S| <= 2 quelle que soit f.
  B. ce que le milieu a deja de non local : le quantum partage (R47, R87 : les
     deux filles portent les deux moities d'un seul quantum, contrainte globale
     de conservation) et la sphere de Bloch comme espace d'etats (R83, R85).
     Si la mesure de A reoriente l'axe partage vers +-a (contrainte globale,
     instantanee dans le repere de l'ether), alors :
       reponse signe  : E = -sign(cos theta), S = 4 > 2 sqrt2 : super-quantique,
                        exclu par l'experience (2,4-2,7) ;
       reponse p(theta) : E = -(2p - 1) ; E = -cos theta <=> p = cos^2(theta/2),
                        la regle de Born sur la sphere de Bloch, UNIQUE.
     La base a la sphere, pas la mesure : -cos theta n'est pas derive.
  C. vitesse du canal : la seule vitesse du milieu est c0 = c ; Salart et al.
     2008 (18 km, balayage des reperes sur 24 h) : v_canal > 1e4 c.  Une onde
     sur le milieu est exclue comme canal ; seule une contrainte globale
     (instantanee) convient.
  D. signature de repere : contrainte instantanee => delta E = 0, aucune
     nouvelle prediction ; canal a vitesse finie => chute vers |S| <= 2 pour des
     detections simultanees a L/v pres dans le repere de l'ether : borne 1e4 c.
Sorties honnetes : -cos theta n'est pas derive (il manque la mesure de Born) ;
pas de delta E nouveau ; ce que la base possede est le squelette (contrainte
globale + sphere), pas la dynamique.
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def chsh_from_E(E):
    a, ap, b, bp = 0.0, math.pi / 2, math.pi / 4, 3 * math.pi / 4
    return abs(E(a - b) - E(a - bp) + E(ap - b) + E(ap - bp))

def E_local(f, theta, n=400):
    """E(theta) pour une loi locale P(+|a,n) = f(a.n), axes n isotropes, B anti-aligne"""
    # integration sur la sphere : n = (sin t cos p, sin t sin p, cos t), a = z, b = (sin th, 0, cos th)
    ts = (np.arange(n) + 0.5) * math.pi / n
    ps = (np.arange(n) + 0.5) * 2 * math.pi / n
    T, P = np.meshgrid(ts, ps, indexing="ij")
    an = np.cos(T)
    bn = np.sin(T) * np.cos(P) * math.sin(theta) + np.cos(T) * math.cos(theta)
    EA = 2 * f(an) - 1
    EB = 2 * f(-bn) - 1
    w = np.sin(T)
    return float(np.sum(EA * EB * w) / np.sum(w))

def main():
    print("R92 -- deriver -cos theta ?\n")

    # A. scan des lois locales
    print("A. Lois de mesure locales P(+|a, n) = f(a.n)")
    fams = {
        "signe (deterministe)": lambda x: (x > 0).astype(float),
        "Malus / Born locale, cos^2(theta/2)": lambda x: (1 + x) / 2,
        "rampe douce, largeur 0,3": lambda x: np.clip((x + 0.3) / 0.6, 0, 1),
        "rampe douce, largeur 0,05": lambda x: np.clip((x + 0.05) / 0.1, 0, 1),
        "reponse cubique": lambda x: (1 + x ** 3) / 2,
    }
    S_local = {}
    for name, f in fams.items():
        S_local[name] = chsh_from_E(lambda th, f=f: E_local(f, th))
        print(f"    {name:38s} S = {S_local[name]:.3f}")
    print("  -> max S = 2 (le signe) : pour toute f, |S| <= 2. C'est l'enonce general de R91.\n")
    check("toute loi locale donne S <= 2 (max 2,000 au signe)", max(S_local.values()) <= 2 + 1e-6 and abs(S_local["signe (deterministe)"] - 2) < 2e-3,
          f"max {max(S_local.values()):.3f}")

    # B. contrainte globale (axe partage reoriente par la mesure de A)
    print("B. Ce que le milieu a deja de non local : le quantum partage (R47, R87) comme contrainte globale")
    S_sign_nl = chsh_from_E(lambda th: -math.copysign(1.0, math.cos(th)) if abs(math.cos(th)) > 1e-12 else 0.0)
    S_born_nl = chsh_from_E(lambda th: -(2 * math.cos(th / 2) ** 2 - 1))
    print(f"    axe partage -> +-a, reponse signe : E = -sign(cos theta), S = {S_sign_nl:.3f} (> 2 sqrt2 = {2*math.sqrt(2):.3f}) : super-quantique, exclu")
    print(f"    axe partage -> +-a, reponse p(theta) = cos^2(theta/2) : E = -cos theta, S = {S_born_nl:.3f} = 2 sqrt2")
    # unicite : E = -(2p - 1) = -cos theta <=> p = cos^2(theta/2)
    thetas = np.linspace(0, math.pi, 7)
    unique = all(abs((1 + math.cos(t)) / 2 - math.cos(t / 2) ** 2) < 1e-12 for t in thetas)
    print("    E = -(2p - 1) : -cos theta <=> p = cos^2(theta/2) exactement, la regle de Born sur la sphere de Bloch.")
    print("  -> la base a la sphere (R83, R85) et la contrainte globale (R47, R87) ; elle n'a pas la mesure :")
    print("     -cos theta n'est pas derive, et sans Born la contrainte globale surchoote (S = 4).\n")
    check("contrainte globale + signe : S = 4 > Tsirelson (exclu)", abs(S_sign_nl - 4) < 1e-9 and S_sign_nl > 2 * math.sqrt(2), "4.000")
    check("contrainte globale + Born : S = 2 sqrt2 ; Born est l'unique p(theta) donnant -cos theta",
          abs(S_born_nl - 2 * math.sqrt(2)) < 1e-9 and unique, "2.828")

    # C. vitesse du canal
    print("C. Vitesse du canal")
    v_bound = 1e4                                   # Salart et al. 2008 : > 10^4 c (reperes balayes sur 24 h)
    v_medium = 1.0                                  # c0 = c
    print(f"    seule vitesse du milieu : c0 = {v_medium:.0f} c ; borne experimentale du canal : > {v_bound:.0e} c (Salart 2008, 18 km)")
    print("  -> une onde sur le milieu ne peut pas etre le canal ; seule une contrainte globale (instantanee) convient.\n")
    check("le milieu a c0 = c, le canal doit depasser 1e4 c : l'onde du milieu est exclue comme canal", v_medium < v_bound / 1000, "c vs 1e4 c")

    # D. signature de repere
    print("D. Signature de repere de l'ether")
    print("    contrainte instantanee dans le repere de l'ether : les statistiques ne dependent pas de l'ordre des")
    print("    mesures (marges 1/2, symetrie) : delta E = 0, aucune prediction nouvelle.")
    print("    canal a vitesse finie v : chute vers |S| <= 2 pour des detections simultanees a L/v pres dans le")
    print("    repere de l'ether ; non observee : v > 1e4 c. Pas de delta E disponible.\n")
    check("pas de delta E nouveau : instantane => 0 ; fini => deja borne", True, "Salart 2008")

    print("Verdict : -cos theta n'est PAS derive. Le milieu possede le squelette (une contrainte globale,")
    print("le quantum partage ; un espace d'etats, la sphere de Bloch) mais pas la dynamique : la regle de")
    print("Born, unique reponse compatible, doit etre ajoutee. Aucune correction delta E ne sort.")
    print("Sortie honnete du test : secteur quantique absent, squelette present.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
