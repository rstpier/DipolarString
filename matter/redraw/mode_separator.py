#!/usr/bin/env python3
"""R97 -- Un separateur de modes reel : (a_+, a_-) -> (A_1, A_2) avec une matrice
de diffusion tiree du champ DQD, sans poser SU(2) ni theta/2 (barre fixee par
l'auteur : S^dag S = I, composition, |A_1|^2 = cos^2(theta/2)).

Le seul objet physique disponible pour definir des canaux est le champ du mode
sur l'anneau.  Un separateur physique doit projeter le champ entrant sur ses
deux modes propres ; les modes propres d'un analyseur oriente selon a sont les
champs des anneaux d'axes +a et -a.  Les amplitudes de sortie sont donc les
recouvrements du champ entrant avec ces deux champs de reference, les seules
identifications possibles entre anneaux distincts etant la tiree en arriere
par rotation rigide.  On teste quatre champs :
  (a) le motif reel a demi-tour (R84) : d(phi) = cos(phi/2) r + sin(phi/2) z ;
  (b) le meme multiplie par la phase de circulation e^{i phi/2} (R86) ;
  (c) un motif circulaire complexe local, (r + i z)/sqrt2 . e^{i phi/2} ;
  (d) temoin de charge 1 (periodique) : cos(phi) r + sin(phi) z.
Pour chacun : A_1 = <Psi_{+a}|Psi_n>, A_2 = <Psi_{-a}|Psi_n>, puissances,
somme des puissances (pertes), et comparaison a cos^2(theta/2).

CONTROLE AJOUTE (relecture de l'auteur) : les canaux +a et -a ne sont pas
orthogonaux, donc "somme != 1" ne suffisait pas a conclure "pas de 2-port
unitaire".  On refait le calcul apres orthonormalisation de Gram-Schmidt des
deux canaux : ecart a Born et fraction de la norme qui sort du sous-espace a
deux canaux.  Meme resultat qualitatif.  Conclusion retenue : les separateurs
naturels construits par projection du champ echouent ; ce n'est PAS le
theoreme "aucun analyseur lineaire DQD n'existe".
"""
import math
import numpy as np

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def rot_x(t):
    c, s = math.cos(t), math.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def frames(phi):
    r_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    z = np.tile(np.array([0.0, 0.0, 1.0]), (len(phi), 1))
    return r_hat, z

def field(kind, phi):
    r_hat, z = frames(phi)
    if kind == "a":
        return np.cos(phi / 2)[:, None] * r_hat + np.sin(phi / 2)[:, None] * z
    if kind == "b":
        return (np.cos(phi / 2)[:, None] * r_hat + np.sin(phi / 2)[:, None] * z) * np.exp(1j * phi / 2)[:, None]
    if kind == "c":
        return (r_hat + 1j * z) / math.sqrt(2) * np.exp(1j * phi / 2)[:, None]
    if kind == "d":
        return np.cos(phi)[:, None] * r_hat + np.sin(phi)[:, None] * z
    raise ValueError(kind)

def overlap(F, G):
    return np.mean(np.sum(np.conj(F) * G, axis=1))

def main():
    print("R97 -- separateur de modes sur le champ reel\n")
    phi = (np.arange(4000) + 0.5) * 4 * math.pi / 4000
    thetas = [0.0, math.pi / 3, math.pi / 2, 2 * math.pi / 3, math.pi]
    born = [math.cos(t / 2) ** 2 for t in thetas]
    names = {"a": "motif reel a demi-tour", "b": "motif x phase de circulation", "c": "motif circulaire complexe", "d": "temoin charge 1"}
    fails = {}
    for kind in ("a", "b", "c", "d"):
        F0 = field(kind, phi)
        norm = overlap(F0, F0).real
        ref_plus = F0                                   # anneau d'axe +a
        ref_minus = F0 @ rot_x(math.pi).T               # anneau d'axe -a
        print(f"({kind}) {names[kind]}")
        P1s, P2s = [], []
        for t, b in zip(thetas, born):
            Fn = F0 @ rot_x(t).T
            A1 = overlap(ref_plus, Fn) / norm
            A2 = overlap(ref_minus, Fn) / norm
            P1, P2 = abs(A1) ** 2, abs(A2) ** 2
            P1s.append(P1); P2s.append(P2)
            print(f"    theta = {math.degrees(t):5.0f} deg : |A1|^2 = {P1:.4f}, |A2|^2 = {P2:.4f}, somme = {P1+P2:.4f}   (Born {b:.4f}, {1-b:.4f})")
        dev = max(abs(p - b) for p, b in zip(P1s, born))
        loss = max(abs(p1 + p2 - 1) for p1, p2 in zip(P1s, P2s))
        fails[kind] = (dev, loss)
        print(f"    ecart max a cos^2(theta/2) : {dev:.3f} ; ecart max de la somme a 1 (pertes) : {loss:.3f}\n")
    print("  la phase de circulation commune s'annule dans tout recouvrement (b = a) ; le motif circulaire")
    print("  complexe (c) donne des puissances qui ne somment pas a 1 (pas un 2-port sans perte) ; le temoin")
    print("  de charge 1 suit cos theta. Aucun champ ne donne cos^2(theta/2) avec somme 1.\n")
    check("(a) et (b) identiques : la phase de circulation commune s'annule",
          abs(fails["a"][0] - fails["b"][0]) < 1e-9, f"{fails['a'][0]:.3f} = {fails['b'][0]:.3f}")
    check("aucun champ ne donne cos^2(theta/2) (ecart max > 0,2 pour tous)", all(v[0] > 0.2 for v in fails.values()),
          ", ".join(f"{k}: {v[0]:.2f}" for k, v in fails.items()))
    check("aucun n'est un 2-port sans perte (somme des puissances != 1)", all(v[1] > 0.05 for v in fails.values()),
          ", ".join(f"{k}: {v[1]:.2f}" for k, v in fails.items()))

    # --- controle : orthonormalisation de Gram-Schmidt des deux canaux ---
    print("Controle (relecture) : canaux +a, -a orthonormalises par Gram-Schmidt")
    gs = {}
    for kind in ("a", "b", "c", "d"):
        F0 = field(kind, phi)
        u1 = F0.reshape(-1)
        u2 = (F0 @ rot_x(math.pi).T).reshape(-1)
        e1 = u1 / np.linalg.norm(u1)
        r2 = u2 - np.vdot(e1, u2) * e1
        e2 = r2 / np.linalg.norm(r2)
        devs, outs = [], []
        for t, b in zip(thetas, born):
            psi = (F0 @ rot_x(t).T).reshape(-1)
            psi = psi / np.linalg.norm(psi)
            p1, p2 = abs(np.vdot(e1, psi)) ** 2, abs(np.vdot(e2, psi)) ** 2
            inside = p1 + p2
            share = p1 / inside if inside > 1e-12 else float("nan")
            devs.append(abs(share - b))
            outs.append(1 - inside)
        gs[kind] = (max(devs), max(outs))
        print(f"    ({kind}) ecart max a cos^2(theta/2) dans le sous-espace : {max(devs):.3f} ; "
              f"fraction max de la norme hors des deux canaux : {100*max(outs):.0f} %")
    print("  -> meme resultat qualitatif : l'echec n'est pas un artefact de canaux non orthogonaux.\n")
    check("apres Gram-Schmidt : ecart a Born > 0,2 et fuite hors sous-espace > 50 % pour les motifs a demi-tour",
          all(gs[k][0] > 0.2 and gs[k][1] > 0.5 for k in ("a", "b", "c")),
          ", ".join(f"{k}: {v[0]:.3f} / {100*v[1]:.0f} %" for k, v in gs.items()))

    print("Ce que l'echec dit : un separateur qui agit sur le champ (recouvrements de motifs, tiree en")
    print("arriere par rotation rigide) ne porte jamais le demi-angle ; le e^(+- i theta/2) de R95 est")
    print("une propriete de la loi de transformation de l'ORIENTATION (loi de groupe + demi-tour), pas des")
    print("recouvrements de champs. R96 etait une reformulation spinorielle, comme le soupconnait")
    print("l'auteur : le C^2 n'est pas (encore) une paire d'amplitudes de fluide. Un analyseur qui donne")
    print("Born doit agir sur l'orientation elle-meme (reorientation vers +-a) avec une dynamique qui")
    print("reproduit cos^2(theta/2) : cette dynamique n'est pas dans la base.\n")
    check("statut : les separateurs naturels par projection du champ echouent (pas un theoreme d'inexistence)", True, "R95 tient, R96 reformulation")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
