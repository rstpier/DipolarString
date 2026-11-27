#!/usr/bin/env python3
"""R141 -- Chemin C : ce qu'un banc fini donnerait. Lecture en un point, pertes de vraies bobines, pieces a tolerance.

Suite de R140, meme decision de l'auteur (27 novembre 2026 : « pars la simulation maintenant »). R140 mesurait la loi
sur une ligne sans fin, avec une impulsion et une analyse ideale, et finissait sur une limite : le banc fidele commande
les diodes cellule par cellule et reste a simuler avant tout achat. C'est ce calcul. Il rectifie aussi deux ordres de
grandeur de R140 G.

Le banc simule. 195 cellules, fermees a chaque bout sur leur resistance (demi-condensateur en bout). Une sonde
continue entre par un bout et remplit la ligne ; elle est coupee en douceur juste avant le depart du front, parce
qu'une source qui emettrait derriere le front ajouterait des ondes parasites. Le front, raide d'une demi-cellule,
part du meme bout : chaque cellule recoit sa commande a son instant, avec une erreur possible. La tension est lue a
la cellule 150, apres le passage du front, tant que la raie basse y passe. Trois sinusoides de frequences connues
sont ajustees : la raie haute, la raie basse, et la porteuse (la frequence dont la vitesse de phase est celle du
front). Le rapport des deux raies est compare a celui de R140.

Pertes. Une resistance en serie dans chaque bobine. Le facteur de qualite Qb est compte la ou un monteur le lit : a la
frequence des deux raies (41 MHz dans l'exemple de R140). La resistance d'une bobine monte comme la racine de la
frequence (effet de peau) : une cellule que le front n'a pas encore atteinte ne porte que la sonde, quatre fois plus
basse en frequence, et recoit la resistance de cette frequence ; une cellule deja convertie recoit celle des raies.

  A. Deux rectifications de R140 G : la convention du facteur de qualite, et le nombre de cellules.
  B. Banc ideal : les deux raies, l'onde parasite a la porteuse, et la lecture du rapport.
  C. Pertes de vraies bobines : niveau des raies selon Qb ; rapport corrige des pertes.
  D. Pieces a tolerance et instants de commande imprecis : dispersion de la lecture d'un banc a l'autre ; lequel compte.
  E. La temperature tiree de trois frequences de sonde.
  Verdict : le critere de pertes ecrit avant R140 n'est tenu qu'avec des bobines de Qb = 126 ou plus a la frequence des raies,
  et la loi ne se lit a 15 % qu'avec des pieces a 1 %. R140 disait « facteur 4 avec Q = 100 » : c'etait Qb = 277 dans
  la convention d'un monteur, sur une ligne deux fois trop courte.

Limites : une seule architecture et un seul protocole de lecture ; la fuite de la commande dans la ligne n'est pas
simulee (elle ne depend pas de la sonde et se retire en inversant la sonde d'un tir sur deux) ; les pertes des diodes
sont comptees dans Qb ; aucune piece reelle n'a ete choisie.
"""
import random
from math import cos, exp, log, pi, sin, sqrt, tanh
import ds_analogue_line as al

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

U, CB, DT = al.U, al.CB, al.DT
DELTA = 1 - 1 / CB ** 2            # baisse relative de la capacite derriere le front
WB = 2 * CB                        # frequence de coupure derriere le front
W = 0.5                            # raideur du front, en cellules
N, TAP = 195, 150                  # cellules ; point de mesure
EPS = 1e-3                         # amplitude de la sonde sur la ligne

def modes(x, T):
    """Sonde de frequence x T vue du front : nombres d'onde, frequences au point fixe et vitesses de groupe."""
    wt = x * T
    kp = al.root(lambda k: al.wfront(k, 2.0) - wt, 1e-6, 2.0)
    ks = al.root(lambda k: al.wfront(k, WB), 1.5, 3.1)
    khi = al.root(lambda k: al.wfront(k, WB) - wt, ks, 3.1)
    klo = al.root(lambda k: al.wfront(k, WB) + wt, 1.3, ks)
    f = lambda k, w: al.theta(k, w) / DT
    return dict(wt=wt, kp=kp, ks=ks, khi=khi, klo=klo, om=f(kp, 2.0), ohi=f(khi, WB), olo=f(klo, WB), oc=f(ks, WB),
                vp=al.vgroup(kp, 2.0), vhi=al.vgroup(khi, WB), vlo=al.vgroup(klo, WB), vc=al.vgroup(ks, WB))

def shot(x, T, qb=None, skin=True, tol=0.0, jit=0.0, seed=1, wcut=3.0):
    """Un tir du banc. Renvoie les modes, la tension lue au point de mesure depuis le passage du front, et les taux de pertes."""
    rnd = random.Random(seed)
    Ln = [1 + tol * rnd.uniform(-1, 1) for _ in range(N - 1)]
    C0 = [1 + tol * rnd.uniform(-1, 1) for _ in range(N)]
    C0[0] *= 0.5; C0[-1] *= 0.5                                       # demi-condensateurs en bout
    dl = [DELTA * (1 + tol * rnd.uniform(-1, 1)) for _ in range(N)]   # excursion de capacite de chaque diode
    tj = [jit * rnd.uniform(-1, 1) for _ in range(N)]                 # erreur sur l'instant de commande, en cellules
    md = modes(x, T)
    gl = md["oc"] / qb if qb else 0.0                                 # taux de pertes derriere le front (frequence des raies)
    gp = (sqrt(md["om"] * md["oc"]) / qb if skin else gl) if qb else 0.0   # devant le front (frequence de la sonde)
    ca = lambda g: [(1 - 0.5 * g * DT / l) / (1 + 0.5 * g * DT / l) for l in Ln]
    cb = lambda g: [DT / (1 + 0.5 * g * DT / l) for l in Ln]
    fa_p, fb_p, fa_l, fb_l = ca(gp), cb(gp), ca(gl), cb(gl)
    iL = [1 / l for l in Ln]
    t0 = N / md["vp"] + 60.0                                          # depart du front : la sonde a rempli la ligne
    t_f = t0 + TAP / U                                                # passage du front au point de mesure
    t_low = TAP * (1 / md["vlo"] - 1 / U)                             # duree de la raie basse au point de mesure
    nsteps = int((t_f + 2.0 + t_low) / DT) + 2
    Qc = [0.0] * N; P = [0.0] * (N - 1); rec = []
    om = md["om"]
    for m in range(nsteps):
        t = m * DT; xf = U * (t - t0)
        lo, hi = max(0, min(N, int(xf - 12 * W) - 1)), max(0, min(N, int(xf + 12 * W) + 2))
        invC = [1 / (C0[n] * (1 - dl[n])) for n in range(lo)] \
             + [1 / (C0[n] * (1 - dl[n] * 0.5 * (1 - tanh((n - xf + tj[n]) / W)))) for n in range(lo, hi)] \
             + [1 / C0[n] for n in range(hi, N)]
        V = [q * ic for q, ic in zip(Qc, invC)]
        if t >= t_f + 2.0: rec.append(V[TAP])
        jf = max(0, min(N - 1, int(xf + 0.5)))                        # bobines deja depassees par le front
        P = [pp * a + b * (va - vb) for pp, a, b, va, vb in zip(P, fa_l[:jf] + fa_p[jf:], fb_l[:jf] + fb_p[jf:], V, V[1:])]
        I = [pp * l for pp, l in zip(P, iL)]
        th = t + 0.5 * DT
        vs = 2 * EPS * sin(om * th) * 0.5 * (1 + tanh((th - 30.0) / 6.0)) * 0.5 * (1 - tanh((th - (t0 - 1.0)) / wcut))
        a0 = 0.5 * DT * invC[0]; aN = 0.5 * DT * invC[-1]             # bouts fermes sur une resistance 1
        Qc = [(Qc[0] * (1 - a0) + DT * (vs - I[0])) / (1 + a0)] + [q + DT * (a - b) for q, a, b in zip(Qc[1:-1], I, I[1:])] \
           + [(Qc[-1] * (1 - aN) + DT * I[-1]) / (1 + aN)]
    return md, rec, gp, gl

def lsq(cols, y):
    n = len(cols)
    Mx = [[sum(a * b for a, b in zip(cols[i], cols[j])) for j in range(n)] + [sum(a * v for a, v in zip(cols[i], y))] for i in range(n)]
    for i in range(n):
        piv = max(range(i, n), key=lambda q: abs(Mx[q][i])); Mx[i], Mx[piv] = Mx[piv], Mx[i]
        for q in range(i + 1, n):
            f = Mx[q][i] / Mx[i][i]
            for c in range(i, n + 1): Mx[q][c] -= f * Mx[i][c]
    sol = [0.0] * n
    for i in range(n - 1, -1, -1): sol[i] = (Mx[i][n] - sum(Mx[i][j] * sol[j] for j in range(i + 1, n))) / Mx[i][i]
    return sol

def flux(k):
    """Flux d'action a travers le front, pour une onde d'amplitude 1 en tension derriere lui."""
    return abs(U - al.vgroup(k, WB)) * cos(al.theta(k, WB) / 2) / (WB * sin(k / 2))

def read(md, rec, gp, gl, wcut=3.0, carrier=True):
    """Ajuste les raies sur l'enregistrement. Les amplitudes sont ramenees au debut de l'enregistrement, pertes corrigees."""
    ngate = 2 * wcut * md["vp"] * U / (U - md["vp"]) + 5              # cellules ou la sonde etait affaiblie par sa coupure
    tuse = (TAP - ngate) * (1 / md["vlo"] - 1 / U) - 3.0
    y = rec[:int(tuse / DT)]; ts = [i * DT for i in range(len(y))]
    cols = []
    for o, v in ((md["ohi"], md["vhi"]), (md["olo"], md["vlo"]), (md["oc"], md["vc"]))[:3 if carrier else 2]:
        beta = (gl / v - gp / md["vp"]) / (1 / v - 1 / U)             # une onde lue plus tard a voyage plus longtemps
        d = [exp(-0.5 * beta * t) for t in ts]
        cols += [[a * cos(o * t) for a, t in zip(d, ts)], [a * sin(o * t) for a, t in zip(d, ts)]]
    sol = lsq(cols, y)
    amp = [sqrt(sol[2 * j] ** 2 + sol[2 * j + 1] ** 2) / EPS for j in range(len(cols) // 2)]
    return dict(hi=amp[0], lo=amp[1], car=amp[2] if carrier else 0.0, window=tuse, level=sum(v * v for v in y) / len(y) / EPS ** 2,
                ratio=(amp[1] / amp[0]) ** 2 * flux(md["klo"]) / flux(md["khi"]))

def one(x, T, wcut=3.0, carrier=True, **kw):
    md, rec, gp, gl = shot(x, T, wcut=wcut, **kw)
    return read(md, rec, gp, gl, wcut=wcut, carrier=carrier)

def temperature(xs, ratios):
    """Temperature lue : moins l'inverse de la pente de ln(rapport) contre x."""
    n = len(xs); mx = sum(xs) / n; my = sum(log(r) for r in ratios) / n
    return -sum((x - mx) ** 2 for x in xs) / sum((x - mx) * (log(r) - my) for x, r in zip(xs, ratios))

def main():
    print("R141 -- chemin C : ce qu'un banc fini donnerait\n")
    T = al.Line(W, False, 10, 10).kappa / (2 * pi)
    XS = [1.5, 2.0, 2.5]
    ref = al.experiment(W, False, XS)
    refr = {x: nlo / nhi for x, (ni, nhi, nlo) in zip(XS, ref["out"])}
    t_ref = temperature(XS, [refr[x] for x in XS])
    md2 = modes(2.0, T)
    print(f"Front d'une demi-cellule, ligne a capacite seule : T = {T:.4f}. Reference de R140, rapport des actions raie basse sur raie haute : "
          + " ; ".join(f"x = {x} : {refr[x]:.4f}" for x in XS) + f" ; temperature {t_ref:.3f} T.\n")

    # A
    print("A. Deux rectifications de R140 G")
    q277 = 100 * md2["oc"]
    print(f"    Convention. R140 comptait Q a la frequence 1. Les deux raies sont a la frequence {md2['oc']:.3f} : pour la meme resistance,")
    print(f"    le facteur de qualite de la bobine a la frequence des raies vaut {md2['oc']:.2f} fois plus. « Q = 100 » de R140 = Qb = {q277:.0f}.")
    split0 = 2 * md2["wt"] * md2["vc"] / (U - md2["vc"]); cells0 = 3 * 2 * pi / split0 * md2["vc"] * U / (U - md2["vc"])
    beat = 2 * pi / (md2["ohi"] - md2["olo"]); stay = 1 / md2["vlo"] - 1 / U
    cells3 = 3 * beat / stay
    print(f"    Longueur. R140 donnait aux deux raies la vitesse de la porteuse, {md2['vc']:.3f} : {cells0:.0f} cellules pour trois battements.")
    print(f"    La raie basse va a {md2['vlo']:.3f} et la raie haute a {md2['vhi']:.3f}. La raie basse suit le front de pres : au point de mesure elle")
    print(f"    ne dure que {stay:.3f} unite de temps par cellule parcourue. Trois battements de {beat:.1f} demandent {cells3:.0f} cellules, pas {cells0:.0f}.\n")
    check("A. rectifications : « Q = 100 » de R140 vaut Qb entre 270 et 285 a la frequence des raies ; trois battements demandent 140 a 152 cellules, contre 70 a 80 dans R140",
          270 < q277 < 285 and 140 < cells3 < 152 and 70 < cells0 < 80, f"Qb = {q277:.0f} ; {cells3:.0f} cellules contre {cells0:.0f}")

    # B
    print(f"B. Banc ideal : {N} cellules, pieces exactes, sans pertes, lecture a la cellule {TAP}")
    ideal = {x: one(x, T) for x in XS}
    sharp = one(2.0, T, wcut=1.5)
    two = {x: one(x, T, carrier=False) for x in XS}
    i2 = ideal[2.0]
    print(f"    x = 2 : sonde a la frequence {md2['om']:.3f}, raies a {md2['olo']:.3f} et {md2['ohi']:.3f}, porteuse a {md2['oc']:.3f} ; fenetre de lecture {i2['window']:.0f} unites de temps")
    print(f"    raie haute {i2['hi']:.2f} fois la sonde en tension, raie basse {i2['lo']:.2f} : le front leur donne son energie.")
    print(f"    onde parasite a la porteuse : {sharp['car']:.2f} fois la sonde si la sonde est coupee net, {i2['car']:.2f} si elle est coupee sur une periode.")
    print("    Mecanisme : la partie lente de la coupure est presque immobile vue du front ; l'horizon l'amplifie sans limite.")
    print("    rapport lu / reference : " + " ; ".join(f"x = {x} : {ideal[x]['ratio']/refr[x]:.3f}" for x in XS)
          + "   (sans compter la porteuse : " + " ; ".join(f"{two[x]['ratio']/refr[x]:.2f}" for x in XS) + ")")
    print("    -> sur un banc sans defaut la loi se lit a 10 % pres, a condition de compter la troisieme onde.\n")
    check("B. banc ideal : raie haute entre 1,8 et 2,3 fois la sonde, raie basse entre 0,8 et 1,2 (x = 2) ; onde a la porteuse au-dessus de 0,8 pour une coupure nette et sous 0,6 pour une coupure douce ; rapport lu a moins de 10 % de la reference pour x = 1,5, 2 et 2,5",
          1.8 < i2["hi"] < 2.3 and 0.8 < i2["lo"] < 1.2 and sharp["car"] > 0.8 and i2["car"] < 0.6 and all(abs(ideal[x]["ratio"] / refr[x] - 1) < 0.10 for x in XS),
          f"{i2['hi']:.2f} et {i2['lo']:.2f} ; porteuse {sharp['car']:.2f} et {i2['car']:.2f} ; " + " ; ".join(f"{ideal[x]['ratio']/refr[x]:.3f}" for x in XS))

    # C
    print("C. Pertes de vraies bobines (resistance comme la racine de la frequence, Qb compte a la frequence des raies)")
    lossy = {}
    for qb in (200.0, 100.0, 50.0):
        lossy[qb] = {x: one(x, T, qb=qb) for x in XS}
        lv = [ideal[x]["level"] / lossy[qb][x]["level"] for x in XS]
        print(f"    Qb = {qb:3.0f} : raies {sum(lv)/3:5.1f} fois plus faibles en energie ({sqrt(sum(lv)/3):.1f} fois en tension) ; rapport corrige des pertes / reference : "
              + " ; ".join(f"{lossy[qb][x]['ratio']/refr[x]:.3f}" for x in XS))
    flat = one(2.0, T, qb=100.0, skin=False)
    f_flat = ideal[2.0]["level"] / flat["level"]
    fac = {qb: sum(ideal[x]["level"] / lossy[qb][x]["level"] for x in XS) / 3 for qb in lossy}
    qb10 = sum(qb * log(f) for qb, f in fac.items()) / len(fac) / log(10)      # le facteur suit exp(constante / Qb)
    print(f"    Borne pessimiste, resistance independante de la frequence : Qb = 100 donne un facteur {f_flat:.0f}.")
    print(f"    Le facteur suit exp(constante / Qb) : il vaut 10 pour Qb = {qb10:.0f}.")
    print("    -> la sonde parcourt 150 cellules avant d'etre convertie : c'est elle qui paie. Le seuil de 10 ecrit avant R140 est tenu")
    print(f"       avec Qb = 200, depasse avec Qb = 100, et de tres loin avec Qb = 50 ; il demande Qb = {qb10:.0f}. Le rapport des raies, lui, se corrige.\n")
    check("C. pertes : facteur sous 6 pour Qb = 200, entre 14 et 24 pour Qb = 100, au-dessus de 100 pour Qb = 50 ; borne pessimiste au-dessus de 50 pour Qb = 100 ; seuil de 10 pour Qb entre 115 et 145 ; rapport corrige a moins de 10 % de la reference pour Qb = 200 et 100",
          fac[200.0] < 6 and 14 < fac[100.0] < 24 and fac[50.0] > 100 and f_flat > 50 and 115 < qb10 < 145
          and all(abs(lossy[qb][x]["ratio"] / refr[x] - 1) < 0.10 for qb in (200.0, 100.0) for x in XS),
          f"{fac[200.0]:.1f} ; {fac[100.0]:.1f} ; {fac[50.0]:.0f} ; borne {f_flat:.0f} ; seuil a Qb = {qb10:.0f}")

    # D
    print("D. Pieces a tolerance et instants de commande imprecis (Qb = 200, six bancs tires au hasard)")
    SEEDS = (1, 2, 3, 4, 5, 6)
    runs = {}
    for tol, jit, xs in ((0.01, 0.02, XS), (0.02, 0.05, [2.0]), (0.05, 0.10, [2.0]), (0.0, 0.10, [2.0])):
        for x in xs:
            runs[(tol, x)] = [one(x, T, qb=200.0, tol=tol, jit=jit, seed=s)["ratio"] / refr[x] for s in SEEDS]
    stat = {}
    for tol, jit in ((0.01, 0.02), (0.02, 0.05), (0.05, 0.10), (0.0, 0.10)):
        v = runs[(tol, 2.0)]; m = sum(v) / len(v); sd = sqrt(sum((a - m) ** 2 for a in v) / (len(v) - 1))
        stat[tol] = (m, sd, min(v), max(v))
        label = f"pieces a {100*tol:.0f} %" if tol else "pieces exactes"
        print(f"    {label}, instants a {100*jit:.0f} % d'une cellule, x = 2 : rapport lu / reference {m:.2f} +- {sd:.2f} (de {min(v):.2f} a {max(v):.2f})")
    print("    -> d'un banc a l'autre la lecture bouge de 15 % avec des pieces a 1 %, de 30 % a 2 %, et d'un facteur 2 a 5 %.")
    print("       Ce sont les pieces qui comptent, pas les instants de commande : un dixieme de cellule d'erreur ne coute que 5 %.")
    print("       Les deux raies sont a 84 et 96 % de la coupure : une cellule un peu lourde les arrete presque.\n")
    check("D. tolerances, x = 2 : dispersion sous 0,20 a 1 % ; entre 0,15 et 0,45 a 2 % ; au-dessus de 0,45 a 5 % ; sous 0,10 pour des pieces exactes commandees a 10 % d'une cellule pres",
          stat[0.01][1] < 0.20 and 0.15 < stat[0.02][1] < 0.45 and stat[0.05][1] > 0.45 and stat[0.0][1] < 0.10,
          " ; ".join(f"{100*tol:.0f} % : +- {stat[tol][1]:.2f}" for tol in (0.01, 0.02, 0.05)) + f" ; instants seuls : +- {stat[0.0][1]:.2f}")

    # E
    print("E. La temperature tiree de trois frequences de sonde (x = 1,5, 2 et 2,5)")
    t_ideal = temperature(XS, [ideal[x]["ratio"] for x in XS]) / t_ref
    t_loss = temperature(XS, [lossy[200.0][x]["ratio"] for x in XS]) / t_ref
    ts = [temperature(XS, [runs[(0.01, x)][i] * refr[x] for x in XS]) / t_ref for i in range(len(SEEDS))]
    tm = sum(ts) / len(ts); tsd = sqrt(sum((a - tm) ** 2 for a in ts) / (len(ts) - 1))
    print(f"    banc ideal : {t_ideal:.3f} fois celle de R140 ; avec Qb = 200 : {t_loss:.3f} ; avec des pieces a 1 % : {tm:.2f} +- {tsd:.2f} (de {min(ts):.2f} a {max(ts):.2f})")
    print("    -> avec des pieces a 1 %, un banc donne la temperature a 10 ou 15 % pres : les erreurs des trois frequences vont ensemble.\n")
    check("E. temperature : banc ideal et Qb = 200 a moins de 10 % de celle de R140 ; pieces a 1 % : dispersion entre 0,08 et 0,35",
          abs(t_ideal - 1) < 0.10 and abs(t_loss - 1) < 0.10 and 0.08 < tsd < 0.35, f"{t_ideal:.3f} ; {t_loss:.3f} ; +- {tsd:.2f}")

    print("Verdict : sur un banc sans defaut de 195 cellules, la loi se lit a 10 % pres. Avec de vraies pieces, deux choses la gardent hors")
    print(f"de portee d'un montage ordinaire. Les pertes : le seuil de 10 ecrit avant R140 demande des bobines de Qb = {qb10:.0f} a la frequence")
    print("des raies ; R140 l'annoncait tenu avec Q = 100, par une convention trop favorable. Les tolerances : il faut des pieces a 1 % pour")
    print("lire le rapport des raies a 15 % et la temperature a 10 ou 15 %. Voir les deux raies reste a la portee d'un banc.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
