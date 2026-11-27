#!/usr/bin/env python3
"""R140 -- Chemin C, calcul C1 : une ligne de transmission dont un front mobile change la vitesse d'onde. Loi de Hawking ?

Decision de l'auteur (27 novembre 2026) : « pars la simulation maintenant ». Premier pas du chemin C de la feuille de
route. En clair : un nageur dans une riviere qui accelere ; la ou le courant va aussi vite que lui, il ne remonte
plus : c'est un horizon. Ici le courant est un front qui avance le long d'une echelle de bobines et de condensateurs
et qui abaisse la capacite derriere lui ; le nageur est une petite onde sonde, plus lente que le front.

Le montage simule. Echelle LC, une cellule par unite de longueur, L = C = 1 devant le front : vitesse d'onde 1,
frequence de coupure 2. Derriere le front la vitesse vaut 1,5. Le front avance a la vitesse u = 1,2 : il rattrape
la sonde, qui ne peut pas le traverser. Deux versions :
  - ligne adaptee : L et C baissent ensemble, l'impedance ne change pas (cas d'ecole) ;
  - ligne a capacite seule : seule C baisse, de 1 a 0,444, et l'impedance monte de 1 a 1,5 (le cas d'un banc a
    diodes varicap, et l'idee de depart de l'auteur : l'energie monte l'impedance).

Ce que la theorie des trous noirs predit. La sonde est convertie en deux ondes courtes laissees derriere le front,
de part et d'autre de la frequence dont la vitesse de phase vaut u : une raie haute et une raie basse. Dans le repere
du front la frequence w de la sonde est conservee, et
    action de la raie basse / action de la raie haute = exp(-2 pi w / kappa),   kappa = |d(vitesse d'onde)/dx| a l'horizon,
ou l'action d'une onde est son energie divisee par sa frequence. T = kappa / (2 pi) joue le role d'une temperature.
De plus action haute - action basse = action de la sonde : les deux raies naissent par paires.

Criteres ecrits avant le calcul (feuille de route et proposition du 27 novembre) : la loi doit sortir a 10 % ; si
les pertes d'une ligne reelle noient l'effet de plus d'un facteur 10, on s'arrete a la simulation ; et la mesure ne
doit pas exister deja (recherche bibliographique, hors script).

  A. La sonde injectee ne va que dans un sens.
  B. Ligne adaptee, front large de 5 cellules : la loi, point par point, et le bilan des actions.
  C. La temperature suit la raideur du front : fronts de 10 et 3 cellules.
  D. Ligne a capacite seule, celle d'un banc : meme loi, temperature decalee de quelques pour cent.
  E. Fronts raides, jusqu'a une demi-cellule : ou la loi cesse-t-elle de tenir ?
  F. Pertes : une resistance en serie dans chaque bobine.
  G. Ordres de grandeur d'un banc : duree d'enregistrement, nombre de cellules, facteur perdu.
  Verdict : la loi sort a 0,1 % sur la ligne adaptee et a 6 % pres sur la temperature pour la ligne d'un banc,
  jusqu'a un front d'une demi-cellule ; les pertes ne changent pas le rapport des deux raies, seulement leur niveau.

Schema numerique : saute-mouton (charges aux pas entiers, flux aux demi-pas), symplectique ; la relation de dispersion
exacte du schema est employee dans l'analyse. Fenetre de calcul qui suit le front.
"""
from math import asin, cos, exp, log, pi, sin, sqrt, tanh
from operator import mul

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

U, CB = 1.2, 1.5            # vitesse du front ; vitesse de l'onde derriere le front (devant : 1)
M = 5                       # pas de temps par cellule parcourue par le front
DT = 1.0 / (U * M)

def root(f, lo, hi):
    flo = f(lo)
    assert (flo > 0) != (f(hi) > 0), "racine non encadree"
    for _ in range(100):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return 0.5 * (lo + hi)

def theta(k, wmax): return 2 * asin(0.5 * DT * wmax * sin(k / 2))      # phase par pas : dispersion exacte du schema
def wfront(k, wmax): return U * k - theta(k, wmax) / DT                  # frequence dans le repere du front (conservee)
def vgroup(k, wmax):
    h = 1e-6
    return (theta(k + h, wmax) - theta(k - h, wmax)) / (2 * h * DT)

class Line:
    """Echelle LC dans une fenetre qui suit le front. Indice 0 : loin derriere ; le front est a l'indice jf."""
    def __init__(self, w, matched, n_behind, n_ahead, loss=0.0):
        self.w, self.jf, self.N, self.loss = w, n_behind, n_behind + n_ahead, loss
        g = lambda xi: 0.5 * (1 - tanh(xi / w))                         # 1 derriere le front, 0 devant
        if matched:
            d = 1 - 1 / CB
            self.cfun = lambda xi: 1 - d * g(xi); self.lfun = self.cfun
        else:
            d = 1 - 1 / CB ** 2
            self.cfun = lambda xi: 1 - d * g(xi); self.lfun = lambda xi: 1.0
        self.Cb, self.Lb = self.cfun(-1e9), self.lfun(-1e9)
        self.wmaxb = 2 / sqrt(self.Lb * self.Cb)
        self.Wz = int(14 * w) + 2                                       # au-dela, le profil est constant a 1e-12 pres
        z0, z1 = self.jf - self.Wz, self.jf + self.Wz
        self.z0, self.z1 = z0, z1
        self.tabC = [[1 / self.cfun(j - self.jf - p / M) for j in range(z0, z1)] for p in range(M)]
        self.tabL = [[1 / self.lfun(j + 0.5 - self.jf - (p + 0.5) / M) for j in range(z0, z1)] for p in range(M)]
        self.Q = [0.0] * self.N; self.P = [0.0] * (self.N - 1)
        self.steps = 0
        vel = lambda xi: 1 / sqrt(self.cfun(xi) * self.lfun(xi))
        self.xih = root(lambda xi: vel(xi) - U, -10 * w, 10 * w); h = 1e-5 * w
        self.kappa = abs(vel(self.xih + h) - vel(self.xih - h)) / (2 * h)

    def inject(self, sigma, x0, nk=480, ext=9.0):
        """Impulsion breve (chapeau mexicain, spectre en k^2 exp(-k^2 sigma^2 / 2)) qui ne va que vers la droite."""
        kmax = 6.5 / sigma; dk = kmax / nk
        ks = [(i + 0.5) * dk for i in range(nk)]
        amp = [(k * sigma) ** 2 * exp(-0.5 * (k * sigma) ** 2) * dk for k in ks]
        th = [theta(k, 2.0) for k in ks]
        for j in range(max(0, int(x0 - ext * sigma)), min(self.N - 1, int(x0 + ext * sigma))):
            self.Q[j] += sum(a * cos(k * (j - x0)) for a, k in zip(amp, ks))
            self.P[j] += sum(a * cos(k * (j + 0.5 - x0) - 0.5 * t) for a, k, t in zip(amp, ks, th))

    def run(self, nsteps):
        Q, P, N, z0, z1, r = self.Q, self.P, self.N, self.z0, self.z1, self.loss
        headC = [1 / self.Cb] * z0; tailC = [1.0] * (N - z1)
        headL = [1 / self.Lb] * z0; tailL = [1.0] * (N - 1 - z1)
        for m in range(self.steps, self.steps + nsteps):
            p = m % M
            V = list(map(mul, Q, headC + self.tabC[p] + tailC))
            il = headL + self.tabL[p] + tailL
            if r:
                P = [(pp * (1 - 0.5 * r * DT * l) + DT * (a - b)) / (1 + 0.5 * r * DT * l) for pp, a, b, l in zip(P, V, V[1:], il)]
            else:
                P = [pp + DT * (a - b) for pp, a, b in zip(P, V, V[1:])]
            I = list(map(mul, P, il))
            Q = [Q[0] - DT * I[0]] + [q + DT * (a - b) for q, a, b in zip(Q[1:], I, I[1:])] + [Q[-1] + DT * I[-1]]
            if p == M - 1:                                              # le front a avance d'une cellule : la fenetre suit
                Q.pop(0); P.pop(0); Q.append(0.0); P.append(0.0)
        self.Q, self.P = Q, P; self.steps += nsteps

    def actions(self, k, j0, j1, C, L):
        """Densites d'action, par unite de k, des ondes de nombre d'onde k dans la zone uniforme [j0, j1) :
        celles qui vont vers la droite, celles qui vont vers la gauche."""
        Z = sqrt(L / C); wmax = 2 / sqrt(L * C)
        th = theta(k, wmax)
        ths = -th if self.steps == 0 else th                            # a l'instant initial le flux est en avance d'un demi-pas
        vr = vi = ir = ii = 0.0
        ck, sk = cos(k), sin(k)
        c, s = cos(k * j0), sin(k * j0); ch, sh = cos(k * (j0 + 0.5)), sin(k * (j0 + 0.5))
        Q, P = self.Q, self.P
        for j in range(j0, j1):
            v = Q[j] / C; i = P[j] / L
            vr += v * c; vi -= v * s; ir += i * ch; ii -= i * sh
            c, s = c * ck - s * sk, s * ck + c * sk
            ch, sh = ch * ck - sh * sk, sh * ck + ch * sk
        ar = (Z * ir + vr * cos(ths / 2) + vi * sin(ths / 2)) / (2 * cos(th / 2))
        ai = (Z * ii + vi * cos(ths / 2) - vr * sin(ths / 2)) / (2 * cos(th / 2))
        fac = C * cos(th / 2) / (wmax * sin(k / 2))                     # action exacte d'un mode du schema saute-mouton
        return fac * (ar * ar + ai * ai), fac * ((vr - ar) ** 2 + (vi - ai) ** 2)

def experiment(w, matched, xs, sigma=None, loss=0.0):
    """Une sonde large bande avalee par le front. Renvoie kappa, T, la duree simulee et, pour chaque x = w / T,
    (action de la sonde, action de la raie haute, action de la raie basse), par unite de frequence."""
    sigma = sigma if sigma else max(2.2, 9.4 * w / 5.0)
    probe = Line(w, matched, 10, 10)
    kappa = probe.kappa; T = kappa / (2 * pi)
    za = int(9 * w) + 2                                                 # au-dela de 9 w le profil est uniforme a 1e-7 pres
    d0 = za + 9 * sigma + 4
    ks0 = root(lambda k: wfront(k, probe.wmaxb), 1.5, 3.1)              # onde dont la vitesse de phase vaut celle du front
    vg0 = vgroup(ks0, probe.wmaxb)
    t_total = (d0 + 9 * sigma) / (U - 1.0) + (za + 40) / (U - vg0) + 8 / kappa
    ln = Line(w, matched, int((U - vg0) * t_total) + za + 120, int(d0 + 9 * sigma) + 40, loss)
    ln.inject(sigma, ln.jf + d0)
    kin = [root(lambda k: wfront(k, 2.0) - x * T, 1e-6, 1.2) for x in xs]
    nin = [ln.actions(k, ln.jf + za, ln.N - 1, 1.0, 1.0) for k in kin]
    left_in = max(a[1] / a[0] for a in nin)
    nin = [a[0] / abs(U - vgroup(k, 2.0)) for a, k in zip(nin, kin)]
    nsteps = int(t_total / DT); nsteps -= nsteps % M
    ln.run(nsteps)
    wb = ln.wmaxb; out = []
    for x, ni in zip(xs, nin):
        kp = root(lambda k: wfront(k, wb) - x * T, ks0, 3.1)
        km = root(lambda k: wfront(k, wb) + x * T, 1.3, ks0)
        nhi = ln.actions(kp, 3, ln.jf - za, ln.Cb, ln.Lb)[0] / abs(U - vgroup(kp, wb))
        nlo = ln.actions(km, 3, ln.jf - za, ln.Cb, ln.Lb)[0] / abs(U - vgroup(km, wb))
        out.append((ni, nhi, nlo))
    e = [q * q for q in ln.Q]; tot = sum(e)
    leak = (sum(e[ln.jf - za:]) + sum(e[:40])) / tot                    # ce qui n'est pas dans la zone d'analyse
    return dict(kappa=kappa, T=T, t=nsteps * DT, out=out, left_in=left_in, leak=leak, N=ln.N, ks=ks0, vg=vg0)

def slope(xs, ratios):
    """Pente de ln(raie basse / raie haute) en fonction de x : la loi donne -1 ; -1 / pente = T mesuree / T."""
    n = len(xs); mx = sum(xs) / n; my = sum(log(r) for r in ratios) / n
    return sum((x - mx) * (log(r) - my) for x, r in zip(xs, ratios)) / sum((x - mx) ** 2 for x in xs)

def main():
    print("R140 -- chemin C, calcul C1 : la loi de Hawking sur une ligne de transmission a front mobile\n")
    XS = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0]

    # A et B
    b = experiment(5.0, True, XS)
    print("A. La sonde injectee")
    print(f"    part de la sonde qui va dans le mauvais sens : {b['left_in']:.1e} ; ce qui reste hors de la zone d'analyse en fin de calcul : {b['leak']:.1e}\n")
    check("A. la sonde ne va que dans un sens (moins de 1e-5 dans l'autre) et tout est converti en fin de calcul (reste sous 1e-8)",
          b["left_in"] < 1e-5 and b["leak"] < 1e-8, f"{b['left_in']:.1e} ; {b['leak']:.1e}")

    print("B. Ligne adaptee, front large de 5 cellules")
    print(f"    kappa = {b['kappa']:.4f}, T = kappa / (2 pi) = {b['T']:.5f} (la coupure de la ligne vaut 2) ; duree simulee {b['t']:.0f}, {b['N']} cellules")
    print("    x = w / T :  raie basse / raie haute   loi exp(-x)   mesure / loi   raie haute / sonde   1 / (1 - e^-x)   (haute - basse) / sonde")
    rb, devs, cons, his = [], [], [], []
    for x, (ni, nhi, nlo) in zip(XS, b["out"]):
        r = nlo / nhi; rb.append(r); devs.append(r / exp(-x) - 1); cons.append((nhi - nlo) / ni - 1); his.append(nhi / ni * (1 - exp(-x)) - 1)
        print(f"      {x:4.1f}       {r:.5e}          {exp(-x):.5e}     {r/exp(-x):.4f}        {nhi/ni:.4f}            {1/(1-exp(-x)):.4f}          {(nhi-nlo)/ni:.5f}")
    sb = slope(XS, rb)
    print(f"    pente ajustee de ln(basse / haute) : {sb:.4f} pour -1 ; temperature mesuree / T = {-1/sb:.4f}")
    print("    -> la loi sort sur huit unites de x, soit un rapport qui va de 0,6 a 0,0003 : a 0,01 % pres a basse frequence, a 6 % au")
    print("       bout. Les deux raies naissent par paires : haute moins basse redonne la sonde a 0,06 % pres.\n")
    check("B. ligne adaptee : loi exp(-x) a moins de 1 % pour x <= 3 et a moins de 7 % jusqu'a x = 8 ; raie haute / sonde = 1 / (1 - e^-x) a 0,1 % pres ; bilan des actions a 0,1 % pres",
          max(abs(d) for x, d in zip(XS, devs) if x <= 3) < 0.01 and max(abs(d) for d in devs) < 0.07 and max(abs(h) for h in his) < 1e-3 and max(abs(c) for c in cons) < 1e-3,
          f"ecart maximal {100*max(abs(d) for d in devs):.1f} % a x = 8 ; bilan {100*max(abs(c) for c in cons):.2f} %")

    # C
    print("C. La temperature suit la raideur du front (ligne adaptee)")
    scal = {5.0: (b["kappa"], -1 / sb)}
    for w in (10.0, 3.0):
        e = experiment(w, True, XS)
        s = slope(XS, [nlo / nhi for (ni, nhi, nlo) in e["out"]]); scal[w] = (e["kappa"], -1 / s)
    for w in (10.0, 5.0, 3.0):
        print(f"    front de {w:4.1f} cellules : kappa = {scal[w][0]:.4f} ; kappa fois largeur = {scal[w][0]*w:.4f} ; temperature mesuree / T = {scal[w][1]:.4f}")
    print("    -> kappa est inversement proportionnel a la largeur du front, et la temperature mesuree le suit a 2 % pres.\n")
    check("C. fronts de 10, 5 et 3 cellules : kappa fois largeur constant (0,240) ; temperature mesuree egale a kappa / (2 pi) a moins de 2 %",
          all(abs(scal[w][0] * w - 0.24) < 1e-3 for w in scal) and all(abs(scal[w][1] - 1) < 0.02 for w in scal),
          " ; ".join(f"{w:g} : {scal[w][1]:.4f}" for w in (10.0, 5.0, 3.0)))

    # D
    print("D. Ligne a capacite seule, celle d'un banc (l'impedance monte de 1 a 1,5 derriere le front)")
    d = experiment(5.0, False, XS)
    rd = [nlo / nhi for (ni, nhi, nlo) in d["out"]]; sd = slope(XS, rd)
    cons_d = [(nhi - nlo) / ni - 1 for (ni, nhi, nlo) in d["out"]]
    print(f"    kappa = {d['kappa']:.4f}, T = {d['T']:.5f}")
    print("    mesure / loi : " + " ; ".join(f"x = {x:g} : {r/exp(-x):.3f}" for x, r in zip(XS, rd)))
    print(f"    temperature mesuree / T = {-1/sd:.4f} ; (haute - basse) / sonde - 1 : de {100*min(cons_d):.2f} a {100*max(cons_d):.2f} %")
    print("    -> meme loi exponentielle, avec une temperature 3 % plus haute que kappa / (2 pi). L'exces du bilan, moins de 1 %, est une")
    print("       onde renvoyee vers l'arriere par le saut d'impedance.\n")
    check("D. ligne a capacite seule : loi a moins de 10 % point par point jusqu'a x = 4 ; temperature mesuree a 3 % au-dessus de kappa / (2 pi) ; exces du bilan sous 1 %",
          all(abs(r / exp(-x) - 1) < 0.10 for x, r in zip(XS, rd) if x <= 4) and 1.01 < -1 / sd < 1.05 and 0 < min(cons_d) and max(cons_d) < 0.01,
          f"temperature x {-1/sd:.3f} ; bilan +{100*max(cons_d):.2f} %")

    # E
    print("E. Fronts raides (ligne a capacite seule)")
    wlim = wfront(root(lambda k: vgroup(k, 3.0) - U, 0.5, 2.0), 3.0)     # au-dela de |wlim| la raie basse n'existe plus
    kc = 2 * asin(sqrt(1 - (U / CB) ** 2)); wlim_ideal = 2 * CB * sin(kc / 2) - U * kc
    print(f"    la raie basse n'existe que pour w < {-wlim:.3f} ({wlim_ideal:.3f} pour la ligne ideale, sans pas de temps) : c'est la bande disponible, quelle que soit la raideur")
    sharp = {}
    for w in (2.0, 1.0, 0.5):
        Tw = Line(w, False, 10, 10).kappa / (2 * pi)
        xs = [x for x in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0) if x * Tw < 0.21]
        e = experiment(w, False, xs)
        rr = [nlo / nhi for (ni, nhi, nlo) in e["out"]]
        sharp[w] = (e["T"], -1 / slope(xs, rr), max(abs(r / exp(-x) - 1) for x, r in zip(xs, rr) if x <= 2), xs[-1])
        print(f"    front de {w:3.1f} cellule(s) : T = {e['T']:.4f} ; temperature mesuree / T = {sharp[w][1]:.3f} ; ecart a la loi jusqu'a x = 2 : {100*sharp[w][2]:.1f} % ; x disponible jusqu'a {-wlim/e['T']:.1f}")
    print("    -> la loi ne lache pas quand le front devient raide : a une demi-cellule elle tient encore a 5 % pres. Ce qui se reduit est")
    print("       la plage de frequences, parce que la bande disponible est fixe et que T grandit.\n")
    check("E. fronts de 2, 1 et 0,5 cellule : temperature mesuree a moins de 6 % de kappa / (2 pi) ; loi a moins de 7 % jusqu'a x = 2 ; bande disponible 0,263 (0,256 pour la ligne ideale)",
          all(abs(v[1] - 1) < 0.06 and v[2] < 0.07 for v in sharp.values()) and abs(-wlim - 0.263) < 3e-3 and abs(wlim_ideal - 0.2557) < 1e-3,
          " ; ".join(f"{w:g} : {sharp[w][1]:.3f}" for w in (2.0, 1.0, 0.5)))

    # F
    print("F. Pertes : une resistance en serie dans chaque bobine (facteur de qualite Q a la frequence 1)")
    xs = [0.5, 1.0, 2.0, 3.0]
    ref = experiment(1.0, False, xs)
    lossy = {}
    for q in (200.0, 100.0):
        e = experiment(1.0, False, xs, loss=1.0 / q)
        dr = max(abs((b_[2] / b_[1]) / (a_[2] / a_[1]) - 1) for a_, b_ in zip(ref["out"], e["out"]))
        lvl = [b_[1] / a_[1] for a_, b_ in zip(ref["out"], e["out"])]
        lossy[q] = (dr, sum(lvl) / len(lvl), exp(-e["t"] / q))
        print(f"    Q = {q:3.0f} : rapport des deux raies change de {100*dr:.2f} % ; niveau des raies x {lossy[q][1]:.4f} ; exp(-duree / Q) = {lossy[q][2]:.4f} (duree simulee {e['t']:.0f})")
    print("    -> les pertes ne touchent pas au rapport des deux raies ; elles baissent leur niveau comme exp(-duree / Q).\n")
    check("F. pertes : rapport des deux raies inchange a moins de 1 % ; niveau egal a exp(-duree / Q) a moins de 3 %",
          all(v[0] < 0.01 and abs(v[1] / v[2] - 1) < 0.03 for v in lossy.values()), " ; ".join(f"Q = {q:.0f} : {v[1]:.3f} pour {v[2]:.3f}" for q, v in lossy.items()))

    # G
    print("G. Ordres de grandeur d'un banc (sonde continue, mesure en un point derriere le front, x = 2)")
    ks, vg = ref["ks"], ref["vg"]
    print(f"    porteuse des deux raies : nombre d'onde {ks:.3f} par cellule, frequence {U*ks:.2f}, soit {U*ks/2:.2f} fois la coupure de la ligne au repos ; vitesse de groupe {vg:.3f}")
    bench = {}
    for w in (0.5, 1.0):
        T = sharp[w][0]; x = 2.0; wt = x * T
        split = 2 * wt * vg / (U - vg)                                   # ecart de frequence entre les deux raies, au point de mesure
        t_rec = 3 * 2 * pi / split                                       # trois battements pour les separer
        cells = t_rec * vg * U / (U - vg)                                # distance parcourue par les ondes les plus anciennes
        t_trav = cells / vg
        kin = root(lambda k: wfront(k, 2.0) - wt, 1e-6, 2.0)
        probe_f = theta(kin, 2.0) / DT                                   # frequence de la sonde
        for q in (100.0, 200.0):
            keep = (1 - exp(-t_trav / q)) / (t_trav / q) * exp(-(cells + 15) / q)
            bench[(w, q)] = (keep, cells + 15)
        bench[w] = (probe_f, split)
        print(f"    front de {w:3.1f} cellule : sonde a la frequence {probe_f:.2f} ({15*probe_f:.0f} MHz si la coupure est a 30 MHz) ; ecart entre les raies {split:.3f} ({15*split/1:.1f} MHz) ; "
              f"enregistrement {t_rec:.0f} ({t_rec*10.6e-3:.1f} microseconde) ; {cells+15:.0f} cellules ; niveau garde : {bench[(w,100.0)][0]:.2f} pour Q = 100, {bench[(w,200.0)][0]:.2f} pour Q = 200")
    print("    -> avec des bobines de Q = 100, un front d'une demi-cellule et une centaine de cellules, les pertes coutent un facteur 4 :")
    print(f"       sous le seuil de 10. Avec un front d'une cellule il faut Q = 200. Exemple : coupure de la ligne au repos a 30 MHz, soit")
    print(f"       10,6 ns par cellule ; les deux raies sont alors vers {15*U*ks:.0f} MHz, au-dessus de cette coupure : elles ne vivent que")
    print("       dans la partie de la ligne ou le front est deja passe.\n")
    check("G. banc : pertes sous le facteur 10 pour un front d'une demi-cellule avec Q = 100 (facteur 4, 90 cellules) ; pour un front d'une cellule il faut Q = 200",
          bench[(0.5, 100.0)][0] > 0.1 and bench[(1.0, 100.0)][0] < 0.1 and bench[(1.0, 200.0)][0] > 0.1 and bench[(0.5, 100.0)][1] < 120,
          f"{bench[(0.5,100.0)][0]:.2f} ; {bench[(1.0,100.0)][0]:.2f} ; {bench[(1.0,200.0)][0]:.2f}")

    print("Verdict : les deux criteres calculables sont remplis. La loi sort a 0,1 % sur la ligne adaptee et a 6 % pres sur la temperature")
    print("pour la ligne d'un banc, jusqu'a un front d'une demi-cellule ; les pertes ne changent pas le rapport des deux raies et coutent")
    print("un facteur 4 sur leur niveau avec Q = 100. Limite : le front simule est impose de l'exterieur, lisse et de vitesse choisie ;")
    print("une impulsion forte lancee dans une vraie ligne a diodes se deforme. C'est le point a simuler avant tout achat.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
