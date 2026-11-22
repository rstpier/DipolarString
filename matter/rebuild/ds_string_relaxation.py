#!/usr/bin/env python3
"""R138 -- La detente du fil appliquee a la corde tournante : explique-t-elle l'ecart entre spectres et reseau ?

Choix fait par delegation le 22 novembre 2026, annonce en R136 et dans la feuille de route : appliquer a la corde
tournante l'aplatissement du potentiel publie par Badalian et Bakker (arXiv:1901.10280, eq. 51 a 54), avec la tension
du reseau. Mecanisme : des paires de quarks legers apparaissent le long d'un fil long, qui tire alors moins fort.

Critere ecrit avant le calcul : l'ecart entre la tension des spectres (423,7 +- 2,5 MeV) et celle du reseau
(445 +- 6,7 MeV) vaut 21 MeV, soit 3,0 ecarts types (R137). La piste tient si la detente publiee le ramene sous
2 ecarts types ; sinon elle est exclue pour ce modele, la corde tournante a extremites massives de R136.

  A. La fonction publiee, relue : sigma_f(r) = sigma (1 - gamma f(r)), f = e^x / (B + e^x), x = sqrt(sigma) (r - R0),
     gamma = 0,40, B = 20, R0 = 6,0 GeV^-1, sigma = 0,182 GeV^2. Ce qu'elle vaut aux longueurs des mesons.
  B. Deux versions de la corde tournante dont la tension depend de la longueur. Naive : la force aux extremites est
     sigma_f(L). Coherente : l'energie d'un fil au repos etant sigma_f(L) L, la force est sa derivee. Controle : une
     rotation stationnaire verifie dE/dJ = vitesse angulaire.
  C. Controles du solveur : sans detente, il redonne le calcul direct de R136 et l'ajustement de R137.
  D. Version coherente : le moment cinetique J(L) passe par un maximum. Au-dela il n'y a plus de rotation stable.
  E. Ajustement des neuf familles au meilleur jeu publie, dans les deux versions : tension demandee par les spectres,
     degradation de l'ajustement, ecart au reseau.
  F. La fourchette publiee : gamma de 0,35 a 0,50, B de 15 a 25, R0 de 5 a 8 GeV^-1.
  G. Les spectres demandent-ils une detente ? Degradation en fonction de gamma, tension libre.
  H. Le critere.
  Verdict : la detente publiee remonte la tension demandee par les spectres de 424 a 430 MeV, mais l'ajustement se
  degrade d'autant : l'ecart complet au reseau reste a 2,9 ecarts types (3,0 sans detente). La piste ne tient pas.

Mesure d'ecart : celle de l'article de 2014, somme des ((M^2 - E^2) / M^2)^2 ; erreur commune s ajustee pour que le
chi2 par degre de liberte vaille 1 sans detente (28 degres de liberte), comme en R136 et R137.
  degradation  = (chi2 minimal du modele - chi2 de reference) / s^2
  ecart complet = racine du minimum sur la tension de (chi2 du modele - chi2 de reference) / s^2 + ((tension - 445) / 6,7)^2
  ecart de tension seul = racine de (ecart complet^2 - degradation)
"""
from bisect import bisect_right
from math import asin, exp, log, pi, sqrt
import ds_b1_reproduction as b1
import ds_pdg2026_update as up

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

HBARC = b1.HBARC
RATE = sqrt(0.182)                      # GeV : racine de la tension de l'article, dans sa fonction f(r)
BEST = (0.40, 20.0, 6.0)                # gamma, B, R0 (GeV^-1) : eq. 54 de l'article
RS_LAT, DRS_LAT = 445.0, 6.7            # MeV : Bulava et al. 2024, erreurs ajoutees en quadrature (R134)
TRAJ = up.TRAJ_2026
PAIRS = sorted({tuple(sorted((q1, q2))) for _, q1, q2 in TRAJ.values()})

def phi(L, par):
    """sigma_f(L) / sigma et sa derivee par rapport a L (GeV)."""
    if par is None: return 1.0, 0.0
    gam, B, R0 = par
    x = RATE * (L - R0)
    if x > 80: return 1 - gam, 0.0
    e = exp(x); f = e / (B + e)
    return 1 - gam * f, -gam * RATE * f * (1 - f)

def solve_L(L, sig, m1, m2, par, mode, guess=None):
    """Corde droite de longueur totale L en rotation stationnaire, extremites de masses m1 et m2.
    Tension uniforme T = sig phi(L). Equilibre de chaque extremite : m w beta gamma = T / gamma + X, avec X = 0
    (version naive, mode "a") ou X = T'(L) fois l'integrale de sqrt(1 - w^2 r^2) le long du fil (version coherente,
    mode "b" : c'est la derivee du lagrangien -T(L) integrale(sqrt(1 - w^2 r^2)) par rapport a la position du bout).
    Inconnues : les vitesses des deux bouts. Renvoie E, J, w, beta1, beta2, ou None s'il n'y a pas de solution."""
    f, df = phi(L, par)
    T = sig * f; Tp = sig * df if mode == "b" else 0.0
    if T + Tp * L <= 0: return None
    if guess is None:
        w = sqrt(T * (1 / m1 + 1 / m2) / L); ba, bb = min(T / (m1 * w), 0.9), min(T / (m2 * w), 0.9)
    else:
        ba, bb = guess
    done = False
    for _ in range(100):
        sa, sb = sqrt(1 - ba * ba), sqrt(1 - bb * bb)
        w = (ba + bb) / L
        G = ba * sa + asin(ba) + bb * sb + asin(bb)                 # 2 w fois l'integrale ; dG/dbeta = 2 s
        fa = m1 * w * ba / sa - T * sa - Tp * G / (2 * w)
        fb = m2 * w * bb / sb - T * sb - Tp * G / (2 * w)
        xa, xb = Tp * (sa / w - G / (2 * w * w * L)), Tp * (sb / w - G / (2 * w * w * L))
        a11 = m1 * ba / (L * sa) + m1 * w / sa ** 3 + T * ba / sa - xa
        a12 = m1 * ba / (L * sa) - xb
        a21 = m2 * bb / (L * sb) - xa
        a22 = m2 * bb / (L * sb) + m2 * w / sb ** 3 + T * bb / sb - xb
        det = a11 * a22 - a12 * a21
        if det == 0: return None
        da, db = (-fa * a22 + fb * a12) / det, (-fb * a11 + fa * a21) / det
        lam = 1.0
        while not (0 < ba + lam * da < 1 and 0 < bb + lam * db < 1):
            lam *= 0.5
            if lam < 1e-12: return None
        ba += lam * da; bb += lam * db
        if abs(da) + abs(db) < 1e-14: done = True; break
    if not done: return None
    w = (ba + bb) / L; E = J = 0.0
    for m, b in ((m1, ba), (m2, bb)):
        s = sqrt(1 - b * b)
        E += m / s + T * asin(b) / w
        J += m * b * b / (s * w) + 0.5 * T * (asin(b) - b * s) / (w * w)
    return E, J, w, ba, bb

class Branch:
    """Table E(J) le long de la branche stable : L croissant jusqu'a J > jmax, ou jusqu'au maximum de J(L)."""
    def __init__(self, sig, m1, m2, par, mode, jmax=12.5, ratio=1.03, lmin=0.1, lmax=18.0):
        self.m0 = m1 + m2
        self.t, self.J, self.E = [], [], []
        L = lmin; guess = None; self.end = "ouverte"
        while True:
            st = solve_L(L, sig, m1, m2, par, mode, guess)
            if st is None: self.end = "sans solution"; break
            if self.J and st[1] <= self.J[-1]:
                self.end = "maximum de J"
                g = guess
                lo = exp(self.t[-2]) if len(self.t) > 1 else exp(self.t[-1])
                top = b1.golden(lambda x: -solve_L(x, sig, m1, m2, par, mode, g)[1], lo, L, 60)
                sx = solve_L(top, sig, m1, m2, par, mode, g)
                if top > exp(self.t[-1]) * (1 + 1e-9) and sx[1] > self.J[-1]:
                    self.t.append(log(top)); self.J.append(sx[1]); self.E.append(sx[0])
                elif len(self.J) > 1 and sx[1] > self.J[-2]:
                    self.t[-1], self.J[-1], self.E[-1] = log(top), sx[1], sx[0]
                break
            self.t.append(log(L)); self.J.append(st[1]); self.E.append(st[0]); guess = (st[3], st[4])
            if st[1] > jmax: break
            L *= ratio
            if L > lmax: self.end = "longueur maximale"; break
        self.jend = self.J[-1]; self.lend = exp(self.t[-1]); self._c = {}
    def _coef(self, lo):
        """Differences divisees de J(t) et E(t) sur les quatre points lo .. lo + 3 (t = ln L)."""
        if lo not in self._c:
            ts = self.t[lo:lo + 4]; out = []
            for ys in (self.J[lo:lo + 4], self.E[lo:lo + 4]):
                d = list(ys)
                for k in range(1, 4):
                    for i in range(3, k - 1, -1): d[i] = (d[i] - d[i - 1]) / (ts[i] - ts[i - k])
                out.append(d)
            self._c[lo] = (ts, out[0], out[1])
        return self._c[lo]
    def at(self, Jc):
        """(E, L) au moment cinetique classique Jc ; None au-dela de la branche. Interpolation cubique en ln L."""
        J = self.J; n = len(J)
        if Jc > self.jend: return None
        if Jc <= J[0]:
            p = log((self.E[1] - self.m0) / (self.E[0] - self.m0)) / log(J[1] / J[0])
            q = (self.t[1] - self.t[0]) / log(J[1] / J[0])
            return self.m0 + (self.E[0] - self.m0) * (Jc / J[0]) ** p, exp(self.t[0] + q * log(Jc / J[0]))
        i = min(bisect_right(J, Jc) - 1, n - 2)
        ts, dj, de = self._coef(max(0, min(i - 1, n - 4)))
        t0, t1, t2 = ts[0], ts[1], ts[2]
        a, b = self.t[i], self.t[i + 1]
        for _ in range(38):
            mid = 0.5 * (a + b)
            if dj[0] + (mid - t0) * (dj[1] + (mid - t1) * (dj[2] + (mid - t2) * dj[3])) > Jc: b = mid
            else: a = mid
        t = 0.5 * (a + b)
        return de[0] + (t - t0) * (de[1] + (t - t1) * (de[2] + (t - t2) * de[3])), exp(t)

BIG = 10.0
def best_a(states, br):
    """Ordonnee optimale d'une famille ; elle ne peut pas descendre sous jtop - jend (etat hors de la branche)."""
    jmin = min(J for J, _ in states); jtop = max(J for J, _ in states)
    lo = max(jmin - 6.0, jtop - br.jend + 1e-9); hi = jmin - 1e-6
    if lo >= hi: return None, BIG, True
    f = lambda a: sum(((M * M - br.at(J - a)[0] ** 2) / (M * M)) ** 2 for J, M in states)
    a = b1.golden(f, lo, hi, 36)
    return a, f(a), (a - lo) < 1e-4 and lo > jmin - 6.0

CACHE = {}
def total(rsig, par, mode):
    """chi2 total des neuf familles a la tension sqrt(sigma0) = rsig (MeV) ; ordonnees optimales."""
    key = (round(rsig, 9), par, mode)
    if key not in CACHE:
        sig = (rsig * 1e-3) ** 2
        brs = {p: Branch(sig, b1.UNIV_MASS[p[0]], b1.UNIV_MASS[p[1]], par, mode) for p in PAIRS}
        tot = 0.0; a = {}; bound = []
        for name, (st, q1, q2) in TRAJ.items():
            a[name], s, atb = best_a(st, brs[tuple(sorted((q1, q2)))]); tot += s
            if atb: bound.append(name)
        CACHE[key] = (tot, a, bound, brs)
    return CACHE[key]

def study(par, mode, ref_tot, s2, n=16):
    """Minimum du modele, degradation, ecart complet et ecart de tension seul."""
    ff = lambda x: total(x, par, mode)[0]
    fq = lambda x: (ff(x) - ref_tot) / s2 + ((x - RS_LAT) / DRS_LAT) ** 2
    xm = b1.golden(ff, 410.0, 470.0, n); xq = b1.golden(fq, 410.0, 470.0, n)
    deg = (ff(xm) - ref_tot) / s2; q = fq(xq)
    return dict(xm=xm, deg=deg, ratio=ff(xm) / ref_tot, q=q, full=sqrt(max(q, 0)), tens=sqrt(max(q - deg, 0)),
                lat=sqrt(max((ff(RS_LAT) - ref_tot) / s2, 0)), bound=total(xm, par, mode)[2], jend=total(xm, par, mode)[3][("u", "u")].jend)

def main():
    print("R138 -- la detente du fil appliquee a la corde tournante\n")
    fm = lambda x: x / HBARC

    # A
    print("A. La fonction publiee, relue a la source (Badalian et Bakker, eq. 51 a 54)")
    gam, B, R0 = BEST
    rows = []
    for r in (0.45, 1.0, R0 * HBARC, 1.33, 1.6, 2.0, 2.5, 3.0, 4.0):
        f, df = phi(fm(r), BEST); rows.append((r, 1 - f, f + df * fm(r)))
        print(f"    r = {r:4.2f} fm : tension abaissee de {100*(1-f):5.2f} % ; force dV/dr = {f + df*fm(r):.3f} sigma")
    rmin = b1.golden(lambda r: phi(fm(r), BEST)[0] + phi(fm(r), BEST)[1] * fm(r), 2.0, 4.0, 50)
    fmin = phi(fm(rmin), BEST)[0] + phi(fm(rmin), BEST)[1] * fm(rmin)
    half = R0 * HBARC + log(B) / RATE * HBARC
    print(f"    la detente commence a R0 = {R0*HBARC:.2f} fm, ou elle vaut gamma / (B + 1) = {100*gam/(B+1):.2f} % ; elle atteint la moitie de 40 % a {half:.2f} fm")
    print(f"    la force passe par un minimum de {fmin:.3f} sigma a {rmin:.2f} fm, puis remonte vers {1-gam:.2f} sigma")
    print("    -> « 40 a 45 % » est la valeur atteinte au-dela de 3 fm. Aux longueurs des mesons ajustes (0,45 a 1,7 fm),")
    print("       la tension ne baisse que de 0,4 a 5 %. R136 G ecrivait que le mecanisme « a la bonne taille » : a verifier.\n")
    at_R0 = 1 - phi(R0, BEST)[0]
    check("A. fonction publiee : detente de gamma / (B + 1) = 1,90 % a R0 = 1,18 fm, moins de 5 % a 1,6 fm, 29 % a 3 fm ; minimum de la force 0,18 sigma a 2,9 fm",
          abs(at_R0 - gam / (B + 1)) < 1e-12 and rows[4][1] < 0.05 and 0.27 < rows[7][1] < 0.30 and 0.15 < fmin < 0.22 and 2.7 < rmin < 3.2,
          f"{100*at_R0:.2f} % a R0 ; {100*rows[4][1]:.1f} % a 1,6 fm ; force minimale {fmin:.3f} sigma a {rmin:.2f} fm")

    # B
    print("B. Deux versions de la corde tournante a tension dependant de la longueur")
    sig_lat = (RS_LAT * 1e-3) ** 2
    dev = {"a": [], "b": []}
    for mode in ("a", "b"):
        for r in (0.5, 1.0, 1.5, 1.9):
            L = fm(r); h = 1e-4
            s0 = solve_L(L, sig_lat, 0.06, 0.06, BEST, mode)
            s1 = solve_L(L * (1 - h), sig_lat, 0.06, 0.06, BEST, mode); s2_ = solve_L(L * (1 + h), sig_lat, 0.06, 0.06, BEST, mode)
            dev[mode].append((s2_[0] - s1[0]) / (s2_[1] - s1[1]) / s0[2] - 1)
        print(f"    version {'naive    ' if mode == 'a' else 'coherente'} : (dE/dJ) / w - 1 a 0,5, 1,0, 1,5 et 1,9 fm = " + ", ".join(f"{d:+.1e}" for d in dev[mode]))
    stat = {}
    for mode in ("a", "b"):
        L = fm(1.5); m = 5000.0
        s = solve_L(L, sig_lat, m, m, BEST, mode)
        f, df = phi(L, BEST)
        stat[mode] = m * s[2] ** 2 * (L / 2) / (sig_lat * (f + df * L))
    print(f"    limite statique (extremites tres lourdes, 1,5 fm) : force centripete / (dV/dL) = {stat['a']:.4f} (naive), {stat['b']:.4f} (coherente)")
    print("    -> la version naive garde la force sigma_f(L) alors que l'energie du fil est sigma_f(L) L : elle tire trop fort et viole")
    print("       dE/dJ = w. La version coherente prend la force dV/dL ; c'est elle qui traduit le potentiel publie.\n")
    check("B. version coherente : dE/dJ = w a mieux que 1e-5 et force statique egale a dV/dL ; version naive : ecart de plus de 1 % et force trop grande de 12 %",
          max(abs(d) for d in dev["b"]) < 1e-5 and max(abs(d) for d in dev["a"]) > 0.01 and abs(stat["b"] - 1) < 1e-3 and 1.08 < stat["a"] < 1.20,
          f"coherente {max(abs(d) for d in dev['b']):.0e} ; naive {max(abs(d) for d in dev['a']):.1e} ; statique {stat['a']:.3f} et {stat['b']:.4f}")

    # C
    print("C. Controles du solveur, sans detente")
    worst = 0.0
    for q1, q2 in PAIRS:
        br = Branch(1 / (2 * pi * 0.8864), b1.UNIV_MASS[q1], b1.UNIV_MASS[q2], None, "a")
        for j in (0.05, 0.3, 1.3, 2.7, 4.4, 7.9, 11.0):
            worst = max(worst, abs(br.at(j)[0] / b1.mass_of(j, 0.8864, b1.UNIV_MASS[q1], b1.UNIV_MASS[q2]) - 1))
    print(f"    table en longueur et interpolation contre le calcul direct de R136, six paires de masses, sept valeurs de J : ecart relatif maximal {worst:.1e}")
    ft = lambda ap: total(1e3 / sqrt(2 * pi * ap), None, "a")[0]
    ap0 = b1.golden(ft, 0.6, 1.2, 30)
    ref_tot = ft(ap0)
    n_states = sum(len(st) for st, _, _ in TRAJ.values()); dof = n_states - len(TRAJ) - 1
    s2 = ref_tot / dof
    curv = (ft(ap0 + 0.01) - 2 * ref_tot + ft(ap0 - 0.01)) / 1e-4
    err0 = sqrt(2 * s2 / curv)
    rs0 = 1e3 / sqrt(2 * pi * ap0); drs0 = rs0 * err0 / (2 * ap0)
    n0 = (RS_LAT - rs0) / sqrt(DRS_LAT ** 2 + drs0 ** 2)
    print(f"    ajustement de reference : alpha' = {ap0:.4f} +- {err0:.4f} GeV^-2, chi2 = {ref_tot/(n_states-1)*1e4:.2f}e-4, sqrt(sigma) = {rs0:.1f} +- {drs0:.1f} MeV ({dof} degres de liberte)")
    print(f"    R137 : 0,8864 +- 0,0106 ; 13,50e-4 ; 423,7 +- 2,5 MeV. Ecart au reseau : {RS_LAT - rs0:.1f} MeV, {n0:.2f} ecarts types\n")
    check("C. sans detente : calcul direct de R136 retrouve a mieux que 1e-6 ; ajustement de R137 retrouve (0,8864 +- 0,0106 ; 13,50e-4 ; 423,7 +- 2,5 MeV)",
          worst < 1e-6 and abs(ap0 - 0.8864) < 2e-4 and abs(err0 - 0.0106) < 2e-4 and abs(ref_tot / (n_states - 1) * 1e4 - 13.50) < 0.02 and abs(rs0 - 423.7) < 0.1,
          f"ecart {worst:.0e} ; alpha' = {ap0:.4f} +- {err0:.4f}")
    ref = study(None, "a", ref_tot, s2)

    # D
    print("D. Version coherente : fin de la branche stable")
    ends = {}
    for q1, q2 in PAIRS:
        br = Branch(sig_lat, b1.UNIV_MASS[q1], b1.UNIV_MASS[q2], BEST, "b"); ends[(q1, q2)] = (br.jend, br.lend * HBARC, br.end)
    print("    tension du reseau, meilleur jeu publie ; moment cinetique classique maximal et longueur ou il est atteint :")
    print("      " + " ; ".join(f"{q1}{q2} : {j:.2f} a {l:.2f} fm" for (q1, q2), (j, l, _) in ends.items()))
    uu = {}
    for label, par in (("gamma 0,35", (0.35, 20.0, 6.0)), ("gamma 0,40", BEST), ("gamma 0,45", (0.45, 20.0, 6.0)), ("coin fort (0,50 ; 15 ; 5)", (0.50, 15.0, 5.0)), ("coin faible (0,35 ; 25 ; 8)", (0.35, 25.0, 8.0))):
        br = Branch(sig_lat, 0.06, 0.06, par, "b"); uu[label] = (br.jend, br.lend * HBARC)
    print("    quarks legers, selon le jeu : " + " ; ".join(f"{k} : {j:.2f} a {l:.2f} fm" for k, (j, l) in uu.items()))
    naive_open = Branch(sig_lat, 0.06, 0.06, BEST, "a").end
    a_ref = total(rs0, None, "a")[1]
    need = 6 - a_ref["rho/a"]
    print(f"    version naive : branche {naive_open} jusqu'a J = 12,5. Les etats J = 6 (a6, f6) demandent un moment classique de {need:.2f} (ordonnee {a_ref['rho/a']:.2f}).")
    print("    -> au-dela de 2,1 fm la force d'extremite tombe plus vite que 1 / L^3 : J decroit quand L croit, la rotation est instable")
    print("       (critere de Rayleigh). Avec le meilleur jeu, un meson leger de la trajectoire dominante ne depasse pas J = 6, le plus haut")
    print("       spin observe ; avec gamma = 0,45 les etats J = 6 ne tiennent qu'en forcant l'ajustement, au coin fort ils n'existent pas.")
    print("       Aucun nombre de la CD n'entre ici, et le profil est phenomenologique.\n")
    check("D. version coherente, tension du reseau, meilleur jeu : J classique maximal 5,70 pour les quarks legers, atteint a 2,10 fm ; six paires entre 2,05 et 2,20 fm ; version naive sans maximum",
          5.6 < ends[("u", "u")][0] < 5.8 and all(2.05 < l < 2.20 and e == "maximum de J" for _, l, e in ends.values()) and naive_open == "ouverte" and uu["coin fort (0,50 ; 15 ; 5)"][0] < 5.0,
          f"uu : {ends[('u','u')][0]:.2f} a {ends[('u','u')][1]:.2f} fm ; coin fort {uu['coin fort (0,50 ; 15 ; 5)'][0]:.2f}")

    # E
    print("E. Ajustement des neuf familles au meilleur jeu publie (gamma = 0,40 ; B = 20 ; R0 = 6,0 GeV^-1)")
    ea, eb = study(BEST, "a", ref_tot, s2), study(BEST, "b", ref_tot, s2)
    for label, r in (("sans detente      ", ref), ("version naive     ", ea), ("version coherente ", eb)):
        print(f"    {label}: tension demandee {r['xm']:.1f} MeV ; chi2 x {r['ratio']:.3f}, degradation {r['deg']:.2f} ({sqrt(max(r['deg'],0)):.2f} ecart type) ; "
              f"a 445 MeV fixe : {r['lat']:.1f} ; ecart de tension seul {r['tens']:.2f} ; ecart complet {r['full']:.2f}")
    print(f"    version coherente : fin de branche a J = {eb['jend']:.2f} ; familles dont l'etat le plus haut est au bout de la branche : {', '.join(eb['bound'])}")
    a4 = {}
    a_lat = total(RS_LAT, None, "a")[1]["rho/a"]
    for label, par, mode in (("sans", None, "a"), ("naive", BEST, "a"), ("coherente", BEST, "b")):
        br = total(RS_LAT, par, mode)[3][("u", "u")]; a4[label] = 1e3 * br.at(4 - a_lat)[0]
    print(f"    a4 a la tension du reseau, ordonnee fixee ({a_lat:.2f}) : {a4['sans']:.0f} MeV sans detente, {a4['naive']:.0f} (naive), {a4['coherente']:.0f} (coherente) ; mesure 1967 MeV")
    print(f"    -> la detente abaisse a4 de {a4['sans']-a4['naive']:.0f} a {a4['sans']-a4['coherente']:.0f} MeV ; il en faudrait {a4['sans']-1967:.0f}. Rectificatif de R136 G : le mecanisme va dans le bon sens,")
    print("       il n'a pas la bonne taille pour la corde tournante.")
    rs_base = 1e3 / sqrt(2 * pi * b1.AP_BASE)
    print(f"    a noter : version naive, tension demandee {ea['xm']:.1f} MeV ; valeur de la base {rs_base:.1f} MeV. Coincidence sans portee : cet ajustement est")
    print("       moins bon que l'ajustement sans detente, et la formule de la base n'a pas de mecanisme (R126).\n")
    check("E. meilleur jeu publie : tension demandee 430,1 MeV (naive) et 434,2 (coherente) ; degradation 3,8 et 6,5 ; ecart complet 2,86 et 3,01 contre 2,97 sans detente",
          429.6 < ea["xm"] < 430.6 and 433.6 < eb["xm"] < 434.8 and 3.5 < ea["deg"] < 4.1 and 6.0 < eb["deg"] < 6.9 and 2.78 < ea["full"] < 2.94 and 2.93 < eb["full"] < 3.09
          and 2.90 < ref["full"] < 3.04 and set(eb["bound"]) == {"rho/a", "omega/f"} and 15 < a4["sans"] - a4["naive"] < 60,
          f"{ea['xm']:.1f} et {eb['xm']:.1f} MeV ; {ea['full']:.2f} et {eb['full']:.2f} ; reference {ref['full']:.2f}")

    # F
    print("F. La fourchette publiee")
    rng = {}
    for label, par in (("gamma 0,35", (0.35, 20.0, 6.0)), ("gamma 0,45", (0.45, 20.0, 6.0)), ("coin faible (0,35 ; 25 ; 8)", (0.35, 25.0, 8.0)), ("coin fort (0,50 ; 15 ; 5)", (0.50, 15.0, 5.0))):
        for mode in ("a", "b"):
            if mode == "b" and par[0] == 0.50:
                print(f"    {label:27s} coherente : branche finie a J = {uu[label][0]:.2f} a la tension du reseau ; les etats J = 5 et 6 n'existent pas")
                continue
            r = study(par, mode, ref_tot, s2, 14); rng[(label, mode)] = r
            print(f"    {label:27s} {'naive    ' if mode == 'a' else 'coherente'} : tension demandee {r['xm']:.1f} MeV ; degradation {r['deg']:5.2f} ; ecart de tension seul {r['tens']:.2f} ; ecart complet {r['full']:.2f}"
                  + (f" ; au bout de la branche : {', '.join(r['bound'])}" if r["bound"] else ""))
    fulls = [ea["full"], eb["full"]] + [r["full"] for r in rng.values()]
    print(f"    -> sur toute la fourchette publiee, l'ecart complet va de {min(fulls):.2f} a {max(fulls):.2f} ; il ne descend jamais sous 2,8.\n")
    check("F. fourchette publiee : ecart complet minimal 2,82, jamais sous 2,7 ; coin fort 4,1 (naive) ; gamma = 0,45 coherent : etats J = 6 au bout de la branche, ecart 7,4",
          2.7 < min(fulls) < 2.9 and 3.9 < rng[("coin fort (0,50 ; 15 ; 5)", "a")]["full"] < 4.3 and rng[("gamma 0,45", "b")]["full"] > 6,
          f"de {min(fulls):.2f} a {max(fulls):.2f}")

    # G
    print("G. Les spectres demandent-ils une detente ? (B = 20, R0 = 6,0 GeV^-1, tension libre)")
    prof = {"a": [(0.0, rs0, 0.0)], "b": [(0.0, rs0, 0.0)]}
    for mode in ("a", "b"):
        for g in (0.10, 0.20, 0.30):
            ff = lambda x: total(x, (g, 20.0, 6.0), mode)[0]
            xm = b1.golden(ff, 415.0, 450.0, 13)
            prof[mode].append((g, xm, (ff(xm) - ref_tot) / s2))
        prof[mode].append((0.35, rng[("gamma 0,35", mode)]["xm"], rng[("gamma 0,35", mode)]["deg"]))
        e = ea if mode == "a" else eb
        prof[mode].append((0.40, e["xm"], e["deg"]))
        prof[mode].append((0.45, rng[("gamma 0,45", mode)]["xm"], rng[("gamma 0,45", mode)]["deg"]))
        print(f"    version {'naive    ' if mode == 'a' else 'coherente'} : " + " ; ".join(f"gamma {g:.2f} : {x:.1f} MeV, {d:+.2f}" for g, x, d in prof[mode]))
    lim = {}
    for mode in ("a", "b"):
        p = prof[mode]
        for (g1, _, d1), (g2, _, d2) in zip(p, p[1:]):
            if d1 < 4 <= d2: lim[mode] = g1 + (g2 - g1) * (4 - d1) / (d2 - d1)
    print(f"    -> la degradation croit avec gamma des la premiere valeur : les spectres preferent l'absence de detente. A deux ecarts types")
    print(f"       ils tolerent gamma jusqu'a {lim['a']:.2f} (naive) ou {lim['b']:.2f} (coherente) : la valeur publiee, 0,40, est au bord.\n")
    mono = all(all(d2 > d1 for (_, _, d1), (_, _, d2) in zip(prof[m], prof[m][1:])) for m in ("a", "b"))
    check("G. degradation croissante avec gamma dans les deux versions, positive des gamma = 0,10 ; borne a deux ecarts types : gamma inferieur a 0,41 (naive), 0,36 (coherente)",
          mono and prof["a"][1][2] > 0.3 and 0.38 < lim["a"] < 0.44 and 0.34 < lim["b"] < 0.40, f"bornes {lim['a']:.2f} et {lim['b']:.2f}")

    # H
    print("H. Le critere ecrit avant le calcul : ecart sous 2 ecarts types")
    print(f"    sans detente : {ref['full']:.2f}. Avec la detente publiee : {ea['full']:.2f} (naive), {eb['full']:.2f} (coherente). Seuil : 2,00.")
    print(f"    lecture la plus favorable, la tension seule sans compter la degradation : {ea['tens']:.2f} (naive), {eb['tens']:.2f} (coherente).")
    print(f"    lecture la plus stricte, tension du reseau imposee sans son erreur : {ref['lat']:.1f} sans detente, {ea['lat']:.1f} et {eb['lat']:.1f} avec.")
    print("    -> la detente deplace le desaccord sans le reduire : ce qu'elle gagne sur la tension (6 a 10 MeV sur 21), elle le perd sur")
    print("       l'ajustement. Le critere n'est atteint que dans la lecture qui ignore cette perte, et dans une seule des deux versions.\n")
    check("H. critere : ecart complet 2,86 et 3,01, au-dessus du seuil de 2 ; tension seule 2,09 et 1,61 ; tension imposee 6,5 et 6,4",
          ea["full"] > 2 and eb["full"] > 2 and 2.0 < ea["tens"] < 2.2 and 1.5 < eb["tens"] < 1.75 and 6.2 < ea["lat"] < 6.8 and 6.1 < eb["lat"] < 6.7,
          f"complet {ea['full']:.2f} et {eb['full']:.2f} ; tension seule {ea['tens']:.2f} et {eb['tens']:.2f}")

    print("Verdict : la piste ne tient pas. La detente publiee, appliquee a la corde tournante, remonte la tension demandee par les spectres")
    print("de 424 a 430 MeV (434 dans la version coherente) et degrade l'ajustement d'autant : l'ecart complet au reseau reste a 2,9 et 3,0")
    print("ecarts types, contre 3,0 sans detente. Les spectres ne demandent pas de detente. B3 reste ouverte.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
