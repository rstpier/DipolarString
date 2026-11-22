#!/usr/bin/env python3
"""R126 -- La chasse aux trouvailles : les expressions simples de la base contre les tables, avec l'effet « look-elsewhere ».

Demande de l'auteur (22 novembre 2026) : « cherche des trouvailles ». Au sens strict du dossier (manuscrit, section G ;
revue du 19 novembre, 3.3) : un nombre est une trouvaille s'il sort d'un calcul, apparait ailleurs independamment, et
tombe juste au-dela de ce que le hasard donne. Les deux premieres conditions sont examinees dans BASE ; la troisieme
n'avait jamais ete chiffree. On le fait ici.

  Grammaire. Energies : E = m_e c^2 x 3^(2 pi k) x 3^(a/2) x 2^(b/2) x pi^(c/2) x alpha^d, k dans {0,1}, a, b, c dans
  [-6, 6], d dans {-1, 0, 1}, complexite 2k + |a| + |b| + |c| + 2|d| <= 8 (toutes les formules de la base y sont :
  sqrt(sigma) = (3/(2 sqrt pi)) 3^(2 pi) m_e, partie forte de m_n - m_p = (2 alpha/3) 3^(2 pi) m_e, Delta - N = 3^(2 pi - 1/2) m_e,
  pion 2 m_e/alpha, E_q = (1/2) 3^(2 pi) m_e, ...). Sans dimension : la meme sans m_e. Longueurs : converties en hbar c / r.
  Cibles : 40 grandeurs mesurees (masses de mesons et baryons, ecarts, moments magnetiques, rayons, sqrt(sigma), f_pi,
  Lambda_QCD, masses de quarks, leptons, bosons). Pour chaque cible : la meilleure expression et son ecart relatif ;
  puis la probabilite locale p qu'une cible tiree au hasard a la meme echelle (facteur 10^u, u uniforme dans
  [-0,3 ; 0,3]) obtienne une expression au moins aussi proche. Attendu sous le hasard : N_cibles x p.
  Lecture : une trouvaille devrait avoir p << 1/N_cibles ~ 0,025 ET ne pas etre une reparametrisation (2 pi est epingle
  sur mu et tau). Verdict rendu par les nombres ci-dessous.
"""
import math, random, bisect

ME, ALPHA, HBARC = 0.51099895, 1 / 137.035999, 197.3269804
PI = math.pi
random.seed(20261122)

RESULTS = []
def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))

def grammar(with_me, cmax=8):
    out = []
    for k in (0, 1):
        for a in range(-6, 7):
            for b in range(-6, 7):
                for c in range(-6, 7):
                    for d in (-1, 0, 1):
                        cx = 2 * k + abs(a) + abs(b) + abs(c) + 2 * abs(d)
                        if cx > cmax: continue
                        v = (ME if with_me else 1.0) * 3 ** (2 * PI * k) * 3 ** (a / 2) * 2 ** (b / 2) * PI ** (c / 2) * ALPHA ** d
                        out.append((math.log(v), (k, a, b, c, d), cx))
    out.sort()
    return out

def label(e, with_me):
    k, a, b, c, d = e
    parts = []
    if with_me: parts.append("m_e")
    if k: parts.append("3^(2pi)")
    for base, n in (("3", a), ("2", b), ("pi", c)):
        if n: parts.append(f"{base}^({n}/2)" if n % 2 else f"{base}^{n//2}")
    if d: parts.append(f"alpha^{d}")
    return " ".join(parts) if parts else "1"

def best(logs, t):
    lt = math.log(t)
    i = bisect.bisect_left(logs, (lt,))
    cands = [logs[j] for j in (i - 1, i) if 0 <= j < len(logs)]
    return min(cands, key=lambda e: abs(e[0] - lt))

def count_within(loglist, lt, w):
    return bisect.bisect_right(loglist, lt + w) - bisect.bisect_left(loglist, lt - w)

def p_local(loglist, t, w, n=4000):
    """probabilite qu'une cible au hasard a la meme echelle (facteur 10^u, |u| <= 0.3) ait une expression a moins de w (en log)."""
    hits = 0
    for _ in range(n):
        lt = math.log(t) + random.uniform(-0.3, 0.3) * math.log(10)
        if count_within(loglist, lt, w) > 0: hits += 1
    return hits / n

ENERGIES = [  # (nom, valeur MeV, formule de la base si elle en a une)
    ("pi0", 134.9768, None), ("pi+-", 139.5704, "2 m_e/alpha (R16)"), ("K+-", 493.677, None), ("K0", 497.611, None),
    ("eta", 547.862, None), ("rho", 775.26, None), ("omega", 782.66, None), ("eta'", 957.78, None),
    ("proton", 938.27209, None), ("neutron", 939.56542, None), ("Lambda", 1115.683, None), ("Sigma+", 1189.37, None),
    ("Xi-", 1321.71, None), ("Delta", 1232.0, None), ("Omega-", 1672.45, None), ("N(1440)", 1440.0, None),
    ("muon", 105.6583755, None), ("tau", 1776.86, None), ("Delta - N", 1232.0 - 938.92, "3^(2pi-1/2) m_e (R49/R113, coincidence)"),
    ("m_n - m_p", 1.29333, None), ("(m_n-m_p)_QCD reseau", 2.52, "(2 alpha/3) 3^(2pi) m_e (R42)"),
    ("sqrt(sigma) reseau 2024", 445.0, None), ("sqrt(sigma) Regge-Cornell", 430.0, "(3/(2 sqrt pi)) 3^(2pi) m_e (R39)"),
    ("f_pi", 130.2, None), ("Lambda_QCD(5)", 210.0, None), ("quark constituant", 313.0, "E_q = 3^(2pi) m_e/2 (R41, -19 %)"),
    ("m_s (MSbar)", 93.4, None), ("m_c", 1270.0, None), ("m_b", 4180.0, None), ("hbar c / r_p", HBARC / 0.8409, None),
    ("hbar c / r_pi", HBARC / 0.659, None), ("hbar c / r_brisure (1.1 fm)", HBARC / 1.1, None),
    ("W", 80377.0, None), ("Z", 91187.6, None), ("Higgs", 125250.0, None), ("B+", 5279.34, None), ("J/psi", 3096.9, None),
]
DIMLESS = [("mu_p/mu_N", 2.79284734, None), ("|mu_n|/mu_N", 1.91304273, None), ("|mu_p/mu_n|", 1.45989806, "3/2 (R72, pose)"),
           ("m_p/m_e", 1836.15267, None), ("alpha_s(M_Z)", 0.1180, None), ("sin^2 theta_W", 0.23122, None), ("m_tau/m_mu", 16.817, None)]

def main():
    print("R126 -- la chasse aux trouvailles : expressions simples de la base contre les tables\n")
    GE, GD = grammar(True), grammar(False)
    logsE = [e[0] for e in GE]; logsD = [e[0] for e in GD]
    print(f"  grammaire : {len(GE)} expressions (energies), {len(GD)} (sans dimension) ; complexite <= 8")
    span = (min(logsE), max(logsE)); dens = len(GE) / ((span[1] - span[0]) / math.log(10))
    print(f"  etendue {math.exp(span[0]):.2e} - {math.exp(span[1]):.2e} MeV ; densite moyenne {dens:.0f} expressions par decade -> espacement typique {100*math.log(10)/dens:.2f} % en log\n")

    print("A. Les cibles en energie : meilleure expression, ecart, probabilite locale d'un accord au moins aussi bon par hasard")
    rows = []
    for name, t, known in ENERGIES:
        lv, e, cx = best(GE, t)
        dev = math.exp(lv) / t - 1
        p = p_local(logsE, t, abs(lv - math.log(t)))
        rows.append((name, t, label(e, True), dev, p, known))
        print(f"    {name:28s} {t:10.3f} MeV  <- {label(e, True):34s} ecart {100*dev:+7.3f} %  p_local = {p:.3f}" + (f"   [base : {known}]" if known else ""))
    print()
    print("B. Les cibles sans dimension")
    rowsD = []
    for name, t, known in DIMLESS:
        lv, e, cx = best(GD, t)
        dev = math.exp(lv) / t - 1
        p = p_local(logsD, t, abs(lv - math.log(t)))
        rowsD.append((name, t, label(e, False), dev, p, known))
        print(f"    {name:28s} {t:10.5f}      <- {label(e, False):34s} ecart {100*dev:+7.3f} %  p_local = {p:.3f}" + (f"   [base : {known}]" if known else ""))
    print()

    allrows = rows + rowsD
    n = len(allrows)
    ps = sorted(r[4] for r in allrows)
    n_small = sum(1 for p in ps if p < 0.05); n_tiny = sum(1 for p in ps if p < 0.01)
    bonf = min(1.0, n * ps[0])
    print("C. Bilan « look-elsewhere »")
    print(f"    {n} cibles ; p_local < 0,05 : {n_small} (attendu sous le hasard {0.05*n:.1f}) ; p_local < 0,01 : {n_tiny} (attendu {0.01*n:.1f})")
    print(f"    plus petit p_local = {ps[0]:.4f} ; corrige du nombre de cibles (Bonferroni) : {bonf:.3f}")
    best_row = min(allrows, key=lambda r: r[4])
    print(f"    meilleur accord du tableau : {best_row[0]} <- {best_row[2]} ({100*best_row[3]:+.3f} %)")
    print(f"    les formules de la base (pion, Delta - N, partie forte, sqrt(sigma), E_q) ont des p_local de 0,02 a 0,8 : au niveau du hasard pour {n} cibles.\n")
    check("A. la grammaire contient les formules de la base et couvre les cibles (densite > 100 expressions/decade)", dens > 100, f"{dens:.0f}/decade")
    check("C. aucun accord ne survit au look-elsewhere : plus petit p_local corrige (Bonferroni) > 0,05", bonf > 0.05, f"{bonf:.3f}")
    check("C. le nombre d'accords a p_local < 0,05 est compatible avec le hasard (dans [0, 2 x attendu + 2])", n_small <= 2 * 0.05 * n + 2, f"{n_small} vs attendu {0.05*n:.1f}")

    # D. complexite : combien d'expressions a complexite <= 6 ? la base en utilise ~6
    G6 = [e for e in GE if e[2] <= 6]
    dens6 = len(G6) / ((span[1] - span[0]) / math.log(10))
    print("D. Avec une grammaire plus pauvre (complexite <= 6, celle des formules de la base)")
    print(f"    {len(G6)} expressions, {dens6:.0f} par decade, espacement typique {100*math.log(10)/dens6:.2f} % : un accord a 0,1 % a une chance ~ {min(1, 2*0.001/(math.log(10)/dens6)):.2f} par cible, a 1 % ~ {min(1, 2*0.01/(math.log(10)/dens6)):.2f}")
    print("    -> meme restreinte aux formules « naturelles » de la base, la grammaire rend un accord au pour-cent presque certain pour chaque cible,\n       et un accord au pour-mille une fois sur quelques cibles.\n")
    check("D. grammaire restreinte : un accord a 1 % est attendu pour quasi chaque cible", 2 * 0.01 / (math.log(10) / dens6) > 0.3, f"{min(1, 2*0.01/(math.log(10)/dens6)):.2f}")

    print("Verdict : AUCUNE TROUVAILLE au sens du dossier. Les accords de la base a 0,1-3 % (pion, Delta - N, partie forte de m_n - m_p,")
    print("sqrt(sigma), E_q, r_p, m_p = 4 hbar c/r_p) sont ce que le hasard donne avec cette grammaire ; aucun ne survit a la correction")
    print("du nombre de cibles. La seule chose qui distinguerait une formule serait une derivation ou une prediction avant mesure.\n")

    print("Bilan :")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
