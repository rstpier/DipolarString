#!/usr/bin/env python3
"""Redraw: does the base's mass law reach the proton and the neutron, and what is the mass table of
the electron's family (the spin-1/2 particles)?

Laws of the base used here (BASE.md, R14-R16):
  - open chain of n strings: m(n) = m_e (n/3)^(2 pi)  (e = 3, mu = 7, tau = 11 at 1 %);
  - ring charge rule q = (n+ - n-) e/3 with n = n+ + n-: n odd -> q in {+-1/3, +-1, +-5/3, ...},
    n even -> q in {0, +-2/3, +-4/3, ...}: the charge fixes the parity of n;
  - closed ring: ~140 MeV = 2 m_e c^2/alpha, nearly independent of n (pions).

Checks:
  A. the parity rule: q = +-1 and +-1/3 need n odd, q = 0 and +-2/3 need n even.
  B. leptons on the ladder at 1 %.
  C. the electron's family (q = -1, n odd) also contains n = 5 (12.7 MeV) and n = 9 (509 MeV),
     charged leptons that do not exist; the family therefore needs n = 3 mod 4 (3, 7, 11, 15, ...),
     and n = 15 (12.6 GeV), 19 (55.6 GeV) are excluded by LEP (m > 102.8 GeV).
  D. proton (q = +1, n odd): the ladder gives 509 MeV (n = 9) or 1795 MeV (n = 11): -46 % / +91 %.
  E. neutron (q = 0, n even): n = 10 gives 986 MeV (+5 %); but the p-n splitting is 0.14 % while
     the ladder's step at n = 10 is 82 %: the pair cannot both sit on the ladder.
  F. closed-ring reading: p/140.05 = 6.70, not an integer (7 rings: +4.5 %).
  G. composite reading (uud from chains n = 4, 4, 5): 18.9 MeV = 2 % of the proton -- 98 % of the
     nucleon mass is binding, for which the base has no mechanism.
  H. quarks on the ladder with the parity rule: deviations 12 % to 169 % -- not a 1 % law.
Exit status is zero only if every check passes.  PDG 2024 masses (MeV); u, d, s in MS-bar at 2 GeV,
c and b at their own scale, t direct.
"""

from __future__ import annotations

import math

ME = 0.51099895
P = 2 * math.pi
OBS = {  # name: (charge in units of e, observed mass MeV)
    "e": (-1, ME), "mu": (-1, 105.6583755), "tau": (-1, 1776.93),
    "u": (2 / 3, 2.16), "d": (-1 / 3, 4.70), "s": (-1 / 3, 93.5),
    "c": (2 / 3, 1273.0), "b": (-1 / 3, 4183.0), "t": (2 / 3, 172570.0),
    "p": (1, 938.27209), "n": (0, 939.56542),
}
LEP_LIMIT = 102.8e3  # MeV, heavy charged lepton (PDG)
RING = 2 * ME / (1 / 137.035999)  # 140.05 MeV
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def m_chain(n: int) -> float:
    return ME * (n / 3) ** P


def parity_for(q: float) -> int:
    """Parity of n required by the ring charge rule: 3q = n+ - n- has the parity of n."""
    return int(round(3 * q)) % 2


def best_n(q: float, m: float, nmax: int = 30):
    par = parity_for(q)
    cands = [n for n in range(1, nmax + 1) if n % 2 == par]
    return min(cands, key=lambda n: abs(math.log(m_chain(n) / m)))


def main() -> int:
    print("Ladder m(n) = m_e (n/3)^(2 pi):")
    for n in range(1, 25):
        cls = "odd: +-1/3, +-1, ..." if n % 2 else "even: 0, +-2/3, ..."
        print(f"  n = {n:2d}  {m_chain(n):12.4f} MeV   charge class {cls}")

    # A. parity rule
    check("A. Ring charge rule q = (n+ - n-) e/3: q = +-1 and +-1/3 need n odd; q = 0 and +-2/3 need n even",
          parity_for(1) == 1 and parity_for(-1 / 3) == 1 and parity_for(0) == 0 and parity_for(2 / 3) == 0,
          "3q = n+ - n- has the parity of n = n+ + n-")

    # B. leptons
    dev = {k: m_chain(n) / OBS[k][1] - 1 for k, n in (("e", 3), ("mu", 7), ("tau", 11))}
    check("B. Leptons on the ladder: e = 3 exact, mu = 7 at -0.8 %, tau = 11 at +1.0 %",
          all(abs(d) < 0.011 for d in dev.values()),
          ", ".join(f"{k}: {d:+.2%}" for k, d in dev.items()))

    # C. the electron's family
    m5, m9, m15, m19, m23 = (m_chain(n) for n in (5, 9, 15, 19, 23))
    check("C. The family q = -1, n odd, also contains n = 5 (12.7 MeV) and n = 9 (509 MeV): charged leptons that are not observed -> the family must be n = 3 mod 4",
          abs(m5 - 12.66) < 0.05 and abs(m9 - 508.7) < 1,
          f"n = 5: {m5:.2f} MeV, n = 9: {m9:.1f} MeV; a 12.7 MeV charged lepton would appear in pi -> l nu and a 509 MeV one in e+e- collisions: neither exists")
    check("C'. Next members n = 15 (12.6 GeV) and n = 19 (55.6 GeV) are excluded by LEP (charged lepton > 102.8 GeV); n = 23 (185 GeV) is the first allowed",
          m15 < LEP_LIMIT and m19 < LEP_LIMIT and m23 > LEP_LIMIT,
          f"n = 15: {m15/1e3:.1f} GeV, n = 19: {m19/1e3:.1f} GeV, n = 23: {m23/1e3:.0f} GeV vs LEP limit {LEP_LIMIT/1e3:.1f} GeV")

    # D. proton
    mp = OBS["p"][1]
    n_cont = 3 * (mp / ME) ** (1 / P)
    d9, d11 = m_chain(9) / mp - 1, m_chain(11) / mp - 1
    check("D. Proton (q = +1 -> n odd): the ladder offers 509 MeV (n = 9, -46 %) or 1795 MeV (n = 11, +91 %): the proton is not on the ladder",
          abs(d9 + 0.458) < 0.01 and abs(d11 - 0.913) < 0.01 and abs(n_cont - 9.92) < 0.01,
          f"continuous n for 938.3 MeV = {n_cont:.2f} (even side); n = 9: {d9:+.1%}, n = 11: {d11:+.1%}")

    # E. neutron and the p-n pair
    mn = OBS["n"][1]
    d10 = m_chain(10) / mn - 1
    step = (11 / 10) ** P - 1
    split = (mn - mp) / mp
    check("E. Neutron (q = 0 -> n even): n = 10 gives 986 MeV (+5 %); but p-n differ by 0.14 % while neighbouring rungs at n = 10 differ by 82 %: the pair cannot both be chains",
          abs(d10 - 0.0495) < 0.002 and abs(step - 0.82) < 0.01 and split < 0.002,
          f"n = 10: {m_chain(10):.1f} MeV ({d10:+.1%}); rung step 10 -> 11: {step:+.0%}; (m_n - m_p)/m_p = {split:.2%}")

    # F. closed rings
    kp, kn = mp / RING, mn / RING
    check("F. Closed-ring reading (140.05 MeV per ring): p = 6.70 rings, n = 6.71 -- not an integer; 7 rings would be +4.5 %",
          abs(kp - 6.70) < 0.01 and abs(7 * RING / mp - 1.045) < 0.003,
          f"p/140.05 = {kp:.3f}, n/140.05 = {kn:.3f}; 7 x 140.05 = {7*RING:.1f} MeV ({7*RING/mp-1:+.1%})")

    # G. composite
    uud = 2 * m_chain(4) + m_chain(5)
    udd = m_chain(4) + 2 * m_chain(5)
    check("G. Composite reading, p = uud from chains n = 4, 4, 5: 18.9 MeV = 2.0 % of the proton (udd: 28.4 MeV = 3.0 % of the neutron): 97-98 % of the nucleon is binding, for which the base has no mechanism",
          abs(uud / mp - 0.0201) < 0.0005 and abs(udd / mn - 0.0303) < 0.0005,
          f"uud = {uud:.2f} MeV ({uud/mp:.1%}), udd = {udd:.2f} MeV ({udd/mn:.1%})")

    # H. quarks
    print("\nFermion table (best n of the parity fixed by the charge):")
    print(f"  {'':5s} {'q':>5s} {'n':>3s} {'ladder (MeV)':>14s} {'observed (MeV)':>15s} {'dev':>8s}")
    qdev = {}
    for name, (q, m) in OBS.items():
        n = best_n(q, m)
        d = m_chain(n) / m - 1
        if name in ("u", "d", "s", "c", "b", "t"):
            qdev[name] = (n, d)
        print(f"  {name:5s} {q:+5.2f} {n:3d} {m_chain(n):14.3f} {m:15.3f} {d:+8.1%}")
    expected = {"u": (4, 0.44), "d": (5, 1.69), "s": (7, 0.12), "c": (10, -0.225), "b": (13, 0.22), "t": (22, -0.19)}
    ok = all(qdev[k][0] == n and abs(qdev[k][1] - d) < 0.01 for k, (n, d) in expected.items())
    check("H. Quarks on the ladder with the parity rule: u n=4 +44 %, d n=5 +169 % (n=3 would be -89 %), s n=7 +12 %, c n=10 -23 %, b n=13 +22 %, t n=22 -19 %: a scatter of 12-169 %, not the 1 % of the leptons",
          ok and abs(m_chain(3) / OBS["d"][1] - 0.109) < 0.002,
          "; ".join(f"{k}: n={n} {d:+.0%}" for k, (n, d) in qdev.items()) + f"; d at n=3: {m_chain(3)/OBS['d'][1]-1:+.0%}")

    print()
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
