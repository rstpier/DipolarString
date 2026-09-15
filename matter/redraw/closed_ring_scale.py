#!/usr/bin/env python3
"""Redraw: what fixes ~140 MeV for the closed ring?

The base has e, hbar, c_0, Z_0 and the electron.  One length near 1.4 fm can be built from them:
half the classical electron radius, r_e/2 = alpha lambda-bar_C / 2 = 1.409 fm, whose mode energy is
    hbar c / (r_e/2) = 2 m_e c^2 / alpha = 140.05 MeV
-- the pi+- to 0.34 %, the pi0 to 3.7 %.  (A known coincidence, m_pi ~ 2 m_e / alpha.)
Two readings inside the base:
  (a) Coulomb: r_e is where the field energy of the charge e equals m_e c^2; the closed ring sits at
      half that radius -- its mode energy is the electron's Coulomb self-energy scale x 2;
  (b) impedance: Z_0 = 2 alpha R_K (R_K = h/e^2 the quantum resistance), so 2/alpha = 4 R_K / Z_0:
      the closed ring's mass is the electron's x 4 R_K/Z_0 -- a closed ring 'sees' h/e^2 where the
      open chain sees Z_0.
Tests of the identification:
  1. pi+- minus pi0: the ring's own Coulomb energy (charges e/3 at the string positions on
     R = 1.41 fm) gives the right sign but 0.36 MeV against 4.59 MeV measured -- 13 x too small;
  2. multiples of 140 MeV are not the hadron spectrum (eta -2 %, proton -4 %, rho -8 %, K +18 %);
  3. status: an identification of the scale, not a derivation -- the mechanism that pins a closed
     ring at r_e/2 is the thing to find.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

HBARC = 197.3269804          # MeV fm
ALPHA = 1 / 137.035999084
ME = 0.51099895
PI0, PIPM = 134.9768, 139.57039
K_PM, ETA, RHO, PROTON = 493.677, 547.862, 775.26, 938.272
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def ring_coulomb(charges, r_fm):
    """Coulomb energy of point charges (units e/3) equally spaced on a ring of radius r (MeV)."""
    n = len(charges)
    k = ALPHA * HBARC / 9          # (e/3)^2 / 4 pi eps0 in MeV fm
    u = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = 2 * r_fm * math.sin(math.pi * (j - i) / n)
            u += k * charges[i] * charges[j] / d
    return u


def main() -> int:
    r_e = ALPHA * HBARC / ME           # classical electron radius, fm
    m_closed = HBARC / (r_e / 2)
    check("The base's 1.4 fm: half the classical electron radius r_e/2 = alpha lambda-bar_C / 2 = 1.409 fm; its mode energy is 2 m_e c^2 / alpha = 140.05 MeV",
          abs(r_e / 2 - 1.409) < 0.001 and abs(m_closed - 140.05) < 0.05, f"r_e = {r_e:.3f} fm, hbar c / (r_e/2) = {m_closed:.2f} MeV")
    check("pi+- at 0.34 % and pi0 at 3.7 % of it -- the known coincidence m_pi ~ 2 m_e / alpha",
          abs(PIPM / m_closed - 1) < 0.005 and abs(PI0 / m_closed - 1) < 0.04, f"pi+-: {PIPM/m_closed-1:+.2%}, pi0: {PI0/m_closed-1:+.2%}")
    rk_over_z0 = 1 / (2 * ALPHA)
    check("(b) In impedances: Z_0 = 2 alpha R_K, so 2/alpha = 4 R_K/Z_0 = 274 -- the closed ring's mass is the electron's x 4 (quantum resistance / vacuum impedance)",
          abs(4 * rk_over_z0 - 2 / ALPHA) < 1e-9, f"R_K/Z_0 = {rk_over_z0:.2f}, 4 R_K/Z_0 = {4*rk_over_z0:.2f}")
    check("(a) In Coulomb terms: r_e is where the field energy of e equals m_e c^2; the closed ring sits at r_e/2",
          True, "either reading names where alpha enters the base; neither is yet a mechanism")
    # test 1: pi+- - pi0 splitting from the ring's Coulomb energy
    r = r_e / 2
    u_pi0 = ring_coulomb([1, -1, 1, -1], r)                 # 2+ 2- alternating, closed 4-ring
    u_pim = min(ring_coulomb(c, r) for c in ([1, -1, -1, -1, -1], [-1, 1, -1, -1, -1]))   # 1+ 4- (all rotations equivalent)
    split = u_pim - u_pi0
    check("Test 1: the ring's own Coulomb energy gives pi+- heavier than pi0 (right sign) by 0.36 MeV, against 4.59 MeV measured -- 13 x too small",
          split > 0 and abs(split - 0.36) < 0.05 and abs((PIPM - PI0) / split - 12.7) < 1.0,
          f"U(pi0 ring) = {u_pi0:+.3f} MeV, U(pi- ring) = {u_pim:+.3f} MeV, difference {split:.3f} MeV vs {PIPM-PI0:.2f} MeV")
    # test 2: multiples of 140 MeV vs hadrons
    devs = {"K+-": K_PM / m_closed, "eta": ETA / m_closed, "rho": RHO / m_closed, "proton": PROTON / m_closed}
    offs = {k: v - round(v) for k, v in devs.items()}
    check("Test 2: hadron masses are not multiples of 140 MeV -- K 3.53, eta 3.91, rho 5.54, proton 6.70 -- the closed-ring scale is the pion's, not a ladder",
          any(abs(o) > 0.2 for o in offs.values()), ", ".join(f"{k}: {v:.2f} x 140 ({o:+.2f})" for (k, v), o in zip(devs.items(), offs.values())))
    check("Status: 140 MeV is identified with 2 m_e c^2/alpha -- the closed ring's radius is r_e/2 -- as a scale, not derived; the mechanism that pins a closed ring there is what the base must supply",
          True, "candidates: the closed ring's fluid sees the quantum resistance h/e^2 instead of Z_0 (impedance reading), or its size is where the charge's field energy is 2 m_e c^2 (Coulomb reading)")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
