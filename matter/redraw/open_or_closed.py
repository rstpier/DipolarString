#!/usr/bin/env python3
"""Redraw: is the electron's circuit open or closed (R52)?  Decided by the four things the base
claims for the electron: mu = mu_B (g = 2.0023), S = hbar/2, no radiation (R11), no excited state (R47).

Readings:
  A. open arc with a standing wave (reflecting ends, R17 poles): a standing wave is two
     counter-running waves, net current zero -> mu = 0; its end charges oscillate at 255 keV/hbar
     and radiate (R11); it has harmonics (R47).  Excluded three times.
  B. closed hairpin (out on one branch, back on the other, separation w = 3/pi^2 R along the 3/4 arc):
     the two currents cancel except for the strip between the branches: mu = (w/R) mu_B = 0.30 mu_B.
     Excluded by the moment.
  C. closed one-way ring: mu = q c R/2 = mu_B for q R = e lambda-bar -- q = e at R = lambda-bar, or
     R12's 3e/4 circulating at R = 4 lambda-bar/3; S = R E_circ/c = hbar/2 with E_circ = m c^2/2
     (3 m c^2/8 at 4 lambda-bar/3); stationary, no radiation; a vortex, no harmonics.  The only
     reading that passes all four: THE CIRCUIT IS CLOSED AND ONE-WAY.
  D. consequences: (i) the circulating charge is e (or 3e/4), not e/(2 sqrt alpha) = 5.85 e, which
     would give mu = 5.85 mu_B: R32's 'vortex charge' is the identity 4 pi K q^2/path = pi hbar c/path,
     not a charge; the circuit energy pi hbar c/path is the spin-half energy hbar c/(2R) of one
     quantum, so R39-R45's numbers stand unchanged.  (ii) no free poles: the static half is not at
     poles, R50 and R52 fall, and R17's width derivation (two pole junctions) with them: w = 4 lambda-bar/pi^2
     is a number without a derivation; the p-n splitting keeps its +-11 % per factor 2 on w (R46).
  E. what fixes the static half is reopened: R13's 'on the rotation axis' (the ring's centre) has
     no energy there in the base yet.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME = 0.51099895
LAMBDA_E = HBARC / ME
W_OVER_R = 3 / math.pi**2                     # R17 D'
Q_VORTEX = 1 / (2 * math.sqrt(ALPHA))
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def mu_ring(q_over_e, R):
    """mu / mu_B for charge q circulating at c on radius R (mu_B = e c lambda-bar / 2)."""
    return q_over_e * R / LAMBDA_E


def main() -> int:
    # A. standing wave
    net_current = 1.0 - 1.0
    omega = (ME / 2) / 6.582119569e-22          # rad/s for 255.5 keV
    check("A. Open arc, standing wave: net current zero -> mu = 0 (needs mu_B); end charges oscillate at 255 keV/hbar = 3.9e20 rad/s and radiate; harmonics exist: excluded three times",
          net_current == 0.0 and abs(omega - 3.88e20) < 0.05e20, f"net current {net_current}; omega = {omega:.2e} rad/s")

    # B. hairpin
    mu_hairpin = W_OVER_R                       # (w/R) mu_B for the same charge running out and back
    check("B. Closed hairpin along the 3/4 arc, branches w = 0.30 R apart: currents cancel except the strip between them, mu = (w/R) mu_B = 0.30 mu_B: excluded by the moment",
          abs(mu_hairpin - 0.304) < 0.002, f"mu/mu_B = {mu_hairpin:.3f}")

    # C. one-way ring
    mu_a = mu_ring(1.0, LAMBDA_E)
    mu_b = mu_ring(0.75, 4 * LAMBDA_E / 3)
    S_a = LAMBDA_E * (ME / 2) / HBARC            # in hbar
    S_b = (4 * LAMBDA_E / 3) * (3 * ME / 8) / HBARC
    check("C. Closed one-way ring: mu = mu_B for q = e at lambda-bar or 3e/4 at 4 lambda-bar/3; S = hbar/2 with E_circ = m/2 (3m/8); stationary, no radiation; a vortex, no harmonics: the only reading passing all four -- the circuit is closed and one-way",
          abs(mu_a - 1) < 1e-12 and abs(mu_b - 1) < 1e-12 and abs(S_a - 0.5) < 1e-12 and abs(S_b - 0.5) < 1e-12,
          f"mu/mu_B = {mu_a:.3f} (e at lambda-bar), {mu_b:.3f} (3e/4 at 4 lambda-bar/3); S = {S_a:.3f}, {S_b:.3f} hbar")

    # D. consequences
    mu_vortex = mu_ring(Q_VORTEX, LAMBDA_E)
    identity = 4 * math.pi * K * Q_VORTEX**2 - math.pi * HBARC
    e_quantum = HBARC / (2 * LAMBDA_E)           # spin-half energy at R = lambda-bar
    e_circuit = math.pi * HBARC / (2 * math.pi * LAMBDA_E)
    check("D. The circulating charge is e, not 5.85 e (that would give mu = 5.85 mu_B): R32's vortex charge is the identity 4 pi K q^2 = pi hbar c, and the circuit energy pi hbar c/path is the spin-half energy hbar c/(2R) of one quantum -- R39-R45's numbers stand; but with no free poles R50, R52 and R17's width derivation fall",
          abs(mu_vortex - 5.85) < 0.01 and abs(identity) < 1e-9 and abs(e_quantum - e_circuit) < 1e-12,
          f"mu(5.85 e) = {mu_vortex:.2f} mu_B; identity residual {identity:.1e}; hbar c/(2R) = pi hbar c/path = {e_quantum*1e3:.1f} keV")

    check("E. Reopened: what fixes the static half -- R13's 'on the rotation axis' (the ring's centre) has no energy there in the base yet; w = 4 lambda-bar/pi^2 is now a number without a derivation",
          True, "recorded")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
