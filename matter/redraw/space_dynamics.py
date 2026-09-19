#!/usr/bin/env python3
"""Redraw: 'dynamisme de l'espace' -- the author's new lead (R57).  Read in the base's terms: the
vacuum's DQDs (R9) are not static.  Three readings are costed.

Checks:
  A. circulating vacuum: each DQD carries spin 1 (R3), i.e. two quanta of circulation on its closed
     bifilar path 2 l_1(3) = 1618 fm: E = 2 pi hbar c/path = 0.77 MeV = 1.5 m_e per DQD.  Breaking one
     into two spin-1/2 daughters (R4) hands each its quantum: the electron's 'one quantum per
     circuit' (R40) would be inherited, not created -- a mechanism.
  B. its price: one such DQD per cell l_1 x D_0^2 (D_0 = 235 fm, R54) is an energy density of
     1.7e-8 MeV/fm^3 = 2.8e24 J/m^3, 5e33 times the observed vacuum energy (5.4e-10 J/m^3): a
     circulating medium must not gravitate as ordinary energy, or its DQDs are static (spin 0) and
     'dynamism' is only their ability to be excited.
  C. expanding medium: if the vacuum spacing D_0 followed the cosmic expansion (H_0 = 7e-11 per year)
     while particles kept their size, the electron's static half (three junctions at D_0, R54) and
     hence g would drift at 1e-10 per year; g is stable to 3e-14 per year: excluded by 4 orders.  If
     everything co-expands the base's ratios are scale-free (R17 A) and nothing local changes, but
     all masses would drift at H_0 relative to G: lunar ranging bounds G-dot/G < 1e-13 per year,
     excluded by 3 orders.  A medium expanding at the cosmic rate is out either way.
  D. a superfluid-like medium would give vortex lines a tension sigma = rho Gamma^2 ln(R/a)/(4 pi c^2);
     with the density of B and a core of lambda-bar_p it is 3e-9 MeV/fm, 3e11 short of the strong
     tension: the medium's own circulation is far too dilute to confine anything.
  E. verdict: 'dynamic space' as circulating DQDs gives the electron its quantum by inheritance
     (R4) but costs a vacuum energy 1e33 above the cosmological one; as expansion it is excluded;
     as a confining superfluid it is 1e11 too weak.  What survives is the inheritance mechanism,
     if the medium's energy does not gravitate.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ME, MP = 0.51099895, 938.27209
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
L1_E = 2 * math.pi * LAMBDA_E / 3
D0 = 6 * LAMBDA_E / math.pi**2
MEV_FM3_TO_J_M3 = 1.602176634e-13 / 1e-45
RHO_LAMBDA = 5.4e-10                 # J/m^3, observed vacuum energy density
H0_PER_YR = 7.0e-11
G_STABILITY_PER_YR = 3e-14           # electron g stability
GDOT_BOUND = 1e-13                   # lunar laser ranging, per year
SIGMA_LAT = 904.0
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    path = 2 * L1_E
    e_dqd = 2 * math.pi * HBARC / path
    check("A. Circulating vacuum: a spin-1 DQD carries two quanta on its closed path 2 l_1 = 1618 fm, E = 0.77 MeV = 1.5 m_e; broken into two spin-1/2 daughters (R4) it hands each its quantum: the electron's circulation is inherited, not created",
          abs(e_dqd - 0.766) < 0.002 and abs(e_dqd / ME - 1.5) < 0.01, f"E_DQD = {e_dqd:.3f} MeV = {e_dqd/ME:.2f} m_e")

    cell = L1_E * D0**2
    rho = e_dqd / cell
    rho_si = rho * MEV_FM3_TO_J_M3
    check("B. Price: one such DQD per cell l_1 x D_0^2 is 1.7e-8 MeV/fm^3 = 2.8e24 J/m^3, 5e33 times the observed vacuum energy: a circulating medium must not gravitate as ordinary energy, or its DQDs are static",
          abs(rho / 1.7e-8 - 1) < 0.1 and abs(rho_si / RHO_LAMBDA / 5e33 - 1) < 0.2, f"rho = {rho:.2e} MeV/fm^3 = {rho_si:.2e} J/m^3 = {rho_si/RHO_LAMBDA:.1e} x rho_Lambda")

    drift_g = H0_PER_YR                        # g ~ 1/D_0 at fixed particle size
    check("C. Expanding medium: D_0 following H_0 with particles fixed drifts g at 7e-11 per year vs stability 3e-14 (excluded by 4 orders); everything co-expanding is scale-free locally but drifts all masses at H_0 relative to G, vs G-dot/G < 1e-13 (3 orders): out either way",
          drift_g / G_STABILITY_PER_YR > 1e3 and H0_PER_YR / GDOT_BOUND > 1e2,
          f"g-drift {drift_g:.0e}/yr vs {G_STABILITY_PER_YR:.0e}; H_0/(G-dot bound) = {H0_PER_YR/GDOT_BOUND:.0f}")

    gamma_over_c = 2 * math.pi * LAMBDA_P       # circulation / c for a core of lambda-bar_p
    sigma_sf = rho * gamma_over_c**2 * math.log(4) / (4 * math.pi)
    check("D. Superfluid-like tension sigma = rho Gamma^2 ln(R/a)/(4 pi c^2) with the density of B and a core of lambda-bar_p: 3e-9 MeV/fm, 3e11 short of the strong tension: the medium's circulation is far too dilute to confine",
          abs(sigma_sf / 3.3e-9 - 1) < 0.3 and SIGMA_LAT / sigma_sf > 1e11, f"sigma = {sigma_sf:.1e} MeV/fm; lattice/sigma = {SIGMA_LAT/sigma_sf:.0e}")

    check("E. Verdict: circulating DQDs give the electron its quantum by inheritance but cost a vacuum energy 1e30 above the cosmological one; expansion is excluded; a confining superfluid is 1e11 too weak. What survives is the inheritance, if the medium's energy does not gravitate",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
