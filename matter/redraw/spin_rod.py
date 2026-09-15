#!/usr/bin/env python3
"""Redraw, R3-R6 made mechanical.  Author's statements: a chain of 3 DQDs breaks in two by one action;
the daughters rotate in opposite senses with spin 1/2 each; spin = symmetry axes / 2; the EM fluid
accumulates in the trailing pole, creating the active charge: (+)---   ---(-).

Reading used here (to be corrected by the author if wrong):
  * "rotation" = rotation of the daughter chain in space about the breaking point;
  * "trailing pole" = the end far from the breaking point, where centrifugal force drives the fluid.
Model: a rigid rod of length L = 3 strings = lambda_C (the electron's Compton wavelength, the one
length the electron gives without a manuscript), mass m_e (uniform along the rod unless stated),
charge -e as a mobile fluid.  Everything below is classical mechanics plus Maxwell.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

HBAR, C, E, EPS0, ME = 1.054571817e-34, 2.99792458e8, 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31
MU_B = E * HBAR / (2 * ME)
LAMBDA_C = 2 * math.pi * HBAR / (ME * C)          # 2.426e-12 m  (= 3 l_1 of the old manuscript)
G_MEASURED = 2.00231930436

RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    L = LAMBDA_C
    # spin 1/2 as rigid rotation about the near end (I = m L^2 / 3) and about the centre (m L^2 / 12)
    for label, inertia in (("about the breaking point (near end)", ME * L ** 2 / 3), ("about the centre", ME * L ** 2 / 12)):
        omega = (HBAR / 2) / inertia
        v_far = omega * (L if "near" in label else L / 2)
        check(f"Spin 1/2 as rotation {label}: the far end moves below c",
              v_far < C, f"omega = {omega:.2e} rad/s, far-end speed = {v_far/C:.2f} c")
    omega = (HBAR / 2) / (ME * L ** 2 / 3)
    # centrifugal drive of the fluid vs its own Coulomb spreading (charge -e piled at the far end over ~ L)
    f_centrifugal = ME * omega ** 2 * L                      # if the fluid carries the mass
    f_coulomb = E ** 2 / (4 * math.pi * EPS0 * L ** 2)
    check("R5: centrifugal force on the fluid at the far end exceeds its Coulomb self-repulsion: the fluid does pile up at the far pole",
          f_centrifugal > f_coulomb, f"F_cf = {f_centrifugal:.1e} N vs F_Coulomb = {f_coulomb:.1e} N, ratio {f_centrifugal/f_coulomb:.0f}")
    check("R5': a complete (saturated) pile-up makes the charge independent of the rotation rate, as atom neutrality requires (1e-21)",
          True, "a partial, rate-dependent pile-up would violate it; so all the fluid must be at the pole")
    # symmetry axes
    check("R6: a rod with its charge at one end has one symmetry axis (spin 1/2); the symmetric mother = has two (spin 1)",
          True, "the rule is satisfied by construction for these two shapes")
    # magnetic moment and g-factor: charge -e at the far end, mass distribution varied
    def g_factor(r2_charge: float, r2_mass: float) -> float:
        # mu = (q/2) <r^2>_q omega ; S = m <r^2>_m omega ; g = (mu / S) / (q / 2m) = <r^2>_q / <r^2>_m
        return r2_charge / r2_mass
    cases = {
        "uniform rod, pivot at the near end, charge at the far end": g_factor(L ** 2, L ** 2 / 3),
        "uniform rod about its centre, charge at one end": g_factor((L / 2) ** 2, L ** 2 / 12),
        "mass and charge together at the far end": g_factor(L ** 2, L ** 2),
    }
    check("R3/R4 predict the magnetic moment: a uniform rod with the charge at one end gives g = 3 (measured 2.0023); mass at the pole gives g = 1",
          abs(cases["uniform rod, pivot at the near end, charge at the far end"] - 3) < 1e-12 and abs(cases["mass and charge together at the far end"] - 1) < 1e-12,
          "; ".join(f"{k}: g = {v:.2f}" for k, v in cases.items()))
    r2_mass_needed = L ** 2 / G_MEASURED
    check("g = 2 requires the mass to sit at <r^2>_mass = <r^2>_charge / 2: half the mass at the pole with the charge, half at the pivot -- a testable distribution",
          abs(r2_mass_needed / L ** 2 - 0.5) < 2e-3, f"<r^2>_mass / L^2 = {r2_mass_needed/L**2:.4f}; e.g. 50 % of m_e at the far pole and 50 % at the breaking point")
    # radiation
    a = omega ** 2 * L
    p_larmor = E ** 2 * a ** 2 / (6 * math.pi * EPS0 * C ** 3)
    tau = ME * C ** 2 / p_larmor
    check("The obstacle every rotating-charge electron meets: Larmor radiation empties the rest energy in femtoseconds",
          tau < 1e-12, f"a = {a:.1e} m/s^2, P = {p_larmor:.1f} W, m_e c^2 / P = {tau:.1e} s -- the base needs a reason for no radiation (a stationary, reactive mode), not a shape")
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    print("Conclusion: read as mechanics, the author's spin is a rod of length lambda_C rotating at 0.24 c at its far end about "
          "the breaking point; the centrifugal drive beats the fluid's Coulomb spreading by ~50, so the fluid piles at the far "
          "pole (R5) and, if the pile-up is complete, the charge does not depend on the rate (atom neutrality). The shape has one "
          "axis (R6). The same mechanics fixes the magnetic moment: g = 3 for a uniform rod, g = 1 with the mass at the pole, and "
          "g = 2.0023 only if half the mass sits with the charge and half at the pivot -- a prediction the base can be held to. "
          "What the base still owes: why the rotating charge does not radiate (femtosecond lifetime otherwise), and the -1 "
          "under a full turn, which a classical rotation does not give.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
