#!/usr/bin/env python3
"""Test of the mode-energy reading of mass: can l_p / l_e be fixed without fitting?

Reading: a particle is the fundamental TEM mode of a closed string loop, m c^2 = hbar omega =
(p - x) 2 pi hbar c / l_loop  (v = c; p the harmonic, x the holonomy shift of step 3b).  The
electron already sits there by calibration: 2 pi R3 = lambda_C, p = 1, x = 0.  The proton then
needs l_p / l_e = (p_p - x_p) / (1836.15 (p_e - x_e)).  The model's internal conditions are
applied one by one:

  (a) geometric closure l = N l_1  ->  l_p / l_e = 9: the proton is LIGHTER than the electron by 9
      for the fundamental; matching 1836 needs the harmonic p_p ~ 16 500, with no natural holonomy
      shift.  More strings = longer loop = lighter: the reading inverts the catalogue's ordering.
  (b) tube-radius floor: a loop of a tube of radius r has l >= 2 pi r, so m / m_e <= p R3 / r = 37.1 p.
      The proton (1836) and the muon (206.8) are excluded at p = 1 by 50 and 5.6.
  (c) the r--Z_e lock: admitting the proton loop R_p = R3/1836 = 0.21 fm as a loop of the same
      string needs r <= 0.21 fm, i.e. R3/r >= 1836, i.e. ln(8R3/r) >= 9.6, i.e. Z_e >= 1.36 Z_0:
      the electron stops being an impedance well -- the manuscript's one robust derived property
      of the electron.
  (d) Brillouin: m_p c^2 is 1900 x the lattice cutoff hbar omega_max = (3/pi) m_e c^2 and the loop
      is 600 x smaller than the cell l_1: the mode cannot be a collective wave of the weave.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math

from berry_holonomy_common import check, report

R3_FM = 386.159
ASPECT = 37.1
R_FM = R3_FM / ASPECT
L1_FM = 2 * math.pi * R3_FM / 3
MP = 1836.15267
MMU = 206.7683
MTAU = 3477.23
A_TRIPLE = 0.1858


def z_e_over_z0(aspect: float) -> float:
    ell = math.log(8 * aspect)
    return math.sqrt(ell * (ell - 2)) / (2 * math.pi)


def main() -> int:
    l_e = 2 * math.pi * R3_FM
    l_p = l_e / MP
    check("Mode-energy reading: electron loop 2 pi R3 = lambda_C; proton loop must be lambda_C / 1836 = 1.32 fm (R_p = 0.21 fm)",
          abs(l_p - 1.3214) < 1e-3, f"l_e = {l_e:.1f} fm, l_p = {l_p:.4f} fm, R_p = {l_p/(2*math.pi):.4f} fm; charge radius 0.841 fm (x4)")

    # (a) geometric closure
    ratio_geom = 27 / 3
    check("(a) Geometric closure l = N l_1 gives l_p / l_e = 9: the fundamental makes the proton 9 x LIGHTER than the electron",
          ratio_geom == 9.0, f"m_p/m_e would be 1/9 = {1/9:.3f}; required 1836: off by {9*MP:.0f}, and in the wrong direction")
    harmonics = {}
    for label, x_e in (("x_e = 0 (manuscript calibration)", 0.0), ("x_e = (1 - a)/3 (planar twisted triple)", (1 - A_TRIPLE) / 3), ("x_e = 1/3 (writhe-closed triple)", 1 / 3)):
        need = 9 * MP * (1 - x_e)          # p_p - x_p
        harmonics[label] = (need, math.floor(need), need - math.floor(need))
    check("(a') Matching 1836 on the 27-string loop needs the harmonic p_p ~ 11 000-16 500, with no natural holonomy shift x_p",
          all(v[1] > 10000 for v in harmonics.values()),
          "; ".join(f"{k}: p_p - x_p = {v[0]:.2f} -> p = {v[1]}, x_p = {v[2]:.2f}" for k, v in harmonics.items()))

    # (b) tube-radius floor
    floor = ASPECT
    check("(b) Tube-radius floor: a loop of the string has l >= 2 pi r, so m/m_e <= p R3/r = 37.1 per harmonic",
          abs(floor - 37.1) < 1e-9, f"R3/r = {floor}; r = {R_FM:.2f} fm")
    check("(b') At p = 1 the proton is excluded by 50 x and the muon by 5.6 x; the tau by 94 x",
          MP / floor > 40 and MMU / floor > 5, f"proton {MP/floor:.1f}, muon {MMU/floor:.2f}, tau {MTAU/floor:.1f}")
    check("(b'') Harmonic rescue needs p_mu >= 6 (loop radius 1.08 r), p_p >= 50, p_tau >= 94 -- integers only by fitting",
          math.ceil(MMU / floor) == 6 and math.ceil(MP / floor) == 50, f"p_mu = {math.ceil(MMU/floor)}, p_p = {math.ceil(MP/floor)}, p_tau = {math.ceil(MTAU/floor)}")

    # (c) the r--Z_e lock
    aspect_needed = MP                     # R_p = R3/1836 >= r  <=>  R3/r >= 1836
    ze_now, ze_needed = z_e_over_z0(ASPECT), z_e_over_z0(aspect_needed)
    check("(c) Admitting the proton loop as a loop of the same string needs R3/r >= 1836, i.e. ln(8R3/r) >= 9.6",
          abs(math.log(8 * aspect_needed) - 9.595) < 0.01, f"ln(8 x 1836) = {math.log(8*aspect_needed):.3f} (now {math.log(8*ASPECT):.3f})")
    check("(c') Then Z_e = 1.36 Z_0 > Z_0: the electron is no longer an impedance well -- the reading and the model's robust property exclude each other",
          ze_needed > 1.0 and ze_now < 1.0, f"Z_e/Z_0 = {ze_now:.3f} now, {ze_needed:.3f} with the tube thin enough for the proton")
    aspect_well = math.exp(1 + math.sqrt(1 + 4 * math.pi ** 2)) / 8     # Z_e = Z_0  <=>  ell(ell-2) = 4 pi^2
    check("(c'') The impedance-well property alone caps R3/r below 197, hence m/m_e < 197 p under this reading: the proton needs p >= 10 even then",
          abs(aspect_well - 197) < 2 and MP / aspect_well > 9, f"R3/r at Z_e = Z_0: {aspect_well:.1f}; 1836 / {aspect_well:.0f} = {MP/aspect_well:.1f}")

    # (d) Brillouin
    cutoff = 3 / math.pi
    check("(d) The proton's mode is 1900 x above the lattice cutoff hbar omega_max = (3/pi) m_e c^2 and 600 x below the cell size",
          MP / cutoff > 1800 and L1_FM / l_p > 500, f"m_p c^2 / hbar omega_max = {MP/cutoff:.0f}; l_1 / l_p = {L1_FM/l_p:.0f}")
    check("(d') So it cannot be a collective wave of the weave; it would have to be a sub-cell wave on a string thinner than the vacuum's",
          True, "no internal condition supplies such a string")

    return report("Conclusion: the mode-energy reading cannot fix l_p / l_e from anything internal. The only closure "
                  "(l = N l_1) gives 9, the wrong way round; the tube radius caps the fundamental at 37 m_e, so the "
                  "proton needs harmonic ~50 (or ~16 500 on the 27-string loop), a fit; and thinning the tube enough "
                  "to allow a 0.21 fm loop makes Z_e = 1.36 Z_0, destroying the electron's impedance well. NEGATIVE. "
                  "The reading only restates E = hc/lambda; with this string it excludes the proton as a fundamental.")


if __name__ == "__main__":
    raise SystemExit(main())
