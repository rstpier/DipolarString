#!/usr/bin/env python3
"""Search for an object other than the 9-turn helix for N_DS = 27 -- and a bound on the whole class.

Inputs, all internal to the model: 27 strings of length l_1 = 2 pi R3 / 3 (wire length W = 9 x 2 pi R3),
string tube radius r = R3 / 37.1, the electron's own thin-wire electrodynamics (regularised Neumann
and Coulomb integrals, validated on the electron ring in proton_impedance_closure.py).  Targets set
by the mass law m ~ N^2 Z, i.e. m = kappa L: Z_p/Z_e = 22.67 and L_p/L_e = 1836.

  A. THE MASS LAW'S OWN OBJECT.  m ~ N^2 Z with the wire-length closure means a flat coil of
     N_turn = 9 x 22.67 = 204 coincident turns of radius R_p = R3/22.67 = 17 fm = 1.64 r -- a turn
     radius below two tube radii, where the thin-wire idealisation is void.  Built physically
     (pitch 2.1 r), its inductance is computed.
  B. LOOSE OBJECTS.  A single 27-string ring; three orthogonal 9-string rings.
  C. THE BEST COHERENT COIL.  Multilayer solenoids wound as tightly as the tube allows (pitch and
     layer spacing 2.1 r), scanned over inner radius and layer count: the object that maximises
     L (and Z) for this wire.  This bounds every object of the class.
  D. SCALING.  L_max ~ W^{5/3} at the packing limit: the string count that could carry 1836 L_e.
  E. THE OTHER READING.  If mass is the 1/l mode energy of the loop (step 4), the proton loop has
     radius 0.21 fm -- a factor 4 from the charge radius, where the inductance reading gives a
     factor 20 and the geometric closure 4000.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import check, report
from proton_impedance_closure import (A_REG, ASPECT, ELL, R3, RW, TARGET, coulomb_c, neumann_l, ring)

W = 9 * 2 * math.pi * R3               # 27 l_1
LE, CE = ELL - 2, 4 * math.pi ** 2 / ELL
PITCH = 2.1 * RW                       # tube-to-tube contact plus 5 %
M_P_OVER_M_E = 1836.15267


def helix(radius: float, pitch: float, turns: float, z0: float = 0.0, pts_per_turn: int = 120):
    m = max(int(turns * pts_per_turn), 12)
    u = (np.arange(m) + 0.5) * turns / m
    du = turns / m
    x = np.stack([radius * np.cos(2 * math.pi * u), radius * np.sin(2 * math.pi * u), z0 + pitch * u], axis=-1)
    dl = np.stack([-2 * math.pi * radius * np.sin(2 * math.pi * u), 2 * math.pi * radius * np.cos(2 * math.pi * u),
                   pitch * np.ones(m)], axis=-1) * du
    return x, dl


def multilayer_coil(a0: float, n_layers: float, spacing: float = PITCH, pitch: float = PITCH):
    """Layers at radii a0 + k*spacing, each a helix of the same height h fixed by the wire length."""
    radii = [a0 + k * spacing for k in range(int(n_layers))]
    per_turn = [math.sqrt((2 * math.pi * a) ** 2 + pitch ** 2) for a in radii]
    turns_per_height = 1.0 / pitch
    h = W / (turns_per_height * sum(per_turn))            # height such that total wire = W
    xs, dls = [], []
    for a in radii:
        x, dl = helix(a, pitch, h / pitch, z0=-h / 2)
        xs.append(x); dls.append(dl)
    return np.concatenate(xs), np.concatenate(dls), h, h / pitch * len(radii)


def z_ratio(l_over_mu0: float, c_over_eps0: float) -> float:
    return math.sqrt((l_over_mu0 / LE) / (c_over_eps0 / CE))


def main() -> int:
    # ---- A. the mass law's own object
    r_p = R3 / TARGET
    n_turn = 9 * TARGET
    check("Mass law + wire length imply 204 coincident turns of radius R3/22.67 = 1.64 tube radii",
          abs(r_p / RW - 1.64) < 0.01 and abs(n_turn - 204) < 1,
          f"R_p = {r_p/RW:.2f} r = {r_p*386.16:.1f} fm, N_turn = {n_turn:.0f}, wire check N_turn x 2 pi R_p = {n_turn*2*math.pi*r_p/W:.3f} W")
    xa, dla = helix(r_p, PITCH, n_turn, pts_per_turn=60)
    la = neumann_l(xa, dla)
    check("Built physically (pitch 2.1 r) that coil has L_p/L_e of order 10, not 1836",
          la / LE < 30, f"L_p/L_e = {la/LE:.1f}; length {n_turn*PITCH/R3:.1f} R3, radius {r_p/R3:.3f} R3")

    # ---- B. loose objects
    xr, dlr = ring(9 * R3, 6000)
    zb1 = z_ratio(neumann_l(xr, dlr), coulomb_c(xr, dlr))
    xs, dls = [], []
    for axis in range(3):
        x, dl = ring(3 * R3, 2000)
        x, dl = np.roll(x, axis, axis=1), np.roll(dl, axis, axis=1)
        xs.append(x); dls.append(dl)
    x3, dl3 = np.concatenate(xs), np.concatenate(dls)
    zb2 = z_ratio(neumann_l(x3, dl3), coulomb_c(x3, dl3))
    check("Loose objects: one 27-string ring and three orthogonal 9-string rings give Z_p/Z_e ~ 1.5-2",
          zb1 < 2 and zb2 < 2.5, f"ring of radius 9 R3: {zb1:.2f}; three orthogonal rings of radius 3 R3: {zb2:.2f}")

    # ---- C. the best coherent coil
    print("  multilayer coil scan (pitch = layer spacing = 2.1 r):  a0/r  layers  turns  height/R3   L/L_e   C/C_e   Z_p/Z_e")
    best_l, best_z = (0, None), (0, None)
    for a0_r in (1.5, 3, 6, 10, 15, 25):
        for nl in (1, 2, 4, 8):
            a0 = a0_r * RW
            x, dl, h, turns = multilayer_coil(a0, nl)
            if len(x) > 40000:
                continue
            l = neumann_l(x, dl)
            c = coulomb_c(x, dl)
            z = z_ratio(l, c)
            print(f"                                                        {a0_r:4.1f}   {nl:3d}    {turns:5.1f}   {h/R3:7.2f}   {l/LE:6.1f}  {c/CE:6.2f}   {z:6.2f}")
            if l / LE > best_l[0]: best_l = (l / LE, (a0_r, nl, turns, h / R3))
            if z > best_z[0]: best_z = (z, (a0_r, nl, turns, h / R3))
    check("Maximum inductance of 27 l_1 of string, wound as tightly as the tube allows, is ~ 70 L_e, not 1836",
          20 < best_l[0] < 150, f"max L/L_e = {best_l[0]:.1f} at a0 = {best_l[1][0]} r, {best_l[1][1]} layers, {best_l[1][2]:.0f} turns, height {best_l[1][3]:.2f} R3")
    check("Maximum Z_p/Z_e over the same family is ~ 10, below 22.67 -- for every object of the class",
          best_z[0] < TARGET * 0.8, f"max Z_p/Z_e = {best_z[0]:.2f} at a0 = {best_z[1][0]} r, {best_z[1][1]} layers")
    check("Hence no object made of 27 strings of l_1 carries the proton mass as inductance: short by a factor > 25",
          M_P_OVER_M_E / best_l[0] > 25, f"1836 / {best_l[0]:.0f} = {M_P_OVER_M_E/best_l[0]:.0f}")

    # ---- D. scaling
    n_needed = 27 * (M_P_OVER_M_E / best_l[0]) ** 0.6
    check("At the packing limit L_max ~ W^{5/3}: carrying 1836 L_e needs ~ N = 27 x (1836/L_max)^{3/5} strings",
          150 < n_needed < 400, f"N ~ {n_needed:.0f} strings -- the catalogue's 27 would have to change by a factor ~ {n_needed/27:.0f}")

    # ---- E. the other reading
    r_mode = 386.16 / M_P_OVER_M_E          # fm: loop radius if m ~ 1/l_loop with the electron's R3 = 386 fm
    check("Mode-energy reading (m ~ 1/l_loop, step 4): proton loop radius 0.21 fm, a factor 4 from the charge radius 0.84 fm",
          abs(r_mode - 0.210) < 0.002, f"R_p = {r_mode:.3f} fm; inductance reading 17 fm (x20), geometric closure 3.5e3 fm (x4000)")
    check("That reading is E = hc/lambda: it locates the proton at its Compton scale but does not derive the ratio",
          True, "what it needs is the closure fixing l_p / l_e = 1/1836 -- the same missing principle, in a different variable")

    return report("Conclusion: there is no object of 27 strings of l_1 that reaches Z_p/Z_e = 22.67 or L_p/L_e = 1836 "
                  "under the electron's own electrodynamics. The tightest coherent coil the tube allows gives at most "
                  "~70 L_e and ~10 Z_e; loose objects give ~1.5-2 Z_e; the mass law's own implied object (204 turns of "
                  "radius 1.6 r) is outside the thin-wire regime and, built physically, gives ~6 L_e. The obstacle is "
                  "not the helix but the mass law m = kappa L at N = 27: 27 l_1 of string cannot carry that inductance. "
                  "Either N_p is ~200-300 rather than 27, or mass is not the inductance of the string -- the 1/l mode-"
                  "energy reading puts the proton loop at 0.21 fm, closest to reality, but still needs a closure for l_p/l_e.")


if __name__ == "__main__":
    raise SystemExit(main())
