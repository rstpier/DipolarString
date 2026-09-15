#!/usr/bin/env python3
"""Proton characteristic impedance: can the 19 % gap on Z_p/Z_e be closed without fitting?

The manuscript's mass law m ~ N^2 Z with m_p/m_e = 1836.15 and N_p/N_e = 9 requires
Z_p/Z_e = 22.67.  Its estimate for the N_DS = 27 soliton -- a closed toroidal helix of 9 turns of
3 strings, "mutual inductance L ~ N^2, series inter-turn capacitance C ~ 1/N" -- gives
N^{3/2} = 27, "within 19 %".  This script computes the same object with the SAME classical
electrodynamics the manuscript uses for the electron, and with only the closures the model
already contains:

  * the electron's own formulas, L_e = mu0 R3 [ln(8R3/r) - 2] and C_e = 4 pi^2 eps0 R3 / ln(8R3/r),
    are reproduced exactly by a Neumann integral and a Coulomb integral regularised at the wire
    radius r (validated below); the same two integrals are then applied to the closed helix;
  * wire length 27 l_1 = 9 x 2 pi R3 (27 strings of l_1 end to end);  N = 9 turns;
  * inter-turn spacing = the model's own matching distance D = 2 cosh(pi) r, or the smallest
    pitch at which the toroidal helix does not cut through itself (horn torus);
  * two readings of the coil's capacitance: the manuscript's series turn-to-turn C_tt/N, and the
    self-capacitance of the whole structure (the reading used for the electron's ring).

No aspect ratio or pitch is adjusted.  R3/r = 37.1 is the manuscript's calibrated value and
enters only through ln(8 R3/r), exactly as for Z_e.  Exit status is zero only if every check
passes; the checks are statements of fact, including negative ones.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import check, report

R3 = 1.0
ASPECT = 37.1
RW = R3 / ASPECT                       # wire (string tube) radius r
D = 2 * math.cosh(math.pi) * RW        # matching distance
ELL = math.log(8 * R3 / RW)            # ln(8 R3 / r) = 5.694
N_TURNS = 9
L_WIRE = N_TURNS * 2 * math.pi * R3    # 27 l_1 = 9 x 2 pi R3
TARGET = 1836.15267 / 81               # Z_p/Z_e required by the mass law
A_REG = RW                             # regulator of the thin-wire integrals


# --------------------------------------------------------------------------- curves

def ring(radius: float, m: int):
    t = (np.arange(m) + 0.5) * 2 * math.pi / m
    x = np.stack([radius * np.cos(t), radius * np.sin(t), np.zeros(m)], axis=-1)
    dl = np.stack([-radius * np.sin(t), radius * np.cos(t), np.zeros(m)], axis=-1) * (2 * math.pi / m)
    return x, dl


def toroidal_helix(rho: float, r_maj: float, n: int, m: int):
    """Closed curve winding n times around a torus of major radius r_maj, minor radius rho."""
    t = (np.arange(m) + 0.5) * 2 * math.pi / m
    dt = 2 * math.pi / m
    c, s, cn, sn = np.cos(t), np.sin(t), np.cos(n * t), np.sin(n * t)
    rad = r_maj + rho * cn
    x = np.stack([rad * c, rad * s, rho * sn], axis=-1)
    drad = -n * rho * sn
    dl = np.stack([drad * c - rad * s, drad * s + rad * c, n * rho * cn], axis=-1) * dt
    return x, dl


def coaxial_stack(rho: float, pitch: float, n: int, m_per: int):
    """n separate coaxial rings (the picture behind 'L ~ N^2'), not a closed curve."""
    xs, dls = [], []
    for k in range(n):
        x, dl = ring(rho, m_per)
        x = x + np.array([0, 0, (k - (n - 1) / 2) * pitch])
        xs.append(x); dls.append(dl)
    return np.concatenate(xs), np.concatenate(dls)


# --------------------------------------------------------------------------- integrals

def neumann_l(x, dl, a: float = A_REG, block: int = 400) -> float:
    """L / mu0 = (1/4pi) sum_ij dl_i . dl_j / sqrt(|x_i - x_j|^2 + a^2)."""
    tot = 0.0
    for i0 in range(0, len(x), block):
        xi, dli = x[i0:i0 + block], dl[i0:i0 + block]
        d2 = np.sum((xi[:, None, :] - x[None, :, :]) ** 2, axis=-1)
        tot += float(np.sum((dli @ dl.T) / np.sqrt(d2 + a * a)))
    return tot / (4 * math.pi)


def coulomb_c(x, dl, a: float = A_REG, block: int = 400) -> float:
    """C / eps0 for a uniformly charged thin wire: Q / <V>, <V> averaged along the wire.
    Exact for a ring; a lower bound (uniform trial charge) for other shapes."""
    ds = np.linalg.norm(dl, axis=-1)
    vsum = 0.0
    for i0 in range(0, len(x), block):
        xi = x[i0:i0 + block]
        d2 = np.sum((xi[:, None, :] - x[None, :, :]) ** 2, axis=-1)
        vsum += float(np.sum(ds[i0:i0 + block] * np.sum(ds[None, :] / np.sqrt(d2 + a * a), axis=1)))
    length = float(np.sum(ds))
    vbar = vsum / length / (4 * math.pi)          # per unit line charge, units 1/eps0
    return length / vbar


def series_c_tt(rho: float, pitch: float, n: int) -> float:
    """Manuscript reading: adjacent turns as a two-wire line, n gaps in series. Units eps0."""
    per_len = math.pi / math.acosh(pitch / (2 * RW))
    return per_len * 2 * math.pi * rho / n


# --------------------------------------------------------------------------- main

def main() -> int:
    # ---- validation on the electron's ring
    xe, dle = ring(R3, 4000)
    le, ce = neumann_l(xe, dle), coulomb_c(xe, dle)
    le_ms, ce_ms = ELL - 2, 4 * math.pi ** 2 / ELL
    check("Regularised Neumann integral reproduces the manuscript's L_e = mu0 R3 [ln(8R3/r) - 2]",
          abs(le / le_ms - 1) < 0.01, f"integral {le:.4f} vs formula {le_ms:.4f} (mu0 R3)")
    check("Regularised Coulomb integral reproduces the manuscript's C_e = 4 pi^2 eps0 R3 / ln(8R3/r)",
          abs(ce / ce_ms - 1) < 0.01, f"integral {ce:.4f} vs formula {ce_ms:.4f} (eps0 R3)")
    ze = math.sqrt(le_ms / ce_ms)
    check("Z_e / Z_0 = 0.730 recovered", abs(ze - 0.730) < 2e-3, f"Z_e/Z_0 = {ze:.4f}")

    # ---- the model's own scaling estimate, restated
    check("Scaling estimate of the manuscript: L ~ N^2, C ~ 1/N gives N^{3/2} = 27, target 22.67 (+19 %)",
          abs(27 / TARGET - 1.191) < 1e-3, f"27 / {TARGET:.3f} = {27/TARGET:.4f}")

    # ---- geometry from the closures
    # (i) wire length: per turn sqrt((2 pi rho)^2 + p^2) = 2 pi R3 ;  (ii) p = D  or  (iii') R_maj = rho
    p_match = D
    rho_match = R3 * math.sqrt(1 - (p_match / (2 * math.pi * R3)) ** 2)
    rmaj_match = N_TURNS * p_match / (2 * math.pi)
    check("At the matching pitch p = D the 9-turn toroidal helix cuts through itself (R_maj < rho)",
          rmaj_match < rho_match, f"R_maj = {rmaj_match:.3f} R3 < rho = {rho_match:.3f} R3 (pitch D = {D:.3f} R3)")
    rho_horn = R3 / math.sqrt(1 + 1 / N_TURNS ** 2)
    rmaj_horn = rho_horn
    p_horn = 2 * math.pi * rmaj_horn / N_TURNS
    check("Smallest self-avoiding closure (horn torus, R_maj = rho) has pitch 1.11 D -- the two closures nearly coincide",
          abs(p_horn / D - 1.11) < 0.02, f"p_horn = {p_horn:.3f} R3 = {p_horn/D:.3f} D, rho = {rho_horn:.4f} R3")

    # ---- nearest-neighbour coupling of two turns at the matching pitch
    x1, dl1 = ring(R3, 3000)
    x2 = x1 + np.array([0, 0, D])
    m12 = neumann_l(np.concatenate([x1, x2]), np.concatenate([dl1, dl1])) - 2 * le
    k1 = (m12 / 2) / le
    check("Two turns at spacing D are weakly coupled: nearest-neighbour k = M/L << 1 (L ~ N^2 assumes k = 1)",
          k1 < 0.3, f"k_1 = {k1:.3f}")

    # ---- the closed helices
    results = {}
    m_hel = 9000
    for name, rho, rmaj, pitch in (("toroidal helix, pitch D", rho_match, rmaj_match, p_match),
                                   ("toroidal helix, horn torus", rho_horn, rmaj_horn, p_horn)):
        x, dl = toroidal_helix(rho, rmaj, N_TURNS, m_hel)
        lp = neumann_l(x, dl)
        cp_self = coulomb_c(x, dl)
        cp_ser = series_c_tt(rho, pitch, N_TURNS)
        results[name] = (lp, cp_self, cp_ser)
    xs, dls = coaxial_stack(R3, D, N_TURNS, 1000)
    lp_stack = neumann_l(xs, dls)
    results["coaxial stack, pitch D (open, the 'N^2' picture)"] = (lp_stack, coulomb_c(xs, dls), series_c_tt(R3, D, N_TURNS))

    print("  geometry                                        L_p/L_e   C_self/C_e  C_series/C_e   Z_p/Z_e(self)  Z_p/Z_e(series)  m_p/m_e = L_p/L_e")
    ratios = []
    for name, (lp, cs, cser) in results.items():
        z_self, z_ser = math.sqrt((lp / le_ms) / (cs / ce_ms)), math.sqrt((lp / le_ms) / (cser / ce_ms))
        ratios += [z_self, z_ser]
        print(f"  {name:48s} {lp/le_ms:8.2f}   {cs/ce_ms:9.3f}   {cser/ce_ms:9.4f}    {z_self:9.2f}       {z_ser:9.2f}        {lp/le_ms:8.1f}")
    print(f"  {'target (mass law)':48s} {'':8s}   {'':9s}   {'':9s}    {TARGET:9.2f}       {TARGET:9.2f}        {1836.15:8.1f}")

    lp_d = results["toroidal helix, pitch D"][0] / le_ms
    check("Actual inductance of the closed 9-turn helix at pitch D is far below 81 L_e",
          lp_d < 30, f"L_p/L_e = {lp_d:.1f} (scaling assumed 81)")
    z_ser_d = math.sqrt(lp_d / (results["toroidal helix, pitch D"][2] / ce_ms))
    check("Series-capacitance reading (the manuscript's) gives Z_p/Z_e ~ 12, not 27",
          8 < z_ser_d < 16, f"Z_p/Z_e = {z_ser_d:.2f}")
    z_self_d = math.sqrt(lp_d / (results["toroidal helix, pitch D"][1] / ce_ms))
    check("Self-capacitance reading (the electron's) gives Z_p/Z_e of order a few",
          z_self_d < 6, f"Z_p/Z_e = {z_self_d:.2f}")
    check("No internal closure reaches the target 22.67 within 5 %",
          all(abs(z / TARGET - 1) > 0.05 for z in ratios), "closest = %.2f" % min(ratios, key=lambda z: abs(z - TARGET)))
    check("Read as inductance = mass (L = kappa m), the same coil gives m_p/m_e ~ 10-20, not 1836",
          lp_d < 100, f"L_p/L_e = {lp_d:.1f}")
    # ---- scan of the one geometric freedom left by the wire-length closure: the pitch
    # The series reading treats adjacent turns as a two-wire line: it applies while the turns
    # are actually side by side, p <= 2 rho (p <= 3.05 D here); beyond that the coil is a wavy ring
    # and only the self reading applies.
    scan = []
    for f in (1.11, 1.5, 2.0, 2.5, 3.0, 5.0, 8.0, 9.0, 9.5):
        pitch = f * D
        rho = R3 * math.sqrt(1 - (pitch / (2 * math.pi * R3)) ** 2)
        x, dl = toroidal_helix(rho, N_TURNS * pitch / (2 * math.pi), N_TURNS, 6000)
        lp = neumann_l(x, dl) / le_ms
        z_self = math.sqrt(lp / (coulomb_c(x, dl) / ce_ms))
        z_ser = math.sqrt(lp / (series_c_tt(rho, pitch, N_TURNS) / ce_ms))
        scan.append((f, pitch / (2 * rho), lp, z_self, z_ser))
    print("  pitch scan (self-avoiding range):  p/D   p/2rho   L_p/L_e   Z_p/Z_e(self)   Z_p/Z_e(series)   series reading applies?")
    for f, adj, lp, zs, zr in scan:
        print(f"                                     {f:4.2f}   {adj:5.2f}   {lp:7.2f}     {zs:7.2f}         {zr:7.2f}          {'yes' if adj <= 1 else 'no'}")
    valid = [row for row in scan if row[1] <= 1.0]
    zs_all = [row[3] for row in scan]
    zr_valid = [row[4] for row in valid]
    i_min = min(range(len(valid)), key=lambda i: valid[i][4])
    check("Series reading has a minimum inside its validity range, ~12 near p = 2D: an extremum exists, at half the target",
          0 < i_min < len(valid) - 1 and 11 < valid[i_min][4] < 13,
          f"min Z_p/Z_e(series) = {valid[i_min][4]:.2f} at p = {valid[i_min][0]:.2f} D; the self reading falls monotonically")
    check("Within the validity range of each reading no pitch reaches 22.67: series spans 12-15, self 1.5-3.2",
          max(zr_valid) < TARGET * 0.95 and max(zs_all) < TARGET * 0.95,
          f"series {min(zr_valid):.1f}-{max(zr_valid):.1f} (p <= 2 rho), self {min(zs_all):.1f}-{max(zs_all):.1f}")
    far = scan[-1]
    check("The series formula crosses 22.67 only where the turns are no longer adjacent (p/2rho ~ 9): outside its validity",
          far[4] > TARGET and far[1] > 5, f"at p = {far[0]} D: Z(series) = {far[4]:.1f}, p/2rho = {far[1]:.1f}")
    check("Coincidences, recorded and NOT used: (N-1)^{3/2} = 22.63 (-0.2 %), 27/2^{1/4} = 22.70 (+0.2 %)",
          abs(8 ** 1.5 / TARGET - 1) < 0.003 and abs(27 / 2 ** 0.25 / TARGET - 1) < 0.003,
          "no mechanism produces either; a match without a mechanism is a fit by another name")
    return report("Conclusion: the '19 %' was the distance between a scaling estimate and the target, not the "
                  "residual of a calculation. Computed with the electron's own electrodynamics and the model's own "
                  "closures, the 9-turn closed helix gives Z_p/Z_e = 14 (series-C reading) or 3 (self-C reading) at the "
                  "matching pitch; within the validity of each reading no pitch reaches 22.67 (series 12-15 with a minimum "
                  "near p = 2D, self 1.5-3.2), and the series formula crosses the target only where the turns are no "
                  "longer adjacent; L_p/L_e = 21 against "
                  "1836 if inductance is mass. The gap is not closable from internal "
                  "conditions: the coil's geometry (turn radius, pitch, major radius) and the meaning of its lumped "
                  "capacitance are not fixed by any axiom, and the toroidal helix self-intersects at the matching "
                  "pitch. What would close it is a principle fixing the coil geometry -- not a better estimate.")


if __name__ == "__main__":
    raise SystemExit(main())
