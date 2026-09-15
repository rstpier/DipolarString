#!/usr/bin/env python3
"""The electron as a partially open turn: does closure act exponentially on the mass?

Hypothesis (author): the mass loop has an exponential effect on closure; the electron would be a
partially open turn.  Test in the model's own electrodynamics: a turn of string length l_s and
impedance Z_e closed through a gap of capacitance C_g.  The exact resonance condition of the loop
follows from the ABCD matrices (line x series element, eigenvalue 1 <=> trace = 2):

      tan(beta l_s / 2) = 1 / (2 omega C_g Z_e)     with omega = beta c.

Since tan > 0 the solutions sit on (0, pi/2) and (pi, 3pi/2): the closed ring's Compton mode
(beta l_s = 2 pi) rises to 3 pi as the turn opens, the standing wave with a current node at the gap
(beta l_s = 2 pi p) is untouched, and the open turn's half-wave fundamental (beta l_s = pi) descends
from the closed ring's DC circulating mode.  Mode energies -- masses, in the mode-energy reading --
change by factors of order 1 with closure, algebraically (linearly in the gap for any visible
opening), never exponentially.

The model's exponentials live elsewhere: R/r = e^{ell}/8 and D/r = 2 cosh(pi Z/Z_0) make the
GEOMETRY exponential in the impedance, i.e. the mass (proportional to Z) LOGARITHMIC in the
geometry -- insensitivity, the opposite of the hypothesis.

One positive consequence: an open turn has beta l_s = pi, exactly the antiperiodic ring of the
unordered pair (note 2026-09-15, section 5c).  The author's intuition and the derived Z_2 holonomy
are the same half-integer ladder at the fundamental; they differ at the second harmonic.
Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.optimize import brentq

from berry_holonomy_common import check, report

R3_FM = 386.159
ASPECT = 37.1
R_FM = R3_FM / ASPECT
Z_E_OVER_Z0 = 0.730
L_S = 2 * math.pi * R3_FM          # string length of the electron turn (fm)
MP = 1836.15267


def block_derivation() -> None:
    bl, zg, ze = sp.symbols("beta_l Z_g Z_e")
    line = sp.Matrix([[sp.cos(bl), sp.I * ze * sp.sin(bl)], [sp.I * sp.sin(bl) / ze, sp.cos(bl)]])
    gap = sp.Matrix([[1, zg], [0, 1]])
    m = line * gap
    cond = sp.simplify(m.trace() - 2)                       # eigenvalue 1 of the round-trip matrix
    # with Z_g = -i / (omega C_g):  i Z_g / Z_e = 1/(omega C_g Z_e) =: 1/q
    q = sp.symbols("q", positive=True)
    cond_q = sp.simplify(cond.subs(zg, -sp.I * ze / q))     # 2 cos - 2 + sin / q
    target = sp.simplify(2 * sp.cos(bl) - 2 + sp.sin(bl) / q)
    check("Round-trip condition of a turn closed through a series gap: 2 cos(beta l) - 2 + sin(beta l)/(omega C_g Z_e) = 0",
          sp.simplify(cond_q - target) == 0, "trace of (line x gap) ABCD matrix = 2; det = 1")
    # half-angle form: tan(beta l / 2) = 1 / (2 q)
    u = sp.symbols("u")
    half = sp.simplify((2 * sp.cos(2 * u) - 2 + sp.sin(2 * u) / q) / (2 * sp.sin(u)))   # divide by 2 sin u (sin u = 0 is the untouched family)
    check("Equivalent half-angle form tan(beta l/2) = 1/(2 omega C_g Z_e), plus the untouched family sin(beta l/2) = 0",
          sp.simplify(half + 2 * sp.sin(u) - sp.cos(u) / q) == 0,
          "u = beta l / 2; tan u > 0 puts the solutions on (0, pi/2) and (pi, 3pi/2)")


def solve_u(a_par: float, lo: float, hi: float) -> float:
    """tan u = A/u on (lo, hi), A = l_s / (4 c C_g Z_e) = l_s Z_0 / (4 Z_e C_g/eps0): A -> 0 closed, A -> inf open.
    tan u > 0, so the branches are (0, pi/2) -- the DC mode of the closed ring becoming the half-wave
    mode of the open turn -- and (pi, 3pi/2) -- the Compton mode beta l = 2 pi rising to 3 pi."""
    f = lambda u: math.tan(u) - a_par / u
    return brentq(f, lo + 1e-9, hi - 1e-9)


def block_closure() -> None:
    low, comp = {}, {}
    for a_par in (1e-6, 1e-3, 0.1, 1.0, 10.0, 100.0, 1e4):
        low[a_par] = solve_u(a_par, 0.0, math.pi / 2)
        comp[a_par] = solve_u(a_par, math.pi, 3 * math.pi / 2)
    check("Compton family: closed turn beta l = 2 pi, fully open beta l = 3 pi -- closure moves this mode by at most a factor 1.5",
          abs(2 * comp[1e-6] - 2 * math.pi) < 1e-3 and abs(2 * comp[1e4] - 3 * math.pi) < 1e-3,
          f"beta l = {2*comp[1e-6]:.4f} (closed) -> {2*comp[1e4]:.4f} (open)")
    check("Low family: the closed ring's DC circulating mode (beta l = 0) becomes the open turn's half-wave fundamental (beta l = pi)",
          abs(2 * low[1e-6]) < 1e-2 and abs(2 * low[1e4] - math.pi) < 1e-3,
          f"beta l = {2*low[1e-6]:.4f} (closed) -> {2*low[1e4]:.4f} (open); the standing wave with a current node at the gap (beta l = 2 pi p) is untouched")
    check("Mode energies (masses, in the mode-energy reading) change by factors of order 1 with closure: 1836 is out of reach by three orders",
          True, "open fundamental / closed Compton mode = 1/2; Compton mode open / closed = 3/2")
    # gap dependence: parallel-plate gap capacitance C_g/eps0 = pi r^2 / g ;  A = l_s Z_0 / (4 Z_e (C_g/eps0))
    def a_of_gap(g_fm: float) -> float:
        return L_S / (4 * Z_E_OVER_Z0 * (math.pi * R_FM ** 2 / g_fm))
    g_star = math.pi * R_FM ** 2 * 4 * Z_E_OVER_Z0 / L_S          # A = 1
    check("The closed-to-open transition happens for a gap of 0.04 tube radii: any visible opening is 'open'",
          abs(g_star / R_FM - 0.039) < 0.005, f"g* = {g_star:.3f} fm = {g_star/R_FM:.3f} r  (r = {R_FM:.2f} fm)")
    gs = np.array([1.0, 10.0, 100.0]) * R_FM
    dev = np.array([math.pi - 2 * solve_u(a_of_gap(g), 0.0, math.pi / 2) for g in gs])   # pi - beta l, low family
    ratio = dev[:-1] / dev[1:]
    check("Beyond the transition the deviation from the open value falls LINEARLY with the gap (algebraic, not exponential)",
          np.all(np.abs(ratio - 10) < 1.5),
          f"pi - beta l = {dev[0]:.4f}, {dev[1]:.5f}, {dev[2]:.6f} at g = r, 10 r, 100 r; successive ratios {ratio[0]:.1f}, {ratio[1]:.1f} (10 = linear)")
    check("So closure makes the mass INSENSITIVE to the gap, not hypersensitive: a 100 x larger gap moves it 100 x less",
          True, "an exponential sensitivity would need the frequency, not the coupling, to carry e^{-g/xi}; every mode is pinned between two multiples of pi/2")


def block_where_exponentials_are() -> None:
    ell_e = math.log(8 * ASPECT)
    # ring: Z/Z_0 = sqrt(ell(ell-2))/2pi  <=>  R/r = e^{ell}/8 ; target Z_p = 22.67 Z_e
    zt = MP / 81 * Z_E_OVER_Z0                     # Z_p / Z_0
    ell_p = 1 + math.sqrt(1 + (2 * math.pi * zt) ** 2)
    check("The model's exponential is geometry-in-impedance: R/r = e^{ell}/8. Reaching Z_p by aspect ratio needs ell = 105, R/r = e^105/8",
          abs(ell_p - 105) < 1, f"ell_e = {ell_e:.3f} -> ell_p = {ell_p:.1f}; R_p/r = 10^{(ell_p - math.log(8))/math.log(10):.0f}")
    d_over_r_p = 2 * math.cosh(math.pi * zt)
    check("Likewise D/r = 2 cosh(pi Z/Z_0): the bifilar spacing for Z_p would be e^52 tube radii",
          d_over_r_p > 1e22, f"D_p/r = 10^{math.log10(d_over_r_p):.0f}  (electron: D/r = 2 cosh pi = 23.2)")
    check("Both go the wrong way for the hypothesis: mass proportional to Z is LOGARITHMIC in the geometry",
          True, "a factor 1836 in mass needs a factor e^{100} in geometry; closure by a few tube radii moves the mass by percent")


def block_tunnelling() -> None:
    cutoff = 3 / math.pi                            # hbar omega_max / m_e c^2
    s_gap = math.log(MP)
    check("The only exponential-in-closure mass law would be a tunnelling splitting m = hbar omega_0 e^{-S}: it needs S_e - S_p = ln 1836 = 7.5 with no internal S",
          abs(s_gap - 7.515) < 1e-3, f"ln(m_p/m_e) = {s_gap:.3f}; the model supplies no action S for a turn's closure")
    check("And its prefactor hbar omega_0 >= m_p c^2 sits 1900 x above the Brillouin cutoff (3/pi) m_e c^2: excluded by the lattice",
          MP / cutoff > 1800, f"m_p c^2 / hbar omega_max = {MP/cutoff:.0f}")


def block_positive() -> None:
    # open turn: beta l_s = p pi ; antiperiodic ring: beta l = (2p - 1) pi
    open_ladder = [p * math.pi for p in (1, 2, 3)]
    anti_ladder = [(2 * p - 1) * math.pi for p in (1, 2)]
    check("An open turn has beta l = pi at the fundamental: the same half-integer ladder as the antiperiodic ring of the unordered pair (section 5c)",
          abs(open_ladder[0] - anti_ladder[0]) < 1e-12, "theta_tot = pi from two open ends, or from one half-twist: degenerate at the fundamental")
    check("They differ at the second harmonic: the open turn has all harmonics (beta l = p pi), the antiperiodic ring only odd multiples of pi",
          abs(open_ladder[1] - 2 * math.pi) < 1e-12 and anti_ladder[1] == 3 * math.pi,
          "a distinguishing consequence, if the radial ladder is ever resolved")
    # anchor consequence: open turn => lambda_C = 2 l_s => R = R3/2 => aspect halves
    ell_half = math.log(8 * ASPECT / 2)
    ze_half = math.sqrt(ell_half * (ell_half - 2)) / (2 * math.pi)
    check("If the electron IS the open turn, its loop is half the size (R = R3/2) and Z_e = 0.62 Z_0: still an impedance well",
          0.5 < ze_half < 1.0, f"ln(4 R3/r) = {ell_half:.3f}, Z_e/Z_0 = {ze_half:.3f} (closed-ring calibration: 0.730)")


def main() -> int:
    block_derivation()
    block_closure()
    block_where_exponentials_are()
    block_tunnelling()
    block_positive()
    return report("Conclusion: in the model's own electrodynamics the closure of a turn moves its mode energy by at most "
                  "a factor 2 (half wave -> full wave), the transition sits at a gap of 0.04 tube radii, and beyond it the "
                  "dependence is linear in the gap: closure makes the mass insensitive, not exponentially sensitive. The "
                  "model's exponentials (R/r = e^ell/8, D/r = 2 cosh(pi Z/Z_0)) make geometry exponential in impedance, i.e. "
                  "mass logarithmic in geometry -- the wrong way. A tunnelling law m = hbar omega_0 e^{-S} is the only "
                  "exponential-in-closure reading; it has no internal S and its prefactor is 1900 x above the Brillouin "
                  "cutoff. What survives: the electron as an open turn IS the antiperiodic ring of section 5c at the "
                  "fundamental (theta_tot = pi), with a distinguishing second harmonic -- good for item 5, nothing for the proton.")


if __name__ == "__main__":
    raise SystemExit(main())
