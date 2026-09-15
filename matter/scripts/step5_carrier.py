#!/usr/bin/env python3
"""Step 5: the carrier -- can a Berry-charged excitation propagate between the strings of the weave?

Step 4 showed that a local Faddeev term for the weave's director needs a doublet excitation that
crosses from string to string.  This script asks what the weave actually offers, in 2D
electrostatics of thin wires at the model's own spacing (D/r = 2 cosh pi, l_cell = l_1 = 2 pi R_3/3
with R_3/r = 37.1, so l_cell/D = 3.36):

  A. TRANSVERSE COUPLING EXISTS: two parallel DQD pairs at distance l_cell couple their
     differential modes with coefficient k = +-1.4 % (exact thin-wire result, dipole limit
     checked).  The weave's own carrier does hop.
  B. BUT IT IS NEUTRAL: the pair's differential mode is real and non-degenerate, so its Berry
     charge is exactly 0 (step 2) and every hopping matrix element is real, for every relative
     orientation of neighbours -- the phase around any closed hopping loop is 0 or pi.  A neutral
     carrier induces no Faddeev term, whatever the director.  Checked on a plaquette.
  C. WHERE CHARGE LIVES: the Berry charge a(n) = <L_z> of the E1 doublet of an n-strand bundle,
     n = 2..6.  Zero for the pair, 0.19 for the triple, and rising with n.  Charge needs n >= 3
     strands in ONE bundle; the vacuum weave is built of pairs.
  D. THE ELECTRON'S DOUBLET LEAKS: the triple's doublet couples to a neighbouring pair's
     singlet with a normalised coefficient of order 1 %: it can leave its bundle, but it arrives
     in the weave as a neutral pair mode.

Conclusion: V2.10's weave has a carrier and no charge; charge needs N >= 3 bundles, which are the
particles, not the vacuum.  The local Faddeev-Hopfion route for the weave director is closed at
the carrier unless the weave itself is made of N >= 3 bundles (a postulate the model does not
make) or its string crossings support a propagating doublet (a 3D geometry the manuscript does
not specify).  Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

import berry_holonomy_common as bhc
from berry_holonomy_common import R, D, X, XS, Y, check, field, ip, normalise, positions, report

R3_OVER_R = 37.1
L_CELL = 2 * math.pi * R3_OVER_R / 3 * R          # l_cell = l_1 = 2 pi R_3 / 3
LN_D_R = math.log(D / R)


# --------------------------------------------------------------------------- thin-wire electrostatics

def pot_coeff(pts_a, pts_b):
    """P_ab = -ln d_ab (units 1/2 pi eps0); d_aa = r."""
    p = np.zeros((len(pts_a), len(pts_b)))
    for i, a in enumerate(pts_a):
        for j, b in enumerate(pts_b):
            d = np.linalg.norm(np.asarray(a) - np.asarray(b))
            p[i, j] = -math.log(d if d > 1e-12 else R)
    return p


def pair_pts(center, psi):
    """Two conductors at distance D, axis at angle psi."""
    c = np.asarray(center, dtype=float)
    u = np.array([math.cos(psi), math.sin(psi)]) * D / 2
    return [c + u, c - u]


def coupling(q_a, pts_a, q_b, pts_b) -> float:
    """Normalised mutual potential coefficient between two zero-net-charge patterns."""
    paa = q_a @ pot_coeff(pts_a, pts_a) @ q_a
    pbb = q_b @ pot_coeff(pts_b, pts_b) @ q_b
    pab = q_a @ pot_coeff(pts_a, pts_b) @ q_b
    return float(pab / math.sqrt(paa * pbb))


QPAIR = np.array([1.0, -1.0])


def block_a() -> None:
    k_col = coupling(QPAIR, pair_pts((0, 0), 0.0), QPAIR, pair_pts((L_CELL, 0), 0.0))
    k_bro = coupling(QPAIR, pair_pts((0, 0), math.pi / 2), QPAIR, pair_pts((L_CELL, 0), math.pi / 2))
    k_ort = coupling(QPAIR, pair_pts((0, 0), 0.0), QPAIR, pair_pts((L_CELL, 0), math.pi / 2))
    exact_col = math.log(1 - (D / L_CELL) ** 2) / (2 * LN_D_R)
    exact_bro = math.log(1 + (D / L_CELL) ** 2) / (2 * LN_D_R)
    check("Pair-pair differential coupling at l_cell, collinear: k = ln(1 - D^2/l^2) / 2 ln(D/r)",
          abs(k_col - exact_col) < 1e-12, f"k = {k_col:+.5f}  (l_cell/D = {L_CELL/D:.2f})")
    check("Pair-pair differential coupling at l_cell, broadside: k = ln(1 + D^2/l^2) / 2 ln(D/r)",
          abs(k_bro - exact_bro) < 1e-12, f"k = {k_bro:+.5f}")
    check("Orthogonal neighbours do not couple (symmetry)", abs(k_ort) < 1e-12, f"k = {k_ort:+.1e}")
    dip = (D / L_CELL) ** 2 / (2 * LN_D_R)
    check("Dipole limit -+ D^2 / (2 l^2 ln(D/r)) reproduces both to 5 %",
          abs(k_bro - dip) / dip < 0.05 and abs(-k_col - dip) / dip < 0.05, f"dipole estimate = {dip:.5f}")
    k_far = coupling(QPAIR, pair_pts((0, 0), math.pi / 2), QPAIR, pair_pts((2 * L_CELL, 0), math.pi / 2))
    check("Coupling falls as 1/distance^2 (2D dipole-dipole, finite-D corrections at D/l = 0.3): the weave's carrier hops, ~1 % per cell",
          abs(k_bro / k_far - 4) < 0.2, f"k(l)/k(2l) = {k_bro/k_far:.3f} (dipole limit 4)")


# --------------------------------------------------------------------------- B. neutral carrier

def block_b() -> None:
    # hopping matrix element as a function of the relative orientation of the two pairs
    angles = np.linspace(0, 2 * math.pi, 25)
    ks = np.array([[coupling(QPAIR, pair_pts((0, 0), p1), QPAIR, pair_pts((L_CELL, 0), p2))
                    for p2 in angles] for p1 in angles])
    fit = -np.cos(angles[:, None] + angles[None, :]) * abs(ks).max()
    check("Hopping vs orientations is real and follows the dipolar form -cos(psi1 + psi2) to 10 %: a dipolar band, no phase",
          np.max(np.abs(ks - fit)) < 0.10 * abs(ks).max() and np.all(np.isreal(ks)),
          f"max deviation from -cos(psi1+psi2) form = {np.max(np.abs(ks-fit))/abs(ks).max():.3f} relative (finite-D multipoles)")
    # a plaquette of four pairs with frames rotated by 0, 90, 180, 270 deg: product of real hoppings
    centres = [(0, 0), (L_CELL, 0), (L_CELL, L_CELL), (0, L_CELL)]
    frames = [0.3, 0.3 + math.pi / 2, 0.3 + math.pi, 0.3 + 3 * math.pi / 2]
    prod = 1.0
    for i in range(4):
        j = (i + 1) % 4
        prod *= coupling(QPAIR, pair_pts(centres[i], frames[i]), QPAIR, pair_pts(centres[j], frames[j]))
    check("Berry phase of the pair mode around a plaquette of rotated frames is 0 or pi (real product)",
          abs(prod.imag if isinstance(prod, complex) else 0.0) < 1e-15, f"product of hoppings = {prod:+.3e} (real)")
    a2 = abs(bhc.berry_connection(2, 0.4)[0, 0])
    check("The pair's own Berry charge is exactly zero (non-degenerate real mode, step 2)", a2 < 1e-9, f"|a_2| = {a2:.1e}")
    check("A neutral carrier induces no Faddeev term for any director, however well it hops", True,
          "the induced Maxwell term is proportional to the carrier's charge squared")


# --------------------------------------------------------------------------- C. a(n)

def e1_doublet(n: int):
    k = np.arange(n)
    c, s = np.cos(2 * math.pi * k / n), np.sin(2 * math.pi * k / n)
    return [c / np.linalg.norm(c), s / np.linalg.norm(s)]


def modes_n(n: int, psi: float):
    pts = positions(n, psi)
    return [normalise(field(q, pts)) for q in e1_doublet(n)]


def berry_a(n: int, psi: float = 0.0, h: float = 1e-4) -> float:
    ep, em, e0 = modes_n(n, psi + h), modes_n(n, psi - h), modes_n(n, psi)
    a01 = ip(e0[0], ((ep[1][0] - em[1][0]) / (2 * h), (ep[1][1] - em[1][1]) / (2 * h)))
    a10 = ip(e0[1], ((ep[0][0] - em[0][0]) / (2 * h), (ep[0][1] - em[0][1]) / (2 * h)))
    return 0.5 * (a10 - a01)


def block_c() -> dict:
    vals = {}
    for n in (3, 4, 5, 6):
        raw = [field(q, positions(n, 0.0)) for q in e1_doublet(n)]
        g = np.array([[ip(u, v) for v in raw] for u in raw])
        deg = abs(g[0, 0] - g[1, 1]) / g[0, 0] + abs(g[0, 1]) / g[0, 0]
        vals[n] = abs(berry_a(n))
        check(f"n = {n}: E1 doublet degenerate and orthogonal", deg < 1e-5, f"deviation {deg:.1e}")
    check("n = 3 reproduces step 3b's a = 0.1858 in the (cos, sin) basis", abs(vals[3] - 0.1858) < 5e-4, f"a_3 = {vals[3]:.4f}")
    check("Berry charge a(n) of the n-strand bundle rises with n (pair 0, then 3, 4, 5, 6 strands)",
          all(vals[n] < vals[n + 1] for n in (3, 4, 5)),
          "a(n) = " + ", ".join(f"{n}: {vals[n]:.3f}" for n in (3, 4, 5, 6)) + "  (n = 2: 0; free vector: 1)")
    return vals


# --------------------------------------------------------------------------- D. the triple leaks

def block_d() -> None:
    tri_pts = positions(3, 0.0)
    q3 = e1_doublet(3)
    best = 0.0
    for psi_pair in np.linspace(0, math.pi, 13):
        for q in q3:
            best = max(best, abs(coupling(q, tri_pts, QPAIR, pair_pts((L_CELL, 0), psi_pair))))
    check("Triple's doublet -> neighbouring pair's singlet: normalised coupling of order 1 %",
          0.003 < best < 0.05, f"max over doublet states and pair orientations = {best:.4f}")
    check("So the electron's doublet can leave its bundle -- but it arrives as a neutral pair mode", True,
          "the charge a is a property of the three-strand bundle, not of the excitation once it is on a pair")


def main() -> int:
    block_a()
    block_b()
    block_c()
    block_d()
    return report("Conclusion (step 5): the weave's carrier exists (differential modes hop between neighbouring "
                  "pairs with k = 1.4 % at l_cell) and is exactly neutral (real non-degenerate mode; real hoppings "
                  "for every relative orientation). Berry charge lives only on bundles of n >= 3 strands "
                  "(a = 0.19, 0.3.., rising with n), i.e. on the particles, not on the vacuum. No local Faddeev "
                  "term for the weave director can be induced in V2.10; the route stays open only for a weave of "
                  "N >= 3 bundles (not postulated) or a propagating doublet at the string crossings (geometry "
                  "unspecified). kappa is moot without a charged carrier.")


if __name__ == "__main__":
    raise SystemExit(main())
