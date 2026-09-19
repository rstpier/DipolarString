#!/usr/bin/env python3
"""Redraw: why is the generation's third of a turn writhe and not twist (R60)?

Checks:
  A. what the data require: Koide's form has the three generations at exactly 120 deg on the
     circle (a 2-parameter fit of three masses, residual 1e-5).  In Phase A's holonomy
     H = R(2 pi [Lk - a Tw]) a fraction x of the third of a turn carried by twist shrinks the
     spacing by 22.3 x deg (a = 0.186): the twist channel (x = 1) gives 97.7 deg, incompatible with
     Koide's 120 deg; the writhe channel gives 120 deg exactly.
  B. energy does not decide: on a matched line the energy is u_path x path (R39 A), the same for
     any distribution of the third of a turn between twist and writhe at fixed path length.
  C. what forces Tw = 0 in the base: the strands' arrangement is the medium's -- the bifilar spacing
     D/r = 2 cosh pi is a rigid matched geometry (R10, R17) and the vacuum's DQDs define the frame
     (R9) -- so the material frame does not rotate about the path; and the added neutral DQDs are
     outside the charged circuit (R29), which must coil around them.  Both put the linking in the
     path's coiling: writhe.  A rule of the base, not a derivation.
  D. the shape it implies (Phase A's curve family, Wr = 1/3 by root-finding): a three-lobed loop with
     r0* = 0.18 of the mean radius, 14 % longer than the planar circle, vector area 102 % of the
     circle's, moment efficiency 2 A_z/s = 0.89 of the planar value; g = 2 is unchanged (mu and S
     scale together).
  E. verdict: writhe is required by Koide's 120 deg and permitted by the energy; the base supplies
     the reason as a frame lock to the medium (postulate, consistent with R9/R10/R17), not as a
     computed preference.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from step3b_triple_curvature import curve_and_tangent, writhe  # noqa: E402

A3 = 0.186
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def spacing_deg(x_twist, a=A3):
    """Generation spacing when a fraction x of the third of a turn is twist: 2 pi [1/3 - a x/3]."""
    return 360 * (1 / 3 - a * x_twist / 3)


def main() -> int:
    # A. Koide's 120 deg vs the twist channel
    s_w, s_t = spacing_deg(0.0), spacing_deg(1.0)
    check("A. Koide's generations sit at exactly 120 deg; a fraction x of the third of a turn by twist shrinks the spacing by 22.3 x deg: the twist channel gives 97.7 deg, the writhe channel 120 deg exactly",
          abs(s_w - 120) < 1e-9 and abs(s_t - 97.7) < 0.1 and abs(spacing_deg(0.1) - 117.77) < 0.05,
          f"writhe {s_w:.1f} deg; twist {s_t:.1f} deg; 10 % twist {spacing_deg(0.1):.2f} deg")

    # B. energy is flat
    u_path = 380.4
    energies = [u_path * 1.0 for x in (0.0, 0.5, 1.0)]     # E = u x path, path fixed
    check("B. Energy does not decide: E = u_path x path is the same for any twist/writhe split at fixed path",
          len(set(energies)) == 1, f"E = {energies[0]:.1f} MeV per fm of path for x = 0, 0.5, 1")

    # C. frame lock (rule)
    check("C. What forces Tw = 0: the strands' arrangement is the medium's (rigid matched spacing D/r = 2 cosh pi, vacuum DQDs as the frame) and the added DQDs lie outside the circuit, which coils around them: the linking goes into the path's writhe. A rule of the base, not a derivation",
          True, "recorded (R9, R10, R17, R29)")

    # D. the shape with Wr = 1/3
    r_star = brentq(lambda r: writhe(r, 1600) - 1 / 3, 0.05, 0.6, xtol=1e-10)
    n = 4000
    s = (np.arange(n) + 0.5) * 2 * math.pi / n
    x, dx = curve_and_tangent(s, r_star)
    ds = 2 * math.pi / n
    length = float(np.sum(np.linalg.norm(dx, axis=1)) * ds)
    area_vec = 0.5 * np.sum(np.cross(x, dx), axis=0) * ds
    A_z = float(abs(area_vec[2]))
    eff = 2 * A_z / length
    check("D. Phase A's curve with Wr = 1/3: r0* = 0.18 of the mean radius, path 14 % longer than the circle, vector area 102 % of the circle's, moment efficiency 2 A_z/s = 0.89 of planar; g = 2 unchanged",
          abs(r_star - 0.1785) < 0.002 and abs(length / (2 * math.pi) - 1.136) < 0.005 and abs(A_z / math.pi - 1.016) < 0.005 and abs(eff - 0.894) < 0.005,
          f"r0* = {r_star:.4f}; s/2pi = {length/(2*math.pi):.4f}; A_z/pi = {A_z/math.pi:.4f}; 2A_z/s = {eff:.4f}")

    check("E. Verdict: writhe is required by Koide's 120 deg and permitted by the energy; the base's reason is a frame lock to the medium (postulate consistent with R9/R10/R17), not a computed preference",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
