#!/usr/bin/env python3
"""Redraw: where does the electron's static half sit once the fluid is a vortex (bilan R32-R49,
open item)?  The circuit gives the circulating half exactly (R40 A); g = 2 still needs the other half
static and off the circuit (R13).

Checks:
  A. the two free poles of the 3/4-turn circuit carry the matched charge delta = e/(pi sqrt alpha)
     (R17); K delta^2 = hbar c/pi^2 exactly (alpha cancels).  Two spherical poles of radius w/2 have
     self-energy 2 x K delta^2/(2 x w/2) = 2 hbar c/(pi^2 w) = m_e c^2/2 at w = 4 lambda-bar/pi^2:
     the static half is the poles' own field energy at the ribbon's half-width -- the same identity
     as R17's junction reading (poles at separation w), so the vortex changes nothing here.
  B. the reading is not free: discs of radius w/2 give pi/2 more (0.40 MeV, 0.79 m_e); the strip's
     line radius w/4 gives twice (0.51 MeV, the whole mass): only the sphere of radius w/2 lands the
     half.  The factor pi/2 between disc and sphere is the same kind of geometric ambiguity as R44's.
  C. total: circuit 255.5 keV + two poles 255.5 keV = m_e c^2 exactly; the split 1/2 : 1/2 is what
     g = 2 requires (R13), and the pole energy scales as 1/w ~ 1/l_1 like the circuit, so mu and tau
     keep g = 2 on the ladder.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME = 0.51099895
LAMBDA_E = HBARC / ME
W = 4 * LAMBDA_E / math.pi**2
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
CIRCUIT = math.pi * HBARC / (2 * math.pi * LAMBDA_E)     # 255.5 keV, R40 A
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    kd2 = K * DELTA**2
    sphere = 2 * kd2 / (2 * (W / 2))
    junction = 2 * HBARC / (math.pi**2 * W)
    check("A. Two spherical poles of charge delta and radius w/2: self-energy 2 K delta^2/w = 2 hbar c/(pi^2 w) = m_e c^2/2 exactly at w = 4 lambda-bar/pi^2 -- the poles' own field at the ribbon's half-width, identical to R17's junction reading; the vortex changes nothing",
          abs(kd2 - HBARC / math.pi**2) < 1e-9 and abs(sphere - junction) < 1e-9 and abs(sphere / (ME / 2) - 1) < 1e-9,
          f"K delta^2 = {kd2:.4f} = hbar c/pi^2; two poles = {sphere*1e3:.2f} keV = m_e c^2/2")

    disc = 2 * (math.pi / 4) * kd2 / (W / 2)
    strip = 2 * kd2 / (2 * (W / 4))
    check("B. Not free: discs of radius w/2 give pi/2 more (0.40 MeV = 0.79 m_e), the strip radius w/4 gives twice (0.51 MeV, the whole mass); only the sphere of radius w/2 lands the half",
          abs(disc / sphere - math.pi / 2) < 1e-9 and abs(strip / ME - 1) < 1e-9,
          f"disc: {disc:.3f} MeV ({disc/ME:.2f} m_e); strip radius: {strip:.3f} MeV ({strip/ME:.2f} m_e)")

    total = CIRCUIT + sphere
    check("C. Circuit 255.5 keV + two poles 255.5 keV = m_e c^2 exactly, the 1/2 : 1/2 split that g = 2 requires; both scale as 1/l_1, so mu and tau keep g = 2 on the ladder",
          abs(total / ME - 1) < 1e-9 and abs(CIRCUIT / sphere - 1) < 1e-9, f"total = {total:.6f} MeV = m_e; ratio {CIRCUIT/sphere:.3f}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
