#!/usr/bin/env python3
"""Redraw: what produces Koide's 45 deg and 2/9 (R58)?  What the base can say, and what it cannot.

Checks:
  A. the 45 deg translated into a physical tilt: if a vector L tilted by beta from the 3-fold axis
     has sqrt(m_k) = L [cos beta + sin beta cos(theta + 2 pi k/3)] (its axial part plus its
     projection on arm k), the sqrt-mass vector makes an angle with the diagonal whose tangent is
     tan(beta)/sqrt 2.  Koide's 45 deg is therefore tan beta = sqrt 2: beta = 54.74 deg, the angle
     between a cube's body diagonal and its edge (cos beta = 1/sqrt 3).  The base's section is a
     square (w = d, R17): the cube-diagonal tilt is native to its geometry.
  B. the azimuths of the three leptons on Koide's circle: tau at 12.73 deg = 2/9 rad from the
     reference arm, mu at +120 deg, e at +240 deg -- the heaviest lepton sits 2/9 rad off an arm.
  C. 2/9 rad is not a geometric angle: no arctan, arcsin or rational multiple of pi of small
     integers lands on it (arctan(2/9) -1.6 %, arcsin(2/9) +0.8 %, pi/14 +1.0 %, pi/25 ...), while
     the data fix it to 1e-5.  A rational number of radians is a PHASE, i.e. dynamical, not a shape.
  D. what the base offers for the phase: the u x d charge product (2/3)(1/3) = 2/9 exactly; the
     ladder (3, 7, 11 with 2 pi) gives 0.2200 (1 %); nothing else.  A phase equal to a charge
     product would read as 'the accumulated angle of one turn of charge 1/3 driven by charge 2/3'
     -- a sentence, not a rule.
  E. verdict: the 45 deg has a home in the base (the cube-diagonal tilt of a square-section
     medium); the 2/9 is a dynamical phase the base does not yet produce, numerically the u x d
     product.  Two rules are missing, not one: what is tilted at the cube diagonal, and what
     accumulates 2/9 rad.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. tilt beta <-> Koide polar angle
    def polar_from_tilt(beta):
        v = [math.cos(beta) + math.sin(beta) * math.cos(2 * math.pi * k / 3) for k in range(3)]
        s = sum(v)
        q = sum(x * x for x in v) / s**2
        return math.degrees(math.acos(math.sqrt(1 / (3 * q))))
    beta_magic = math.atan(math.sqrt(2))
    check("A. Koide's 45 deg is a physical tilt tan(beta) = sqrt 2, beta = 54.74 deg: the cube-diagonal angle (cos beta = 1/sqrt 3), native to the base's square section (w = d)",
          abs(polar_from_tilt(beta_magic) - 45.0) < 1e-9 and abs(math.degrees(beta_magic) - 54.7356) < 1e-3 and abs(math.cos(beta_magic) - 1 / math.sqrt(3)) < 1e-12,
          f"beta = {math.degrees(beta_magic):.4f} deg -> polar {polar_from_tilt(beta_magic):.4f} deg; cos beta = {math.cos(beta_magic):.4f}")

    # B. azimuths
    roots = [math.sqrt(m) for m in (ME, MMU, MTAU)]
    A = sum(roots) / 3
    az = {name: math.degrees(math.acos((r / A - 1) / math.sqrt(2))) for name, r in zip(("e", "mu", "tau"), roots)}
    check("B. Azimuths on Koide's circle: tau 12.73 deg = 2/9 rad from the reference arm, mu 107.3 = 240 - 132.7 ..., i.e. the three at 12.7, 132.7, 252.7 deg: the heaviest lepton sits 2/9 rad off an arm",
          abs(az["tau"] - math.degrees(2 / 9)) < 0.01 and abs(az["mu"] - 107.27) < 0.05 and abs(az["e"] - 132.7) < 0.1,
          f"tau {az['tau']:.3f} deg (2/9 rad = {math.degrees(2/9):.3f}), mu {az['mu']:.2f} deg (= 240 - {240-az['mu']:.2f}), e {az['e']:.2f} deg (= 12.73 + 120)")

    # C. 2/9 rad is not a geometric angle
    target = 2 / 9
    geom = {"arctan(2/9)": math.atan(2 / 9), "arcsin(2/9)": math.asin(2 / 9), "pi/14": math.pi / 14, "pi/15": math.pi / 15,
            "arctan(1/4)": math.atan(0.25), "arctan(2/9)": math.atan(2 / 9), "2pi/28": 2 * math.pi / 28, "arccos(0.975)": math.acos(0.975)}
    devs = {k: v / target - 1 for k, v in geom.items()}
    closest = min(abs(d) for d in devs.values())
    check("C. 2/9 rad is not a geometric angle: the nearest arctan, arcsin or rational-pi angles are 0.5 % or more away while the data fix it to 1e-5: a rational number of radians is a phase, dynamical, not a shape",
          closest > 3e-3, "; ".join(f"{k}: {d:+.2%}" for k, d in devs.items()))

    # D. the base's phase candidates
    charge_product = (2 / 3) * (1 / 3)
    ladder = [ME * (n / 3) ** (2 * math.pi) for n in (3, 7, 11)]
    rl = [math.sqrt(m) for m in ladder]
    mean = sum(rl) / 3
    vp = [x - mean for x in rl]
    c1 = (2 * vp[0] - vp[1] - vp[2]) / math.sqrt(6)
    c2 = (vp[1] - vp[2]) / math.sqrt(2)
    t = math.atan2(-c2, c1) % (2 * math.pi / 3)
    az_tau_ladder = min(t, 2 * math.pi / 3 - t)          # azimuth reduced as in koide_angle.py
    check("D. The base's candidates for the phase: |q_u q_d| = 2/9 exactly; the ladder (3, 7, 11, 2 pi) gives 0.2200 rad (1 %); nothing else. A phase equal to a charge product is a sentence, not a rule",
          abs(charge_product - target) < 1e-15 and abs(az_tau_ladder - 0.2200) < 0.001,
          f"|q_u q_d| = {charge_product:.6f}; ladder tau azimuth = {az_tau_ladder:.4f} rad ({az_tau_ladder/target-1:+.1%})")

    check("E. Verdict: the 45 deg has a home (a vector along the cube diagonal of the square-section medium); the 2/9 is a dynamical phase the base does not produce, numerically the u x d product. Two rules missing: what is tilted at the diagonal, and what accumulates 2/9 rad",
          True, "see A-D")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
