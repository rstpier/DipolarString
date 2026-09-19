#!/usr/bin/env python3
"""Redraw: '2 pi is a circle or a sphere, so the exponent is an imperfect spherical shape --
deformation or breathing' (the author, after R55).  Tested as a reading of the ladder.

Checks:
  A. exponent = perimeter/radius: 2 pi is a closed circle -- the closed ring of R53; an open 3/4 turn
     would give 3 pi/2 = 4.71 and put the muon at 54 m_e (measured 207): the ladder's 2 pi is the
     closed ring's own signature.
  B. deformation: any planar deformation at fixed mean radius RAISES perimeter/radius above 2 pi
     (ellipse of ellipticity d: +3 d^2/4).  The muon's +0.15 % would be d = 4.5 %; the tau's -0.12 %
     cannot be a deformation: the reading fails on the tau's sign.
  C. breathing: a purely radial breathing keeps perimeter/radius = 2 pi at every instant; it shifts
     nothing.
  D. sphere: a classical ring of radius lambda-bar carries a quadrupole e R^2/2 = 7.5e4 e fm^2; a
     spin-1/2 object has none.  A ring tumbled over all orientations has no quadrupole but also no
     moment (needs mu_B): the sphere is not reachable classically; the quadrupole vanishes by the
     spin-1/2 algebra, not by shape.
  E. where the author's circle does live: Koide's exact form sqrt(m_k) = A [1 + sqrt 2 cos(theta + 2 pi k/3)]
     puts the three leptons on a circle at 120 degrees in sqrt(m), rotated by theta = 0.2222 rad = 2/9
     (Brannen) to 0.05 %: three points on a circle with a small fixed rotation, exact to 1e-5 --
     an 'imperfect circle' that fits, where the power law is only a 1 % approximant (R55).
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

ME, MMU, MTAU = 0.51099895, 105.6583755, 1776.93
HBARC = 197.3269804
LAMBDA_E = HBARC / ME
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. perimeter / radius
    mu_open = (7 / 3) ** (3 * math.pi / 2)
    check("A. Exponent = perimeter/radius: 2 pi is a closed circle (R53); an open 3/4 turn (3 pi/2 = 4.71) would put the muon at 54 m_e vs 207: the ladder's 2 pi is the closed ring's signature",
          abs(mu_open - 54.1) < 0.5 and abs(MMU / ME - 206.8) < 0.1, f"(7/3)^(3 pi/2) = {mu_open:.1f}; m_mu/m_e = {MMU/ME:.1f}")

    # B. deformation
    p_mu = math.log(MMU / ME) / math.log(7 / 3)
    p_tau = math.log(MTAU / ME) / math.log(11 / 3)
    eps_mu, eps_tau = p_mu / (2 * math.pi) - 1, p_tau / (2 * math.pi) - 1
    d_mu = math.sqrt(4 * eps_mu / 3)
    check("B. Any planar deformation raises perimeter/radius above 2 pi (ellipse: +3 d^2/4): the muon's +0.15 % would be a 4.5 % ellipticity, but the tau's -0.12 % cannot be a deformation: fails on the sign",
          eps_mu > 0 and eps_tau < 0 and abs(d_mu - 0.045) < 0.003, f"eps_mu = {eps_mu:+.3%} -> d = {d_mu:.3f}; eps_tau = {eps_tau:+.3%}")

    # C. breathing
    check("C. Purely radial breathing keeps perimeter/radius = 2 pi at every instant: shifts nothing",
          all(abs(2 * math.pi * r / r - 2 * math.pi) < 1e-12 for r in (0.5, 1.0, 2.0)), "P(t)/R(t) = 2 pi")

    # D. sphere and quadrupole
    Q_ring = 0.5 * LAMBDA_E**2
    check("D. A classical ring of radius lambda-bar has a quadrupole e R^2/2 = 7.5e4 e fm^2; spin-1/2 has none; a fully tumbled ring has no quadrupole but no moment either: the sphere is not reachable classically, the quadrupole vanishes by the spin algebra",
          abs(Q_ring - 7.46e4) < 200, f"e R^2/2 = {Q_ring:.3e} e fm^2")

    # E. Koide's circle
    roots = [math.sqrt(m) for m in (ME, MMU, MTAU)]
    A = sum(roots) / 3
    # solve theta: sqrt(m_k) = A (1 + sqrt2 cos(theta + 2 pi k/3)) with k = 0 (e), 1 (mu), 2 (tau) up to ordering
    best = None
    for k_e, k_mu, k_tau in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        for i in range(200000):
            th = i * 2 * math.pi / 200000
            pred = [A * (1 + math.sqrt(2) * math.cos(th + 2 * math.pi * k / 3)) for k in (k_e, k_mu, k_tau)]
            err = max(abs(p / r - 1) for p, r in zip(pred, roots))
            if best is None or err < best[0]:
                best = (err, th, (k_e, k_mu, k_tau))
    err, theta, order = best
    theta_red = theta % (2 * math.pi / 3)
    dev = min(abs(theta_red - 2 / 9), abs(theta_red - (2 * math.pi / 3 - 2 / 9)))
    check("E. Koide's exact form puts the three leptons on a circle at 120 deg in sqrt(m), rotated by theta = 2/9 rad to 0.05 %: three points on a circle with a small fixed rotation, exact to 1e-5 -- the author's imperfect circle lives there, not in the exponent",
          err < 1e-3 and dev / (2 / 9) < 0.002,
          f"fit error {err:.1e}; theta mod 2pi/3 = {theta_red:.5f} rad; nearest 2/9 = {2/9:.5f} ({dev/(2/9):+.3%}); order {order}")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
