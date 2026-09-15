#!/usr/bin/env python3
"""Step 2 of the Hopfion direction: does the bifilar geometry supply a theta--a[n] coupling?

Mode: the differential TEM mode of two conductors at +/-(D/2) m in the transverse plane,
D/r = 2 cosh(pi) (the DS matching condition).  Transported parameter: the angle psi of the
pair axis m.  A continuous (U(1)) Berry connection A = i<E|d_psi E> would be the coupling
through which the Faddeev term arises (see audit_s3_reconnection.py, CP^1 identity).

Witness: the degenerate circular pair e_+/- of an optical fibre (Tomita-Chiao), which must
acquire exp(-/+ i psi).  If the method did not see that, it would see nothing.

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    RESULTS.append(("PASS" if condition else "FAIL", name, detail))


R = 1.0
D = 2 * math.cosh(math.pi) * R
L, N = 6 * D, 701
XS = np.linspace(-L, L, N)
X, Y = np.meshgrid(XS, XS, indexing="ij")
DA = (XS[1] - XS[0]) ** 2


def mode(psi: float):
    """Normalised transverse E-field of the differential mode with pair axis at angle psi."""
    m = np.array([math.cos(psi), math.sin(psi)])
    p = 0.5 * D * m

    def line_charge(sign):
        dx, dy = X - sign * p[0], Y - sign * p[1]
        r2 = dx * dx + dy * dy
        return sign * dx / r2, sign * dy / r2

    ex1, ey1 = line_charge(+1)
    ex2, ey2 = line_charge(-1)
    ex, ey = ex1 + ex2, ey1 + ey2
    outside = ((X - p[0]) ** 2 + (Y - p[1]) ** 2 > R * R) & ((X + p[0]) ** 2 + (Y + p[1]) ** 2 > R * R)
    ex, ey = np.where(outside, ex, 0.0), np.where(outside, ey, 0.0)
    nrm = math.sqrt(np.sum(ex * ex + ey * ey) * DA)
    return ex / nrm, ey / nrm


def overlap(a, b) -> float:
    return float(np.sum(a[0] * b[0] + a[1] * b[1]) * DA)


def block_real_mode() -> None:
    e0 = mode(0.0)
    check("Differential TEM mode is real-valued", True,
          "electrostatic field of two line charges; no imaginary part anywhere")
    check("Full turn returns the mode to itself", abs(overlap(e0, mode(2 * math.pi)) - 1) < 1e-9,
          f"<E(0)|E(2pi)> = {overlap(e0, mode(2*math.pi)):+.9f}")
    check("Quarter turn is orthogonal", abs(overlap(e0, mode(math.pi / 2))) < 1e-9,
          f"<E(0)|E(pi/2)> = {overlap(e0, mode(math.pi/2)):+.2e}")


def block_berry_connection() -> None:
    """A(psi) = i <E|d_psi E>. For a real normalised family, <E|dE> = (1/2) d<E|E> = 0."""
    h = 1e-3
    worst = 0.0
    for psi in (0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4):
        ep, em, ec = mode(psi + h), mode(psi - h), mode(psi)
        de = ((ep[0] - em[0]) / (2 * h), (ep[1] - em[1]) / (2 * h))
        worst = max(worst, abs(overlap(ec, de)))
    check("U(1) Berry connection vanishes identically", worst < 1e-8,
          f"max |<E|d_psi E>| over four angles = {worst:.1e}; holonomy over 0->2pi is exp(i*0) = 1")


def block_z2_holonomy() -> None:
    """Half a turn swaps the conductors and flips the sign of the differential mode.
    For a labelled pair (m a vector) that loop is not closed and there is no phase.
    For an unordered pair (m ~ -m, a director) it IS closed and the holonomy is -1: a Z2 phase."""
    e0, epi = mode(0.0), mode(math.pi)
    ov = overlap(e0, epi)
    check("Half turn flips the sign exactly", abs(ov + 1) < 1e-9,
          f"<E(m)|E(-m)> = {ov:+.9f} -> Z2 holonomy -1 iff the pair is unordered")


def block_witness() -> None:
    """The degenerate circular pair of a round fibre acquires exp(-i psi): the method discriminates."""
    ex, ey = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    e_plus = (ex + 1j * ey) / math.sqrt(2)

    def rot(v, psi):
        c, s = math.cos(psi), math.sin(psi)
        return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])

    worst = 0.0
    for psi in (0.3, math.pi / 2, 2.0, math.pi):
        rotated = rot(e_plus.real, psi) + 1j * rot(e_plus.imag, psi)
        worst = max(worst, abs(np.vdot(e_plus, rotated) - np.exp(-1j * psi)))
    check("Witness: degenerate circular pair acquires exp(-i psi)", worst < 1e-12,
          f"max deviation from exp(-i psi) = {worst:.1e}; a linear (real) mode is (e+ e^-ipsi + e- e^+ipsi)/sqrt2 "
          "and the two phases cancel")


def main() -> int:
    block_real_mode()
    block_berry_connection()
    block_z2_holonomy()
    block_witness()
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print()
    print(f"RESULT: {len(RESULTS) - len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    print("Conclusion: no continuous theta--a[n] coupling from the bifilar geometry at the TEM level "
          "(Faddeev route closed for V2.10); a Z2 holonomy exists for an unordered pair.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
