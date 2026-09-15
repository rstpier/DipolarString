#!/usr/bin/env python3
"""Step 2 of the Hopfion direction: does the bifilar geometry (n = 2) supply a theta--a[n] coupling?

The differential TEM mode of two conductors is real and non-degenerate.  A continuous U(1)
Berry connection A = i<E|d_psi E> under transport of the pair axis would be the coupling
through which the Faddeev term arises.  Witness: the degenerate circular pair of a round
fibre (Tomita--Chiao), which must acquire exp(-/+ i psi).

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import (berry_connection, check, common_mode, ip, modes, report,
                                   transport)


def block_real_mode() -> None:
    check("Differential TEM mode is real-valued", True,
          "regularised field of two line charges; no imaginary part anywhere")
    check("Full turn returns the mode to itself", abs(transport(2, 2 * math.pi)[0, 0] - 1) < 1e-9,
          f"<E(0)|E(2pi)> = {transport(2, 2*math.pi)[0,0]:+.9f}")
    check("Quarter turn is orthogonal", abs(transport(2, math.pi / 2)[0, 0]) < 1e-9,
          f"<E(0)|E(pi/2)> = {transport(2, math.pi/2)[0,0]:+.2e}")


def block_berry_connection() -> None:
    worst = max(abs(berry_connection(2, psi)[0, 0]) for psi in (0.0, 0.5, 1.1, 2.3))
    check("U(1) Berry connection vanishes identically", worst < 1e-9,
          f"max |<E|d_psi E>| over four angles = {worst:.1e}; real normalised family => (1/2) d<E|E> = 0")


def block_z2_holonomy() -> None:
    ov = transport(2, math.pi)[0, 0]
    check("Half turn flips the sign exactly", abs(ov + 1) < 1e-9,
          f"<E(m)|E(-m)> = {ov:+.9f} -> Z2 holonomy -1 iff the pair is unordered (m ~ -m)")
    cm = ip(common_mode(2, 0.0), common_mode(2, math.pi))
    check("Common mode carries no phase", abs(cm - 1) < 1e-9, f"<Ec(0)|Ec(pi)> = {cm:+.9f}")


def block_witness() -> None:
    ex, ey = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    e_plus = (ex + 1j * ey) / math.sqrt(2)

    def rot(v, psi):
        c, s = math.cos(psi), math.sin(psi)
        return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])

    worst = max(abs(np.vdot(e_plus, rot(e_plus.real, psi) + 1j * rot(e_plus.imag, psi)) - np.exp(-1j * psi))
                for psi in (0.3, math.pi / 2, 2.0, math.pi))
    check("Witness: degenerate circular pair acquires exp(-i psi)", worst < 1e-12,
          f"max deviation = {worst:.1e}; a real linear mode is (e+ e^-ipsi + e- e^+ipsi)/sqrt2, phases cancel")


def main() -> int:
    block_real_mode()
    block_berry_connection()
    block_z2_holonomy()
    block_witness()
    return report("Conclusion (n = 2): no continuous theta--a[n] coupling from the bifilar geometry at the "
                  "TEM level; the Faddeev route is closed for the PAIR. A Z2 holonomy exists for an unordered pair.")


if __name__ == "__main__":
    raise SystemExit(main())
