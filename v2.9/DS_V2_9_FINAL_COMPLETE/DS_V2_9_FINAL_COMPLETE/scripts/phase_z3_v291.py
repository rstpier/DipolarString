#!/usr/bin/env python3
"""Conditional compact-phase and Z3 coefficient checks for DS V2.9.1.

This script verifies algebraic normalization only.  It does not supply the
fractional monodromy, core action, screening length, or a map from phase-slip
fugacity to Newton's constant.
"""

from scipy.constants import alpha as alpha_codata
from scipy.constants import e, epsilon_0, hbar, mu_0, pi


def main() -> None:
    z_0 = (mu_0 / epsilon_0) ** 0.5
    alpha_impedance = e**2 * z_0 / (4.0 * pi * hbar)

    q = e / 3.0
    winding = 1.0 / 3.0
    strings = 3.0

    # Coefficient multiplying ln(L/a) in S_slip/hbar.
    one_string = pi * hbar * winding**2 / (q**2 * z_0)
    z3_coefficient = strings * one_string
    expected = 3.0 / (4.0 * alpha_impedance)

    print(f"Z0                         = {z_0:.12g} ohm")
    print(f"alpha from impedance       = {alpha_impedance:.12g}")
    print(f"CODATA alpha               = {alpha_codata:.12g}")
    print(f"one-string coefficient     = {one_string:.12g}")
    print(f"three-string coefficient   = {z3_coefficient:.12g}")
    print(f"3/(4 alpha_impedance)      = {expected:.12g}")
    print(f"algebraic relative residual= {abs(z3_coefficient / expected - 1.0):.3e}")
    print("STATUS: conditional scaling coefficient, not a complete instanton action or a derivation of G")


if __name__ == "__main__":
    main()
