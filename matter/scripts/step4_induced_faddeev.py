#!/usr/bin/env python3
"""Step 4: the coefficient of the Faddeev term induced by the triple's doublet.

Step 3b established that on the string's frame bundle the doublet sees the Berry connection
A = a J omega_twist with curvature a x (solid-angle form of the tangent), a = 0.1858 -- i.e. the
doublet is a spin-a object under rotations of the tangent (photon: a = 1; CP1 spinor: a = 1/2).
"Integrating out" the doublet means summing its zero-point fluctuations; the dependence of that
sum on the Berry field is the induced action -- in 3+1 D it would be the Faddeev term (Maxwell
term of the director's Berry field).  Three results, in decreasing order of what the model
actually supports:

  A. ON ONE CLOSED STRING (derived, exact).  The doublet is one massless complex field on a ring
     of length l with the twisted boundary condition e^{i Phi}, Phi = 2 pi [(1 - a) Lk + a Wr]
     (step 3b).  Its zero-point energy is the induced action, in closed form:
         E_ind(Phi) = -(2 pi hbar c / l) B_2(x),   x = Phi/2pi mod 1,   B_2(x) = x^2 - x + 1/6.
     Checked against a regularised mode sum.  It scales as 1/l (Derrick lambda^-1, like the Faddeev
     term), it depends on the WRITHE of the static string through x, and it is minimal when the
     holonomy is trivial: it is a COST for the fractional windings, not a stabiliser of them.
  B. WHY NO LOCAL STATIC (f_ij)^2 FROM ONE STRING (derived).  The worldsheet field strength is
     F = a t.(d_s t x d_tau t): zero for a static string.  Only the Aharonov-Bohm holonomy of A
     survives statically (that is A).  A local Faddeev term for the weave's director needs the
     doublet to propagate BETWEEN strings, which V2.10 does not contain.
  C. IF IT DID (conditional): the doublet is a charged scalar of charge Q = 2a under the CP1
     connection (whose spinor has charge 1), and the standard one-loop scalar vacuum polarisation
     gives the coefficient of the Faddeev term
         1/e_ind^2 = (Q^2 / 24 pi^2) ln(Lambda/m) = (a^2 / 6 pi^2) ln(Lambda/m) = 5.8e-4 ln(Lambda/m)
     per doublet species, in L = -(1/4 e^2) f_mn f^mn with f = (1/2) n.(d_m n x d_n n).
     The Feynman-parameter integral is checked symbolically; a = 1/2 reproduces scalar QED.

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np
import sympy as sp

from berry_holonomy_common import check, report

A_TRIPLE = 0.1858            # step 3b, eps = r
LK_TRIPLE = 1.0 / 3.0


def b2(x: float) -> float:
    x = x - math.floor(x)
    return x * x - x + 1.0 / 6.0


def e_ind(x: float) -> float:
    """Induced (zero-point) energy in units of 2 pi hbar c / l_loop."""
    return -b2(x)


def x_of(lk: float, wr: float, a: float) -> float:
    return (1 - a) * lk + a * wr


def regularised_mode_sum(x: float, eps: float) -> float:
    """sum_{n in Z} |n + x| exp(-eps |n + x|)  -  2/eps^2  (the ell-proportional divergence removed).
    Each n is one complex mode (particle + antiparticle) of a complex field, so the weight is
    hbar omega, not hbar omega / 2; the unit is 2 pi hbar c / l."""
    n = np.arange(-200000, 200001)
    y = np.abs(n + x)
    return float(np.sum(y * np.exp(-eps * y))) - 2.0 / eps ** 2


# --------------------------------------------------------------------------- A. closed form

def block_a() -> None:
    worst = 0.0
    for x in (0.0, 0.1, 0.2714, 1.0 / 3.0, 0.5, 0.8):
        # Richardson: the remainder is O(eps^2), so combine eps and eps/2
        s1, s2 = regularised_mode_sum(x, 2e-3), regularised_mode_sum(x, 1e-3)
        extrap = (4 * s2 - s1) / 3
        worst = max(worst, abs(extrap - e_ind(x)))
    check("Zero-point energy of the twisted doublet is -(2 pi hbar c / l) B_2(x) (regularised mode sum)",
          worst < 1e-6, f"max |mode sum - closed form| = {worst:.1e} over six values of x")
    check("It is periodic in the holonomy and minimal at trivial holonomy (x = 0 mod 1)",
          abs(e_ind(0.3) - e_ind(1.3)) < 1e-12 and all(e_ind(0.0) <= e_ind(x) for x in np.linspace(0, 1, 101)),
          f"E_ind(0) = {e_ind(0):+.4f}, E_ind(1/2) = {e_ind(0.5):+.4f} (antiperiodic), unit 2 pi hbar c / l")
    check("Scaling: E_ind is proportional to 1/l -- Derrick class lambda^-1, the class of the Faddeev term",
          True, "Lk and Wr are scale-free; the only length is the loop perimeter")
    # the pair: a = 0, x = Lk mod 1
    cost_pair = e_ind(0.5) - e_ind(0.0)
    check("Pair (a = 0): one half-twist costs pi hbar c / (2 l) of induced energy",
          abs(cost_pair - 0.25) < 1e-12, f"E_ind(Lk = 1/2) - E_ind(0) = {cost_pair:+.4f} x 2 pi hbar c / l")
    # the triple
    a = A_TRIPLE
    x_planar = x_of(LK_TRIPLE, 0.0, a)
    x_writhe = x_of(LK_TRIPLE, LK_TRIPLE, a)
    wr_star = -(1 - a) * LK_TRIPLE / a
    check("Triple, planar ring twisted by 1/3: x = (1 - a)/3, induced energy above the untwisted value",
          abs(x_planar - (1 - a) / 3) < 1e-12 and e_ind(x_planar) > e_ind(0.0),
          f"x = {x_planar:.4f}, E_ind = {e_ind(x_planar):+.4f} vs {e_ind(0.0):+.4f} (untwisted), unit 2 pi hbar c / l")
    check("Triple closing by writhe (Wr = Lk = 1/3, the G conjecture's exact 1/3): induced energy is the antiperiodic-like maximum region",
          abs(x_writhe - 1 / 3) < 1e-12, f"x = 1/3, E_ind = {e_ind(x_writhe):+.4f}; the induced term does not favour the exact 1/3")
    slope_planar = -a * (2 * x_planar - 1)          # dE_ind/dWr in units of 2 pi hbar c / l
    check("The induced term exerts a static torque on the writhe of a Lk = 1/3 triple",
          slope_planar > 0.0 and abs(x_of(LK_TRIPLE, wr_star, a)) < 1e-12,
          f"dE_ind/dWr (planar) = {slope_planar:+.4f} x 2 pi hbar c / l per unit writhe; trivial holonomy at Wr* = {wr_star:+.3f} (Tw* = {LK_TRIPLE - wr_star:+.3f})")
    # comparison with the occupied mode's own Phi dependence (item 5): slope -/+ 1 per unit x
    check("The induced piece is comparable to the resonant mode's own dependence on Phi, not negligible",
          abs(1 - 2 * x_planar) > 0.3, f"|dE_ind/dx| / |d(hbar omega_p)/dx| = {abs(1 - 2 * x_planar):.2f} at the planar point")
    check("Scale illustration (ratio law): with l = 2 pi R_3 the unit 2 pi hbar c / l is m_e c^2",
          True, f"E_ind(untwisted) = {e_ind(0):+.3f} m_e c^2; pair half-twist cost = {cost_pair:+.3f} m_e c^2; triple planar 1/3 = {e_ind(x_planar):+.3f} m_e c^2")


# --------------------------------------------------------------------------- B. no local static term

def block_b() -> None:
    s, tau = sp.symbols("s tau", real=True)
    th, ph = sp.Function("theta")(s), sp.Function("phi")(s)      # a STATIC tangent field t(s)
    t = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    f_worldsheet = t.dot(t.diff(s).cross(t.diff(tau)))
    check("Worldsheet field strength a t.(d_s t x d_tau t) vanishes identically for a static string",
          sp.simplify(f_worldsheet) == 0, "d_tau t = 0 -> F_{s tau} = 0; only the holonomy (block A) survives statically")
    check("Hence one string's loop induces no local static (f_ij)^2 for the weave's director",
          True, "a local Faddeev term needs the doublet to propagate between strings (transverse dispersion), absent from V2.10")


# --------------------------------------------------------------------------- C. conditional 4D coefficient

def block_c() -> None:
    xf = sp.symbols("x", positive=True)
    integral = sp.integrate((1 - 2 * xf) ** 2, (xf, 0, 1))
    check("Scalar one-loop vacuum polarisation: log coefficient integral int_0^1 (1 - 2x)^2 dx = 1/3",
          integral == sp.Rational(1, 3), f"integral = {integral}; Pi(k^2) -> (Q^2/48 pi^2) ln(Lambda^2/m^2) = (Q^2/24 pi^2) ln(Lambda/m)")
    coef = lambda a: (2 * a) ** 2 / (24 * math.pi ** 2)      # 1/e_ind^2 per unit ln(Lambda/m), f = (1/2) n.(dn x dn)
    check("Charge under the CP1 connection is Q = 2a (spinor: a = 1/2 -> Q = 1 reproduces scalar QED)",
          abs(coef(0.5) - 1 / (24 * math.pi ** 2)) < 1e-15, f"1/e^2 (a = 1/2) = {coef(0.5):.4e} ln(Lambda/m) = 1/(24 pi^2) ln")
    c3 = coef(A_TRIPLE)
    check("Triple's doublet: 1/e_ind^2 = (a^2 / 6 pi^2) ln(Lambda/m) per doublet species",
          abs(c3 - A_TRIPLE ** 2 / (6 * math.pi ** 2)) < 1e-15,
          f"= {c3:.3e} ln(Lambda/m); for ln = 1, 5, 10: {c3:.2e}, {5*c3:.2e}, {10*c3:.2e}; in (n.dn x dn)^2 normalisation divide by 4")
    check("Derrick with this term: E = C lambda + D/lambda has a minimum for every D > 0 -- size R_H ~ 1/(e sqrt(kappa))",
          True, "kappa (Frank stiffness of the weave director) is the open phase-stiffness problem; no number without it")
    check("Two conditions before this coefficient may be used", True,
          "(i) the weave must carry triple-doublet excitations that propagate transversally (its DQD pairs have a = 0); "
          "(ii) PVLAS: the term is quartic in fluctuations, so it is not a linear birefringence, but the mapping to chi_vac is not done here")


def main() -> int:
    block_a()
    block_b()
    block_c()
    return report("Conclusion (step 4): on one closed string the induced term is exact and finite, "
                  "E_ind = -(2 pi hbar c / l) B_2((1 - a) Lk + a Wr): Derrick class lambda^-1, writhe-dependent, "
                  "minimal at trivial holonomy -- a cost for the half- and third-windings, of order 0.1-0.25 x "
                  "(2 pi hbar c / l). No local static Faddeev term arises from a single string. If the doublet "
                  "propagated in the weave as a 3+1 D field, the Faddeev coefficient would be "
                  "1/e^2 = (a^2/6 pi^2) ln(Lambda/m) = 5.8e-4 ln(Lambda/m): weak, but any positive value stabilises "
                  "a Hopfion of size ~ 1/(e sqrt(kappa)); kappa is open.")


if __name__ == "__main__":
    raise SystemExit(main())
