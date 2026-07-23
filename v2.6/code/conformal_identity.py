"""
Symbolic verification of the conformal-representative clause (paper Section 8.2).

Two claims are verified:

(1) Matrix identity: the isotropic representative
        ds^2 = -c0^2 e^{-chi} dt^2 + e^{+chi} (dr^2 + r^2 dOmega^2)
    is conformally related to the optical representative
        g_opt = diag(-c_eff^2, 1, 1, 1),  c_eff = c0 e^{-chi},
    by  g_iso = e^{chi} g_opt  (exact matrix identity).

(2) Factor-of-2 resolution: the static geodesic radial acceleration
        d^2 r / dtau^2 = -Gamma^r_tt (dt/dtau)^2,
    with the exterior solution chi = 2GM/(r c0^2), equals EXACTLY
        -2GM/r^2            in the OPTICAL representative (twice Newton), and
        -GM/r^2 * e^{-chi}  in the ISOTROPIC representative (Newton, weak field),
    reproducing the paper's Section 8.2 verbatim. Matter must therefore couple
    to the isotropic representative; the static geodesic of (31) coincides with
    the force law (30).
"""
import sympy as sp

r, G, M, c0 = sp.symbols('r G M c0', positive=True)
chi_s = sp.Function('chi', positive=True)(r)      # generic profile for (1)
chi_e = 2*G*M/(r*c0**2)                           # exterior solution for (2)

# ---------- (1) Matrix identity g_iso = e^chi g_opt ----------
c_eff = c0*sp.exp(-chi_s)
g_opt = sp.diag(-c_eff**2, 1, 1, 1)
g_iso = sp.diag(-c0**2*sp.exp(-chi_s), sp.exp(chi_s), sp.exp(chi_s), sp.exp(chi_s))
res1 = sp.simplify(g_iso - sp.exp(chi_s)*g_opt)
print("Matrix residual  g_iso - e^chi g_opt  =", res1.is_zero_matrix)
assert res1 == sp.zeros(4, 4), "Conformal matrix identity FAILED"
print("(1) verified: g_iso = e^chi g_opt exactly.\n")

# ---------- (2) Static geodesic radial acceleration ----------
def static_geodesic_acc(g_tt, g_rr):
    """d^2 r/dtau^2 = -Gamma^r_tt (dt/dtau)^2 for a momentarily static particle.
    Gamma^r_tt = -(1/2) g^rr d_r g_tt ;  (dt/dtau)^2 = -c0^2/g_tt."""
    Gamma_r_tt = -sp.Rational(1, 2)*(1/g_rr)*sp.diff(g_tt, r)
    dt_dtau2 = -c0**2/g_tt
    return sp.simplify(-Gamma_r_tt*dt_dtau2)

# Optical representative: g_tt = -c0^2 e^{-2chi}, g_rr = 1
a_opt = static_geodesic_acc(-c0**2*sp.exp(-2*chi_e), sp.Integer(1))
# Isotropic representative: g_tt = -c0^2 e^{-chi}, g_rr = e^{+chi}
a_iso = static_geodesic_acc(-c0**2*sp.exp(-chi_e), sp.exp(chi_e))

print("a_optical   =", a_opt)          # expected: -2*G*M/r**2 exactly
print("a_isotropic =", a_iso)          # expected: -G*M*exp(-chi)/r**2 exactly

assert sp.simplify(a_opt + 2*G*M/r**2) == 0, \
    "Optical representative does not give exactly 2GM/r^2"
assert sp.simplify(a_iso + (G*M/r**2)*sp.exp(-chi_e)) == 0, \
    "Isotropic representative does not give exactly GM/r^2 e^{-chi}"
ratio = sp.simplify(sp.series(a_opt/a_iso, G, 0, 1).removeO())
print("weak-field ratio optical/isotropic =", ratio)
assert ratio == 2, "Factor-of-2 between representatives not reproduced"
print("(2) verified: optical gives exactly -2GM/r^2 (twice Newton); "
      "isotropic gives exactly -GM/r^2 e^{-chi} (Newton in the weak field); "
      "ratio = 2 at leading order -- Section 8.2 reproduced verbatim.")
