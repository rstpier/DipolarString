"""
Symbolic verification of Theorem 3 (exact covariance of the scalar sector).

Claim (paper Eq. 26-27): with chi = ln(Z/Z0), the Lagrangian density

    L_DQD = -(c0^4 / 32 pi G) [ (grad Z)^2 / Z^2  -  (dZ/dt)^2 / (Z0^2 c0^2) ]

is IDENTICALLY equal to the covariant scalar-field form

    L_cov = -(c0^4 / 32 pi G) (e^chi / c0) sqrt(-g) g^{mu nu} d_mu chi d_nu chi,

with g = diag(-c_eff^2, 1, 1, 1) and c_eff = c0 e^{-chi} (optical representative).

The script computes both sides for a generic field Z(t,x,y,z) and asserts a null
symbolic residual.
"""
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
c0, Z0, G = sp.symbols('c0 Z0 G', positive=True)
Z = sp.Function('Z', positive=True)(t, x, y, z)

pref = -c0**4 / (32*sp.pi*G)

# --- Left-hand side: the Lagrangian as written in Eq. (26) ---
gradZ2 = sp.diff(Z, x)**2 + sp.diff(Z, y)**2 + sp.diff(Z, z)**2
L_dqd = pref * (gradZ2 / Z**2 - sp.diff(Z, t)**2 / (Z0**2 * c0**2))

# --- Right-hand side: covariant form on the optical effective metric ---
chi   = sp.log(Z / Z0)
c_eff = c0 * sp.exp(-chi)

# g_{mu nu} = diag(-c_eff^2, 1, 1, 1)  =>  sqrt(-det g) = c_eff
sqrt_mg = c_eff
# g^{mu nu} d_mu chi d_nu chi = -(d_t chi)^2 / c_eff^2 + (grad chi)^2
gradchi2 = sp.diff(chi, x)**2 + sp.diff(chi, y)**2 + sp.diff(chi, z)**2
kin = -sp.diff(chi, t)**2 / c_eff**2 + gradchi2
L_cov = pref * (sp.exp(chi) / c0) * sqrt_mg * kin

# --- Residual ---
residual = sp.simplify(L_dqd - L_cov)
print("Residual  L_DQD - L_cov  =", residual)
assert residual == 0, "Theorem 3 identity FAILED"
print("Theorem 3 verified: the scalar sector is identically covariant "
      "on the effective metric (null symbolic residual).")
