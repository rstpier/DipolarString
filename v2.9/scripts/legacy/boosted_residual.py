"""
Symbolic verification of the boosted-solution residual (paper Sections 8 (Theorem 3) and 11).

Claim: the convective (Heaviside-ellipsoid) profile

    chi_boost(t,x,y,z) = k / sqrt( gamma^2 (x - v t)^2 + y^2 + z^2 ),
    k = 2GM/c0^2,   gamma^2 = 1 / (1 - v^2/c0^2),

is the exact boosted image of the static harmonic solution chi = k/r, and the
residual of the wave operator of the linearized scalar sector,

    Box chi = laplacian(chi) - (1/c0^2) d^2 chi / dt^2,

applied to it is IDENTICALLY null away from the source (Figure: 'Boosted solution
of the scalar sector'). This is the statement that the Heaviside ellipsoid is not
an additional postulate but a consequence of the covariance of the scalar sector
(Theorem 3).
"""
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
G, M, c0 = sp.symbols('G M c0', positive=True)
v = sp.symbols('v', positive=True)

k = 2*G*M/c0**2
gamma2 = 1/(1 - v**2/c0**2)

# --- Static solution: harmonic, laplacian = 0 away from r = 0 ---
r = sp.sqrt(x**2 + y**2 + z**2)
chi_static = k/r
lap_static = sum(sp.diff(chi_static, s, 2) for s in (x, y, z))
lap_static = sp.simplify(lap_static)
print("laplacian(chi_static) =", lap_static)
assert lap_static == 0, "Static solution is not harmonic"

# --- Boosted profile: x -> gamma (x - v t) ---
rstar = sp.sqrt(gamma2*(x - v*t)**2 + y**2 + z**2)
chi_boost = k/rstar

box = (sum(sp.diff(chi_boost, s, 2) for s in (x, y, z))
       - sp.diff(chi_boost, t, 2)/c0**2)
residual = sp.simplify(box)
print("Box(chi_boost) =", residual)
assert residual == 0, "Boosted-solution wave residual is NOT null"
print("Verified: the Heaviside-ellipsoid profile is an exact solution of the "
      "wave operator (null symbolic residual) -- the boosted image of the "
      "static solution, as required by Theorem 3.")
