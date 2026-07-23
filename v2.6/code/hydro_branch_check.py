"""
Second, independent cross-check of the substrate-star hydrostatic branch (DS model).

This code is deliberately built on a DIFFERENT footing from hydro_branch.py, so that
agreement between the two is a genuine cross-check rather than a re-run:

  * Dimensionless formulation. The structure system (paper Eq. 39) is reduced to a
    parameter-free ODE by the scaling
        r = a x,   M = b m,   u = u0 theta,
        a = c0^2 / (4 sqrt(pi G u0)),   b = a c0^2 / (4 G) = c0^4 / (16 G sqrt(pi G u0)).
    The single physical mass scale b (propto c0^4 / sqrt(G^3 u0) propto l1^2) is thus
    ANALYTIC, not obtained from a numerical l1-scan -- an independent confirmation of
    M_crit propto l1^2.

  * Reduced variables (m, theta) instead of (M, u).
        dm/dx     =  x^2 theta
        dtheta/dx = -theta m / x^2
    Compactness follows in closed form:  r_s/R = 2GM/(R c0^2) = m_surf / (2 x_surf).

  * Hand-coded fixed-step classical RK4 (no SciPy solver), surface located by linear
    interpolation between the two bracketing steps where theta crosses 1.

If this reproduces M_crit, the critical central density, the compactness, the branch
checkpoints and the LRD distribution statistics obtained by hydro_branch.py, the two
independent implementations agree.

NOTE (methodological): both implementations originate from the same human-AI
collaboration; this is a cross-check, not independent third-party verification.
"""
import numpy as np

# --- Constants (SI) ---
G    = 6.67430e-11
c0   = 2.99792458e8
h    = 6.62607015e-34
Msun = 1.98892e30
l1   = 8.09e-13                      # nominal string length (m)

u0 = h * c0 / l1**4                  # native vacuum energy density
print(f"u0 = {u0:.4e} J/m^3   (paper: ~4.6e23)")

# --- Dimensionless scales ---
a = c0**2 / (4.0 * np.sqrt(np.pi * G * u0))   # length scale (m)
b = a * c0**2 / (4.0 * G)                      # mass scale (kg)
print(f"mass scale b = {b/Msun:.4e} Msun   (propto l1^2, analytic)")


def integrate_star(theta_c, dx=2.0e-4, x_max=60.0):
    """Hand-coded fixed-step RK4 on the reduced (m, theta) system.
    Returns (M_solar, R_m, rs_over_R) or None if no surface (theta=1) is reached."""

    def rhs(x, y):
        m, th = y
        return np.array([x * x * th, -th * m / (x * x)])

    # Regular series start away from x=0:  m ~ (theta_c/3) x^3,  theta ~ theta_c - (theta_c^2/6) x^2
    x = 1.0e-3
    m = (theta_c / 3.0) * x**3
    th = theta_c - (theta_c**2 / 6.0) * x**2
    y = np.array([m, th])

    while x < x_max:
        th_prev, x_prev, y_prev = y[1], x, y.copy()
        k1 = rhs(x, y)
        k2 = rhs(x + 0.5 * dx, y + 0.5 * dx * k1)
        k3 = rhs(x + 0.5 * dx, y + 0.5 * dx * k2)
        k4 = rhs(x + dx, y + dx * k3)
        y = y + (dx / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        x = x + dx
        if y[1] <= 1.0:                       # theta crossed the surface value
            # linear interpolation for the crossing point
            f = (th_prev - 1.0) / (th_prev - y[1])
            x_s = x_prev + f * dx
            m_s = y_prev[0] + f * (y[0] - y_prev[0])
            M = b * m_s
            R = a * x_s
            rs_over_R = m_s / (2.0 * x_s)     # = 2GM/(R c0^2), closed form
            return M / Msun, R, rs_over_R
    return None


# --- Scan the sequence of equilibria (theta_c = u_c/u0) ---
theta_grid = np.logspace(np.log10(1.01), 6, 400)
Ms, Rs, comps, ucs = [], [], [], []
for tc in theta_grid:
    out = integrate_star(tc)
    if out is None:
        continue
    M, R, comp = out
    ucs.append(tc); Ms.append(M); Rs.append(R); comps.append(comp)

ucs = np.array(ucs); Ms = np.array(Ms); Rs = np.array(Rs); comps = np.array(comps)

imax = int(np.argmax(Ms))
print(f"\nM_crit = {Ms[imax]:.3e} Msun at u_c/u0 = {ucs[imax]:.2f}, r_s/R = {comps[imax]:.3f}")
print(f"log10 M_crit = {np.log10(Ms[imax]):.3f}  (paper S15.2: 6.21)")

print("\nBranch checkpoints (compare paper Table, Sec. 8.2):")
for target in [1.05, 1.9, 3.3, 5.9, 14.0]:
    i = np.argmin(np.abs(ucs - target))
    print(f"  u_c/u0 = {ucs[i]:7.2f} : M = {Ms[i]:.3e} Msun , r_s/R = {comps[i]:.3f}")

i_hi = np.argmin(np.abs(ucs - 1e4))
print(f"High-density attractor (u_c/u0 ~ 1e4): M = {Ms[i_hi]:.3e} Msun, r_s/R = {comps[i_hi]:.3f}")

# --- LRD distribution statistics from THIS code's own branch ---
stable = slice(0, imax + 1)
logM = np.log10(Ms[stable])
order = np.argsort(logM); v = logM[order]
w = np.ones_like(v); w = w / w.sum()
cw = np.cumsum(w) - 0.5 * w
med = np.interp(0.50, cw, v)
q10 = np.interp(0.10, cw, v)
logMc = np.log10(Ms[imax])
f_half = np.mean((logM >= logMc - 0.5) & (logM <= logMc))
print("\nLRD distribution statistics (log-uniform prior in u_c):")
print(f"  predicted median log M = {med:.4f}  (paper S15.2: 6.10 [6.096])")
print(f"  10% quantile   log M = {q10:.4f}  (paper: 5.32)")
print(f"  half-decade fraction  = {f_half:.3f}  (paper: ~0.80)")
