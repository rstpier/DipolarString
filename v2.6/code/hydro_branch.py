"""
Independent reproduction of the substrate-star hydrostatic branch (DS model, Eq. 39).

System:
    dM/dr   = 4 pi r^2 u / c0^2
    dchi/dr = 2 G M / (r^2 c0^2)
    dP/dr   = -(u+P)/c0^2 * G M / r^2 ,   P = u/3  =>  du/dr = -4 u G M / (r^2 c0^2)
Surface: P(R) = u0/3  i.e.  u(R) = u0, with u0 = h c0 / l1^4.
Scan central densities u_c/u0, find M(u_c), compactness r_s/R = 2GM/(R c0^2),
stability by first-maximum criterion.
"""
import numpy as np
from scipy.integrate import solve_ivp

# --- Constants (SI) ---
G   = 6.67430e-11
c0  = 2.99792458e8
h   = 6.62607015e-34
Msun = 1.98892e30
l1  = 8.09e-13                    # nominal string length (m)

u0 = h * c0 / l1**4               # native vacuum energy density
print(f"u0 = {u0:.4e} J/m^3   (paper: ~4.6e23)")

def integrate_star(uc_ratio, r_max_factor=1e6, rtol=1e-9):
    """Integrate from center outward; return (M_solar, R, rs_over_R) or None."""
    uc = uc_ratio * u0
    # Characteristic length scale of the problem:
    L_star = c0**2 / np.sqrt(G * uc)          # ~ gravitational length for density uc
    r0 = 1e-8 * L_star                         # tiny seed radius (series start)
    # Series start: M0 = (4pi/3) r0^3 uc / c0^2 ; u ~ uc
    M0 = (4*np.pi/3) * r0**3 * uc / c0**2
    def rhs(r, y):
        M, u = y
        dM = 4*np.pi * r**2 * u / c0**2
        du = -4.0 * u * G * M / (r**2 * c0**2)
        return [dM, du]
    def hit_surface(r, y):
        return y[1] - u0
    hit_surface.terminal = True
    hit_surface.direction = -1
    sol = solve_ivp(rhs, [r0, r_max_factor * L_star], [M0, uc],
                    events=hit_surface, rtol=rtol, atol=[1e-20, u0*1e-14],
                    method='RK45', dense_output=False, max_step=np.inf)
    if sol.t_events[0].size == 0:
        return None
    R = sol.t_events[0][0]
    M = sol.y_events[0][0][0]
    rs = 2*G*M/c0**2
    return M/Msun, R, rs/R

# --- Scan the sequence of equilibria ---
uc_ratios = np.logspace(np.log10(1.01), 6, 400)
Ms, Rs, comps, ucs = [], [], [], []
for q in uc_ratios:
    out = integrate_star(q)
    if out is None:
        continue
    M, R, c = out
    ucs.append(q); Ms.append(M); Rs.append(R); comps.append(c)

ucs = np.array(ucs); Ms = np.array(Ms); Rs = np.array(Rs); comps = np.array(comps)

# --- First maximum (M_crit) ---
imax = np.argmax(Ms)
print(f"\nM_crit = {Ms[imax]:.3e} Msun at u_c/u0 = {ucs[imax]:.2f}, r_s/R = {comps[imax]:.3f}")

# Compare with paper checkpoints
for target in [1.05, 1.9, 3.3, 5.9, 14.0]:
    i = np.argmin(np.abs(ucs - target))
    print(f"u_c/u0 = {ucs[i]:7.2f} : M = {Ms[i]:.3e} Msun , r_s/R = {comps[i]:.3f}")

# Attractor at high density
i_hi = np.argmin(np.abs(ucs - 1e4))
print(f"High-density attractor (u_c/u0 ~ 1e4): M = {Ms[i_hi]:.3e} Msun, r_s/R = {comps[i_hi]:.3f}")

np.savez('branch_data.npz', ucs=ucs, Ms=Ms, Rs=Rs, comps=comps,
         imax=imax, u0=u0, l1=l1)
print("\nSaved branch_data.npz")
