"""V2.4 supplementary: prior sensitivity + analytic fraction + M_crit(l1) scaling scan."""
# (Consolidated from the analysis session; requires branch_data.npz from hydro_branch.py)
import numpy as np
from scipy.integrate import solve_ivp
G, c0, h, Msun = 6.67430e-11, 2.99792458e8, 6.62607015e-34, 1.98892e30

# --- Analytic fraction and prior sensitivity ---
d = np.load('branch_data.npz')
ucs, Ms, imax = d['ucs'], d['Ms'], int(d['imax'])
u, M = ucs[:imax+1], Ms[:imax+1]
logM = np.log10(M); logMc = logM[-1]
u_star = np.interp(logMc - 0.5, logM, u)
print("Analytic f_1/2 =", np.log(u[-1]/u_star)/np.log(u[-1]/u[0]))
N = 400000
for name, s in [
  ("log-uniform u_c", np.interp(np.random.uniform(np.log(u[0]), np.log(u[-1]), N), np.log(u), logM)),
  ("flat u_c",        np.interp(np.random.uniform(u[0], u[-1], N), u, logM)),
  ("flat M",          np.log10(np.random.uniform(M[0], M[-1], N))),
  ("flat log M",      np.random.uniform(logM[0], logMc, N))]:
    print(f"{name:16s}: f_1/2 = {np.mean(s > logMc - 0.5)*100:.0f}%")

# --- M_crit(l1) scaling scan ---
def Mcrit(l1):
    u0 = h*c0/l1**4
    def star(q):
        uc = q*u0; L = c0**2/np.sqrt(G*uc); r0 = 1e-8*L
        rhs = lambda r,y: [4*np.pi*r**2*y[1]/c0**2, -4*y[1]*G*y[0]/(r**2*c0**2)]
        ev = lambda r,y: y[1]-u0; ev.terminal=True; ev.direction=-1
        s = solve_ivp(rhs,[r0,1e6*L],[(4*np.pi/3)*r0**3*uc/c0**2,uc],events=ev,rtol=1e-8,atol=[1e-20,u0*1e-12])
        return s.y_events[0][0][0]/Msun if s.t_events[0].size else np.nan
    qs = np.logspace(np.log10(2), np.log10(60), 40)
    return np.nanmax([star(q) for q in qs])
for f in [0.5, 1.0, 1.63, 4.0]:
    Mc = Mcrit(f*8.09e-13)
    print(f"l1 x {f}: M_crit = {Mc:.3e} Msun ; M_crit/f^2 = {Mc/f**2:.3e}")
