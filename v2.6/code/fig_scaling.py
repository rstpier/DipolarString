"""Figure: M_crit vs l1 quadratic scaling (paper fig:scaling). Recomputes the 4-point scan."""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 11, 'figure.dpi': 150})

G, c0, h, Msun = 6.67430e-11, 2.99792458e8, 6.62607015e-34, 1.98892e30
l1_nom = 8.09e-13

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

factors = np.array([0.5, 1.0, 1.63, 4.0])
Mc = np.array([Mcrit(f*l1_nom) for f in factors])
for f, m in zip(factors, Mc):
    print(f"l1 x {f}: M_crit = {m:.3e} Msun ; M_crit/f^2 = {m/f**2:.3e}")

fig, ax = plt.subplots(figsize=(6.4, 4.8))
ff = np.logspace(np.log10(0.4), np.log10(5), 100)
ax.plot(ff, Mc[1]*ff**2, 'k--', lw=1.4, label=r'$M_{\rm crit} \propto \ell_1^2$')
ax.plot(factors, Mc, 'o', ms=8, color='#1a5fb4', label='Full numerical integrations')
ax.axhline(4.297e6, color='#e66100', lw=1.2, ls=':',)
ax.annotate(r'Sgr A$^*$ (GRAVITY): $4.30\times10^6\,M_\odot$' + '\n' + r'$\Leftrightarrow \ell_1\times1.63$',
            xy=(1.63, Mc[2]), xytext=(0.55, 8e6),
            arrowprops=dict(arrowstyle='->', lw=1), fontsize=9.5)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$\ell_1 / \ell_{1,\rm nom}$')
ax.set_ylabel(r'$M_{\rm crit}\ \ [M_\odot]$')
ax.set_title(r'Quadratic scaling of the critical mass (Eq. 39)')
ax.legend(loc='upper left', fontsize=9.5)
ax.grid(alpha=0.25, which='both')
fig.tight_layout()
fig.savefig('fig_scaling.pdf'); fig.savefig('fig_scaling.png')
print("Scaling figure saved.")
