"""Figures for V2.3: (1) equilibrium branch M(u_c); (2) predicted log M distribution vs LRDs."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'figure.dpi': 150})

d = np.load('branch_data.npz')
ucs, Ms, comps, imax = d['ucs'], d['Ms'], d['comps'], int(d['imax'])

# ============ Statistics of the predicted distribution ============
# Log-uniform prior in central density over the stable branch [~1, u_crit]
mask_stable = np.arange(len(ucs)) <= imax
lu = np.log(ucs[mask_stable]); logM = np.log10(Ms[mask_stable])
# Sample log-uniform in u_c:
n_samp = 400000
lu_s = np.random.uniform(lu.min(), lu.max(), n_samp)
logM_s = np.interp(lu_s, lu, logM)
med = np.median(logM_s); q10 = np.quantile(logM_s, 0.10)
frac_half_decade = np.mean((logM_s > np.log10(Ms[imax]) - 0.5) & (logM_s <= np.log10(Ms[imax])))
print(f"Predicted median log M = {med:.2f}  (paper: 6.10)")
print(f"10% quantile log M    = {q10:.2f}  (paper: 5.32)")
print(f"Fraction in half-decade below M_crit = {frac_half_decade*100:.0f}%  (paper: ~80%)")
print(f"log M_crit = {np.log10(Ms[imax]):.2f} (paper: 6.21)")

# ============ FIGURE 1 : Equilibrium branch ============
fig, ax = plt.subplots(figsize=(7.0, 5.0))
ax.plot(ucs[mask_stable], Ms[mask_stable], color='#1a5fb4', lw=2.2, label='Stable branch')
ax.plot(ucs[~mask_stable], Ms[~mask_stable], color='#c01c28', lw=2.2, ls='--',
        label=r'Unstable ($\to$ frozen collapse)')
ax.axhline(Ms[imax], color='gray', lw=0.8, ls=':')
ax.plot(ucs[imax], Ms[imax], 'k*', ms=16, zorder=5)
ax.annotate(r'$M_{\rm crit} \approx 1.62\times10^{6}\,M_\odot$' + '\n' +
            r'($u_c/u_0 \approx 14$, $r_s/R \approx 1.22$)',
            xy=(ucs[imax], Ms[imax]), xytext=(30, 6.0e5),
            arrowprops=dict(arrowstyle='->', color='k', lw=1),
            fontsize=10.5)
# Shade the "dark" region r_s/R > 1 along the stable branch
i_dark = np.where(comps[mask_stable] >= 1.0)[0]
if i_dark.size:
    ax.axvspan(ucs[mask_stable][i_dark[0]], ucs[imax], color='k', alpha=0.10,
               label=r'Dark, horizonless ($r_s/R>1$)')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'Central density $u_c/u_0$')
ax.set_ylabel(r'Total mass $M\ \ [M_\odot]$')
ax.set_title('Substrate-star equilibrium sequence (independent reproduction, Eq. 39)')
ax.legend(loc='lower right', fontsize=9.5, framealpha=0.9)
ax.grid(alpha=0.25, which='both')
fig.tight_layout()
fig.savefig('fig_branch.pdf'); fig.savefig('fig_branch.png')

# ============ FIGURE 2 : Predicted mass distribution vs LRDs ============
fig2, ax2 = plt.subplots(figsize=(7.0, 5.0))
bins = np.arange(4.2, 6.5, 0.1)
ax2.hist(logM_s, bins=bins, density=True, color='#1a5fb4', alpha=0.55,
         label='DS prediction (stable branch,\nlog-uniform prior in $u_c$)')
ax2.axvline(np.log10(Ms[imax]), color='k', lw=1.5, label=r'$\log M_{\rm crit} = 6.21$ (cutoff)')
ax2.axvline(med, color='#1a5fb4', lw=1.5, ls='--', label=f'Predicted median = {med:.2f}')
# Observations (Rusakov et al. 2026, electron-scattering-corrected):
ax2.axvline(6.10, color='#e66100', lw=2.2, label='LRD median, 7 firm masses = 6.10\n(Rusakov et al. 2026)')
ax2.axvspan(5.3, 5.9, color='#e66100', alpha=0.15,
            label=r'18 fainter LRDs: median $\log M \lesssim 5.6$')
# Registered prediction window
ax2.axvspan(np.log10(Ms[imax])-0.5, np.log10(Ms[imax]), ymin=0, ymax=1,
            facecolor='none', edgecolor='k', hatch='///', lw=0.0, alpha=0.35,
            label=f'Registered prediction: ~{frac_half_decade*100:.0f}% of objects\nin half-decade below cutoff')
ax2.set_xlabel(r'$\log_{10}(M/M_\odot)$')
ax2.set_ylabel('Probability density')
ax2.set_title('Predicted substrate-star mass function vs. Little Red Dot masses')
ax2.legend(loc='upper left', fontsize=8.4, framealpha=0.95)
ax2.set_xlim(4.2, 6.6)
ax2.grid(alpha=0.25)
fig2.tight_layout()
fig2.savefig('fig_lrd.pdf'); fig2.savefig('fig_lrd.png')
print("Figures saved.")
