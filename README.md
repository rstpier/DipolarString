# Dipolar Strings (DS) Model

Manuscript and supplementary material for *The Universe as an RLC Network: Analogue Gravity and Quantum Coherence via Impedance-Matched Transmission Lines — The Dipolar Strings Model* (R. St-Pierre, 2026).

Archived on Zenodo (time-stamped record of the registered prediction, Section 14.1):
https://doi.org/10.5281/zenodo.21196099

## Contents

**Manuscript**
- `CD_journal_EN.tex` / `CD_journal_EN.pdf` — consolidated manuscript (V1.8 → V2.4 review trajectory documented in Appendix A)

**Numerical integration (Section 8)**
- `supplementary/hydro_branch.py` — hydrostatic integration of the structure system (Eq. 32), SciPy adaptive RK45 with event-based surface detection; reproduces M_crit = 1.616e6 Msun at u_c/u_0 = 14.0, r_s/R = 1.22; writes `branch_data.npz`
- `supplementary/branch_data.npz` — precomputed equilibrium-sequence data
- `supplementary/supplementary_v24.py` — analytic accumulation fraction f_1/2, prior-sensitivity analysis (log-uniform / flat u_c / flat M / flat log M), and M_crit ∝ ℓ₁² scaling scan; requires `branch_data.npz`

**Symbolic verification (Section 7)**
- `supplementary/ppn_2pn.py` — SymPy 1PN/2PN expansions of the isotropic exponential metric vs. Schwarzschild (γ = β = 1; Δg₀₀ = −U³/6, Δg_ij = +U²δ_ij/2)

**Figures**
- `supplementary/make_figures.py`, `supplementary/fig_heaviside.py` — figure generation from `branch_data.npz`
- `fig_branch.pdf`, `fig_lrd.pdf`, `fig_scaling.pdf`, `fig_heaviside.pdf` — precompiled figures (Figs. 2–5 of the manuscript)

## Reproduce

```
pip install numpy scipy matplotlib sympy
cd supplementary
python hydro_branch.py        # regenerates branch_data.npz and the M(u_c) sequence
python supplementary_v24.py   # f_1/2, prior sensitivity, l1-scaling scan
python ppn_2pn.py             # symbolic 1PN/2PN verification
python make_figures.py        # regenerates the figures
```

All numerical results are reproducible from the structure equations (Eq. 32 of the manuscript), which fully specify the integration.

## Status

Exploratory working hypothesis — not a validated contribution. Both numerical implementations originate from the same collaboration (human–AI, disclosed in the Acknowledgments); independent third-party verification has not yet been performed and is invited. The registered prediction and its frozen statistical criteria are those of the Zenodo deposit above.

## License

Manuscript: CC-BY 4.0. Code: MIT (or adapt as preferred).
