# Dipolar Strings (DS) Model

Manuscript and supplementary material for *Dipolar String Model: an Exploratory Framework for the Vacuum from Transmission-Line Electrodynamics* (R. St-Pierre, 2026, V2.6).

Archived on Zenodo:
- Concept DOI (always latest version): https://doi.org/10.5281/zenodo.21196098
- v1.0 version DOI (time-stamped record of the registered prediction, Section 15.2): https://doi.org/10.5281/zenodo.21196099

## Contents

**Manuscript**
- `DS_model_V2_6.tex` / `DS_model_V2_6_final.pdf` — consolidated manuscript (V1.8 → V2.5 review trajectory documented in Appendix A)

**Numerical integration (Section 9)**
- `hydro_branch.py` — hydrostatic integration of the structure system (Eq. 39), SciPy adaptive RK45 with event-based surface detection; reproduces M_crit = 1.616e6 Msun at u_c/u_0 = 14.0, r_s/R = 1.22; writes `branch_data.npz`
- `hydro_branch_check.py` — independent re-implementation: dimensionless reduction (hand-coded fixed-step RK4, linear surface interpolation); confirms M_crit, log M_crit = 6.21, median log M = 6.10, and M_crit ∝ ℓ₁² analytically
- `supplementary_v24.py` — analytic accumulation fraction f_1/2 (Eq. 65), prior-sensitivity analysis (log-uniform / flat u_c / flat M / flat log M), and ℓ₁-scaling scan; requires `branch_data.npz`

**Symbolic verification (Sections 7, 8, 11)**
- `theorem3_covariance.py` — Theorem 3: L(Z) identically covariant on the effective metric (Eq. 26–27), null residual
- `conformal_identity.py` — conformal clause g_iso = e^χ g_opt and factor-of-2 resolution (Section 8.2, Eq. 31)
- `ppn_2pn.py` — SymPy 1PN/2PN expansions of the isotropic exponential metric vs. Schwarzschild (γ = β = 1; Δg₀₀ = −U³/6, Δg_ij = +U²δ_ij/2; Section 8.3, Eq. 32–34)
- `boosted_residual.py` — Heaviside ellipsoid as exact boosted static solution (Section 11, null wave residual)

**Full verification suite**
- `verify_all.py` — 40 automated checks spanning Sections 5–15 (Ze, D/r = 2cosh π, Theorem 4 / Johns-node uniqueness and parity, Theorems 1–3, PPN/2PN, photon sphere, Bell bound, M_crit, caustic fraction)

**Figures**
- `make_figures.py`, `fig_heaviside.py`, `fig_scaling.py` — figure generation (Figs. 2–5 of the manuscript; `make_figures.py` requires `branch_data.npz`)

## Reproduce

```
pip install numpy scipy matplotlib sympy
python hydro_branch.py        # regenerates branch_data.npz and the M(u_c) sequence
python hydro_branch_check.py  # independent dimensionless re-integration
python supplementary_v24.py   # f_1/2, prior sensitivity, l1-scaling scan
python verify_all.py          # full suite: 40/40 checks
python make_figures.py        # regenerates the figures
```

All numerical results are reproducible from the structure equations (Eq. 39 of the manuscript), which fully specify the integration.

## Status

Exploratory working hypothesis — not a validated contribution. Both numerical implementations originate from the same collaboration (human–AI, disclosed in the Acknowledgments); independent third-party verification has not yet been performed and is invited. The registered prediction and its frozen statistical criteria are those of the v1.0 Zenodo deposit above. Open problems are inventoried in Section 16 of the manuscript.

## License

Manuscript: CC-BY 4.0. Code: MIT.
