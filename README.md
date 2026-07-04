# Dipolar Strings (DS) Model — V2.3
Manuscript and supplementary material.
Archived on Zenodo: https://doi.org/10.5281/zenodo.21196099

## Contents
- `CD_V2.3_EN.tex/pdf` — manuscript
- `hydro_branch.py` — hydrostatic integration (Eq. 24), reproduces M_crit = 1.616e6 Msun
- `ppn_2pn.py` — symbolic 1PN/2PN verification (gamma = beta = 1)
- `make_figures.py`, `fig_heaviside.py` — figure generation

## Reproduce
pip install numpy scipy matplotlib sympy
python hydro_branch.py && python make_figures.py && python ppn_2pn.py
