# Dipolar String Model — V2.11 (erratum edition, incremental package)

Author: R. St-Pierre
Edition date: 21 November 2026 (package frozen 22 November 2026)
Manuscript version: **V2.11 (erratum edition)**
Zenodo: **not deposited** (author's decision, 22 November 2026). The Zenodo record remains V2.10 (concept DOI
10.5281/zenodo.21196098); this edition is distributed through the repository only.
Supersedes: **V2.10**, which remains archived unchanged in [`../v2.10/`](../v2.10/)

## What this package is

V2.11 is the frozen V2.10 manuscript with an **erratum** applied in place. Twelve items, every one
marked *Erratum E*n in the text and listed in the manuscript's Document History; all of them reduce the
model's claims, and two of them (E6, E7) perform computations V2.10 had listed as not done — both come
out against the model. The itemized erratum, with the V2.10 statement, the correction and the evidence
for each item, is [`ERRATUM.md`](ERRATUM.md).

In one paragraph: charge by chirality and mass by inductance are withdrawn (E1); the electron as an
extended soliton is withdrawn (E2); the proton-impedance item is blocked in principle (E3); Γ_pole = 1/3
is a reflection coefficient, its identification with the charge is an identification (E4); A5 is
reworded (E5); the preferred-frame parameters of the scalar sector are computed, α₁ = −8, excluded by
binary pulsars by five orders of magnitude (E6); with the general-relativistic factors restored the
substrate-star branch is not dark (M_crit = 3.9×10⁵ M☉, r_s/R ≤ 0.31), the frozen collapses lose their
mechanism, and falsification items 1–3 and the registered little-red-dot prediction are withdrawn (E7);
under the amplitude axiom adopted by the author the photon carries the lattice dispersion and LHAASO
excludes ℓ₁ as the lattice spacing (E8); a spin-½ excitation of the cubic weave travels at 2c₀/3 (E9);
the zero-point energy of a lattice at ℓ₁ is 10³² times the dark-energy density (E10); the status
ledger is re-filed (E11); the author's adoption of the amplitude axiom A9 is recorded (E12).

The companion note **`negative_results/DS_negative_results.pdf`** (6 pages) collects the fourteen
structural negative results established by the programme, each with its proof sketch, its classical source where one exists, and the script
that reproduces it. They hold for the whole class of "vacuum as a network of matched lines" models.

## Verify

```bash
python verification/verify_v2_11.py     # exit 0
```

The script checks that each erratum item is present in the V2.11 source (and that the V2.10 hygiene
checks still pass), that the companion note carries its fourteen results, and re-runs the seven
evidence scripts of the erratum (R119–R125, in `../matter/redraw/` and `../matter/rebuild/`), requiring
a `RESULT` line with zero failures from each. The symbolic script of E6 needs `sympy`; the V2.9.2
numerical suite (unchanged, validating the equations as integrated in V2.10) is run when `numpy`,
`scipy` and `matplotlib` are available and skipped otherwise:

```bash
python -m pip install sympy             # for E6
python -m pip install -r ../v2.9/requirements.txt   # optional, for the delegated suite
```

A passing run means the listed consequences follow from the encoded equations. It is not third-party
validation, and it is not independent verification: the scripts and this edition were produced by the
same human–AI collaboration as the manuscript.

## Build

```bash
cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error DS_model_V2_11_erratum.tex
cd ../negative_results && latexmk -pdf -interaction=nonstopmode -halt-on-error DS_negative_results.tex
```

71 and 6 pages, zero undefined references. `figures/` holds the four figures the manuscript includes
(unchanged from V2.10; Figs. branch and lrd are superseded by E7 and kept as the V2.10 record).

## Contents

```
ERRATUM.md          the itemized erratum (E1–E12), with evidence
manuscript/         V2.11 LaTeX source and PDF
figures/            the four figures included by the manuscript (unchanged)
negative_results/   companion note, LaTeX source and PDF
verification/       verify_v2_11.py
```

`FILE_SHA256SUMS.txt` and `PACKAGE_CONTENTS.txt` cover this directory only.

## Licences

Manuscript, note and figures: CC BY 4.0. Code: MIT. See `../v2.9/LICENSE.md`.
