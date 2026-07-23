# Dipolar Strings (DS / CD) model

Reproducible material for *Dipolar String Model: an Exploratory Framework for the Vacuum from
Transmission-Line Electrodynamics* (R. St-Pierre, 2026) — a deterministic framework in which the
physical vacuum is modelled as a gas of dipolar cells behaving as bifilar transmission lines, and
gravity as refraction of an impedance field.

**Status: exploratory working paper. No human peer review has taken place.** Both numerical
implementations originate from the same human–AI collaboration, disclosed in the manuscript's
Acknowledgments; independent third-party verification is invited and has not yet been performed.

## Current version — V2.7

| | |
|---|---|
| Manuscript | [`v2.7/DS_model_V2_7_final.pdf`](v2.7/DS_model_V2_7_final.pdf) (54 pp.) · [LaTeX source](v2.7/DS_model_V2_7_final.tex) |
| Package | [`v2.7/DS_V2.7_package.zip`](v2.7/DS_V2.7_package.zip) — manuscript, scripts, verification logs, frozen environment |
| Verification | `v2.7/verification/verify_all.py` — 40 symbolic/numerical checks, expected **40/40 PASS** |
| Website | https://dipolarstrings.org |

New in V2.7: an appendix on the emergent relativistic dynamics of autonomous sub-gap breathers
(proof of mechanism on a discrete sine-Gordon lattice, with measured second-order convergence),
plus the master script reproducing it.

Previous version: [`v2.6/`](v2.6/) — archived on Zenodo as
[10.5281/zenodo.21434867](https://doi.org/10.5281/zenodo.21434867).

## Archive

- Concept DOI (always resolves to the latest version): [10.5281/zenodo.21196098](https://doi.org/10.5281/zenodo.21196098)
- v1.0 version DOI — time-stamped record of the registered LRD prediction (Section 15.2): [10.5281/zenodo.21196099](https://doi.org/10.5281/zenodo.21196099)

## Reproduce

```bash
pip install -r v2.7/requirements.txt
cd v2.7/verification && python3 verify_all.py     # 40/40 PASS, exit 0
python3 ../scripts/p2_master.py rest              # breather appendix, first table
```

`p2_master.py all` includes the h = 0.25 runs and takes several hours.

## Falsification targets

The model is meant to be shot at. The five open confrontations — black-hole shadow diameter
(+4.63 %, ngEHT), Little Red Dot mass function, exponentially spaced gravitational-wave echoes,
absence of any primordial-black-hole evaporation burst, and the 2PN double-pulsar residuals — are
stated in the Discussion of the manuscript, and the known failures are inventoried in its Open
Problems section. Contributions are organised by
[discussion category](https://github.com/rstpier/DipolarString/discussions); see
[contribute.html](https://dipolarstrings.org/contribute.html).

## Licence

Manuscript and figures: CC BY 4.0 ([`LICENSE-MANUSCRIPT.md`](LICENSE-MANUSCRIPT.md)).
Code: MIT ([`LICENSE`](LICENSE)).
