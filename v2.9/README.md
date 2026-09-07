# Dipolar String Model — final V2.9 branch publication package

Author: R. St-Pierre  
Package date: 7 September 2026  
Manuscript version: **V2.9.2**  
Concept DOI: 10.5281/zenodo.21196098

This archive is the final consolidated package for the V2.9 development branch. The
manuscript is V2.9.2, the version encoded in the PDF, the LaTeX source and the citation
metadata; V2.9.1 was the 31 August 2026 Zenodo edition, superseded by the two corrections
recorded under "V2.9.2 (September 2026)" in the Document History. The package name uses
"V2.9 FINAL" only as a branch label; it does not rewrite the manuscript's version history.

The bundle combines the current manuscript and figures, current verification scripts,
A8-consistent substrate calculations, Johns-node reconstruction, Appendix D breather
support, historical V2.7/V2.8 archives, and a transparent reconstruction of the five
missing Appendix E programs.

## Quick verification

```bash
python -m pip install -r requirements.txt
python verification/verify_v291.py
python scripts/johns/scn_from_dqd.py
python scripts/substrate/substrate_stars_v291.py
python scripts/schrodinger/calc1.py
python scripts/schrodinger/calc1b.py
python scripts/schrodinger/calc3.py
python scripts/schrodinger/calc2.py --localized
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error DS_model_V2_9_2_Zenodo.tex
```

Recorded final-package check on 31 August 2026: **40/40 current-claim verification tests
PASS**. This means the listed consequences follow from the encoded equations and current
numerical conventions. It is not third-party experimental validation of the physical
model.

## Directory map

- `manuscript/` — V2.9.1 PDF, LaTeX source and build products.
- `figures/` — publication figures; regenerated V2.9.1 figures use `_v291`.
- `data/` — regenerated A8-consistent branch data and summary.
- `scripts/substrate/` — current hydrostatic integration and figure generation.
- `scripts/johns/` — independent Johns-node reconstruction and band calculations.
- `scripts/schrodinger/` — transparent reconstruction of the five missing Appendix E
  programs and the reconstructed bound-band cache.
- `scripts/breather/` — Appendix D breather master calculation.
- `scripts/legacy/` — retained diagnostics with narrowed or historical scope.
- `verification/` — current-claim verification and recorded logs.
- `archive/` — immutable historical packages and source snapshots, including a superseded
  alternate Appendix E reconstruction kept only for provenance.
- `APPENDIX_E_AUDIT.md` — exact evidential status of the reconstructed Appendix E layer.
- `FINAL_PACKAGE_NOTES.md` — package-level scientific and provenance boundary.
- `MISSING_ARTIFACTS.md` — list of historical originals that remain unavailable.
- `BUILD_PROVENANCE.md` — SHA-256 hashes of the three supplied input archives.

## Appendix E status

The historical `calc1.py`, `calc1b.py`, `calc1_common.py`, `calc2.py` and `calc3.py`
were not recovered. The versions in `scripts/schrodinger/` are reconstructions dated
31 August 2026. They reproduce nearly all printed Appendix E numerical anchors, but that
agreement must be described as **reconstruction consistency**, not recovery of the
original source or independent confirmation of the historical execution.

Three specific caveats remain important:

1. `calc3.py` does not recover the historical finite-difference convention that produced
   the reported mean `m_F * omega'' = 0.975`. The exact reconstructed single-band identity
   gives 1.000000; the historical 0.975 value and its finite-force convention remain unresolved.
2. The `calc2.py` "on-band" value is a projection onto the reconstructed 30-mode continued
   bound-band subspace used by this reconstruction. It is a proxy for the historical
   off-band diagnostic, not proof that the exact historical projector has been recovered.
3. The two-packet interference cross term follows algebraically from linear superposition.
   Its numerical agreement checks packet phase and Fourier-bin conventions; it is not, by
   itself, evidence for nonclassical measurement physics or the Born rule.

See `APPENDIX_E_AUDIT.md` for the full wording boundary.

## Current scientific status

The package verifies consequences of the model's equations, not the physical truth of
the equations. In particular:

- the A8-consistent substrate-star branch is reproducible;
- the Johns constraint scan is reproducible: 4096 sign motifs -> 8 symmetric/unitary candidates -> 2 proper-cubic candidates -> 1 node passing the deterministic 50-direction isotropy scan, equal to `S_johns.npy`;
- the coefficient `3/(4 alpha)` is a conditional logarithmic phase-slip coefficient;
- no complete instanton action or equation connecting instanton fugacity to `G` is supplied;
- entanglement, the Born rule and measurement are not derived;
- scalar-profile PPN differences are diagnostics, not exterior predictions after the
  General Relativity bootstrap.

## Historical integrity

Files in `archive/` are preserved as provenance. Their successful historical tests must
not automatically be interpreted as V2.9.1 validation. See `archive/ARCHIVE_STATUS.md`.

## Licences

Manuscript and figures: CC BY 4.0. Code: MIT. See `LICENSE.md`.
