# Manifest — final DS V2.9 branch package (manuscript V2.9.1)

## Publication record

| Path | Function | Status |
|---|---|---|
| `manuscript/DS_model_V2_9_1_Zenodo.pdf` | Publication manuscript | Current |
| `manuscript/DS_model_V2_9_1_Zenodo.tex` | Complete source | Current |
| `figures/fig_branch_v291.png` | A8-corrected equilibrium branch | Regenerated |
| `figures/fig_scaling_v291.png` | Quadratic length scaling | Regenerated |
| `figures/fig_lrd_v291.png` | A8-corrected mass distribution | Regenerated |
| `figures/fig_heaviside.png` | Linearized Heaviside profile | Current with narrowed caption |

## Reproduction code

| Path | Function | Scope |
|---|---|---|
| `verification/verify_v291.py` | Current-claim verification | Authoritative for V2.9.1 |
| `scripts/substrate/substrate_stars_v291.py` | A8 branch, data and three figures | Current |
| `scripts/phase_z3_v291.py` | Compact-phase algebra | Conditional result only |
| `scripts/johns/scn_from_dqd.py` | 2^12 Johns sign scan, 24 proper-cubic rotations, 50-direction Bloch isotropy scan | Reconstructed supplementary artifact named by Theorem 4 |
| `scripts/johns/johns_scn.py`, `bloch.py`, `bands.py`, `select.py` | Johns reconstruction components and band diagnostics | Independent V2.8 implementation |
| `scripts/schrodinger/calc1_common.py` | D3 band continuation and full 3-D map | Reconstructed from archived V2.8 engine |
| `scripts/schrodinger/calc1.py` | Free envelope and truncation hierarchy | Reconstructed; reproduces printed anchors |
| `scripts/schrodinger/calc1b.py` | Two-packet interference | Reconstructed; reproduces printed anchors |
| `scripts/schrodinger/calc3.py` | Action count and gauge-force test | Reconstructed; finite-force mass estimator remains non-identical |
| `scripts/schrodinger/calc2.py` | Local potential and band leakage | Reconstructed; reproduces printed anchors |
| `scripts/breather/p2_master.py` | Appendix D calculations | Historical/current appendix support |
| `scripts/legacy/theorem3_covariance.py` | Scalar-action identity | Mathematical identity retained |
| `scripts/legacy/conformal_identity.py` | Conformal representatives | Mathematical identity retained |
| `scripts/legacy/boosted_residual.py` | Flat wave-operator residual | Linearized scope only |
| `scripts/legacy/ppn_2pn.py` | Exponential/Schwarzschild comparison | Diagnostic only |
| `scripts/legacy/supplementary_v24.py` | Old hydrostatic source | Historical; superseded |

## Historical packages

The `archive/` directory contains eight untouched ZIP archives: the V2.7 reproducibility
package, six Johns sessions, and the closure archive. Checksums are recorded in
`archive/SHA256SUMS.txt`.

## Verification contract

The current suite must exit with status zero. Optional SymPy checks are skipped when
SymPy is unavailable; installing `requirements.txt` enables them. A successful run means
only that the listed consequences follow from the encoded equations.

The Appendix E filenames are the historical names, but their present contents are
reconstructions dated 31 August 2026. See `scripts/schrodinger/RECONSTRUCTION_NOTES.md`.

## Final-package provenance documents

| Path | Function | Status |
|---|---|---|
| `FINAL_PACKAGE_NOTES.md` | Defines final-package/version boundary | Current |
| `APPENDIX_E_AUDIT.md` | Scientific/evidential status of Appendix E reconstruction | Current |
| `BUILD_PROVENANCE.md` | Input archive SHA-256 record | Current |
| `archive/source_snapshots/` | Exact supplied source snapshots | Provenance only |

## Final freeze verification

`verification/FINAL_PACKAGE_VERIFICATION.txt` records the final 40/40 current-claim run
plus the quick Appendix E reconstruction checks executed after package consolidation.
`verification/scn_from_dqd_log.txt` separately records the reconstructed Theorem 4
supplement: 4096 -> 8 -> 2 -> 1, with the selected node isotropic to better than 10^-6
over 50 deterministic random directions.
