# Dipolar Strings (DS / CD)

Repository for the **Dipolar Strings (DS)** exploratory framework by R. St-Pierre.

The project investigates a deterministic description of the physical vacuum based on
dipolar cells modelled as transmission-line structures. The repository contains the
scientific manuscript, numerical simulations, verification scripts, and archived
releases.

> **Status:** Exploratory research project. The model has not undergone independent
> peer review. Independent verification and critical evaluation are encouraged.

---

## Latest Stable Release

**Version:** V2.9 — manuscript revision V2.9.2, package frozen 7 September 2026

- Manuscript: `v2.9/manuscript/DS_model_V2_9_2_Zenodo.pdf`
- LaTeX source: `v2.9/manuscript/DS_model_V2_9_2_Zenodo.tex`
- Verification: `v2.9/verification/verify_v291.py` — 40/40 current-claim tests PASS
- Package documentation: `v2.9/README.md`, `v2.9/MANIFEST.md`
- Website: https://dipolarstrings.org

A passing verification run means only that the listed consequences follow from the
encoded equations and the current numerical conventions. It is not third-party
experimental validation of the physical model.

---

## Repository Structure

```
v2.9/      Latest stable release — manuscript, scripts, verification, historical archives
v2.8/      Johns-node microstructure sessions (DQD development)
v2.7/      Archived release
v2.6/      Archived release
```

The `main` branch always contains the latest stable version; development is performed on
dedicated Git branches.

---

## Reproducibility

```bash
python -m pip install -r v2.9/requirements.txt
cd v2.9
python verification/verify_v291.py      # 40/40 PASS, exit 0
```

The release also ships the substrate, Johns-node, Schrödinger and breather calculations
under `v2.9/scripts/`, the regenerated figures under `v2.9/figures/`, and the untouched
historical V2.7 and V2.8 packages under `v2.9/archive/`.

**Appendix E caveat:** the five Appendix E programs in `v2.9/scripts/schrodinger/` are
reconstructions dated 31 August 2026, not the recovered historical originals. Their
agreement with the printed anchors is reconstruction consistency, not recovery. See
`v2.9/APPENDIX_E_AUDIT.md` and `v2.9/MISSING_ARTIFACTS.md`.

---

## Releases

Concept DOI (always points to the latest release):

- https://doi.org/10.5281/zenodo.21196098

Archived version records:

- V2.6 — https://doi.org/10.5281/zenodo.21434867
- v1.0 — https://doi.org/10.5281/zenodo.21196099

---

## Contributing

The project welcomes independent verification, reproduction attempts, criticism,
and discussions. See [contribute.html](https://dipolarstrings.org/contribute.html).

---

## License

- Manuscript and figures: CC BY 4.0
- Code: MIT
