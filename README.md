# Dipolar Strings (DS / CD)

Repository for the **Dipolar Strings (DS)** exploratory framework by R. St-Pierre.

The project investigates a deterministic description of the physical vacuum based on
dipolar cells modelled as transmission-line structures. The repository contains the
scientific manuscript, numerical simulations, verification scripts, and archived
releases.

> **Status:** Exploratory research project. The model has not undergone independent
> peer review. Independent verification and critical evaluation are encouraged.
>
> An internal referee-style scientific review of the whole programme (in French,
> 19 November 2026) is in [`REVUE_SCIENTIFIQUE_2026-11-19.md`](REVUE_SCIENTIFIQUE_2026-11-19.md).
> It was written by the project's AI collaborator, not by an independent reviewer.
>
> What can follow the programme — string-theory study, the effective hadronic string, analogue
> gravity on transmission lines — is laid out in a roadmap (in French, 22 November 2026):
> [`FEUILLE_DE_ROUTE_2026-11-22.md`](FEUILLE_DE_ROUTE_2026-11-22.md).
>
> **22/11/2026 — Dipolar Strings changes direction.** The original programme — the vacuum as a network of
> impedance-matched transmission lines — is closed by the V2.11 erratum. From this date the project turns to string
> theory: the dipolar filament is re-read as a relativistic string, an open string with opposite charges at its ends.
> A fundamental string in four dimensions being inconsistent (R127), the author has taken ten dimensions: the
> framework then becomes a known construction, the type IIA superstring on three 2-tori with intersecting D6-branes,
> in which its motifs have an exact place but the gauge group, the three generations and the hadronic tension are
> chosen, not derived (R128). The author then specified nine space dimensions in three groups of three, "probably
> not tori": in a product, three-dimensional groups must be flat tori or their quotients, hidden spheres and infinite
> groups are excluded, and the remaining option is a Calabi–Yau space (R129). The author then chose a Calabi–Yau
> space (R130): a consistent choice that costs nothing in the visible space, but on which only integers can be
> computed, the number of generations remains a chosen integer, and 52 to 410 massless fields must be fixed by a
> further mechanism. Finally the filament is taken as a chain of links (R131), the string-bit formulation: an exact
> place for the original dipolar filament, no new number, and a link of finite size excluded by Lorentz invariance.
> No new claim is made. See
> [`matter/rebuild/TDC.md`](matter/rebuild/TDC.md) (in French).

---

## Latest Release

**Version:** V2.11, erratum edition of 21 November 2026 (package frozen 22 November 2026) — the
frozen V2.10 manuscript with a twelve-item erratum applied in place, every item reducing the
model's claims; plus a companion note of fourteen negative results. Not deposited on Zenodo
(author's decision, 22 November 2026): it is distributed through this repository only.

- Erratum: `v2.11/ERRATUM.md` — manuscript: `v2.11/manuscript/DS_model_V2_11_erratum.pdf`
- Negative results: `v2.11/negative_results/DS_negative_results.pdf`
- Verification: `v2.11/verification/verify_v2_11.py` — source checks + the seven evidence scripts

**Zenodo record:** V2.10, frozen 7 September 2026 — an incremental package over V2.9.2

- Manuscript: `v2.10/manuscript/DS_model_V2_10_Zenodo.pdf`
- LaTeX source: `v2.10/manuscript/DS_model_V2_10_Zenodo.tex`
- Verification: `v2.10/verification/verify_v2_10.py` — 40/40 delegated + 9 source checks PASS
- Package documentation: `v2.10/README.md`; scripts, data and archives in `v2.9/`
- Website: https://dipolarstrings.org

A passing verification run means only that the listed consequences follow from the
encoded equations and the current numerical conventions. It is not third-party
experimental validation of the physical model.

---

## Repository Structure

```
v2.11/     Erratum edition (21 November 2026) — manuscript with erratum E1–E12, note of negative results
v2.10/     Zenodo record — corrected manuscript and its figures (incremental)
v2.9/      V2.9.2 package — scripts, data, verification suite, historical archives
v2.8/      Johns-node microstructure sessions (DQD development)
v2.7/      Archived release
v2.6/      Archived release

matter/    Working directory for the matter sector — not part of any release
```

V2.10 changes three sentences of the manuscript and nothing else, so it does not duplicate
the V2.9.2 code and data; `v2.9/` stays frozen and its checksums remain valid.

`matter/rebuild/` holds the November 2026 reconstruction attempt and its stopping criteria (all
three fired on 21 November 2026: `matter/rebuild/A3.md`).

`matter/` is a research front, not a release. The framework is a theory of the medium,
not of matter; that directory is where an attempt at the matter sector is built, and it
opens with the constraints that make several obvious routes dead ends before they start
([`matter/CONSTRAINTS.md`](matter/CONSTRAINTS.md)).

The `main` branch always contains the latest stable version; development is performed on
dedicated Git branches.

---

## Reproducibility

```bash
python -m pip install -r v2.9/requirements.txt
cd v2.10
python verification/verify_v2_10.py     # exit 0
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

- Latest deposit — https://doi.org/10.5281/zenodo.22559109
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
