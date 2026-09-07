# Dipolar String Model — V2.10 (incremental package)

Author: R. St-Pierre
Package date: 7 September 2026
Manuscript version: **V2.10**
Concept DOI: 10.5281/zenodo.21196098
Supersedes: **V2.9.2**, which remains archived unchanged in [`../v2.9/`](../v2.9/)

## What this package is

V2.10 changes **three sentences** of the manuscript and nothing else. No equation,
no numerical value, no figure, no script is modified. This package is therefore
incremental: it carries the corrected manuscript and the four figures the
manuscript includes. Everything else — scripts, data, verification logs, the
historical V2.7 and V2.8 archives, the Appendix E reconstruction and its audit —
is unchanged from V2.9.2 and lives in `../v2.9/`, whose checksums remain valid.

## The three corrections

All three remove an internal inconsistency in favour of the more conservative
reading. None of them changes a result.

1. **Primary falsification target.** A leftover sentence in the little-red-dot
   section still named the shadow-diameter deviation as the model's primary
   target, contradicting its withdrawal in the same revision and the revised
   falsification hierarchy. It now names the permanent absence of PBH evaporation
   bursts, the first entry of that hierarchy.
2. **Count of maturity levels.** The epistemological-status list announced six
   levels of maturity while presenting seven headings, the seventh being the
   record of what the revision withdraws rather than a level. The announcement is
   corrected, not the list.
3. **Weave-chirality entry.** It was titled *resolved, with a residual question*,
   which reads as a closed item. What is dissolved is the chiral closure of
   Theorem 4; the converse question — by what mechanism a chirality that is mute
   in the electromagnetic sector could acquire physical relevance — remains open,
   and the title now says so.

They are recorded in the manuscript's own Document History under
*V2.10 (September 2026)*.

## Verify

```bash
python -m pip install -r ../v2.9/requirements.txt
python verification/verify_v2_10.py     # exit 0
```

The script runs the **unchanged V2.9.2 numerical suite** from `../v2.9/` (40/40
PASS) and adds nine source-level checks: the four manuscript-hygiene checks
re-run against the V2.10 source, and five checks that the three corrections and
the version stamp are actually present. It exits non-zero if either half fails.

A passing run means the listed consequences follow from the encoded equations and
the current numerical conventions. It is not third-party experimental validation,
and it is not independent verification: both implementations originate from the
same collaboration, and the first has not been re-run on the A8-consistent source.

## Build the manuscript

```bash
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error DS_model_V2_10_Zenodo.tex
```

67 pages, zero undefined references. `figures/` holds the four figures the
manuscript includes; the regenerated `_v291` variants and the full figure set
remain in `../v2.9/figures/`.

## Contents

```
manuscript/   V2.10 PDF, LaTeX source and build products
figures/      the four figures included by the manuscript
verification/ verify_v2_10.py — delegates to v2.9, adds the V2.10 source checks
```

`FILE_SHA256SUMS.txt` and `PACKAGE_CONTENTS.txt` cover **this directory only**.
For the V2.9.2 package integrity, use `../v2.9/FILE_SHA256SUMS.txt`.

## Licences

Manuscript and figures: CC BY 4.0. Code: MIT. See `../v2.9/LICENSE.md`.
