# Appendix E code reconstruction

## Provenance

The original files `calc1.py`, `calc1b.py`, `calc1_common.py`, `calc2.py`, and
`calc3.py` were not found in the supplied archives, the public `main` branch, the
`v2.8-foundations` branch, or the visible GitHub commit history. The files in this
directory were reconstructed on 31 August 2026. They are **not** presented as recovered
historical source code.

The reconstruction uses only:

- the Johns scattering matrix `../johns/S_johns.npy`;
- the D3 braid map in archived V2.8 `defect.py` and `track.py`;
- the protocol and numerical anchors printed in Appendix E of the V2.9 manuscript.

## Mapping

| File | Reconstructed role |
|---|---|
| `calc1_common.py` | Sparse D3 Bloch operator, 30-point eigenvector continuation, packet and full-map utilities |
| `calc1.py` | Free-packet dispersion and quadratic/cubic truncation diagnostics |
| `calc1b.py` | Two-packet interference and ring-bin test |
| `calc3.py` | Single action-count closure and single-band inertial response |
| `calc2.py` | Projected local-potential map; optional full 185,856-component evolution |

## Reproduction

From the package root:

```bash
python scripts/schrodinger/calc1_common.py
python scripts/schrodinger/calc1.py
python scripts/schrodinger/calc1b.py
python scripts/schrodinger/calc3.py
python scripts/schrodinger/calc2.py --localized
python scripts/schrodinger/calc2.py --height 0.015 --full-map
```

The first command creates the compressed 30-point band cache. The full-map option is
the expensive path. The quick potential calculation is the exact projection of one
map step onto the measured 30-mode band.

## Numerical cross-check against Appendix E

| Observable | Appendix E | Reconstruction |
|---|---:|---:|
| Free variance at 480 steps | 23.455 | 23.454608 |
| Exact dispersive predictor | 23.452 | 23.454328 |
| Quadratic prediction | 22.87 | 22.861924 |
| Quadratic / cubic phase residual | 0.084 / 0.016 | 0.08166 / 0.01460 |
| Modal amplitude / phase error | `1.9e-8` / `1.5e-9` | `1.65e-9` / `1.09e-9` |
| Interference peak | 0.392699 | 0.392699 |
| Interference correlation | 0.99962 | 0.999826 |
| `E/N` | 0.36713 | 0.3671338 |
| `P/N` | 1.5024 | 1.5023599 |
| Gauge velocity error | at most 1.5% | 1.425% |
| Positive-potential full/projected delay | 9.2 / 9.3 | 9.176 / 9.135 |
| Negative-potential full/projected delay | 8.0 / 5.1 | 8.188 / 5.018 |
| Negative-potential off-band fraction | 0.78% | 0.738% |
| Localized frequency | 0.41819 | 0.418189 |
| Localized weight in stated window | 88% | 88.24% (`|z-z_c| <= 11`) |

No manuscript anchor was used as a fitted continuous parameter. The discrete packet
placement and the well-weight integration window are stated in the code because those
conventions were not printed in Appendix E.

## Scientific boundary

Agreement with a printed anchor validates the reconstructed implementation, not the
physical premise of the model. Disagreement must be reported, not tuned away. In
particular, the historical finite-force leakage correction behind the reported
`m_F omega'' = 0.975` was not recoverable from the surviving material; `calc3.py`
therefore reports the exact single-band value and prints the historical full-map value
only as a comparison target.

## Final-package audit note

For publication, the `calc2.py` on-band number should be described as projection onto the
reconstructed 30-mode continued bound-band subspace. The historical full-band projector
was not recovered. The interference calculation is a useful implementation-consistency
check, but the existence of the cross term follows from linear superposition. See
`../../APPENDIX_E_AUDIT.md`.
