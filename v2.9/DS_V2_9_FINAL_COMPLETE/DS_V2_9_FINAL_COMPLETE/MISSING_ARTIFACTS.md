# Missing historical originals and reconstructed replacements

The final V2.9 branch package contains the current publication material and all supplied
source snapshots, but it is **not** a byte-for-byte recovery of every historical program.

## Appendix E originals not recovered

The following historical files named by Appendix E were not available in the supplied
archives or recovered source history:

- `calc1.py`
- `calc1b.py`
- `calc1_common.py`
- `calc2.py`
- `calc3.py`

Transparent replacements with the same functional roles are included under
`scripts/schrodinger/`. They are explicitly marked `RECONSTRUCTED` in their output and
must not be called recovered originals.

The reconstruction recovers the principal printed anchors for free-packet dispersion,
phase-truncation hierarchy, interference, the action-count closure, potential delays,
the reconstructed leakage proxy and the localized state. The historical finite-force
mass convention behind `m_F * omega'' = 0.975` remains unresolved.

The historical definition of the Appendix E "off-band" projector is also unavailable.
The reconstructed implementation measures projection onto a 30-mode continued
bound-band subspace. Its numerical agreement with the printed leakage anchor is useful,
but does not establish that the original projector was identical.

## Other historical gaps

A second independently written dimensional re-integration of the current A8 source
`u-u0` was not supplied. The package therefore uses the dimensionless implementation
and does not present the superseded V2.7 dimensional code as independent confirmation.

The early historical `scn_from_dqd.py` source remains lost. A current reconstructed
replacement is now included at `scripts/johns/scn_from_dqd.py`. It performs the
exhaustive 2^12 sign scan, the 24 proper-cubic rotations, and a deterministic
50-direction Bloch isotropy scan, and verifies that the unique selected node equals
`scripts/johns/S_johns.npy`. It must not be described as recovery of the historical bytes.

## Superseded alternate reconstruction

The separately supplied `annexeE_reconstruction(1).zip` is preserved under
`archive/source_snapshots/` for provenance. Its two-defect splitting/chirality analysis
and several potential diagnostics do not reproduce the final Appendix E reconstruction
and are **not** current verification inputs.
