# Appendix E reconstruction audit — final V2.9 branch package

Date: 31 August 2026

## Verdict

The Appendix E computational layer is **substantially reconstructible but not historically
recovered**. The current scripts reproduce most printed numerical anchors from a plausible
implementation constrained by the manuscript and surviving V2.8 machinery. That supports
internal consistency of the published calculation protocol. It does not prove that the
same historical code produced the published numbers.

## What is reproducible in the included reconstruction

- bound-band continuation and the quoted local band anchor;
- free-packet variance and the exact/quadratic dispersive comparison;
- quadratic-to-cubic phase-residual hierarchy;
- two-packet interference peak and correlation;
- static `E/N` and `P/N` action-count consistency;
- gauge-sweep velocity relation to the measured band derivative;
- positive and negative local-potential traversal delays;
- the localized projected state frequency and stated integration-window weight.

The final current-claim suite reports 40/40 PASS for its declared scope.

## Limits that must remain visible

### 1. Missing historical source

The five Appendix E programs are reconstructions, not recovered originals. Numerical
agreement therefore tests whether the manuscript is reproducible under the reconstructed
protocol, not whether the historical run has been independently authenticated.

### 2. Inertial estimator

The historical mean `m_F * omega'' = 0.975` and its quoted range are not reproduced by
the reconstructed finite-difference estimator. The exact single-band identity gives 1,
while the full-map numerical estimate depends on derivative-window conventions. The
package keeps this discrepancy open rather than tuning it away.

### 3. "On-band" / leakage metric

`calc2.py` projects the full state onto the reconstructed 30-mode continued bound-band
subspace. Because the historical projector is unavailable, this metric should be called
**reconstructed subspace leakage** or **30-mode bound-band projection leakage** when
precision matters. Its agreement with the printed 0.78% anchor is not proof that every
state outside those 30 modes is physically outside the complete bound band.

### 4. Interference test

The cross term of two linearly superposed complex fields is algebraically expected. The
numerical value is still useful as a check that the lattice reconstruction, packet phases
and Fourier-bin conventions are mutually consistent, but it is not by itself evidence of
nonclassical quantum measurement physics.

### 5. Action-quantum interpretation

The `E/N`, `P/N` and `m_F omega''` closures establish consistency of one normalization
within the reconstructed single-band dynamics. They do not independently derive the
physical value of Planck's constant or the Born rule.

## Superseded reconstruction

The archived alternate `annexeE_reconstruction(1).zip` contains a different defect model
whose two-defect splitting is non-monotonic and whose exponential/power fits have poor
R^2. Its printed statement that `1/d` is excluded is not supported by that dataset. Its
orientation comparison also must not be presented as a verified reconstruction of the
published opposite-braid-handedness factor. It is retained only as provenance of the
reconstruction process.

## Safe publication wording

Recommended description:

> "The original Appendix E source files were unavailable at final packaging. We provide a
> transparent reconstruction based on the published protocol and surviving lattice code.
> The reconstruction recovers nearly all printed numerical anchors, with the historical
> finite-force mass estimator remaining unresolved. These files establish reproducibility
> of a compatible implementation, not byte-for-byte recovery of the original computation."
