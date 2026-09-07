# Matter sector — working directory

**Status: empty research front. Nothing here is part of any released version.**

The Dipolar Strings model is, in its own words, *a theory of the medium, not of
matter*. The vacuum, light and collective gravitation are in scope; the weak and
strong sectors, neutrino masses, cosmological evolution, the Born rule and a
microscopic account of entanglement are not derived. This directory is where an
attempt at the matter sector is built, and it is deliberately kept outside the
versioned release folders (`v2.6/` … `v2.10/`) until it produces something that
survives `CONSTRAINTS.md`.

## Read this first

[`CONSTRAINTS.md`](CONSTRAINTS.md) — **the hard boundary**. The manuscript does
not merely leave the matter sector undone; on several fronts it establishes that
particular routes are closed, one of them by more than thirty orders of
magnitude. Any work started here that violates a constraint in that file is
wasted before it begins. Read it before writing a line.

[`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md) — the matter-sector entries of the
manuscript's own inventory, with their status and what would close each one.
This is the work list.

## What the model already fixes

| | |
|---|---|
| Anchor | R₃ = ℏ/m_e c₀, hence ℓ₁ = 2πR₃/3 — the single experimental anchor |
| Electron impedance | Z_e ≈ 0.73 Z₀ ≈ 275 Ω — *Calibrated*, not *Derived* |
| Mass law | m ∝ N_DS² · Z_{N_DS} |
| Resonance radius | R_{N_DS} = (p·ℓ₁/2π)(Z₀/Z_{N_DS}) |
| States in scope | photon (N_DS = 1), DQD (N_DS = 2), electron (N_DS = 3) |

The muon, tau and the weak bosons were tabulated in earlier versions as
reverse-engineered mass fits without predictive content, and are **withdrawn**
from the catalogue. The surviving table is a *reverse-engineered
parametrization, not a set of predictions* — the manuscript labels it so itself.

## Layout

```
CONSTRAINTS.md    what is closed, and why — read before starting
OPEN_PROBLEMS.md  the work list, from the manuscript's own inventory
notes/            derivations in progress, one file per attempt
scripts/          code; anything numerical ships with the prose
verification/     checks a matter-sector claim must pass to be merged
paper/            article skeleton for deposit — compiles clean, awaits content
```

## Rules

The repository's ground rules apply here unchanged, plus one:

- **Declared maturity level.** Every statement is filed at one of the paper's six
  levels — *Derived*, *Conditionally derived*, *Calibrated*, *Postulated*,
  *Reparametrized*, *Closed*. A ratio obtained by fitting a measured mass is
  *Reparametrized*, not *Derived*. This is the failure mode that produced the
  withdrawn muon/tau/boson entries.
- **Ratios, not masses.** See `CONSTRAINTS.md` §1. An absolute mass is not a
  target in this framework; a spectrum of ratios given one scale is.
- **Quantified or nothing.** A numerical claim ships with the code that produces
  it and a check under `verification/`.

## Licence

Manuscript and figures: CC BY 4.0. Code: MIT.
