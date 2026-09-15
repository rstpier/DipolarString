# Proton characteristic impedance — can the 19 % gap close without fitting?

**Date:** 15 September 2026
**Status:** NEGATIVE, derived. The "19 %" was the distance between a *scaling estimate* and the
target, not the residual of a calculation. Computed with the electron's own electrodynamics and
the model's own closures, the 9-turn closed helix gives `Z_p/Z_e = 14` (the manuscript's
series-capacitance reading) or `3` (the electron's self-capacitance reading) at the matching pitch;
within the validity of each reading no pitch reaches 22.67 — the series reading spans 12–15 with
a minimum near `p = 2.5 D`, the self reading 1.5–3.2. The derivation stays open; what it needs is
a principle that fixes the coil geometry, not a better estimate.
**Script:** `../scripts/proton_impedance_closure.py`, 16/16 PASS.

## 1. The target and the estimate

The mass law `m ∝ N²Z` with `m_p/m_e = 1836.15` and `N_p/N_e = 9` requires `Z_p/Z_e = 22.67`.
The manuscript's estimate for the `N_DS = 27` soliton — a closed toroidal helix of 9 turns of 3
strings, "mutual inductance `L ∝ N²`, series inter-turn capacitance `C ∝ 1/N`" — gives
`N^{3/2} = 27`, quoted as "within 19 %, a residual of logarithmic order comparable to the
`ln(8R₃/r)` corrections of the electron derivation".

## 2. Method — the electron's electrodynamics, nothing else

- The electron's `L_e = μ₀R₃[ln(8R₃/r) − 2]` and `C_e = 4π²ε₀R₃/ln(8R₃/r)` are reproduced to
  0.02 % by a Neumann integral and a Coulomb integral regularised at the wire radius `r`. The same
  two integrals are then applied to the closed helix. `R₃/r = 37.1` enters only through
  `ln(8R₃/r)`, exactly as for `Z_e`.
- Closures, all internal: wire length `27ℓ₁ = 9 × 2πR₃` (27 strings of `ℓ₁` end to end); `N = 9`
  turns; inter-turn pitch = the model's own matching distance `D = 2cosh(π) r`, or the smallest
  pitch at which the toroidal helix does not cut through itself (horn torus, `R_maj = ρ`, which is
  `1.11 D`). Nothing adjusted.
- Two readings of the coil's lumped capacitance: the manuscript's series turn-to-turn `C_tt/N`
  (adjacent turns as a two-wire line — meaningful while the turns are actually side by side,
  `p ≤ 2ρ`), and the self-capacitance of the whole coil — the reading the manuscript itself uses
  for the electron's ring.

## 3. Results

```
  geometry                                        L_p/L_e   C_self/C_e  C_series/C_e   Z_p/Z_e(self)  Z_p/Z_e(series)  m_p/m_e = L_p/L_e
  toroidal helix, pitch D                             20.66       1.953      0.1002         3.25           14.36            20.7
  toroidal helix, horn torus                          20.89       2.045      0.0968         3.20           14.69            20.9
  coaxial stack, pitch D (open, the 'N^2' picture)    13.88       2.882      0.1007         2.19           11.74            13.9
  target (mass law)                                                                        22.67           22.67          1836.2
```

Pitch scan over the self-avoiding range (the one geometric freedom left by the wire-length
closure; `ρ` shrinks as the pitch grows and the coil degenerates into a wavy ring):

```
  pitch scan (self-avoiding range):  p/D   p/2rho   L_p/L_e   Z_p/Z_e(self)   Z_p/Z_e(series)   series reading applies?
                                     1.11    0.35     20.90        3.20           14.69          yes
                                     1.50    0.47     13.49        2.24           12.37          yes
                                     2.00    0.64     11.43        1.89           11.90          yes
                                     2.50    0.81     10.66        1.72           11.88          yes
                                     3.00    0.98     10.35        1.62           12.06          yes
                                     5.00    1.80     10.68        1.46           13.60          no
                                     8.00    4.13     12.64        1.45           18.56          no
                                     9.00    6.31     13.45        1.47           22.57          no
                                     9.50    9.06     13.87        1.48           26.88          no
```

1. **The toroidal helix self-intersects at the matching pitch:** `R_maj = 0.895 R₃ < ρ = 0.995 R₃`.
   The smallest self-avoiding closure has pitch `1.11 D`; the two closures nearly coincide and
   give the same answer.
2. **`L ∝ N²` assumes perfect coupling; two turns at spacing `D` have `k₁ = 0.19`.** The actual
   inductance of the closed helix is `L_p/L_e = 20.7`, not 81.
3. **`Z_p/Z_e = 14.4` (series reading) or `3.25` (self reading)** at the matching pitch. Within the
   validity of each reading the series value spans `11.9–14.7`, with a **minimum of 11.9 near
   `p = 2.5 D`** — an internal extremum exists, at half the target — and the self value falls
   monotonically from `3.2` to `1.5`. **No admissible configuration reaches 22.67.** The series
   formula does cross the target, near `p ≈ 9 D`, but there `p/2ρ ≈ 6`: the turns are no longer
   adjacent and the two-wire formula does not apply.
4. Read as inductance = mass (`L = κm`, the manuscript's own equivalence), the same coil gives
   `m_p/m_e = 21`, against 1836.
5. Coincidences, recorded and **not used**: `(N−1)^{3/2} = 22.63` (−0.2 %) and `27/2^{1/4} = 22.70`
   (+0.2 %). No mechanism produces either; a match without a mechanism is a fit by another name.

## 4. What this changes

- Item 1 of `OPEN_PROBLEMS.md` is not a "19 % residual of logarithmic order". The calculation lands
  a factor 1.6 (series) to 7 (self) *below* the target, where the estimate landed 19 % *above* it.
  The residual is not logarithmic and its sign is the other way.
- The obstacle is not a coefficient but the object: the coil's turn radius, pitch, major radius and
  the meaning of its lumped capacitance are fixed by no axiom — and the manuscript's two size
  closures (17 fm from the resonance-radius law, 3.5×10³ fm from `2πR_N = Nℓ₁`) already disagree by
  a factor 200. The proton stays *Reparametrized*.
- A4 (relaxation toward matching, `|Γ| → 0`) cannot help: it pushes any structure toward `Z₀`, i.e.
  toward `Z_p/Z_e ≈ 1.4`, the opposite of a high-impedance proton. An extremum principle on `Z(p)`
  would select the series minimum, `≈ 12`, not 22.67.
- What would close it: a principle that fixes the coil geometry and the capacitance reading, or a
  different object for `N_DS = 27` than a 9-turn helix. Once such a principle exists, the two
  integrals of §2 give `Z_p` with no further input; the machinery is in place.

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §4 (a fit is not a derivation).
