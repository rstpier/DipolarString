# The electron as a partially open turn — does closure act exponentially on the mass?

**Date:** 15 September 2026
**Origin:** the author's long-held hypothesis: *the mass loop has an exponential effect on
closure; the electron would be a partially open turn.*
**Status:** NEGATIVE for the proton, derived — closure moves a mode's energy by factors of order 1,
algebraically; the model's exponentials relate geometry to impedance, which makes mass
*logarithmic* in geometry. POSITIVE for item 5 — the open turn and the antiperiodic ring of the
unordered pair are the same half-integer ladder at the fundamental, and differ at the second
harmonic.
**Script:** `../scripts/open_turn_closure.py`, 16/16 PASS.

## 1. The object and its exact resonance condition

A turn of string length `ℓ_s` and impedance `Z_e`, closed through a gap of capacitance `C_g`.
The round-trip ABCD matrix (line × series element) must have eigenvalue 1, i.e. trace 2:

    2 cos(βℓ_s) − 2 + sin(βℓ_s) / (ωC_gZ_e) = 0,   ω = βc,

equivalently `tan(βℓ_s/2) = 1/(2ωC_gZ_e)` together with the untouched family `sin(βℓ_s/2) = 0`
(standing waves with a current node at the gap, which a series element cannot see). Verified
symbolically. Since the tangent is positive, the solutions live on `(0, π/2)` and `(π, 3π/2)`:

| family | closed turn | fully open turn |
|---|---|---|
| Compton mode | `βℓ_s = 2π` | `3π` (×1.5) |
| current node at the gap | `2πp` | `2πp` (untouched) |
| DC circulating current | `0` | `π` — the open turn's half-wave fundamental |

**Closure moves mode energies by factors of order 1.** The open fundamental is half the closed
Compton mode; the Compton mode itself rises by 3/2 on opening. A factor 1836 is three orders of
magnitude out of reach.

## 2. The gap dependence is algebraic, and it saturates

With `C_g = ε₀πr²/g` (parallel plates), the closed-to-open transition sits at `g* = 0.04 r`
(0.4 fm): any visible opening is "open". Beyond it the deviation from the open value falls
*linearly* with the gap (ratios 9.6 and 10 per decade). Closure makes the mass **insensitive** to
the gap, not hypersensitive: a gap 100× wider moves the mode 100× less. An exponential sensitivity
would need the *frequency* to carry `e^{−g/ξ}`; here every mode is pinned between two multiples
of `π/2`, and only the coupling across the gap (a Q factor, a splitting) can be exponential.

## 3. Where the model's exponentials actually are

`R/r = e^{ℓ}/8` with `ℓ(ℓ−2) = (2πZ/Z₀)²`, and `D/r = 2cosh(πZ/Z₀)`: the **geometry is exponential
in the impedance**, so a mass proportional to `Z` is **logarithmic in the geometry**. That is the
intuition's exponential, read the right way round — and it goes against the hypothesis: reaching
`Z_p = 22.67 Z_e` by aspect ratio needs `ℓ = 105`, `R/r = 10⁴⁵`; by bifilar spacing, `D/r = 10²³`.
A factor 1836 in mass needs a factor `e^{100}` in geometry; closing a turn by a few tube radii
moves the mass by percent.

## 4. The only exponential-in-closure mass law, and why it is excluded

A tunnelling splitting `m = ℏω₀ e^{−S}` is the one reading in which closure (an action barrier
`S`) enters an exponent. It needs `S_e − S_p = ln 1836 = 7.5` with no internal `S`, and a prefactor
`ℏω₀ ≥ m_pc²`, 1900× above the Brillouin cutoff `(3/π)m_ec²`. Excluded by the lattice, and a new
postulate anyway.

## 5. What survives — the open turn *is* the antiperiodic ring

The open turn's fundamental, `βℓ_s = π`, is exactly the antiperiodic ring of the unordered pair
(`2026-09-15-hopfion-direction.md` §5c): `θ_tot = π` from two open ends, or from one
half-twist — **degenerate at the fundamental**. The author's intuition and the derived ℤ₂ holonomy
are the same half-integer ladder. They differ at the second harmonic: the open turn has all
harmonics (`βℓ_s = pπ`), the antiperiodic ring only odd multiples of `π` — a distinguishing
consequence if the radial ladder is ever resolved. If the electron *is* the open turn, its loop is
half the size (`R = R₃/2`, `λ_C = 2ℓ_s`), `ln(4R₃/r) = 5.00` and `Z_e = 0.62 Z₀`: still an impedance
well. Good for item 5; nothing for the proton.

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §1c (the proton is not a loop of the
vacuum's string under any reading).
