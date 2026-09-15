# Proton characteristic impedance — can the 19 % gap close without fitting?

**Date:** 15 September 2026
**Status:** NEGATIVE, derived — extended in §5 to the whole class (no object of 27 strings
reaches the targets) and in §6 to the mode-energy reading (its only closure gives the wrong
direction; the tube radius caps a fundamental at `37 m_e`; admitting the proton loop destroys the
electron's impedance well). The proton is not a loop of the vacuum's string under any reading. The "19 %" was the
distance between a *scaling estimate* and the
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

## 5. Another object than the helix? — a bound on the whole class

Script: `../scripts/proton_object_search.py`, 9/9 PASS. Same inputs, same integrals; the question
is now whether *any* closed structure made of 27 strings of `ℓ₁` can reach the targets set by the
mass law `m = κL`, i.e. `Z_p/Z_e = 22.67` and `L_p/L_e = 1836`.

**A. The mass law's own object.** `m ∝ N²Z` with the wire-length closure is the statement "a flat
coil of `N_turn = 9 × 22.67 = 204` coincident turns of radius `R_p = R₃/22.67 = 17.0 fm`" — and
`17.0 fm = 1.64 r`: a turn radius below two tube radii, where the thin-wire idealisation the law
rests on is void. Built physically (pitch `2.1 r`, a spring 11.5 R₃ long), that coil has
`L_p/L_e = 5.8`.

**B. Loose objects.** One 27-string ring (radius `9R₃`): `Z_p/Z_e = 1.49`. Three orthogonal
9-string rings: `1.69`.

**C. The best coherent coil — the bound.** Multilayer solenoids wound as tightly as the tube allows
(pitch = layer spacing = `2.1 r`), scanned over inner radius and layer count:

```
  multilayer coil scan (pitch = layer spacing = 2.1 r):  a0/r  layers  turns  height/R3   L/L_e   C/C_e   Z_p/Z_e
                                                         1.5     1    217.3     12.30      5.6    2.19     1.60
                                                         1.5     2    129.6      3.67     12.7    1.00     3.57
                                                         1.5     4     71.5      1.01     37.1    0.58     7.98
                                                         1.5     8     37.7      0.27     66.4    0.60    10.49
                                                         3.0     1    110.6      6.26      8.9    1.45     2.47
                                                         3.0     2     82.1      2.32     20.2    0.81     4.99
                                                         3.0     4     54.2      0.77     47.1    0.57     9.13
                                                         3.0     8     32.2      0.23     66.8    0.64    10.25
                                                         6.0     1     55.6      3.15     18.1    1.08     4.10
                                                         6.0     2     47.3      1.34     36.0    0.71     7.09
                                                         6.0     4     36.5      0.52     62.1    0.59    10.24
                                                         6.0     8     25.0      0.18     67.2    0.71     9.72
                                                        10.0     1     33.4      1.89     29.3    0.97     5.49
                                                        10.0     2     30.2      0.85     51.1    0.73     8.37
                                                        10.0     4     25.4      0.36     70.8    0.67    10.26
                                                        10.0     8     19.2      0.14     67.2    0.81     9.09
                                                        15.0     1     22.3      1.26     39.3    0.99     6.31
                                                        15.0     2     20.8      0.59     60.0    0.81     8.60
                                                        15.0     4     18.4      0.26     71.6    0.79     9.52
                                                        15.0     8     14.9      0.11     61.5    0.95     8.03
                                                        25.0     1     13.4      0.76     47.3    1.15     6.42
                                                        25.0     2     12.8      0.36     61.2    1.02     7.75
                                                        25.0     4     11.9      0.17     63.5    1.03     7.84
                                                        25.0     8     10.3      0.07     57.8    1.18     6.99
```

Maximum `L_p/L_e = 72` (4 layers, 18 turns, a squat coil 0.26 R₃ high); maximum `Z_p/Z_e = 10.5`.
Since a coherent tight coil maximises inductance for a given wire, and the uniform-charge
capacitance is a lower bound (so `Z` is over-estimated), **every object of the class falls short:
by a factor 26 in `L` and 2.2 in `Z`.**

**D. Scaling.** At the packing limit `L_max ∝ W^{5/3}`: carrying `1836 L_e` as inductance needs
`N ≈ 27 × 26^{3/5} ≈ 190` strings, seven times the catalogue's 27.

**E. The other reading.** If mass is the `1/ℓ` mode energy of the loop (§5f of the Hopfion note),
the proton loop has radius `R₃/1836 = 0.21 fm` — a factor 4 from the charge radius `0.84 fm`,
where the inductance reading gives 17 fm (×20) and the geometric closure 3.5×10³ fm (×4000). It is
the least inconsistent reading, but it is `E = hc/λ`: it locates the proton at its Compton scale
and derives nothing until a closure fixes `ℓ_p/ℓ_e = 1/1836`.

**Conclusion.** The obstacle is not the helix. Under the electron's own electrodynamics, **27
strings of `ℓ₁` cannot carry the proton's mass as inductance, whatever their shape** — the mass law
`m = κL` and the index `N_p = 27` are jointly unrealizable. Either `N_p` is of order 200, or mass
is not the inductance of the string. Both are changes to the catalogue, not to the object.
**Status:** NEGATIVE, derived (bound by explicit construction; thin-wire integrals accurate to
~10 % at contact, far inside the margin).

## 6. The mode-energy reading, tested — can `ℓ_p/ℓ_e` be fixed without fitting?

Script: `../scripts/proton_mode_energy_reading.py`, 11/11 PASS.

The reading: a particle is the fundamental TEM mode of a closed loop, `mc² = (p − x)·2πℏc/ℓ_loop`.
The electron already sits there by calibration (`2πR₃ = λ_C`, `p = 1`). The proton then needs
`ℓ_p = λ_C/1836 = 1.32 fm`, `R_p = 0.21 fm`. The model's internal conditions, one by one:

- **(a) Geometric closure `ℓ = Nℓ₁`** gives `ℓ_p/ℓ_e = 9`: the fundamental makes the proton *nine
  times lighter* than the electron — the reading inverts the catalogue's ordering (more strings,
  longer loop, lighter). Matching 1836 on the 27-string loop needs the harmonic `p_p ≈ 16 500`
  (`11 000–16 500` with the holonomy shifts `x_e ∈ {0, 0.27, ⅓}`), with no natural `x_p`: a fit.
- **(b) Tube-radius floor.** A loop of a tube of radius `r` has `ℓ ≥ 2πr`, so `m/m_e ≤ p·R₃/r = 37.1 p`.
  At `p = 1` the proton is excluded by 50, the muon by 5.6, the tau by 94. The harmonic rescue
  needs `p_μ ≥ 6`, `p_p ≥ 50`, `p_τ ≥ 94` — integers only by fitting.
- **(c) The `r`–`Z_e` lock.** Admitting `R_p = 0.21 fm` as a loop of the same string needs
  `r ≤ 0.21 fm`, i.e. `R₃/r ≥ 1836`, i.e. `ln(8R₃/r) ≥ 9.6`, i.e. **`Z_e = 1.36 Z₀ > Z₀`**: the
  electron stops being an impedance well, the manuscript's one robust derived property of the
  electron. The impedance-well property alone caps `R₃/r` below 197 (where `Z_e = Z₀`), hence
  `m/m_e < 197 p`: even the loosest admissible tube needs `p_p ≥ 10`.
- **(d) Brillouin.** `m_pc²` is 1900× the lattice cutoff `ℏω_max = (3/π)m_ec²` and the loop is
  600× smaller than the cell `ℓ₁`: the mode cannot be a collective wave of the weave; it would have
  to be a sub-cell wave on a string thinner than the vacuum's, which nothing internal supplies.

**Conclusion.** NEGATIVE, derived. The mode-energy reading fixes nothing: its only closure gives
the wrong direction, its floor excludes the proton as a fundamental, and thinning the string to
admit the proton destroys the electron's impedance well. It is `E = hc/λ`, no more — and with this
string it *excludes* the proton (and the muon) as fundamentals of a loop of the vacuum's string.

**Where the proton problem now stands.** Three readings of "mass", three failures for `N_p = 27`
with the vacuum's string: inductance (max `72 L_e`, §5), charge energy (0.4 % of the mass,
`CONSTRAINTS.md` §1b), mode energy (floor `37 m_e` per harmonic, §6). The common cause is the
string itself — `r = R₃/37.1`, `ℓ₁ = 2πR₃/3` — which is sized for the electron and cannot be
sized for the proton without breaking `Z_e < Z₀`. **The proton is not a loop of the vacuum's
string under any reading the model offers.** What would change that is a second string species
or a second scale, neither of which the axioms contain.

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §4 (a fit is not a derivation).
