# Matter-sector work list

Drawn from the Open Problems inventory of the V2.10 manuscript, restricted to the
matter sector. Each entry gives the manuscript's status, what would close it, and
the constraint it must not violate.

Discussion groups for all of these exist at
[dipolarstrings.org/contribute.html](https://dipolarstrings.org/contribute.html).

---

## Live targets

### 1. Proton characteristic impedance
*Manuscript status: open derivation. Group: `proton-impedance`.*

The mass law `m ∝ N²Z` with `m_p/m_e = 1836` and `N_p/N_e = 9` requires
`Z_p/Z_e ≈ 22.7`. The manuscript's scaling estimate for the `N_DS = 27` soliton
(closed toroidal helix, 9 turns of 3 strings, `L ∝ N²`, series `C ∝ 1/N`) gives
`N^{3/2} = 27`, quoted as "within 19 %".

**Attempt (15 September 2026), negative:** `notes/2026-09-15-proton-impedance-closure.md`,
`scripts/proton_impedance_closure.py`, 16/16. Computed with the electron's own
electrodynamics (the regularised Neumann and Coulomb integrals reproduce `L_e`, `C_e`
to 0.02 %) and only the model's closures (wire length `27ℓ₁`, 9 turns, pitch = the
matching distance `D` or the smallest self-avoiding pitch `1.11 D`), the closed helix
gives `L_p/L_e = 20.7` (not 81: two turns at `D` have coupling 0.19) and
`Z_p/Z_e = 14` (series-capacitance reading) or `3` (self-capacitance reading, the
electron's). Within the validity of each reading **no pitch reaches 22.67**: series
`12–15` with a minimum of 12 near `p = 2.5 D`, self `1.5–3.2`; the series formula crosses
the target only where the turns are no longer adjacent. The "19 %" was the distance
between an estimate and the target, not the residual of a calculation; the
calculation lands below the target by a factor 1.6 to 7. The toroidal helix also
self-intersects at the matching pitch. A4 cannot select a pitch: it pushes toward
`Z₀`, the wrong direction. Coincidences recorded and not used: `(N−1)^{3/2} = 22.63`,
`27/2^{1/4} = 22.70`.

**Search for another object (15 September), negative for the whole class:**
`scripts/proton_object_search.py`, 9/9; note §5. The mass law's own object is a flat
coil of 204 turns of radius `R₃/22.67 = 1.64 r` — below two tube radii, outside the
thin-wire regime; built physically it gives `L_p/L_e = 5.8`. Loose objects (27-string
ring, three orthogonal rings) give `Z_p/Z_e = 1.5–1.7`. The tightest coherent
multilayer coil the tube allows — the inductance-maximising object for this wire —
gives at most `L_p/L_e = 72` and `Z_p/Z_e = 10.5`: **every object of 27 strings of
`ℓ₁` falls short by a factor 26 in `L` and 2.2 in `Z`.** Carrying `1836 L_e` as
inductance needs `N ≈ 190` strings (`L_max ∝ W^{5/3}`). The `1/ℓ` mode-energy reading
puts the proton loop at 0.21 fm (a factor 4 from the charge radius, the least
inconsistent reading) but is `E = hc/λ` until a closure fixes `ℓ_p/ℓ_e`.

**Mode-energy reading tested (15 September), negative:** `scripts/proton_mode_energy_reading.py`,
11/11; note §6. Its only closure, `ℓ = Nℓ₁`, makes the proton nine times *lighter* than
the electron (harmonic `~16 500` needed: a fit); the tube radius caps a fundamental at
`m/m_e ≤ R₃/r = 37.1` (proton excluded by 50, muon by 5.6); and admitting the 0.21 fm
proton loop needs `R₃/r ≥ 1836`, which makes `Z_e = 1.36 Z₀` — the electron stops being
an impedance well. Three readings of mass (inductance, charge energy, mode energy), three
failures: **the proton is not a loop of the vacuum's string under any reading the model
offers**; the string is sized for the electron and cannot be resized without breaking
`Z_e < Z₀`.

**Closes when:** the axioms supply a second string species or a second scale — the
catalogue's `N_p = 27` on the vacuum's string is unrealizable under every reading. Until
then the proton stays *Reparametrized*, and this item is no longer "the most
contribution-shaped": it is blocked at the level of the string.
**Must respect:** `CONSTRAINTS.md` §1 — the target is the ratio `Z_p/Z_e`, never
an absolute `m_p`; and §4 — a fit is not a derivation.
**Formerly the most contribution-shaped item in the inventory;** now blocked in principle
(see below), like the neutrino sector.

### 2. Equation of state of the `e/3` energy fluid
*Manuscript status: unspecified. Group: `charge-chirality`.*

Charge is not carried by the connectivity of the network (`CONSTRAINTS.md` §2).
If it exists here it resides in the constitutive content of the string, whose
equation of state the framework does not give.

**Step 1 of the matter topology (15 September):** `notes/2026-09-15-string-end-charge.md`,
`scripts/string_end_charge.py`, 11/11. The manuscript's string is a conductor carrying a
charged fluid (`±e/3` per branch), not a flux line; the divergence source the model *does*
have is the **unpaired branch** — a broken DQD — Coulomb `e/3` beyond `ℓ₁`. Complete DQDs
and their A7 chains carry only bound charge. The three structural consequences are
unpairing (0.47 keV, the manuscript's `E_coh`), re-pairing, and the exchange of forward and
return branches under inversion. What stays open is the core: *why* the fluid carries
`e/3` — `Γ_pole = 1/3` is a reflection coefficient and does not fix a charge.

**Closes when:** an equation of state is written for the fluid that fixes its charge at
`e/3` per branch and holds it there (the unpaired branch already supplies the divergence
source and the three structural consequences).
**Must respect:** §2. Also carries the declared baryon-asymmetry debt: strict pair
creation implies exact matter–antimatter symmetry at formation.

### 3. Stabilization of three-dimensional solitons — the Hopfion route
*Manuscript status: absent, tied to the muteness of the weave. Group: `chiral-closure`.
Direction recorded 15 September 2026: `notes/2026-09-15-hopfion-direction.md`. Closed for the
pair; for `N ≥ 3` the Faddeev field strength exists (step 3b) and the induced term is computed
(step 4); the local route is closed at the carrier (step 5): the weave has no charged excitation.*

Without a Derrick-evading term the energy `E(λ) = Aλ³ + Bλ² + Cλ` is monotonic and a localized
configuration collapses. The Dzyaloshinskii–Moriya carrier is **excluded by more than thirty
orders of magnitude** (`CONSTRAINTS.md` §3). The class that works without threshold is a term
quartic in first derivatives — the Faddeev term `(∂ᵢn × ∂ⱼn)²` — and **string bending does not
supply it** (it gives Frank `K₃`, which scales as `λ¹`; `CONSTRAINTS.md` §3a).

The known theory for a unit-vector field with that term is **Faddeev–Niemi**: stable knotted
solitons in 3+1 D classified by the Hopf invariant `π₃(S²) = ℤ` — a *linking number of closed
loops* — with fermionic quantization available (Krusch–Speight 2006). DS already has two fields
with target `RP²`: the DQD pole axis `n` (the nematic order of A7) and the string tangent `t`.

**Step 1 done (15 September), negative for V2.10 as written:** quadratic rod elasticity
coarse-grains to Frank terms only; the Faddeev term is the Maxwell term of the director's CP¹
Berry gauge field (verified identity); it arises when a U(1) phase coupled to that gauge field is
integrated out; and **the V2.10 telegrapher action contains no coupling between its phase `θ` and
any orientation variable** — its only variables are `Φ` and `θ`. The direction is located, not
closed.
**Step 2 done (15 September), negative for the pair:** the differential TEM mode of two
conductors is real and non-degenerate, its U(1) Berry connection is identically zero, and the
bifilar geometry supplies no continuous `θ`–`a[n]` coupling. **The Hopfion-via-Faddeev route is
closed for the pair.** Witness (fibre circular pair) confirms the method. Script:
`scripts/step2_bifilar_berry_phase.py`, 7/7. What the pair *does* supply is a **ℤ₂ holonomy**,
derived: half a turn flips the differential mode. With the author's decision that the pair is
unordered (15 September), this is physical: a closed loop with an odd number of half-twists makes
the mode antiperiodic, contributes `π` per half-twist to `θ_tot` (item 5) and carries
`L_z = n + ½` about its axis (item 4).
**Step 3a done (15 September), reopens the question for `N ≥ 3`:** three indistinguishable
conductors on an equilateral triangle carry a **degenerate doublet** of differential modes (E
irrep of C₃v, degenerate to `10⁻⁷`) with a **non-zero, constant, antisymmetric Berry connection**
`A_ab = ∓0.1858` — exactly the structure the pair lacked. Transport over the closed loop of the
unordered triple (`ψ: 0 → 2π/3`) is a rotation by 120.000° of the doublet, eigenvalues
`e^{±2πi/3}` — the **permutation part** of the holonomy. Script:
`scripts/step3_triple_holonomy.py`, 12/12.
**Step 3b done (15 September):** the connection and the permutation composed, and the curvature
computed. `a = 0.186` is the doublet's mean angular momentum `⟨L_z⟩` about the string axis
(photon: 1; pair: 0; grid-calibrated, `0.15–0.25` for `ε/r = 0.5–2`). Tilting the string has no
matrix element in the doublet, so on the string's frame bundle the connection is exactly
`a J ω_twist` and the **curvature is `a ×` the solid-angle form of the tangent — the photon's
Berry curvature scaled by `a`, i.e. `F = 2a f_μν[t]`: the Faddeev field strength of the tangent
field exists for the triple.** Closed-string holonomy, verified against Fuller/Călugăreanu on an
explicit curve: `H = R(2π[(1 − a)Lk + aWr])`, `Lk ∈ ℤ/3`. A planar ring twisted by ⅓ gives
`0.271`, not ⅓; **the G conjecture's ⅓ is reached exactly only when the triple closes by writhe
with zero twist** (such a curve exists). Script: `scripts/step3b_triple_curvature.py`, 17/17.
**Step 4 done (15 September), the coefficient:** on one closed string the induced term is exact
— the doublet's zero-point energy with the holonomy as twisted boundary condition,
`E_ind = −(πℏc/ℓ) B₂((1 − a)Lk + aWr)`, `B₂(x) = x² − x + ⅙`, in the manuscript's `ℓ` and the
unit of its `E_conf`. Derrick class `λ⁻¹`, writhe-dependent, **minimal at trivial holonomy**: a
cost of `πℏc/(4ℓ)` for the pair's half-twist and `+0.03–0.06` units for the triple's third, and a
static torque driving a `Lk = ⅓` triple toward `Wr ≈ −1.5`. One string induces **no local static
Faddeev density** (its worldsheet field strength vanishes for a static string). If the doublet
propagated between strings as a 3+1 D field, the coefficient would be
`1/e² = (a²/6π²) ln(Λ/m) = 5.8×10⁻⁴ ln(Λ/m)` per species (charge `Q = 2a` under the CP¹
connection; scalar one-loop), enough for Derrick — a Hopfion of size `~1/(e√κ)`. Script:
`scripts/step4_induced_faddeev.py`, 16/16.
**Step 5 done (15 September), the carrier:** the weave's differential modes do hop between
neighbouring pairs (`k = ∓1.4 %` at `ℓ_cell`, exact thin-wire result; the manuscript's own channel
is the node scattering at the orthogonal 3-DQD crossings) — but the carrier is **exactly neutral**:
a real non-degenerate mode, real hoppings for every relative orientation, phase 0 or π around any
closed loop. Berry charge `a(n) = ⟨L_z⟩` lives only on bundles of `n ≥ 3` strands (`0, 0.19, 0.33,
0.43, 0.51` for `n = 2…6`), i.e. on the particles, not on the vacuum; the electron's doublet can
leak into a neighbouring pair (`1.7 %`) but arrives neutral. **The local Faddeev–Hopfion route for
the weave director is closed at the carrier for V2.10 as written.** Script:
`scripts/step5_carrier.py`, 17/17.
**What is now open — and it is no longer a Hopfion:** (i) the *single closed bundle* — its
stability is tension against the `1/ℓ` mode energy, with the induced writhe torque (step 4) acting
against elastic coefficients the model does not give; (ii) the two postulates that would reopen
the local route, which V2.10 does not make: a weave of `N ≥ 3` bundles, or a physical (not merely
regularizing) crossing frame carrying a propagating doublet; (iii) `κ` and PVLAS, moot until one
of those is made. The
`L_z ∈ ℤ + ⅓` caveat is reduced, not closed: neither part of the phase is an SO(3) angular
momentum (permutation label + Berry phase; the twisted ring's conserved generator is a screw
combination), so there is no contradiction — and no `j = ½` either.
**Must respect:** §3, §3a, §9 — a *topologically protected* label cannot also be mixed, so a
protected Hopf charge is a quantum number, not a mass mechanism. Any coefficient obtained must
be confronted with the PVLAS bound before use.

### 4. Spin ½, `g = 2`, fermionic statistics
*Manuscript status: declared half-success. No dedicated group.*

The semi-classical magneton `μ = μ_B` is obtained; spin ½, `g = 2` and fermionic
statistics are not.

**Derived lead (15 September, conditional on the pair being unordered — decided):** an unordered
bifilar loop with an odd number of half-twists carries an antiperiodic differential mode,
`ψ(φ+2π) = −ψ(φ)`, hence `L_z = n + ½` about the loop axis and `e^{2πiL_z} = −1` — half-integer
angular momentum and double-valuedness under `2π`, the two properties spin ½ is made of
(`notes/2026-09-15-hopfion-direction.md` §5c). **Not yet spin ½:** a half-integer `L_z` about one
axis is not a `j = ½` representation of SO(3); `g = 2` and statistics are untouched. Compare
Williamson–van der Mark (1997) before going further. For the triple (`N_DS = 3`) the analogous
holonomy is `2π[(1 − a)Lk + aWr]` with `a = 0.186` (item 3, step 3b): a permutation label plus a
Berry phase, neither an SO(3) angular momentum, so the `L_z ∈ ℤ + ⅓` caveat is not an
obstruction — and not a derivation of `j = ½` either.

**Closes when:** fermionic statistics emerge from the medium rather than being
imposed. Until then the matter sector has produced a magnetic moment, not a
fermion.

---

## Conditional on inputs that are themselves open

### 5. Aspect ratio `R₃/r` and the loop phase `θ_tot`
*Group: `aspect-ratio`.* Two undetermined dimensionless inputs. `Z_e` depends
logarithmically on the first; the second fixes the offset of the radial mode
ladder. Neither follows from the matching condition, which constrains `D/r` only.
Any particle-sector result inherits their status.

**`θ_tot` partly derived (15 September, given the unordered-pair decision):** the ℤ₂ holonomy of
the nematic pair is one of the "topological phase shifts accumulated at the internal interfaces"
that `θ_tot` collects — a closed loop on which the pair axis makes `n_half` half-turns contributes
`π·n_half`. For odd `n_half`, `2βℓ + θ_tot = 2pπ` gives `βℓ = (p − ½)π` and the resonance-radius
law runs on `p − ½`. Still undetermined: the interface part of `θ_tot`, and `n_half` for a given
particle. For `N_DS = 3` the contribution is `±2π[(1 − a)Lk + aWr]` on the two circular doublet
states, `a = 0.186`, `Lk ∈ ℤ/3` (item 3, step 3b): it equals `2π/3` per third of linking only when
the linking is carried by writhe, not twist.
**A second, equivalent origin (15 September):** a *partially open turn* has `βℓ_s = π` at the
fundamental — two open ends give `θ_tot = π` exactly as one half-twist does
(`notes/2026-09-15-open-turn-closure.md`). The two readings are degenerate at the fundamental and
differ at the second harmonic (all harmonics vs odd multiples of `π`). Closure itself moves mode
energies by factors of order 1, algebraically; it is not a mass mechanism.

### 6. Local cell spacing `ℓ_cell`
*Group: `cell-spacing`.* Postulated equal to `ℓ₁` at rest — the manuscript calls
its derivation *the missing brick*. `Π_sat`, `ω_c`, the saturation coefficient
and hence `a = 1/Π_sat²` are conditional on it.

---

## Blocked in principle

### 7. Neutrino sector
*Manuscript status: **major**, first in the inventory. Group: `neutrino-sector`.*

Not merely open. The obstacle is a **sign in the exponent** (`CONSTRAINTS.md` §7):
every mechanism available here makes neutrino masses *directly* proportional to
the available scale, where the Standard-Model mechanism makes them *inversely*
proportional, and no inversion mechanism exists in the framework. Work here needs
an inversion mechanism first; without one, the three requirements (quantization
to exactly three levels, an opening→mass law reproducing the `Δm²` ratio ≈ 33, a
mixing mechanism) cannot be met.

### 8. Absolute mass scale
*Manuscript status: **closed / conceded**. Group: `mass-scale`.*

Closed by the conformal obstruction (`CONSTRAINTS.md` §1). Recorded here so that
it is not reopened by accident. `m_e` is an input.

---

## Out of scope

The weak and strong sectors, and any contact with a regime where P and C are
separately violated, lie outside the model's scope as stated. Cosmological
evolution likewise — the baryon asymmetry is recorded as a debt, not as a target.
