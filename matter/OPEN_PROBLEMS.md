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
`Z_p/Z_e ≈ 22.7`. A closed toroidal-helix geometry for the `N_DS = 27` soliton
(9 turns of 3 strings, mutual inductance `L ∝ N²`, series inter-turn capacitance
`C ∝ 1/N`) yields `Z_p/Z_e = N^{3/2} = 27` — **within 19 %**, a residual of
logarithmic order comparable to the `ln(8R₃/r)` corrections of the electron
derivation.

**Closes when:** the gap is closed from an internal matching condition, *without
fitting*. That would render the proton mass geometric.
**Must respect:** `CONSTRAINTS.md` §1 — the target is the ratio `Z_p/Z_e`, never
an absolute `m_p`; and §4 — a fit is not a derivation.
**Most contribution-shaped item in the inventory:** self-contained and quantified.

### 2. Equation of state of the `e/3` energy fluid
*Manuscript status: unspecified. Group: `charge-chirality`.*

Charge is not carried by the connectivity of the network (`CONSTRAINTS.md` §2).
If it exists here it resides in the constitutive content of the string, whose
equation of state the framework does not give.

**Closes when:** an equation of state is written that supplies a genuine
divergence source, and reproduces the three structural consequences already
established (mirror-pair creation, mechanical annihilation, P acting as C).
**Must respect:** §2. Also carries the declared baryon-asymmetry debt: strict pair
creation implies exact matter–antimatter symmetry at formation.

### 3. Stabilization of three-dimensional solitons
*Manuscript status: absent, tied to the muteness of the weave. Group:
`chiral-closure`.*

Without a chiral term the energy `E(λ) = Aλ³ + Bλ² + Cλ` is monotonic and a
localized configuration collapses (Derrick). The Dzyaloshinskii–Moriya carrier is
**excluded by more than thirty orders of magnitude** (`CONSTRAINTS.md` §3).

**Closes when:** a stabilization mechanism of a different kind is exhibited, or
the absence of stable 3D solitons is accepted as structural and its consequences
for the electron-as-quantum reading (§1) are drawn.
**Must respect:** §3. Do not restart from Dzyaloshinskii–Moriya.

### 4. Spin ½, `g = 2`, fermionic statistics
*Manuscript status: declared half-success. No dedicated group.*

The semi-classical magneton `μ = μ_B` is obtained; spin ½, `g = 2` and fermionic
statistics are not.

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
