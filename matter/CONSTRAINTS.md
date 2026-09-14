# Constraints on any matter-sector attempt

Every item here is established, not assumed — in the V2.10 manuscript, or in an
audit recorded under `notes/` with a script under `scripts/` that reproduces it.
Each names what establishes it. A line of work that violates one of these is not
merely unpromising — it is refuted before it starts.

---

## 1. There is no intrinsic mass scale. Absolute masses are not a target.

A vacuum that is purely electromagnetic and impedance-matched is conformally
invariant and therefore possesses no intrinsic scale. **Four independent
analyses** converge on this:

1. the equation of state `P = u/3`;
2. the single-velocity character of the linear lattice;
3. the zero-point sum over the braided defect — a pure ultraviolet term
   proportional to the cutoff;
4. at the non-linear level, the soliton-to-quantum mass ratio of a transmission
   line is fixed by its characteristic impedance alone.

The fourth is the one that is not blind to non-perturbative structure: for a line
of impedance `Z` the relevant dimensionless coupling is `Z e²/2ℏ`, so a matched
substrate (`Z = Z₀`) gives `2πα`, and **any** soliton of that line — whatever
non-linear constitutive relation is assumed — weighs `O(1/α)` times its own
quantum. Strong coupling would require `Z ≈ 174 Z₀`, incompatible with the
transparency of the vacuum.

> `m_e` is not "not yet derived" in this framework. It is structurally an
> **input**, exactly as the corresponding Yukawa coupling is an input of the
> Standard Model; and the electron is to be read as the **quantum** of the line
> rather than as its soliton.

**What remains attainable: a spectrum of *ratios* given one scale.**

### 1a. Coiling a line does not raise the coupling

The characteristic impedance of a **line** (`√(L/C)` of a conductor geometry) and
the impedance of the **medium** (`√(μ/ε)`, which sets `α = e²Z₀/2h`) are distinct
quantities. Raising the former by coiling does not raise the coupling, so **no
confinement-type strong-coupling regime is reachable by winding a line** — only
by changing the medium itself.

---

## 2. Charge is not carried by the connectivity of the network

A stationary state bearing a net flux through a closed surface requires a genuine
divergence source within it. The braided defect — with either handedness, and at
equal cost — provides none.

If charge exists in this framework it must reside in the **constitutive content
of the string** (the `e/3` energy fluid), **whose equation of state is
unspecified**. That equation of state is the object to construct; a topological
argument alone will not produce charge.

Three structural consequences do hold, and are usable:

- the vacuum node being achiral (Theorem 4, `O_h`-invariant scattering), net
  chirality is conserved and zero, so chiral solitons can only be created in
  **mirror pairs**;
- annihilation is mechanical: two mirror windings untwist into the achiral
  vacuum, releasing achiral `N_DS = 1` solitons;
- spatial inversion exchanges the two winding senses, so **parity acts as charge
  conjugation** — the mirror image of an electron is a positron. This is
  invisible to the emergent electromagnetic sector, which respects P and C
  separately; any contact with the weak sector, where P and C are separately
  violated, lies outside the model's scope.

**Declared cosmological debt:** strict pair creation implies exact
matter–antimatter symmetry at formation, making the observed baryon asymmetry an
open problem of the framework, as it is of the Standard Model.

---

## 3. No stable three-dimensional solitons — and the obvious repair is excluded

The geometric chirality of the weave is real but **mute** in the electromagnetic
sector, confirmed on two independent fronts: the achirality of the scattering
matrix, and the degeneracy of the two mirror weaves in the electrostatics of a
charged braided defect.

A chiral constitutive term of the **Dzyaloshinskii–Moriya** type would be the
natural carrier, and it would simultaneously evade Derrick's theorem: for an
energy functional `E(λ) = Aλ³ + Bλ² + Cλ` in three dimensions, a finite
equilibrium size exists **iff** `|B| ≥ √(3AC)`; without the chiral term the
energy is monotonic and the localized configuration collapses.

**But** stabilizing an object at the Compton scale that way requires a vacuum
chirality of order unity at that scale, **excluded by bounds on the rotation of
astrophysical polarization by more than thirty orders of magnitude**.

> The muteness of the weave and the absence of stable three-dimensional solitons
> in this framework are **one constraint, not two**.

Do not start from the Dzyaloshinskii–Moriya route. If a stabilization mechanism
exists, it is a different one.

### 3a. Which class of term could work, and what the model does not have

Derrick scaling in three dimensions, where a term `∫(∂ᵏφ)² d³x` goes as `λ^(3−2k)`:

| term | scaling | stabilises? |
|---|---|---|
| potential | `λ³` | no |
| chiral, one derivative (DM) | `λ²` | **only if `\|B\| ≥ √(3AC)`** — a threshold |
| gradient | `λ¹` | no |
| **curvature, four derivatives** | **`λ⁻¹`** | **yes, for every `D > 0` — no threshold** |

`E = Aλ³ + Cλ + D/λ` gives `λ² = (−C + √(C²+12AD))/6A`, positive for every `D > 0`. A bending
rigidity would therefore evade Derrick *more robustly* than the excluded chiral term, and it is
not of DM type.

**But V2.10 contains no such term:** *bending* and *Skyrme* each appear zero times. What exists
is three different things, which must not be conflated with a string bending stiffness:

- **A5**, topological rigidity of the DQD — *"possesses no continuous deformation mode"*:
  discrete, not elastic, and it serves as the short-distance regularization;
- the **shear modulus of A7** — a collective, `S`-dependent property of the *medium*;
- the **phase stiffness `ρ_φ`** — undetermined, a stated target of the phase-lock problem.

Two cautions before pursuing this. A stiffness calculation was already rejected once, for
displacing `D` about `D₀`, as illegitimate under A5 — and it gave `χ_vac ~ 2×10⁻⁴`, eighteen
orders above the PVLAS bound. Any bending coefficient must be confronted with that bound *before*
being used to stabilise anything. And a free consistency check exists: the Genesis conjecture
already assumes closed loops with a *minimum circumference* `~ ℓ₁`, which is what a bending
rigidity would produce. A derived stiffness must reproduce it, or one of the two is wrong.

---

## 4. The particle catalogue is a parametrization, not a prediction

The surviving table retains only the three states within scope — the photon
(`N_DS = 1`), the DQD constituent of the vacuum (`N_DS = 2`) and the electron
(`N_DS = 3`), the model's unique calibration anchor. The manuscript labels the
table a *reverse-engineered parametrization, not a set of predictions*.

The muon, tau and the weak bosons (`W`, `Z`, `H`) were tabulated in earlier
versions and are **withdrawn**: reverse-engineered mass fits without predictive
content. **Reproducing a measured mass by choosing an index is that same failure
mode.** A new entry earns its place only if its index follows from something
other than the mass it is meant to explain.

---

## 5. Two dimensionless inputs of the particle sector are undetermined

- the aspect ratio `R₃/r ≈ 37.1`, on which `Z_e` depends logarithmically — this
  is why `Z_e` is filed *Calibrated* rather than *Derived*, and why the model
  counts **two** inputs and not one;
- the accumulated interface phase `θ_tot` of the single-valuedness condition,
  which fixes the offset of the radial mode ladder.

Neither follows from the matching condition, which constrains `D/r` only. Any
result that depends on either is conditional until they are fixed.

---

## 6. The spin sector is a declared half-success

The semi-classical magneton `μ = μ_B` is obtained, but **spin ½, `g = 2` and
fermionic statistics are absent**. A matter sector that does not produce
fermionic statistics has not produced matter.

---

## 7. The neutrino sector is blocked in principle, not merely open

Three readings all fail against table data, and the obstacle is structural: the
mechanism favoured in the Standard-Model literature makes neutrino masses small
by being **inversely** proportional to a large scale, whereas every mechanism
available in this framework makes them **directly** proportional to the available
scale. The mismatch is one of **sign in the exponent, not of magnitude**, and no
inversion mechanism is available here.

---

## 8. The connectivity sector cannot supply a charge-conserving doublet

*Established by audit of the S₃ reconnection attempt, 14 September 2026. See
[`notes/2026-09-14-s3-reconnection-audit.md`](notes/2026-09-14-s3-reconnection-audit.md);
reproduce with `scripts/audit_s3_reconnection.py`.*

Let `G` be the relabelling group whose action on a closure is gauge. For a three-string closure
the physical states of the connected sector are the `G`-orbits of the two 3-cycles, and the two
cases were computed:

- **`G` contains an odd permutation** → the two 3-cycles merge into one physical state. No
  doublet. (The manuscript puts us here: the Johns node is `O_h`-invariant, *reflections
  included*, so odd relabellings are gauge.)
- **`G` contains only even permutations** → the two stay distinct, and what distinguishes them is
  the **sense of traversal of the loop**. Under the framework's charge-as-chirality assignment
  that is the charge, since parity acts as charge conjugation (§2). Then

      [Q̂,H] = 0  ⟹  (q₊−q₋)⟨−|H|+⟩ = 0  ⟹  ⟨−|H|+⟩ = 0

  and any mass-mixing term between them vanishes identically.

**Either way there is no charge-conserving two-state mass doublet in the connectivity sector.** A
doublet requires a degree of freedom carried by neither the pairing nor the winding sense.

---

## 9. Topological protection and mediated mixing are mutually exclusive

*Same audit.*

A recurring proposal is to supply the missing two-state degree of freedom `η = ±` as a
*topologically protected* twist, framing or self-linking parity. That is self-defeating:

- if `η` is topologically protected, no **local** operator changes it, so it commutes with any
  local reconnection operator, the Hamiltonian is block-diagonal in `η`, and the mixing amplitude
  `⟨−|H_eff|+⟩` vanishes **at every order**;
- if the mixing amplitude is non-zero, `η` is changed by a local operator and is therefore **not**
  protected — and an unprotected two-state degree of freedom has no symmetry reason to be
  degenerate at zeroth order, so the diagonal term generically survives and destroys the
  Dirac-like form.

The mechanism needs `η` protected and unprotected at once. Only an approximate symmetry with
explicitly controlled breaking could thread this, and it has to be exhibited, not assumed.

Combined with §3, this closes the chirality-like candidates entirely: the weave chirality is mute
in the electromagnetic sector, so it does not couple; and anything that does couple is not
protected. **Any `η` is a new field.**

---

## What is *not* closed

The obstruction above bites on absolute scales and on the routes named. It does
not forbid:

- a spectrum of **mass ratios** given the one anchor;
- closing `Z_p/Z_e` from an internal matching condition — the toroidal-helix
  geometry gives `N^{3/2} = 27` against `22.7` required, **19 % away**, a
  residual of logarithmic order comparable to the `ln(8R₃/r)` corrections of the
  electron derivation. Closing it *without fitting* would render the proton mass
  geometric;
- constructing the equation of state of the `e/3` energy fluid (§2), which is
  the declared home of charge;
- a stabilization mechanism for three-dimensional solitons that is **not** of
  Dzyaloshinskii–Moriya type (§3). The scaling in §3a says which class works: a four-derivative
  curvature term, which stabilises without a threshold. Deriving a **bending energy for a bifilar
  line at fixed `D`** from the telegrapher parameters is the one route that would supply it
  without a new field — bending the axis changes `L` and `C` at order `(D/R_c)²`, and the
  geometry is fixed, `D/r = 2 cosh π` being derived. It must be argued explicitly that bending
  the axis at fixed `D` is not the internal deformation A5 forbids.

These four are the live targets. `OPEN_PROBLEMS.md` states each one and what
would close it.
