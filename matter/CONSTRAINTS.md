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

### 1b. The mass law is a ratio law, not an energy budget

"Mass is the macroscopic equivalent of the loop inductance" reads as if `m_e c²` were stored in
the loop. It is not, and the manuscript's own numbers say so: the inductive self-energy of the
electron loop at its own current is `½L_eI_e² = 2.19 keV` — **0.43 %** of `m_e c² = 511 keV` —
and the quantum of the loop's LC mode is `ℏ/√(L_eC_e) ≈ 101 keV`. Neither is the rest energy.
`m ∝ N²Z` fixes *ratios* between particles; the 511 keV is the anchor, an input stored nowhere
in the model. This is consistent with §1 — and the headline phrase should eventually be
corrected to say so.

### 1c. Inductance = mass cannot reach the proton with 27 strings

The mass law reads `m = κL` with `L` the inductance of the coiled string. Under the electron's
own thin-wire electrodynamics (the regularised Neumann and Coulomb integrals that reproduce `L_e`
and `C_e`), 27 strings of `ℓ₁` — wire length `9 × 2πR₃`, tube radius `r = R₃/37.1` — cannot be
wound into anything with `L > ~72 L_e` or `Z > ~10.5 Z_e`, against `1836` and `22.67` required
(`scripts/proton_object_search.py`, 9/9; `notes/2026-09-15-proton-impedance-closure.md` §5). The
object the law itself implies — 204 coincident turns of radius `1.64 r` — is outside the thin-wire
regime and gives `5.8 L_e` when built. **`m = κL` and `N_p = 27` are jointly unrealizable**; the
index would have to be `~190`, or mass is not inductance. Do not look for a cleverer shape.

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

Derrick scaling in three dimensions. Under `x → λx`, a term homogeneous of degree `k` in
**first** derivatives of a field scales as `λ^(3−k)`; a term quadratic in **second**
derivatives scales as `λ^(3−4) = λ⁻¹`:

| term | form | scaling | stabilises? |
|---|---|---|---|
| potential | `V(φ)` | `λ³` | no |
| chiral, one derivative (DM) | `φ·(∇×φ)` | `λ²` | **only if `\|B\| ≥ √(3AC)`** — a threshold |
| gradient, **including Frank elasticity** | `(∂n)²` | `λ¹` | **no** |
| **quartic in first derivatives (Faddeev)** | `(∂ᵢn × ∂ⱼn)²` | **`λ⁻¹`** | **yes, for every `D > 0` — no threshold** |
| quadratic in second derivatives | `(∂²φ)²` | `λ⁻¹` | yes |

`E = Aλ³ + Cλ + D/λ` gives `λ² = (−C + √(C²+12AD))/6A`, positive for every `D > 0`. A
four-derivative term evades Derrick *more robustly* than the excluded chiral term, and it is not
of DM type.

**Correction (15 September 2026).** An earlier version of this section said a *bending rigidity
of the string* would supply such a term. **It would not.** The bending energy of a string is
`κ² = |(n·∇)n|²` — quadratic in *first* derivatives of the tangent field. Coarse-grained over a
gas of strings it gives a Frank *bend* constant `K₃`, which scales as `λ¹` and does not
stabilise. The Faddeev term `(∂ᵢn × ∂ⱼn)²` is quartic in `∂n`, and string bending does not
produce it. What could, inside DS, is an open question with two candidates:

- the **twist** of the bifilar *ribbon* — a two-conductor line has a torsional degree of freedom
  in addition to bending, and Călugăreanu–White–Fuller (`Lk = Tw + Wr`) ties twist to the writhe
  of the centreline, a global geometric quantity; whether A5 forbids that twist as a "continuous
  deformation mode" must be settled first;
- **crossing interactions** between distinct strings, which are non-local in the string picture
  and could coarse-grain to a quartic gradient term.

**Two distinct vector fields.** The nematic order `Q_ij = S(n_in_j − δ_ij/3)` is built on the
**DQD pole axis** `n` — transverse to the string — not on the string tangent `t`. Both have
target `RP²` and `π₃(RP²) = π₃(S²) = ℤ`, so both admit Hopf solitons in principle; but Frank
elasticity of `n` and bending of `t` are different energies. The question "does DS contain a
Faddeev term" must name its field.

**V2.10 contains no such term, and its action cannot generate one:** *bending*, *Skyrme*,
*Faddeev*, *Hopf*, *knot* and *linking* each appear zero times, and the telegrapher action
`S = ½∫dt dx[C'(∂ₜΦ)² − (1/L')(∂ₓΦ)²]` is a 1+1 D scalar action whose only variables are `Φ`
and `θ` — it carries no orientation of the string. The Faddeev term is the Maxwell term of the
director's Berry gauge field (verified identity, `scripts/audit_s3_reconnection.py`); it would
require a coupling between `θ` and that connection, which the action does not contain — **and
which the bifilar geometry does not supply**: the differential TEM mode is real, its U(1) Berry
connection vanishes identically, and the holonomy over a full turn of the pair axis is 1
(`scripts/step2_bifilar_berry_phase.py`, 7/7). **The Faddeev route is closed for the pair.** What
the pair *does* supply is a ℤ₂ holonomy — a half-turn flips the mode — physical because the pair
is unordered (author's decision, 15 September); see `notes/2026-09-15-hopfion-direction.md` §5c.

**For `N ≥ 3` the Faddeev field strength exists; the Faddeev term is not yet derived**
(`scripts/step3_triple_holonomy.py`, 12/12; `scripts/step3b_triple_curvature.py`, 17/17). Three
indistinguishable conductors carry a degenerate doublet of differential modes with a non-zero
Berry connection `a J` per unit twist, `a = ⟨L_z⟩ = 0.186` (grid-calibrated; photon 1, pair 0).
Tilting the string has no matrix element in the doublet, so on the string's frame bundle the
curvature is exactly `a ×` the solid-angle form of the tangent: `F = 2a f_μν[t]`, the CP¹ field
strength of the identity above, on the *tangent* field. Integrating the doublet out
(`scripts/step4_induced_faddeev.py`, 16/16) gives, on one closed string, an exact Casimir-type
energy `−(πℏc/ℓ) B₂((1 − a)Lk + aWr)` — Derrick class `λ⁻¹`, writhe-dependent, and *minimal at
trivial holonomy*, i.e. a cost for the fractional windings — but **no local static Faddeev
density**: the doublet's worldsheet field strength vanishes for a static string. A local term
needs the doublet to propagate *between* strings, which V2.10 does not contain; if it did, the
coefficient would be `1/e² = (a²/6π²) ln(Λ/m) ≈ 5.8×10⁻⁴ ln(Λ/m)` per species, positive, hence
Derrick-sufficient for a Hopfion of size `~1/(e√κ)` with `κ` the open director stiffness. PVLAS
(below) is not evaluated for it. **And V2.10 has no such carrier** (`scripts/step5_carrier.py`,
17/17): the weave's differential modes hop between neighbouring pairs (`∓1.4 %` at `ℓ_cell`; the
model's own channel is the node scattering at the orthogonal crossings) but are exactly neutral —
real, non-degenerate, real hoppings, phase 0 or π on any closed loop. Berry charge exists only on
bundles of `n ≥ 3` strands (`a = 0.19, 0.33, 0.43, 0.51` for `n = 3…6`), which are the particles.
**The local Faddeev route is closed at the carrier for V2.10 as written**; it reopens only under a
postulate the model does not make (a weave of `N ≥ 3` bundles, or a physical crossing frame — which
the manuscript itself calls a statistical homogenization device). The closed-string holonomy is `2π[(1 − a)Lk + aWr]`, `Lk ∈ ℤ/3`; it is the G
conjecture's ⅓ exactly only for a triple that closes by writhe with zero twist. This section's
conclusion — no Derrick-evading term in V2.10 — stands; what changed is that the model's own
`N = 3` geometry contains the *field strength* such a term is built from.

**A5, as worded, is too strong for this.** "No continuous deformation mode" would forbid the
frame rotation the ℤ₂ and ℤ₃ holonomies require, and already forbids the axis bending the loop
sector requires. Its own justification (pole–pole contact regularisation, `ε₀` as an orientational
response without internal deformation) only needs **`D` fixed / no pole fusion**; that is the
reading adopted here, and the wording should be amended in the next revision. With the pair
unordered, the `±` poles are labels of the differential *mode*, not static properties of the
conductors. Against this reading, what exists in V2.10 is three different things, none of them a
quartic gradient term:

- **A5**, topological rigidity of the DQD — *"possesses no continuous deformation mode"*:
  discrete, not elastic, and it serves as the short-distance regularization;
- the **shear modulus of A7** — a collective, `S`-dependent property of the *medium*;
- the **phase stiffness `ρ_φ`** — undetermined, a stated target of the phase-lock problem.

Two cautions before pursuing this. A stiffness calculation was already rejected once, for
displacing `D` about `D₀`, as illegitimate under A5 — and it gave `χ_vac ~ 2×10⁻⁴`, eighteen
orders above the PVLAS bound. Any new elastic coefficient must be confronted with that bound
*before* being used to stabilise anything. And a free consistency check exists: the Genesis
conjecture already assumes closed loops with a *minimum circumference* `~ ℓ₁`. A derived
elasticity must reproduce it, or one of the two is wrong.

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
result that depends on either is conditional until they are fixed. Since 15 September one
contribution to `θ_tot` is derived, given the decision that the pair is unordered: `π` per
half-twist of the pair axis on a closed loop (`OPEN_PROBLEMS.md` item 5). The interface part
and the number of half-twists per particle remain inputs.

---

## 6. The spin sector is a declared half-success

The semi-classical magneton `μ = μ_B` is obtained, but **spin ½, `g = 2` and
fermionic statistics are absent**. A matter sector that does not produce
fermionic statistics has not produced matter. A derived lead exists since 15 September — an
unordered bifilar loop with an odd number of half-twists carries `L_z = n + ½` and is
double-valued under `2π` — but a half-integer `L_z` about one axis is not a `j = ½`
representation of SO(3); see `OPEN_PROBLEMS.md` item 4. This section stands.

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
  *estimate* gives `N^{3/2} = 27` against `22.7` required; the *calculation* with the
  electron's own electrodynamics and the model's closures gives `14` (series-capacitance
  reading) or `3` (self-capacitance reading), and within the validity of each reading no
  pitch reaches the target (`notes/2026-09-15-proton-impedance-closure.md`). The "19 %"
  was an estimate's distance, not a residual — and no other shape of 27 strings does
  better (§1c: at most `10.5 Z_e`, `72 L_e`). Closing it means changing the catalogue
  (`N_p`, or mass ≠ inductance), not the object;
- constructing the equation of state of the `e/3` energy fluid (§2), which is
  the declared home of charge;
- a stabilization mechanism for three-dimensional solitons that is **not** of
  Dzyaloshinskii–Moriya type (§3). §3a gives the scaling: the class that works is a term quartic
  in first derivatives (Faddeev) or quadratic in second derivatives. **String bending does not
  supply it** (it gives Frank `K₃`, `λ¹`). The DS-internal candidates are the twist of the
  bifilar ribbon and string-crossing interactions, and the field must be named — DQD pole axis
  `n` or string tangent `t`. **Closed for the pair** (no Berry connection, ℤ₂ holonomy only).
  **For `N ≥ 3` the field strength exists** on the tangent field — the doublet's Berry curvature
  is `0.19 ×` the photon's (step 3b); the induced coefficient is computed (step 4): exact and
  writhe-dependent on one string, `1/e² = (a²/6π²) ln(Λ/m)` in 3+1 D *if* a charged doublet
  crossed between strings — and the weave has none (step 5): its carrier is neutral, charge lives
  only on `N ≥ 3` bundles. **The local Faddeev–Hopfion route is closed at the carrier for V2.10 as
  written.** What stays open is the single closed bundle's stability (tension against the `1/ℓ`
  mode energy, elastic coefficients unknown). See `notes/2026-09-15-hopfion-direction.md` §5c–5g.

These four are the live targets. `OPEN_PROBLEMS.md` states each one and what
would close it.
