# Direction — Hopf solitons of the director field

**Date:** 15 September 2026
**Status (15 September, after step 3b):** the Faddeev route is CLOSED for the *pair* (steps 1–2,
derived). For `N ≥ 3` the Faddeev *field strength* EXISTS on the tangent field of the string:
the doublet's Berry curvature is `a ×` (solid-angle form of the tangent), `a = 0.19` (step 3b,
derived; `a` is grid-calibrated). The coefficient of an induced Faddeev *term* is OPEN. With the
author's decision that the pair is unordered: the pair's holonomy contributes `π` per half-twist
to `θ_tot` (DERIVED); the triple's holonomy is `2π[(1 − a)Lk + aWr]` (DERIVED) and equals the G
conjecture's ⅓ exactly only when the triple closes by writhe with zero twist. §5d's "exact ℤ₃
holonomy" was the permutation part alone (corrected in §5e). Spin, charge and 3D stability remain
OPEN. The rest is POSTULATED.
**Origin:** a discussion thread following the S₃ audit, starting from the working postulate
*there is no static mass value* and ending on a structure the literature already knows.

---

## 1. The thread, in the order it ran

Each step is a fact about the V2.10 manuscript or about known physics, not a new claim.

1. **No static mass** is not against the model: it already defines mass as loop inductance (a
   response to `dI/dt`), reads the electron as *the quantum of the line, not its soliton*, and
   measures in Appendix "breathers" the rest energy `E₀(ω)` of an object whose mass *is* a
   function of its internal frequency.
2. **Derrick's theorem applies to static solutions only.** Breathers and Q-balls evade it by
   time dependence. The postulate is the known loophole to blocker #1 (3D localization).
3. **L and C are the two faces of one object.** The DQD is a bifilar line: a distributed
   capacitor between its two conductors and a distributed inductor around the loop. The real
   mass law is `m ∝ N²√(L/C)`. Which one "is" mass is the choice of electromechanical analogy
   (impedance: `L`; mobility: `C`), and it dissolves in the dynamic picture where energy
   alternates between them at `ω = 1/√(LC)`.
4. **"C ↔ particle, L ↔ wave" is exact as number–phase conjugacy**: `[Q̂, Φ̂] = iℏ`. Charge
   number is countable (particle); flux phase interferes (wave). The manuscript already has the
   wave half — `θ = qΦ/ℏ` as compact canonical phase. The capacitive face of a *bifilar* pair
   carries `+Q/−Q`: a **dipole**, net zero — which is why connectivity gives no Gauss source.
5. **Two morphologies** — a polarized cloud when delocalized, structureless when confined — are
   half in the manuscript verbatim (*"a resonant cavity mode spread at the Bohr scale, its
   field Π carries the extension and the exponential decay of the orbital mode"*) and half
   forced by the closed size channel (386 fm loop vs `< 10⁻³ fm` measured: the mode's profile
   vs the quantum's interaction point). The transition between them is the measurement
   problem, open problem 3. The known analogue is the **polaron**; `χ_vac ≪ 1` would make the
   dressing thin, i.e. the transparency axiom predicts a near-point-like electron.
6. **Diffuse charge and volume.** Far field: Gauss, volume irrelevant, `1/r` for free — which
   the negative result requires. Self-energy: `E₀r₀ = e²/8πε₀`, *"an exact Coulomb identity in
   which α cancels"*; energy and size are locked in product. Consequences recorded in
   `CONSTRAINTS.md` §1b: **neither the electrostatic self-energy nor the loop's inductive energy
   (`½L_eI_e² = 2.19 keV`, 0.43 % of `m_ec²`) is the rest mass.** And the Rydberg constant is
   missing exactly one number — a radius, `r₀ ∈ [0.27, 0.85]R₃` — which is the volume of the
   diffuse charge. A Hartree-type self-consistency (`ρ = e·|Π-cloud|²`, Poisson, iterate) would
   derive it. Hohenberg–Kohn guarantees the 3D density carries the whole ground state.
7. **Quantum numbers.** `n`: present (`n*`, structural). `l, m`: absent but free from a 3D cavity
   mode; **the accidental `l`-degeneracy of hydrogen is an immediate falsifiable test** of whether
   the DS potential is truly `1/r`. `s`: **cannot emerge from the present field content** — `Z`
   (scalar), `Ω` (vector), `Π`/`Q_ij` (vector/tensor) — because spin ½ is not a mode number but a
   property of how the field transforms; it needs a spinor field **or a topological double cover**
   (Finkelstein–Rubinstein). Pauli is the same missing piece.
8. **A fluid more fundamental than matter** has a theorem (Madelung 1927: Schrödinger *is*
   hydrodynamics), a specific signature (quantized circulation from single-valued phase — the
   manuscript's own quantization mechanism), a lineage the manuscript cites five times
   (Volovik), and a proof that a *bosonic* field can make fermions: **Skyrme**, quantized by
   Finkelstein–Rubinstein and Witten, stable in 3D thanks to a four-derivative term. The
   cost: for `N > 1` the Madelung fluid lives in `3N` dimensions unless the field is
   fundamental and particles are its quanta — i.e. quantum field theory with a mechanical
   substrate, which has a rest frame (open problem 6).
9. **Extra dimensions: the minimum is zero.** The dimensions that do the work are *internal*
   (the field's target space), and DS has two with the right topology: the compact phase (a
   circle — charge by winding, Kaluza–Klein's mechanism without the spatial dimension) and the
   nematic director (`RP²`, `π₃ = ℤ`). "Opening a dimension locally" is, in known physics, an
   internal space whose size is a field — and DS already has such fields, `S(x)` with
   `S = S(T)`, `T ∝ ⟨Ω²⟩`, and `ρ_φ(x)`. Where they vanish, the internal space opens: **that is
   the core of a topological defect.**

---

## 2. The direction

**The DS director field admits knotted solitons.** `π₃(RP²) = π₃(S²) = ℤ`. The known theory is
Faddeev–Niemi (1997): unit-vector field, gradient term plus a term *quartic in first
derivatives*, stable solitons in 3+1 D (Battye–Sutcliffe 1998) classified by the Hopf
invariant — the **linking number of preimage loops** — with **fermionic quantization possible**
(Krusch–Speight 2006). Experimentally, stable Hopfions exist in chiral nematic liquid crystals
(Ackerman–Smalyukh 2017).

That is: closed linked loops, stable in 3D, with a topological charge that is an integer and a
route to spin ½, in a nematic field — with no extra spacetime dimension and no spinor field.
Every one of those is something the matter sector lacks, and they arrive together.

**The manuscript cites none of it.** Zero occurrences of Madelung, Skyrme, Faddeev, Hopf,
knot, linking, Finkelstein.

---

## 3. The correction that has to be made before building on this

The first statement of this idea (14 September, audit §10) said a *bending rigidity of the
string* would be the stabilising term. **That was wrong.** String bending is
`κ² = |(n·∇)n|²`, quadratic in first derivatives of the tangent; coarse-grained it is a Frank
bend constant `K₃`, scaling as `λ¹` — it does not stabilise. The Faddeev term
`(∂ᵢn × ∂ⱼn)²` is quartic in `∂n`. String bending does not produce it. Corrected in
`CONSTRAINTS.md` §3a, in the audit note, and in the verification script.

---

## 4. Two vector fields, not one

The nematic order `Q_ij = S(n_in_j − δ_ij/3)` is built on the **DQD pole axis** `n` —
*"statistical alignment of DQD poles"* — which for a bifilar line is **transverse** to the
string. The **string tangent** `t` is a different field. Both have target `RP²` and admit Hopf
solitons in principle; but Frank elasticity of `n` and bending of `t` are different energies,
and a quartic term for one is not a quartic term for the other. Any calculation must name its
field.

---

## 5. The decisive calculation — step 1 done, result negative for V2.10 as written

**(a) Quadratic rod elasticity gives Frank only.** A bifilar string with Kirchhoff bending and
twist has an energy quadratic in first derivatives of its tangent and frame. Any local
coarse-graining of a quadratic energy is quadratic in gradients at leading order: the Frank
constants `K₁, K₂, K₃`, all `λ¹`. No quartic term. (Homogeneity; degree check in the script.)

**(b) Where a quartic term comes from — verified.** Write `n = z†σz` with `z ∈ CP¹`. The Berry
connection `a_μ = −i z†∂_μz` has curvature `f_μν = ½ n·(∂_μn × ∂_νn)`, so the Faddeev term is
**exactly** `4 f_μν²`: the Maxwell term of the director's emergent gauge field. Verified
symbolically (`../scripts/audit_s3_reconnection.py`, "Faddeev density = 2 × Berry curvature"). A
Maxwell term for that gauge field is generated when a **U(1)-charged field coupled to `a[n]`**
is integrated out — the standard route in CP^{N−1} models.

**(c) DS has the U(1).** The compact phase `θ = qΦ/ℏ`. The mechanism therefore needs one thing:
a coupling between `θ` and the Berry connection `a[n]` — physically, the flux phase of the line
must advance when the ribbon's frame rotates. That is what the twist of a framed curve *is*
(Călugăreanu–White–Fuller: `Lk = Tw + Wr`, the twist being a U(1) connection along the curve).

**(d) The V2.10 action has no such coupling — and no variable that could carry it.** The
telegrapher action is

    S = ½ ∫ dt dx [ C'(∂ₜΦ)² − (1/L')(∂ₓΦ)² ]

a 1+1 D scalar action along the string coordinate. Its variables are `Φ` and `θ` only: no
tangent `t`, no director `n`, no `Π`, no `Ω`. Verified by inspection of the action block. **The
action does not know how the string sits in space.** (The words *framing*, *twist* and
*Kirchhoff* appear once each in the manuscript, all in unrelated senses.)

**Conclusion.** The quartic term cannot arise from the V2.10 action — not because A5 forbids
the twist, but because the action has no variable for it. The direction is not closed; it is
**located**: it requires one addition to the action, a coupling between the existing phase `θ`
and the existing director `n`. Not a new field — a new coupling between two fields already in
the model, whose minimal form is dictated by the geometry of a framed curve rather than chosen.
It is nonetheless an addition, and adding it is a physics decision, not a derivation.

**What the coupling would buy, if justified.** Integrating out `θ` gives the Faddeev term (3D
stability) as the Maxwell term of `a[n]`; the Hopf charge appears as linking; and, because the
*charge* phase is what couples to the geometry, a first bridge between the charge sector and
the topology — which nothing in DS currently provides.

**Step 2 — done. Negative for the Faddeev route; one unexpected positive.**
Script: `../scripts/step2_bifilar_berry_phase.py`, 7/7 PASS (regularised smooth fields, shared
machinery in `../scripts/berry_holonomy_common.py`).

The differential TEM mode of two conductors at `±(D/2)m`, `D/r = 2 cosh π`, is a **real**
field profile. For a real normalised family, `⟨E|∂_ψE⟩ = ½ ∂_ψ⟨E|E⟩ = 0`, so the U(1) Berry
connection `A = i⟨E|∂_ψE⟩` is **identically zero** (numerically `< 10⁻¹³` at four angles) and
the holonomy over a full turn of the pair axis is `1`. **The bifilar geometry supplies no
continuous `θ`–`a[n]` coupling at the TEM level. The Faddeev term cannot be generated this
way. The Hopfion-via-Faddeev route is closed for V2.10.**

The method is not blind: the witness — the degenerate circular pair `e± = (eₓ ± i e_y)/√2` of
a round optical fibre (Tomita–Chiao) — acquires `e^{∓iψ}` exactly. A linearly polarised real
mode is `(e₊e^{−iψ} + e₋e^{+iψ})/√2`: the two Berry phases cancel. That is *why* a mode locked
to a material frame has no anholonomy, and it is a theorem, not a feature of this geometry.

The only loophole is **beyond TEM**: higher angular modes `e^{±ilφ}` around the pair come in
degenerate `±` pairs and would carry a U(1) Berry phase like the fibre. They live at
`ω ~ c/D`, near the Brillouin cutoff, outside the telegrapher description. Whether DS reaches
them is a question about extending the model, not about V2.10.

**The unexpected positive: a ℤ₂ holonomy.** `⟨E(m)|E(−m)⟩ = −1` exactly: half a turn of the
ribbon swaps the conductors and flips the sign of the differential mode. If the two conductors
are *labelled* (`m` a vector, parameter space `S¹`) the half-turn is not a closed loop and
nothing happens. If they are **unordered** (`m ~ −m`, a director, parameter space `RP¹`) the
half-turn *is* a closed loop and its holonomy is **−1** — a ℤ₂ Berry phase, the Möbius bundle
over `RP¹`. The manuscript's own nematic order is exactly of this kind: *"compensated
head-to-tail alignment, null net polarization"*. Consequences, if the DQD pair is unordered:

- on a closed string loop with an **odd number of half-twists**, the differential mode is
  **antiperiodic**, `Φ(s+L) = −Φ(s)`, and the single-valuedness quantization shifts by a
  half-integer — the "fractional winding" the S₃ dossier reached for, now with a geometric
  origin;
- a mode that returns to minus itself around the minimal closed loop of its parameter space is
  the signature of **double-valuedness** — the structure spin ½ needs (Finkelstein–Rubinstein)
  — though the loop here is a twist of the *internal* frame, not a `2π` rotation of *space*;
  identifying the two would require `Lk = Tw + Wr` on a closed loop and is not done.

This is not the Faddeev term (which needs continuous Berry *curvature*), and it is not spin.
It is a derived geometric fact about a bifilar ribbon, conditional on the pair being unordered.
**Status: DERIVED (unordered pair) / OPEN (whether the DQD pair is unordered, and the link to
spatial rotation).**

## 5c. Decision (15 September): the pair is unordered, the string is nematic

Taken by the author. It makes the ℤ₂ holonomy physical, and the following are derived from it.

**A contribution to `θ_tot`, one of the two undetermined inputs.** The manuscript's
quantization condition is `2βℓ + θ_tot = 2pπ`, where `θ_tot` "collects the topological phase
shifts accumulated at the internal interfaces" and is listed (open problem 7) as not following
from the matching condition. The ℤ₂ holonomy *is* such a shift: a closed loop on which the
nematic pair axis makes `n_half` half-turns contributes `π·n_half` to `θ_tot`. For odd
`n_half`, `βℓ = (p − ½)π`: **half-integer winding**, and the resonance-radius law
`R = (p·ℓ₁/2π)(Z₀/Z)` runs on `p − ½`. Derived, given the decision. What remains undetermined
in `θ_tot` is the interface part and the value of `n_half` for a given particle.

**The ring chain.** A mode antiperiodic around a ring, `ψ(φ+2π) = −ψ(φ)`, is
`e^{i(n+½)φ}`: its angular momentum about the ring axis is `L_z = n + ½`, and it changes sign
under a `2π` rotation about that axis, `e^{2πiL_z} = −1`. Verified symbolically. So an
unordered bifilar loop with an odd number of half-twists carries **half-integer angular
momentum about its own axis and is double-valued under `2π`** — the two properties spin ½ is
made of. **Not yet spin ½:** a single half-integer `L_z` about one axis is not a spin-½
representation of SO(3) (rotations about other axes tilt the ring and must reproduce the full
`j = ½` spectrum); `g = 2` and fermionic statistics are untouched. The literature has this
picture — Williamson–van der Mark (1997), the electron as a photon on a Möbius-like double loop
— and it is worth reading before going further.

**A5 must be re-read.** As written, "no continuous deformation mode" forbids the frame
rotation the decision requires — and already forbids the axis bending the loop sector requires.
Its own justification (pole–pole contact regularisation, `ε₀` as an orientational response
"without internal deformation") only needs **`D` fixed / no pole fusion**. That is the reading to
adopt; the wording should be amended in the next revision.

**The poles become mode labels.** If the two conductors are indistinguishable, the `+`/`−` of
the DQD are not static labels of the conductors but the sign structure of the *differential
mode* — an excitation, not a property. That is consistent with the working postulate that
nothing here is static, and it reclassifies the `±` poles of A5 and of the S₃ dossier as
properties of the mode, not of the string.

**Mismatch, then resolved (§5d):** the G conjecture uses "fractional winding 1/3"; the ℤ₂
holonomy gives ½. See below.

## 5d. Step 3a — the triple (N_DS = 3) reopens what the pair closed

Script: `../scripts/step3_triple_holonomy.py`, 12/12 PASS.

The electron is `N_DS = 3`. Three indistinguishable conductors on an equilateral triangle have
**two** differential modes, and they form a **degenerate doublet** (the E irrep of C₃v) —
degenerate to `10⁻⁷`, orthogonal to `10⁻⁷`, with the common mode `2.96×` higher. A degenerate
doublet is exactly what the pair lacked and what the fibre has. Results:

- the non-abelian Berry connection `A_ij = ⟨E_i|∂_ψE_j⟩` is **antisymmetric, non-zero and
  constant** — `A_ab = ∓0.1858` per unit angle, an so(2) generator — where the pair's was
  identically zero;
- the overlap of the doublet with itself after a third of a turn, `ψ: 0 → 2π/3`, is a
  **rotation by 120.000°**, `det = 1`, eigenvalues `e^{±2πi/3}` to `10⁻¹⁵` — **the permutation
  part of the holonomy** (the three conductors exchanged). *Correction, step 3b:* this is not the
  full holonomy, which composes it with parallel transport by the connection; see §5e;
- the full turn is the identity; the common mode carries no phase.

**The ½ / ⅓ mismatch dissolves — for the permutation part:** `n` indistinguishable conductors on
an `n`-gon give a `ℤ_n` permutation and fractional winding `1/n`. The pair gives ½ (and, having
`a = 0`, nothing else), the triple gives ⅓ *plus* a geometric part weighted by `a` (§5e). **The G
conjecture takes "three strings each with fractional winding 1/3" as an input; the holonomy of the
indistinguishable triple reproduces it exactly when the triple closes by writhe with zero twist.**
Whether it is the *same* quantity the conjecture means must be checked against that section
before it is claimed.

**The Faddeev question is reopened for `N ≥ 3`, not answered.** A non-zero connection is
necessary for Berry *curvature*, not sufficient. What Faddeev needs is `f_μν ≠ 0` over a
two-parameter family — the tangent direction `t ∈ S²` of the string *and* the frame angle `ψ`
— i.e. the transport of the doublet along a curve whose tangent moves on the sphere, with the
material frame following the conductors. That is where Călugăreanu enters: the holonomy on a
closed loop would be `2π(Lk − Wr)` with `Lk ∈ ℤ/3`, and the writhe `Wr` is the continuous
geometric quantity whose variation is the curvature. **Step 3b: compute that curvature.** Done,
§5e.

**A caveat that must not be skipped.** A phase `e^{2πi/3}` on a *propagating* mode around a
ring would give `L_z ∈ ℤ + ⅓`, which is not a representation of SO(3) — only integer and
half-integer angular momenta exist in three dimensions. So either the ℤ₃ label is *internal*
(which of the two doublet states — an isospin-like charge, not an angular momentum), or a
triple with a third-turn twist is not a rotation eigenstate, or the twist is compensated by
writhe. In two dimensions a ℤ₃ exchange phase is anyonic; in three it needs resolving. **OPEN.**

Companion tests, cheaper and independent: the `l`-degeneracy of the Bohr-scale cavity mode
(§1.7); the Rydberg radius by Hartree self-consistency (§1.6).

## 5e. Step 3b — the curvature, and the honest holonomy

Script: `../scripts/step3b_triple_curvature.py`, 17/17 PASS.

**Correction to §5d.** Step 3a measured two things and never composed them: the overlap of the
doublet with itself after a third of a turn (a rotation by 120°, which is just the permutation of
the three conductors) and the Berry connection `A = aJ` per unit frame angle, `a = 0.1858`. The
holonomy of a closed string is the *product* of parallel transport by the connection and that
identification. Composed, a planar ring whose triple is twisted by one third of a turn has
holonomy `(1 − a) × 120° = 97.70°`, not 120°. The "exact ℤ₃ holonomy" of §5d was the permutation
part alone.

**What `a` is.** `a = ⟨L_z⟩` of the circular doublet state — its mean angular momentum about the
string axis — obtained identically by finite difference of the rigidly rotated pattern and by
applying the rotation generator to the pattern. A pure vector (the photon's polarisation) has
`a = 1`; a mode fully locked to the conductors has `a = 0`. The triple's doublet sits at
`a = 0.186`: 81 % locked to the conductors, 19 % free. It is a *model* number, not a topological
one: `a = 0.148, 0.186, 0.248` for `ε/r = 0.5, 1, 2`.

**The curvature — the object step 3b was defined to compute.** On the frame bundle of the string
(tangent `t ∈ S²` × material frame angle) the tilt generators `L_x, L_y` have *no* matrix element
inside the doublet (`4×10⁻¹⁸` against `0.186` for the twist; exactly zero by the symmetry of the
slab), so the connection is exactly `A = a J ω_z` with `ω_z` the twist 1-form, and by
Maurer–Cartan the curvature is

    F = −a J ω_x ∧ ω_y = a × (area form of S², pulled back by the tangent t).

That is the photon's Berry curvature (Tomita–Chiao) scaled by `a`. In the language of the CP¹
identity of the audit, `F = 2a f_μν[t]`: **the Faddeev field strength of the tangent field
exists for the triple, with weight `a`** — and does not exist for the pair (`a = 0`, a real
non-degenerate mode). The Faddeev *term* `(f_μν)²` would be induced by integrating out a charged
excitation of the doublet coupled to `F`; its coefficient needs the doublet's propagator and is
**not computed here**.

**The holonomy on a closed string.** With Călugăreanu `Lk = Tw + Wr` — verified numerically:
the parallel-transported frame returns rotated by `2πWr`, equal to the solid angle swept by the
tangent, both to `10⁻⁶` mod `2π` — the holonomy of the doublet around a closed unordered triple
is

    H = R( 2π [ (1 − a) Lk + a Wr ] ),   Lk = Tw + Wr ∈ ℤ/3,

checked on an explicit closed curve (chain 148.051° = formula 148.051°). It interpolates between
the pure permutation `2πLk` (`a = 0`, the ℤ₃ of §5d) and the photon `2πWr` (`a = 1`).
Consequences:

- **the G conjecture's "fractional winding ⅓" is reached exactly only when the triple closes by
  writhe with zero twist** (`Wr = Lk = ⅓`; such a curve exists — found by root-finding on a
  wobbling-ring family, holonomy 120.000°). A planar ring twisted by ⅓ gives `0.271`, not `⅓`.
  So if the conjecture's input is this holonomy, the electron's triple is a genuinely
  three-dimensional curve, not a twisted planar ring — which is also what a Hopf soliton needs;
- the pair is unchanged: `a = 0` exactly (real mode), holonomy `2πLk`, `Lk ∈ ℤ/2`: `π` per
  half-twist regardless of how it is split between twist and writhe;
- the contribution to `θ_tot` (§5c) for the triple is `±2π[(1 − a)Lk + aWr]` on `e±`.

**The `L_z ∈ ℤ + ⅓` caveat, reduced.** The phase has two parts and neither is an SO(3) angular
momentum of the ring: the permutation part is the internal label "which conductor" (81 %), the
geometric part is a Berry phase (19 %). A twisted ring is not rotationally symmetric about its
axis; the conserved generator is the screw combination (rotation about the axis + rotation of
the material frame), whose eigenvalues are not constrained to `ℤ/2`. No contradiction with SO(3)
— but also no `j = ½` and no `g = 2`; those remain OPEN.

**Status.** Curvature of the Faddeev form for the triple's tangent field: DERIVED (form and
weight `a`, grid-calibrated). Holonomy formula: DERIVED. Coefficient of an induced Faddeev term,
the PVLAS bound on it, and whether the vacuum weave's excitations are triple-doublet modes: OPEN.

## 6. Guard-rails

- `CONSTRAINTS.md` §1: ratios only. A Hopfion spectrum gives ratios with no adjustable index —
  falsifiable against measured mass ratios, the opposite of the withdrawn table. Do not promise
  it fits.
- §3: not Dzyaloshinskii–Moriya. The Faddeev term is not chiral; that is the point.
- §9: a topologically protected Hopf charge is a **quantum number**, not a mass mechanism. It
  is conserved, not mixed.
- Any elastic coefficient must be confronted with the PVLAS bound before use (the last
  stiffness calculation gave eighteen orders above it).
- Genesis already assumes closed loops with minimum circumference `~ ℓ₁`; a derived elasticity
  must reproduce that or contradict it.
- Charge is still not here. A Hopfion has a linking charge, not an electric one.

---

## References to add to the manuscript

- E. Madelung, *Quantentheorie in hydrodynamischer Form*, Z. Phys. 40 (1927) 322.
- T. H. R. Skyrme, *A non-linear field theory*, Proc. R. Soc. A 260 (1961) 127.
- D. Finkelstein, J. Rubinstein, *Connection between spin, statistics, and kinks*, J. Math.
  Phys. 9 (1968) 1762.
- L. Faddeev, A. J. Niemi, *Stable knot-like structures in classical field theory*, Nature 387
  (1997) 58.
- R. A. Battye, P. M. Sutcliffe, *Knots as stable soliton solutions in a three-dimensional
  classical field theory*, Phys. Rev. Lett. 81 (1998) 4798.
- S. Krusch, J. M. Speight, *Fermionic quantization of Hopf solitons*, Commun. Math. Phys. 264
  (2006) 391.
- P. J. Ackerman, I. I. Smalyukh, *Static three-dimensional topological solitons in fluid
  chiral ferromagnets and colloids*, Nature Materials 16 (2017) 426.
- P. Hohenberg, W. Kohn, *Inhomogeneous electron gas*, Phys. Rev. 136 (1964) B864.
