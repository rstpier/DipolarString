# Direction — Hopf solitons of the director field

**Date:** 15 September 2026
**Status:** the Faddeev route is CLOSED for V2.10 (steps 1–2, both derived, 15 September); a ℤ₂
holonomy of the bifilar mode is DERIVED for an unordered pair and recorded as a lead. The rest is
POSTULATED.
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
Script: `../scripts/step2_bifilar_berry_phase.py`, 6/6 PASS.

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

Companion tests, cheaper and independent: the `l`-degeneracy of the Bohr-scale cavity mode
(§1.7); the Rydberg radius by Hartree self-consistency (§1.6).

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
