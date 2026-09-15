# Direction — Hopf solitons of the director field

**Date:** 15 September 2026
**Status:** POSTULATED. A research direction, not a result. Nothing here is derived.
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

## 5. The decisive calculation (not done)

Coarse-grain the energy of the bifilar string gas into a director-field energy density and read
off what appears at order `(∂n)⁴`. At `(∂n)²` the Frank constants `K₁, K₂, K₃` will appear;
that is expected and does not stabilise. The question is the quartic coefficient. Two candidate
sources inside DS:

- **the twist of the bifilar ribbon.** A two-conductor line has a torsional degree of freedom
  as well as bending. Călugăreanu–White–Fuller, `Lk = Tw + Wr`, ties local twist to the
  writhe of the centreline — a global geometric quantity of exactly the self-linking kind the
  S₃ dossier reached for. **First question, before anything is spent: does A5 ("no continuous
  deformation mode") forbid the ribbon twist?** The manuscript's own particle sector needs bent
  strings, so A5 cannot mean "no curvature of the axis"; whether it means "no twist" is
  undecided;
- **crossing interactions** between distinct strings, non-local in the string picture,
  possibly local at order `(∂n)⁴` after coarse-graining.

If the quartic coefficient is non-zero: Hopfions exist in DS with no new term; 3D stability,
an integer linking charge and a route to spin ½ open together. If it is zero: the direction
closes cleanly and the framework has no internal stabilizer, which is itself a result.

Companion tests, cheaper and independent: the `l`-degeneracy of the Bohr-scale cavity mode
(§1.7); the Rydberg radius by Hartree self-consistency (§1.6).

---

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
