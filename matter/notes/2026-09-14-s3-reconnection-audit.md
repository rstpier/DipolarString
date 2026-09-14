# Hostile audit — three-string S₃ closure and unitary reconnection

**Date:** 14 September 2026
**Audited:** `DS_matter_S3_reconnection_audit.md` + `verify_s3_reconnection.py` (matter-sector
attempt, produced in a separate AI collaboration)
**Reproduce:** `python3 ../scripts/audit_s3_reconnection.py` — 16/16 PASS, exit 0. Nothing below
is taken from the dossier; every quantity was recomputed from the definitions.
**Verdict:** (B) exploratory note only, as framed. (A) available under a rewritten *negative* claim.

---

## 1. What survives independent recomputation

All of it. The dossier contains no mathematical error.

| Claim | Independent result |
|---|---|
| 2 × [3], 3 × [2+1], 1 × [1+1+1] | confirmed |
| Cay(S₃, transpositions) = K₃,₃ | confirmed complete bipartite; spectrum {−3,0,0,0,0,+3} |
| One-step 3-cycle → 3-cycle forbidden | confirmed (parity flips each move) |
| [2+1] the unique intermediate | confirmed: exactly 3 two-step paths, all of type [2+1] |
| Class adjacency `[[0,√6,0],[√6,0,√3],[0,√3,0]]` | **re-derived from scratch** by projection onto the normalised class basis; identical. Spectrum {−3,0,3} |
| `E³−(Λ₃+Λ_μ)E²+(Λ₃Λ_μ−9g²)E+6Λ₃g²` | exact; symbolic difference = 0 |
| `E_low = −6g²/Λ_μ + O(g⁴)`, Λ₃-independent | exact; the g² coefficient factors as Λ₃(Λ_μa+6), so Λ₃ cancels |
| `Δ = (√(Λ_μ²+24g²)−Λ_μ)/4` | exact for the truncated block |
| `Δ = 3g²/Λ_μ − 18g⁴/Λ_μ³ + …` | exact |
| K_rec Hermitian, U_rec unitary | trivially correct (R² = I, R† = R) |

The dossier is honest about its own limits. What follows is what it does not say.

---

## 2. NO-GO A — the parity/charge dichotomy

*Stronger than the dossier's §7, which states only half of it.*

Let G ≤ S₃ be the relabelling group whose action σ ↦ πσπ⁻¹ is gauge. Physical [3]-sector states
are the G-orbits inside C₃ = {(123),(132)}. Computed both cases:

- **G contains an odd permutation** → 3 orbits, the two 3-cycles merge → one physical [3] state,
  no doublet.
- **G ⊆ A₃ ≅ C₃** → 4 orbits, the two 3-cycles stay distinct → doublet exists.

The dossier asserts the first without justifying G = S₃ rather than A₃. **Its §6 argument is
incomplete** — but the conclusion is rescuable, and the manuscript supplies the missing step: the
Johns node is O_h-invariant, *"reflections included (verified in the supplementary script)"*, so
odd relabellings are gauge.

Now the hard part. In the second case what distinguishes the two orbits is the **sense of traversal
of the loop**. The manuscript establishes: charge is carried by winding chirality; spatial
inversion exchanges the two senses; *"parity therefore acts as charge conjugation — the mirror
image of an electron is a positron."* So the doublet, if it exists, is charge-conjugate:
Q̂|±⟩ = ±q|±⟩. Then

    [Q̂,H] = 0  ⟹  (q₊−q₋)⟨−|H|+⟩ = 0  ⟹  ⟨−|H|+⟩ = 0  ⟹  Δ ≡ 0.

**For every choice of G, the connectivity sector yields either no doublet, or a doublet whose
mixing is forbidden by charge conservation.** η can be neither the pairing σ nor the winding sense.

---

## 3. NO-GO B — protection versus mixing

*Not in the dossier at all. Kills candidates 1–2 of its §7.*

If η is **topologically protected** — which is what "protected" means, and what would let it
survive the quotient — then no local operator changes it. K_rec is local by construction (a
two-pole exchange at a contact). Hence [η̂, K_rec] = 0, H_eff is block-diagonal in η, and
⟨−|H_eff|+⟩ = 0 **at every order**. Δ ≡ 0.

Conversely, if ⟨−|H_eff|+⟩ ≠ 0 then η is changed by a local operator and is **not** protected. An
unprotected two-state internal degree of freedom has no symmetry reason to be degenerate at zeroth
order; generically h_z ≠ 0, which the dossier's own §10 shows destroys the Dirac-like form.

**The mechanism needs η simultaneously protected and unprotected.** A framed/ribbon embedding or a
topologically protected twist cannot supply it. The only escape is an approximate symmetry with
controlled breaking, which would have to be exhibited.

---

## 4. NO-GO C — η cannot come from existing DS microstructure

Framing, twist and self-linking parity are all chirality-like. `CONSTRAINTS.md` §3 already
establishes that the weave chirality is real but **mute** in the electromagnetic sector — confirmed
on two independent fronts — and that its only named carrier, a Dzyaloshinskii–Moriya term, is
excluded by more than thirty orders of magnitude.

Mute ⟹ no coupling ⟹ Δ = 0 (NO-GO B again). Not mute ⟹ the excluded carrier.

**Any η is a new field.** This answers the dossier's own audit question 4 in the negative.

---

## 5. Findings the dossier misses

**F1 — the sign of Λ_μ is unconstrained.** §12 defines Λ_μ = E₂₁ − E₃ without fixing its sign;
§8 then treats it as a positive detuning. If Λ_μ < 0, the [2+1] class is the ground state and the
connected three-string object is not the ground state at all. → the cheapest decisive test; see §8.

**F2 — [1+1+1] is a different component-number sector.** The identity permutation is three
disconnected loops. Including it at amplitude g makes the low state a superposition of one
three-string loop and three one-string loops: a "particle" with no definite topology.

**F3 — the finite-state model does not evade the conformal obstruction, it relocates it.**
Δ = 3g²/Λ_μ has energy dimensions only because g and Λ_μ do, and a conformally invariant medium
fixes neither without an external scale. **No absolute mass can ever come out of this
construction**, however well g and Λ_μ are computed. The best attainable result is the
dimensionless Δ/E₁ as a function of Λ_μ/E₁ — which is exactly the shape of §13, and should be
advertised as the target rather than as a diagnostic.

**F4 — direct conflict with "the electron is the quantum of the line, not its soliton".** A
three-string closure is a soliton-like extended object. The manuscript's conformal-obstruction
section resolves the electron the other way. The dossier's question 9 raises this and never
answers it.

**F5 — §6 and §8 use different physical assumptions.** §6 quotients [2+1] to one class state,
forced by indistinguishability; §8 re-expands it into three |μ_a⟩ with two "dark" combinations.
Numerically consistent (√6 recovered) but logically not: if the quotient is a genuine gauge
identification the dark states do not exist; if they do, the quotient was too strong.

**F6 — the §13 diagnostic is worse than "not quantitatively justified".** Computed: the
perturbative formula overestimates Δ by **9.5×** at Λ_μ = 0.12 MeV and **2.8×** at 0.58 MeV. An
order of magnitude at one end. Good that it is flagged; the number should be given.

**F7 — §13 is a tautology.** Δ = m_ec² and E₁ = 3m_ec² both come from the anchor, so Δ/E₁ = 1/3 is
m_e/m_e: zero physical content. Its only content is the resulting g/Λ_μ, a statement about the old
*rejected* micro-loop numbers, not about the construction.

**F8 — the three-string restriction is never declared.** R_ij is right multiplication (exchange of
outgoing destinations); left multiplication would exchange incoming poles — isomorphic graph, but a
different operator, and a physical contact does one or the other. More importantly the model admits
**no reconnection that changes the number of strings**, which is what real string reconnection does.

---

## 6. Classification

**DERIVED** — S₃ representation; connected ⟺ 3-cycle; K₃,₃; one-step transition forbidden; [2+1]
the unique intermediate; class adjacency and its spectrum; characteristic polynomial; E_low and its
Λ₃-independence; exact and perturbative Δ; Hermiticity/unitarity of K_rec, U_rec.

**CONDITIONAL** — H_eff = (E₃ − 3g²/Λ_μ)I − (3g²/Λ_μ)τ_x and the exact truncated Δ, *given* a
physical doublet in a fixed charge sector, η-blind C₃-symmetric couplings, and Λ_μ > 0.

**POSTULATED** — existence of η; Λ_μ > 0; η-blindness of the couplings; the frozen three-string
sector.

**REJECTED** — η as the winding sense or the charge-conjugate pair (NO-GO A); η as a topologically
protected twist or framing (NO-GO B); η from existing DS microstructure (NO-GO C); e/3 from
Γ_pole = 1/3 or from spatial C₃; any absolute mass from this construction (F3).

**OPEN** — g; Λ_μ and its sign; spin ½, statistics; localization; charge; reconciliation with
"electron = quantum of the line".

---

## 7. Does V2.10 determine g or Λ_μ?

**No.** The word *reconnection* appears exactly once in the whole manuscript, in the Genesis
conjecture, which states that the reconnection dynamics is not derived. Missing, precisely:

1. a contact action between two bifilar strings during pole re-partnering — the manuscript has the
   telegrapher action for an isolated line and the Johns scattering matrix for the vacuum cell,
   neither of which describes two strings exchanging partners;
2. the 2×2 unitary that follows, giving γ = gτ_c/ℏ and hence g once τ_c is defined;
3. **an energy functional on closure topologies** giving E₂₁ − E₃, sign included. The manuscript
   assigns no energy to a closure as a function of its cycle type.

The dossier under-weights (3): it treats Λ_μ as merely uncalculated, but its sign decides whether
the object exists.

---

## 8. Corrected priority for the next calculation

The dossier's §18 puts g first. **Do Λ_μ first, with its sign.** If Λ_μ < 0 the branch falls
without a single matrix element being computed (F1). Cheapest and most decisive test available.

---

## 9. Publication verdict — (B)

The positive content is coursework: Cay(S₃, transpositions) = K₃,₃ is a standard exercise, and the
Schrieffer–Wolff elimination is a textbook technique correctly applied. The real content is
**negative**, and this audit hardens it: §8's conditional extension is not awaiting a derivation of
η, it is **blocked on every candidate η the dossier names**. Presenting it as a live route, as the
§20 claim does, is already an overclaim.

**(A) becomes available if and only if the claim is rewritten as a no-go:**

> Labelled closures of three oriented strings form S₃; connected closures are the 3-cycle class;
> elementary partner exchanges generate Cay(S₃, transpositions) = K₃,₃, so no elementary
> reconnection links two 3-cycles and every transition passes through the [2+1] class. The
> connectivity sector **cannot** supply a charge-conserving mass doublet: if the relabelling group
> contains an odd permutation the two 3-cycles are one physical state; if it does not, they differ
> by winding sense, hence — under the framework's charge-as-chirality assignment — by charge, and
> [Q̂,H] = 0 annihilates the mixing term. A doublet therefore requires a degree of freedom carried
> by neither the pairing nor the winding sense. No such degree of freedom exists in the present
> microstructure: its only chirality-like quantity is established mute in the electromagnetic
> sector, and its natural carrier is excluded by more than thirty orders of magnitude.

That closes a route cleanly with a four-line proof and directs the next work. It is not matter, and
the article must not use the word.

---

## 10. One opening, found while auditing

A **bending rigidity** would be a Derrick-evading term that is *not* of Dzyaloshinskii–Moriya type
— the escape NO-GO C leaves open. Computed scaling in 3D, where ∫(∂ᵏφ)²d³x goes as λ^(3−2k):

| term | scaling | stabilises? |
|---|---|---|
| potential | λ³ | no |
| chiral, one derivative (DM) | λ² | **only if \|B\| ≥ √(3AC)** — threshold |
| gradient | λ¹ | no |
| **curvature, four derivatives** | **λ⁻¹** | **yes, for every D > 0 — no threshold** |

`E = Aλ³ + Cλ + D/λ` gives `λ² = (−C + √(C²+12AD))/6A`, always positive. Strictly better than the
excluded chiral route.

But V2.10 has **no** such term: *bending* and *Skyrme* both appear zero times. What exists is A5
(topological rigidity of the DQD — discrete, not elastic), the medium's S-dependent shear modulus
(A7, collective, not a property of one string), and the undetermined phase stiffness ρ_φ. Three
different things, none of them a string bending stiffness. See `../CONSTRAINTS.md` §3.
