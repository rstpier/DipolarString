# The author's spin mechanism — spin 1 broken in two

**Date:** 15 September 2026
**Stated by the author:** every particle comes from an assembly of spin-1 DQDs (symmetric); a chain
of 3 DQDs breaking in two produces a positron and an electron, with opposite rotation from one
breaking action, spin ½ each — spin 1 broken in two.
**Status:** CONSISTENT as wave arithmetic, with one testable consequence (the electron ring has
radius `R₃/2`); the fermion itself is still not derived.
**Script:** `../scripts/spin_from_breaking.py`, 8/8 PASS.

## 1. The arithmetic, and it is exact

Mother = closed ring of 3 DQDs = 6 strings, circumference `6ℓ₁`, travelling mode of winding 1
(spin 1). Daughters = two rings of 3 strings (three `−` branches → electron, three `+` → positron),
circumference `3ℓ₁`. At the mother's frequency **half a wavelength fits a daughter**: the daughter's
mode has winding ½ — the antiperiodic ring of the Hopfion note §5c, `L_z = ½` about its axis.
"Spin 1 → ½ + ½" is wavelength inheritance, and `ℏω` is identical on both sides (0.2555 MeV with
the manuscript's `ℓ₁`).

## 2. Two consequences

- **Energy.** The mother's quantum is one quarter of the pair's rest energy: at least 0.77 MeV
  must come from outside, as in real pair creation — and that photon carries its own angular
  momentum. The ½ + ½ is inherited, not forced by conservation.
- **Size.** With the manuscript's `ℓ₁` the daughter's winding-½ mode sits at 0.256 MeV, *half* the
  electron's mass, while the anchor (`2πR₃ = λ_C`) is the winding-1 mode. The arithmetic closes
  only if the electron ring has **radius `R₃/2`** (strings of `ℓ₁/2`): then mother = `m_ec²`
  (winding 1), each daughter = `m_ec²` (winding ½), pair threshold = mother + one quantum. The
  open-turn note reached the same `R₃/2` independently (`Z_e = 0.62 Z₀`, still an impedance well).
  This halves `ℓ₁`, the cell, the Brillouin cutoff (`3/π → 6/π m_ec²`) and `E_conf`: a
  recalibration of the manuscript, testable against everything that depends on `ℓ₁`.

## 3. What the mechanism needs and the model does not supply

1. **A ℤ₂ structure on a ring of three single branches.** An antiperiodic travelling wave needs a
   sign flip around the ring. The derived source (§5c) is the half-twist of the unordered *pair*;
   three single branches joined end to end have no pair, and their junctions have equal impedance
   (`Γ = 0`): no flip. Either the electron keeps a pair structure — a bundle, which §6 of the
   step-1 note found disfavoured against the end-to-end ring — or a new rule flips the sign at
   pole junctions.
2. **A rule that does not sort by parity.** If every head-to-tail junction flipped the sign, odd
   rings would be antiperiodic (fermions: `e` 3, `d` 5) and even rings periodic (bosons: DQD 2 —
   the author's premise — but also `u` 6 and `ν` 6, which is wrong). A junction rule cannot be the
   whole story.
3. **The fermion.** `L_z = ½` about one axis is a mode property. The `−1` under *any* `2π`
   rotation and the exchange sign live in a quantum superposition of the medium's configurations
   (`2026-09-15-wen-link.md`); the mechanism gives the value ½ of one component, inherited from
   the mother's wavelength, not the statistics.

## 4. What stays

The mechanism unifies three things the notes had separately: the mother as a neutral DQD
assembly (pair creation = unpairing), the daughters' opposite circulation (charge sign), and the
half-integer ladder of the daughters (the antiperiodic ring, the open turn). Its one hard
prediction — the electron ring at `R₃/2` — is worth carrying into the next manuscript revision as
a question, not a change.

## 5. The author's correction — an open-chain mother, `(+)———   ———(−)`, and the axis rule

Script: `../scripts/open_turn_stability.py`, 5/5 PASS. Corrections stated by the author: the mother
is not a closed ring but an **open chain of DQDs**, `=:=:=:=:=:=:=` (three DQDs are `===`); the
chain splits lengthwise into a `+` row and a `−` row; the fluids keep opposite momenta and pile
up at **opposite extremities**, `(+)———   ———(−)`; and **spin = number of symmetry axes / 2**:
the string `—` and the DQD `=` have two axes (spin 1), a chain with its charge at one end has one
axis (spin ½).

- **Arithmetic unchanged.** An open bifilar chain of 3 DQDs (length `3ℓ₁`, open ends) has the
  half-wave fundamental at `ℏc/2R₃ = 0.256 MeV`, and so has each daughter chain of 3 strings:
  half the electron mass, strings of `ℓ₁/2` (`R₃/2`) needed — as in §2.
- **The axis rule** reproduces photon 1, DQD 1, one-ended chain ½ — and gives 3/2 for a triangle
  of three strings, 2 for a square, infinity for a closed ring. It is a shape heuristic; spin is
  the double cover of the rotation group, and the rule says nothing about the `−1` under `2π` or
  the exchange sign.
- **A static charge at one end of a conductor cannot stay there.** In the model's own
  electrodynamics the fluid of an open chain either spreads along it (equipotential — then the
  chain has two axes and the rule gives spin 1) or oscillates as the half-wave mode, whose
  antinodes sit at *both* ends with opposite signs and average to a uniform net charge. A
  persistent one-way flow is impossible on an open chain. The picture `———(−)` needs a **trap at
  the pole** — the "reactive barrier of `Γ_pole` transferring momentum onto the walls" the
  manuscript mentions (l. 315), a nonlinearity it does not write down: the fluid's equation of
  state, item 2 once more.
- **Stability: the chain closes.** The free ends of a chain carry opposite poles and attract;
  the interaction energy falls monotonically as the gap closes — 3.3 keV released at `δ = e/3`,
  0.6 MeV at the mode's `δ ≈ 4e` — and every mode family loses energy on closing (Compton
  `3π → 2π`: 0.77 → 0.51 MeV; low `π → 0`: 0.26 → 0). Nothing in the model holds the gap open:
  the like net charges lose to the poles by the chaining margin (~45), and the classical gap mode
  is continuous (no topological obstruction). The one-axis object decays into the closed ring,
  the assembly the fluid prefers (step-1 note §6). That 0.6 MeV is of the order of the electron
  mass: either `δ ≈ 4e` is too large to be physical, or the closed ring is the deep ground state.

**What the corrected mechanism needs, in one line:** something that traps the fluid at a pole
and keeps a chain from closing — a nonlinearity of the fluid at the pole. V2.10 names it
(`Γ_pole` as a wall, `Π_sat`) and does not give it. Until it does, the model's electron is the
closed ring, with the consequences already recorded (`L_z = n + ½` only through a pair's
half-twist, no exchange sign).

## 6. Two further statements, tested

**"The EM fluid accumulates in the trailing pole, creating the active charge"**
(`../scripts/trailing_pole.py`, 6/6). As a motion effect it is excluded: a charge depending on
velocity is bounded by atom neutrality (`|q_p + q_e| < 10⁻²¹ e`, electron at `v/c = α`) to
`κ < 2×10⁻¹⁷` in `q(1 + κv²/c²)`. As an internal-flow effect it is the standing wave of the pole
reflections (`Γ_pole = 1/3`, SWR 2): symmetric at both poles, oscillating, zero time average — no
net charge at one pole. A net pile needs a one-way, nonlinear barrier at the pole.

**"Each junction produces a binding mass generating a displacement curvature"**
(`../scripts/junction_curvature.py`, 5/5, poles held in contact at `2r` while the angle varies).
A pole-to-pole junction of two like-charge strings binds by 6.3 keV (1.2 % of `m_e`, the
manuscript's `E_coh` scale) and is angle-blind to 5 %: what dependence exists prefers the straight
chain, with a stiffness of 0.05 keV/rad². The binding mass is real; it produces no curvature. The
author states the curvature is a simple mechanical consequence; what mass, what motion and what
force he means was not yet stated, so this computation may not be the intended one.

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §6 (the spin sector is a declared
half-success), §9.
