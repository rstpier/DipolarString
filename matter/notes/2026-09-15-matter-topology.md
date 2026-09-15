# A matter topology for the DS medium — drawn from the theories that do more

**Date:** 15 September 2026
**Status:** CLOSED AS A MECHANISM (15 September, `2026-09-15-wen-link.md`): the DS string is a
conductor, not a flux line, and the vacuum is a classical network, not a string condensate — Wen's
mechanism has nothing to live in. Kept as vocabulary and as the record of what was tested. Below,
the original design; everything in it is POSTULATED unless marked otherwise. It reuses derived
pieces (steps 2–5, the open-turn resonance) and borrows its structure from two published
constructions that obtain fermions from a string medium. **Step 1 done (15 September,
`2026-09-15-string-end-charge.md`):** the DS string is a conductor carrying a charged fluid, not a
flux line; the charge is the *unpaired branch*, not a string end; pair creation is unpairing
(0.5 keV); the field sector confines nothing. Level 1 stands with that substitution, Level 3's
confinement falls. Figure:
`../figures/matter_topology.png` (`../scripts/draw_matter_topology.py`).

## 1. Which theories do more, and what they use

| theory | what it gets from a string/lattice medium | what it lacks |
|---|---|---|
| **Levin–Wen** (string-net condensation, PRB 2003, 2005) | Maxwell *and* fermions from a bosonic lattice: the string-net's wave is the photon, the **ends of open strings are the charges**, and with the "fermionic" string rule those ends are fermions (spin-statistics from the string's framing) | no masses; abstract lattice |
| **Bilson-Thompson** (helon model 2005; with Markopoulou–Smolin 2007) | the whole first generation as **braids of three ribbons**, each ribbon twisted `+1, −1, 0` ↔ charge `+e/3, −e/3, 0`; charge = sum of twists / 3; colour = position of the odd ribbon; chirality = handedness of the braid; **antiparticle = mirror braid** | no dynamics, no masses, no generations |
| Skyrme / Faddeev–Niemi | 3D stability from a quartic term; fermionic quantization from topology (Finkelstein–Rubinstein) | needs a field the DS weave does not carry (step 5) |

The DS medium already contains, by its own axioms, exactly the ingredients these two need: a
network of strings (the weave), a wave on the string (the photon, Theorem 4), **three strings
per electron**, **`e/3` per string** (`Γ_pole = 1/3`), the sign of charge from the **winding
chirality** of the string (the catalogue's methodological note), and the structural triple
"mirror-pair creation, mechanical annihilation, P acting as C" (`OPEN_PROBLEMS.md` item 2) —
which is word for word Bilson-Thompson's antiparticle rule.

## 2. The dictionary

| DS object | topological role | source |
|---|---|---|
| string (`N = 1`) | a helon: a ribbon-like line carrying chirality `+`, `−` or `0` | manuscript catalogue |
| winding chirality of a string | the helon's twist: charge `±e/3` or 0 | manuscript (postulated), `Γ_pole = 1/3` (derived) |
| unordered pair, half-twist | the ℤ₂ framing of a ribbon: `π` in `θ_tot`, antiperiodic mode | step 2, §5c (derived) |
| triple bundle, closed cyclically | the closure of a 3-braid with one cyclic permutation — one loop of `3ℓ₁`; transverse doublet, ℤ₃ permutation part, holonomy `2π[(1−a)Lk + aWr]` | steps 3a–3b (derived) |
| **unpaired branch** (a broken DQD) | **a charge `±e/3` — the divergence source; the end of a string is nothing** | step 1 (derived): `2026-09-15-string-end-charge.md` |
| the weave (closed, endless) | the neutral vacuum: no ends, no charge | Theorem 4; `CONSTRAINTS.md` §2, §8 |
| open turn at the end of a braid | the massive resonance: `βℓ_s = π` | `2026-09-15-open-turn-closure.md` (derived) |

## 3. The topology, drawn

**Level 0 — vacuum.** A closed string-net of DQDs. No ends, hence no charge: the vacuum is
neutral by construction, and `CONSTRAINTS.md` §2/§8 (connectivity carries no charge) becomes a
*feature*: closed networks cannot be charged, only string ends can.

**Level 1 — charge** (corrected by step 1). Unpairing a DQD frees two branches of opposite
charge: **mirror-pair creation** (0.5 keV); re-pairing them is **annihilation**; inversion
exchanges the forward and return branches: **P acts as C**. The three structural consequences of
item 2 are the three moves on branches. Charge magnitude `e/3` per branch is the fluid's content,
asserted to follow from `Γ_pole = 1/3` — the open core of item 2.

**Level 2 — fermions.** A **braid of three string ends**, its far end attached to the weave.
Charge = (sum of chiralities)/3; colour = position of the odd string; chirality of the particle =
handedness of the braid; antiparticle = mirror braid. First generation:

| particle | strands | charge |
|---|---|---|
| `ν` | `(0, 0, 0)` | 0 |
| `e⁻` | `(−, −, −)` | −1 |
| `u` | `(+, +, 0)` × 3 positions | +2/3, three colours |
| `d` | `(−, 0, 0)` × 3 positions | −1/3, three colours |

The DS electron of the manuscript — three strings end to end in one ring — **is** the closure of
the `e⁻` braid with a cyclic permutation: the same object, seen as a closed 3-braid. Its transverse
structure is the triple of steps 3a–3b.

**Level 3 — mass.** The free end of the braid is folded into a turn of size `R` — the author's
*partially open turn*. The mass is the mode energy of that end-resonance, `∝ 1/ℓ`, with
`βℓ_s = π` at the fundamental (derived). The electron's anchor fixes `ℓ_s`. Quark braids are open
at both ends; the note first read this as "an open braid must terminate on another braid" —
topological confinement — and a hadron as a **Θ-graph** of three braids (proton: `u, u, d`,
charge +1). **Step 1 removes the support for that:** unpaired branches are free charges in a
linear medium beyond half a cell; nothing in the field sector confines them. The Θ-graph stays as
a drawing, not a mechanism. Their masses are *not* in this topology: the model has no strong sector and cannot reach
hadron masses with the vacuum's string (`CONSTRAINTS.md` §1a, §1c). That limit is inherited, not
created, by the topology.

**Level 4 — spin and statistics.** Derived so far: the half-twisted (antiperiodic) turn carries
`L_z = n + ½` and changes sign under `2π` about its axis (§5c). Not derived: the full `j = ½`
representation, `g = 2`, and the exchange statistics of two ends. Levin–Wen's mechanism is the
reference: two string ends exchanged drag their strings, and the string-net rule decides the sign.
The belt-trick version (a `2π` rotation of the end puts a full twist in the attached ribbon, phase
`−1` for a spinor ribbon) needs a carrier with Berry charge `a = ½`; DS's bundles have
`a = 0, 0.19, 0.33, 0.43, 0.51` for `n = 2…6` — the ℤ₂ of the antiperiodic mode, not the belt
trick, is the derived source of the sign.

**Level 5 — photon.** The wave on a string with no end — DS native (Theorem 4). Bilson-Thompson's
boson triples are not needed.

## 4. What the topology explains that the manuscript only lists

- why charge comes in `e/3` and why the electron has three strings (a braid of three);
- the three structural consequences of charge (pair creation, annihilation, `P = C`) as the three
  moves on string ends;
- why the vacuum is neutral and why connectivity carries no charge (no ends);
- quark colour (position in the braid), quark confinement (open braids must join), the neutrino
  as the untwisted braid, antiparticles as mirrors;
- the same half-integer ladder from two readings (half-twist, open end) — item 5;
- where the weave's mute chirality could become audible: not in the vacuum but in the
  *handedness of a braid*, a property of the particle.

## 5. What it demands — the calculable steps, in order

1. **Is a DS string a flux line or a conductor?** The end of a *bifilar differential mode* is
   neutral (`+λ` and `−λ` terminate together); the end of a *single string carrying net line
   charge* is a charge. Levin–Wen's ends are charges because their strings are electric-flux
   lines of the emergent gauge field. The topology needs the DS string (`N = 1`, the "unfolded
   unit string") to be a flux line of the Theorem-4 Maxwell field with flux `e/3`; `Γ_pole = 1/3`
   says as much. **To derive:** the end field of a semi-infinite DS string in the weave, and
   whether Gauss gives `e/3`. This is item 2 (charge) restated as a computation.
2. **The end-resonance mass.** A braid terminated by a turn: stub resonance with the ABCD
   machinery (`open_turn_closure.py` extended to a stub on a line). Gives the ratio of the
   electron's `ℓ_s` to `R₃`, and the harmonic content that distinguishes the open turn from the
   half-twisted ring.
3. **Exchange statistics of two ends.** The DS node scattering matrix (Theorem 4) carries `±½`
   signs from flux conservation; Levin–Wen's fermionic string-net rule is a sign of the same
   kind. **To check:** whether the node's rule is the bosonic or the fermionic one. If fermionic,
   the ends are fermions by Levin–Wen's theorem — spin-statistics *derived*, the matter sector's
   declared half-success closed.
4. **`j = ½`.** From the antiperiodic mode (one axis) to the full rotation group: needs the
   response of the turn to tilts — the tangent-frame machinery of step 3b applies.
5. **What it cannot do:** hadron masses, generations, weak interactions — out of the model's
   declared scope, and the strong sector is not reachable by winding.

## 6. Status, honestly

A topology is a bookkeeping of what is conserved and how it moves; this one makes the manuscript's
charge rule, its three structural consequences, its electron and its neutrality of the vacuum into
one picture, borrowed from constructions that are known to yield fermions. It is POSTULATED. It
becomes physics at step 1 (are string ends charges?) and step 3 (is the node rule fermionic?) —
both computable inside V2.10, neither done. The proton's mass stays where §1c left it.

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §1c (no hadron mass from the vacuum's
string), §2 (a divergence source is required — here, the string end), §9 (a protected label is
not a mass mechanism).

## 7. Rings of mixed signs — the author's proposal for the quarks (n = 5, 6)

Script: `../scripts/ring_catalogue.py`, 8/8 PASS. Proposal: quarks as rings of `n = 5` (charge
`1/3`) and `n = 6` (charge `2/3`), each string carrying `±e/3`, joined end to end by their poles
(§6 of `2026-09-15-string-end-charge.md`).

- **Charge rule, derived:** a ring of `n` strings has `q = (n₊ − n₋) e/3`, so **the parity of `n`
  fixes the charge class**: odd rings carry odd thirds (`±1/3`, `±1`, …: `d`, `e`), even rings even
  thirds (`0`, `±2/3`, …: `u`, neutral states). The proposal `d = 5`, `u = 6` obeys it; so would
  `d = 3`, `u = 2` or `d = 7`, `u = 4`. Parity is derived, the value of `n` is not — it needs another
  criterion.
- **Stability:** every mixed ring is bound at every pole junction (arcs on `R_n = nℓ₁/2π`,
  `δ = q`), whatever the order of the signs; the lowest-energy order alternates the signs
  (`−+−+…`) as much as possible.
- **Arrangements:** distinct sign orders up to rotation and reflection — `u (6: 4+2−)`: 3,
  `d (5: 2+3−)`: 2, `d (3: 1+2−)`: 1. Three states for `u` only: the count is not colour.
- **Masses:** the inductance scaling `(n/3)² Z_e` gives `u (6) = 2.0 MeV` (PDG 2.2) and
  `d (5) = 1.4 MeV` (PDG 4.7), with `d` lighter than `u` — the wrong order; the mode-energy reading
  makes both lighter than the electron. Current-quark masses are scheme-dependent and dominated
  by binding: a weak target either way.
- **Counting:** proton `uud` = 6 + 6 + 5 = 17 strings in three rings, neutron `udd` = 16 — not the
  27 of the reparametrized catalogue, which §1c had already excluded as a mass carrier.

**Status:** the parity rule is DERIVED from "`e/3` per string, rings end to end"; the assignment
`5/6` is POSTULATED (parity-compatible, not selected). Colour and the quark masses are not
reproduced.

## References
- M. Levin, X.-G. Wen, *Fermions, strings, and gauge fields in lattice spin models*, Phys. Rev. B
  67, 245316 (2003); *String-net condensation*, Phys. Rev. B 71, 045110 (2005).
- S. O. Bilson-Thompson, *A topological model of composite preons*, hep-ph/0503213 (2005);
  S. O. Bilson-Thompson, F. Markopoulou, L. Smolin, *Quantum gravity and the Standard Model*,
  Class. Quantum Grav. 24, 3975 (2007). (Boson assignments and the twist convention — full or
  half twist per helon — to be read from the papers before use.)
- T. H. R. Skyrme (1961); D. Finkelstein, J. Rubinstein (1968); L. H. Kauffman, *Knots and
  Physics* (belt trick).
