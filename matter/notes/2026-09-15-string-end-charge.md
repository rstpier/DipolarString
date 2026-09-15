# Step 1 of the matter topology — flux line or conductor? What is a charge in V2.10?

**Date:** 15 September 2026
**Status:** DERIVED from the manuscript's own definitions plus thin-line electrostatics. Answer:
**a conductor carrying a charged fluid, not a flux line.** The model's free charge is the
*unpaired branch* of a DQD, not the end of a string. Levin–Wen's mechanism survives with that
substitution; the topology's confinement does not.
**Script:** `../scripts/string_end_charge.py`, 11/11 PASS.

## 1. What the manuscript says the string is

- "a string containing an energy fluid of fractional charge `e/3`, possessing an intrinsic negative
  pole and positive pole" (l. 91); "the charge is produced by the accumulation of the energy fluid
  contained in the strings, and the polarity of the charge by the winding sense of the loop" (l. 248);
- the DQD is two anti-parallel branches with `Q_A = +e/3`, `Q_B = −e/3` at spacing `D₀`, attracting
  with `(e/3)²/(4πε₀D₀²)` (l. 286–289); each cell carries the maximal dipole `p_max = (e/3)D₀`
  (l. 962); the matching distance `D₀ = 2cosh(π) r` comes from balancing that attraction against
  the native linear power density (l. 339);
- `Γ_pole = 1/3` is the impedance reflection coefficient at the half-DQD interface,
  `(Z_bif − Z_half)/(Z_bif + Z_half) = (2 − 1)/(2 + 1)` (l. 302–306). Its identification with the
  charge `e/3` (l. 261: "the magnitude is fixed by the interface coefficient") is asserted: a
  dimensionless reflection coefficient does not fix a charge.

The branch is a transmission-line conductor: its TEM field is transverse, `E_∥ = 0` on a conductor.
**It carries no longitudinal flux.** So the string is not a Faraday tube; its charge is a *substance*
distributed along it.

## 2. Gauss, object by object (thin-line electrostatics, `e²/4πε₀ = 1.44 MeV·fm`)

| object | monopole | far field | free charge? |
|---|---|---|---|
| complete DQD (two anti-parallel branches of `ℓ₁`, `±e/3`, spacing `D₀`) | 0 | dipole `(e/3)D₀ = p_max` | no — its ends carry nothing |
| **unpaired branch** (`ℓ₁`, `−e/3`) | `−e/3` | Coulomb to 2 %, 0.3 %, 0.02 % at `R = 2, 5, 20 ℓ₁` | **yes** — the whole branch is the source |
| A7 head-to-tail chain of complete DQDs (pitch `ℓ_cell`) | 0 | bound charges `∓P` at the ends, `P = (e/3)D₀/ℓ_cell = 0.30 e/3`, reached as `1/R` (deviations 19, 11, 5, 3 % at 2, 4, 8, 16 cells) | no — a polarisation line; any sphere of complete cells encloses 0 |

**The end of a string is nothing.** The Levin–Wen "end of an open string" is realised in DS by the
*unpaired branch* — a broken DQD — whose whole length is the divergence source `CONSTRAINTS.md`
§2 requires. This is the reading under which the manuscript's electron (three strings, `−e`) is
Gauss-consistent: three unpaired negative branches.

## 3. Pair creation is unpairing, and it is cheap

Unpairing a DQD into two free `±e/3` branches costs its Coulomb binding: `U_bind = −0.47 keV`
(two anti-parallel lines of `ℓ₁` at `D₀`), against the manuscript's `E_coh = (2/9)αℏc/ℓ₁ = 0.40 keV`
(Theorem B). That is `10⁻³` of `m_ec²`: the creation threshold of a charge is not the electron
mass — consistent with `CONSTRAINTS.md` §1b (the mass law is a ratio law, not an energy budget).
The three moves of item 2 — mirror-pair creation, mechanical annihilation, `P` acting as `C` — are
unpairing, re-pairing, and the exchange of the forward and return branches under inversion.

## 4. No confining string between charges

Around a free `e/3` the DQD polarisation saturates (`p_max` per cell) only within
`R_sat = ℓ_cell √(ℓ_cell/4πD₀) = 0.52 ℓ_cell`; beyond it the medium is linear (`ε₀`) and the
potential is Coulomb. A saturated tube would need a cross-section of `3.4` cells² and carry
`0.46 eV/fm` — negligible and confined to half a cell. **The field sector gives no flux-tube
between unpaired branches**; the A7 chains are polarisation lines, not strings of the Levin–Wen
kind connecting charges. The topology's Level 3 ("an open braid must end on another braid",
topological confinement) has no support here; if quarks are unpaired branches, nothing in V2.10
confines them.

## 5. What this does to the topology note

- Dictionary: "end of a string = charge" → **"unpaired branch = charge"**; "cut a string" → "unpair
  a DQD"; the helon "0" is a *paired* DQD, the helons "±" are unpaired branches. String counts
  become `e: 3`, `ν: 6`, `u: 4`, `d: 5` (a different catalogue from the reparametrized one).
- Level 1 stands (charge, pairs, annihilation, `P = C`), at the price still of the assertion
  `Γ_pole = 1/3 ↔ e/3` — item 2 keeps its open core: *why* the fluid carries `e/3`.
- Level 3's confinement falls. Level 4 (statistics of two unpaired branches exchanged) is
  untouched and is the next computable step.

**Must respect:** `CONSTRAINTS.md` §2 (a divergence source is required: here the unpaired
branch), §1b, §9.
