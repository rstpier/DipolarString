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

## 6. How polarised strings assemble — end to end, not side by side

Script: `../scripts/string_assembly.py`, 7/7 PASS. Question raised by the author: if the fluid is
polarised (`+` and `−` poles, manuscript l. 91), do strings assemble along their length? Model:
a string = net charge `q = −e/3` spread over `ℓ₁` plus pole charges `±δ` at its ends; poles touch at
`2r`, side-by-side spacing `D₀`.

```
  binding energies (keV; negative = bound), same net charge q = -e/3 unless noted:
   delta/q   end-to-end (pole contact)   side-by-side, opposite q (DQD)   side-by-side anti-par., same q   side-by-side parallel, same q
    0.00          +0.25                        -0.47                          +0.47                            +0.47
    0.25          +0.04                        -0.52                          +0.41                            +0.52
    0.50          -0.61                        -0.70                          +0.23                            +0.70
    0.75          -1.68                        -0.99                          -0.06                            +0.99
    1.00          -3.19                        -1.40                          -0.47                            +1.40
```

- **End to end, pole against pole**, like-charge strings bind as soon as `δ/q > 0.28`; at `δ = q` a
  junction binds by 3.2 keV. This is how chains and rings form.
- **Side by side**, only opposite charges bind without poles (the DQD, −0.47 keV). Like charges
  need `δ/q > 0.71` and bind ten times more weakly; a *parallel* like-charge pair never binds, so
  three strings on a triangle are frustrated.
- **Three like-charge polar strings joined end to end close into a ring** on `R₃` — bound at
  `δ = q`, unbound without poles: the manuscript's electron, three strings end to end, is the
  assembly the fluid prefers. A side-by-side bundle of three is not.

So yes: strings assemble along their length — into chains and closed rings by their poles — and
into a side-by-side pair only when the charges are opposite (the DQD). What the model does not
give is the pole strength `δ`: "intrinsic poles" are stated, not quantified; `δ/q ≥ 0.3` is the
condition for like-charge strings to chain at all.

**Must respect:** `CONSTRAINTS.md` §2 (a divergence source is required: here the unpaired
branch), §1b, §9.
