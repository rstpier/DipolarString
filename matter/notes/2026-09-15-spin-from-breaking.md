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

**Must respect:** `CONSTRAINTS.md` §1 (ratios only), §6 (the spin sector is a declared
half-success), §9.
