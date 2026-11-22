# Erratum to the Dipolar String Model manuscript — V2.11 (erratum edition, 21 November 2026)

Author: R. St-Pierre. Applies to: **V2.10** (frozen 7 September 2026, Zenodo concept DOI
10.5281/zenodo.21196098). The corrected text is `manuscript/DS_model_V2_11_erratum.tex` (and its
PDF); every change is marked *Erratum E*n in place and listed in the manuscript's Document History.
The frozen V2.10 remains the record of what was claimed before the erratum. Every item below is
reproduced by a script of this repository that prints a `RESULT` line with zero failures.

**All corrections reduce the model's claims.** None changes an equation as integrated in V2.10; two
of them (E6, E7) add computations the manuscript had listed as not performed, and both come out
against the model.

| # | Where (V2.10) | V2.10 statement | Correction | Evidence |
|---|---|---|---|---|
| E1 | Abstract, *Methods* | "Electric charge is dictated topologically by the winding chirality of the closed loop; mass follows from closed-loop inductance." | Withdrawn. Charge is not carried by the connectivity of the network (the manuscript's own open problem on charge); the inductive self-energy of the electron loop at its own current is 2.19 keV, 0.43 % of m_e c². Charge and mass are inputs. | `matter/CONSTRAINTS.md` §1b, §2 |
| E2 | §Electron impedance; "Closure of the size channel"; particle table | The electron is the first closed soliton, loose winding R₃/r ≈ 37.1, loop of 386 fm. | Withdrawn as an extended object. A charge e on a ring of radius R₃ has form factor j₀(qR₃), excluded by Compton, e–p and LEP data (G_E = 1 to 10⁻³ fm); a Bohr magneton on a charge loop needs R ≥ R₃; the loop's LC quantum (101 keV) fits neither m_e with its zero point nor the absence of excited electrons. The electromagnetic electron is Dirac's; Z_e is retained as the electrodynamics of the withdrawn object. | R100–R103, R121: `matter/redraw/ring_form_factor.py`, `internal_lc_oscillator.py` |
| E3 | Open problems, "Proton characteristic impedance — open derivation" | Toroidal-helix estimate within 19 % of Z_p/Z_e = 22.7; derivation open. | Blocked in principle: with the electron's own thin-wire electrodynamics, 27 strings of ℓ₁ reach at most ≈ 72 L_e and 10.5 Z_e against 1836 and 22.7 required; m = κL and N_p = 27 are jointly unrealizable. | `matter/CONSTRAINTS.md` §1c; `matter/scripts/proton_object_search.py` (9/9) |
| E4 | Status ledger, "Derived (Predictive)" | Γ_pole = 1/3 listed as derived and predictive. | Derived as a reflection coefficient of the pole geometry; its identification with the charge fraction e/3 is an identification, not a derivation. | `matter/CONSTRAINTS.md` §2 |
| E5 | Axiom A5 | "The DQD possesses no continuous deformation mode." | Too strong: it forbids the frame rotation that the ℤ₂/ℤ₃ holonomies require and the axis bending of the loop sector. Reading adopted: "D fixed, no pole fusion", which is all its justification needs. | `matter/CONSTRAINTS.md` §3a |
| E6 | Open problems, "Preferred-frame parameters not computed"; Theorem 3 discussion; abstract | α₁, α₂ not computed; standing against Einstein-aether bounds undetermined. | Computed. The scalar sector as written is a diagonal metric in the substrate frame; a moving mass carries no gravitomagnetic term, while the PPN form requires g₀ᵢ = −½(4γ+3+α₁−α₂)V_i − ½(1+α₂)W_i. Hence **α₁ = −(4γ+4) = −8** without a free parameter; α₂ = −1 (instantaneous potential) or 0 (Heaviside profile); α₃ = −3 or −4 with a rest-mass source. Two independent readings (substrate frame; comoving frame after an exact Lorentz boost); GR reproduced as control. Bound: \|α₁\| < 4×10⁻⁵ (binary pulsars): excluded by 2×10⁵. The only exit, the tensor bootstrap to GR, makes the frame unobservable. Open problem closed, negatively. The value is not new: it is the classical one for the standard "stratified" theories of gravity of the early 1970s, ruled out on this ground by preferred-frame effects (Ni 1972; Nordtvedt & Will 1972); the computation recovers it for the DS scalar sector. | R124: `matter/rebuild/alpha12_preferred_frame.py` (SymPy, 5/5) |
| E7 | §Substrate stars (table, Fig. branch, Fig. lrd), §Frozen collapses, Falsification hierarchy, §LRD "Registered prediction", abstract, conclusion | Stable branch with M_crit = 1.29×10⁶ M☉, r_s/R = 0.81; stable dark objects (r_s/R > 2/3); frozen collapses; falsification items 1–3; registered LRD pile-up prediction. | The structure equations omit the two TOV factors (M + 4πr³P_g/c²) and (1 − 2GM/rc²)⁻¹. Restored with every other convention kept (A8 source u − u₀, inertia u + P, surface P = u₀/3): stable branch with **M_crit = 3.9×10⁵ M☉ (×0.30)**, **maximal stable compactness 0.31** (each factor alone 0.40, 0.53), below 2/3 and below the radiation-sphere bound 0.49; the Newtonian 0.81 exceeded that bound. No stable dark configuration; surface redshift 1.20 for the most compact star; frozen collapses lose their mechanism (c_eff → 0 belonged to the scalar dynamics excluded by E6). Falsification items 1–3 and the registered prediction are withdrawn; M_crit ∝ ℓ₁² survives as a dimensional statement. | R125: `matter/rebuild/limit_states_tov.py` (5/5; control reproduces the V2.10 table) |
| E8 | §Emergent Lorentz, Brillouin cutoff | Photons above ~1 eV propagate as the N = 1 soliton, non-dispersive by construction. | Under the amplitude axiom (E12) the photon is the quantum of the lattice wave and carries its dispersion; the escape is withdrawn. LHAASO (GRB 221009A, 2026): E_QG,2 > 6×10⁻⁸ E_Pl = 7×10¹¹ GeV, against √12 ℏc/ℓ₁ = 0.85 MeV for the lattice: ℓ₁ excluded as lattice spacing by 10¹⁵; a photon-carrying lattice needs a spacing < 10⁻²⁷ m. ℓ₁ survives as the electron's Compton scale. | R119, R123: `matter/redraw/veille_21nov2026.py`, `a3_adoption_consequences.py` |
| E9 | (new structural limitation; electron section) | — | On the cubic weave, and on every pulse network of its type, a helicity-locked spin-½ excitation propagates at 1/3 of the link speed against 1/2 for the photon (ratio 2/3, excluded at 10⁻¹¹); the one lattice carrying both at one speed (Weyl automaton) needs a non-reciprocal T_h node with a physical quarter-wave, outside C1–C5 and against A5. | R105–R112 (`matter/redraw/BASE.md`) |
| E10 | (new; A8) | A8 postulated as non-gravitating native energy. | Under the amplitude axiom the zero-point energy of a quantized lattice at ℓ₁ is ≈ 9×10²² J/m³ = 1.7×10³² × the dark-energy density: A8 becomes mandatory, and it is the cosmological-constant problem under another name. | R123 |
| E11 | Status ledger | — | Re-filed in place: Γ_pole (E4); substrate-star branch and t_freeze ≃ Δt_echo/2 (E7). | — |
| E12 | Document History | — | Record of the author's decision (21 November 2026) to adopt an amplitude axiom, A9: a complex amplitude over the medium's configurations, linear evolution, Born rule. Its first consequences are E8–E10; its programme and stopping criteria are in `matter/rebuild/A3.md`. | `matter/redraw/BASE.md` (decision note), `matter/rebuild/A3.md` |

## What is unchanged

The V2.9.2 numerical suite (40/40) validates the equations as integrated in V2.9.2/V2.10; the erratum
does not alter those equations, it adds the TOV re-integration (E7) as a separate script and reads the
results differently. Theorem 4 (uniqueness of the Johns node under C1–C5, in the narrowed scope of
V2.9.2), Theorem 3 (covariance identity), γ = β = 1 for a system at rest in the substrate frame, and the
withdrawals already recorded in V2.9–V2.10 stand.

## What survives after the erratum

Identities about the representation (D/r = 2cosh π; Maxwell in one dimension rewritten in circuit
variables and the uniqueness of the node; the covariance identity; γ = β = 1; M_crit ∝ ℓ₁² as a
dimensional statement); combinatorics once the rules are posed (the charge spectrum and a colour label
from three strands; Koide's relation with a posed phase); two parameter-free hadronic post-dictions
resting on two posed integers (the strong part of m_n − m_p, 2.47 MeV vs 2.52 ± 0.29 lattice; √σ =
430 MeV, the Regge–Cornell value, 2–7σ below lattice determinations with dynamical quarks); and a set of
structural negative results, collected in the companion note `negative_results/DS_negative_results.pdf`.
No measured quantity is derived from inputs fixed elsewhere.

## Provenance

The computations behind E6–E10 were carried out between 19 and 21 November 2026 and are recorded step by
step in `matter/redraw/BASE.md` (R119–R125) and `matter/rebuild/`. They were written by the project's AI
collaborator working from the author's statements and decisions; they are not independent verification.
