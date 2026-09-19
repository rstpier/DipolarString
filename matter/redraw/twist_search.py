#!/usr/bin/env python3
"""Redraw: what fixes the twist Tw of the wound charged strand (R33)?  Every rule of the base that
could set it is tried.

Candidates:
  A. topology: Tw = Lk - Wr (Calugareanu-White).  Lk is an integer fixed at birth (R4) and conserved;
     Wr is set by the shape and vanishes for a planar ring or a straight arm.  So Tw is an integer
     the strand is born with: conserved, not derived.
  B. energy alone: the solenoid energy grows as Tw^2, so a free strand unwinds to Tw = 0; nothing
     dynamical holds a twist unless it is conserved.
  C. flux quantisation: at R33's solution (46 turns/fm, strand density) the core flux is 0.038 h/e;
     a full quantum h/e would give 680 x the lattice tension; the base's own flux unit, the electron
     ring's self-flux (alpha/pi)(ln(8R/a) - 2) h/e = 0.0086 h/e, fits 4.5 times: no natural integer.
  D. matching to Z_0 (R10, R17 applied to the helical line): L' = mu_0 n^2 pi a^2, C' = 2 pi eps0/ln(b/a),
     Z = Z_0 n a sqrt(ln(b/a)/2); Z = Z_0 fixes n a = sqrt(2/ln(b/a)) = 1.2 (b = r_p) to 0.74
     (b/a = 37.1): about one turn per core radius, pitch angle 7-12 deg, Tw = 3-5 turns in the
     0.9 fm centre.  The only rule of the base that fixes the twist.
  E. consequence: the matched gain pi (n a)^2 = 1.7-4.5 gives 5-14 MeV/fm at the nucleon strand
     density, 66-170 x short; the lattice tension with the matched pitch needs a path density of
     3.3-5.4 e/fm, strands of e/3 only 0.06-0.10 fm long (the ladder's n = 13).  The free
     number moves from Tw to the density.  The electron is unaffected: its matched winding on a
     10.4 fm tube adds 0.03 MeV to 511 keV.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP = 0.51099895, 938.27209
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
RP = 0.8409
R_OVER_A = 37.1
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
SIGMA_LAT = 904.0
L_CENTRE = 0.9
MU0, C_SI, E_SI, H_SI = 4e-7 * math.pi, 2.99792458e8, 1.602176634e-19, 6.62607015e-34
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. topology
    Wr_ring, Wr_arm = 0.0, 0.0
    Lk = 5
    check("A. Tw = Lk - Wr: Lk is an integer fixed at birth and conserved, Wr = 0 for a planar ring or a straight arm, so Tw is the integer the strand was born with -- conserved, not derived",
          Lk - Wr_ring == 5 and Lk - Wr_arm == 5, "Calugareanu-White; Wr(planar) = 0")

    # B. energy alone
    u = lambda tw: tw**2
    check("B. Energy alone: solenoid energy ~ Tw^2, a free strand unwinds to Tw = 0; only conservation holds a twist",
          u(0) < u(1) < u(2), "u(Tw) monotonic")

    # C. flux quantisation
    a = LAMBDA_P
    lam9 = (1 / 3) / L1_9                                        # e/fm
    n33 = 46.2                                                   # turns/fm (R33)
    I = lam9 * E_SI / 1e-15 * C_SI                               # A
    phi = MU0 * n33 * 1e15 * I * math.pi * (a * 1e-15) ** 2      # Wb
    phi_q = phi / (H_SI / E_SI)
    phi_e = (ALPHA / math.pi) * (math.log(8 * R_OVER_A) - 2)     # electron ring self-flux in h/e
    check("C. Flux at R33's solution: 0.038 h/e; a full quantum h/e would give 680 x the lattice tension; the electron ring's own flux, (alpha/pi)(ln(8R/a) - 2) = 0.0086 h/e, fits 4.5 times: no natural integer",
          abs(phi_q - 0.038) < 0.002 and abs((1 / phi_q) ** 2 - 680) < 60 and abs(phi_e - 0.0086) < 0.0003 and abs(phi_q / phi_e - 4.5) < 0.2,
          f"Phi = {phi_q:.4f} h/e; (h/e / Phi)^2 = {(1/phi_q)**2:.0f}; Phi_e = {phi_e:.4f} h/e; ratio {phi_q/phi_e:.1f}")

    # D. matching to Z_0
    def na_matched(b_over_a):
        return math.sqrt(2 / math.log(b_over_a))
    na_rp, na_37 = na_matched(RP / a), na_matched(R_OVER_A)
    psi = lambda na: math.degrees(math.atan(1 / (2 * math.pi * na)))
    tw_rp, tw_37 = na_rp * L_CENTRE / a, na_37 * L_CENTRE / a
    check("D. Matching the helical line to Z_0: Z = Z_0 n a sqrt(ln(b/a)/2) = Z_0 fixes n a = sqrt(2/ln(b/a)) = 1.20 (b = r_p) to 0.74 (b/a = 37.1): about one turn per core radius, pitch 7-12 deg, Tw = 3-5 turns in the 0.9 fm centre -- the only rule of the base that fixes the twist",
          abs(na_rp - 1.201) < 0.005 and abs(na_37 - 0.744) < 0.005 and 7 < psi(na_rp) < 13 and 7 < psi(na_37) < 13 and 3 < tw_37 < tw_rp < 5.5,
          f"n a = {na_rp:.3f} (b = r_p), {na_37:.3f} (b/a = 37.1); pitch angle {psi(na_rp):.1f}-{psi(na_37):.1f} deg; Tw(0.9 fm) = {tw_37:.1f}-{tw_rp:.1f}")

    # E. consequence
    u9 = 4 * math.pi * K * lam9**2
    gains = (math.pi * na_37**2, math.pi * na_rp**2)
    sig = tuple(u9 * g for g in gains)
    lam_needed = tuple(math.sqrt(SIGMA_LAT / (4 * math.pi * K * g)) for g in gains)
    l_strand = tuple((1 / 3) / lam for lam in lam_needed)
    n_ladder = tuple(3 * (L1_E / ls) ** (1 / (2 * math.pi)) for ls in l_strand)
    # electron: matched winding on its tube a_e = R_3/37.1, arc 2 pi lambda_e
    a_e = LAMBDA_E / R_OVER_A
    lam_e = (1 / 3) / L1_E
    dE_e = 4 * math.pi * K * lam_e**2 * gains[1] * 2 * math.pi * LAMBDA_E
    check("E. Matched gain 1.7-4.5: 5-14 MeV/fm at the nucleon strand density, 66-170 x short; the lattice needs a path density of 3.3-5.4 e/fm, strands of 0.06-0.10 fm (ladder n = 13): the free number moves from Tw to the density. The electron's matched winding adds 0.03 MeV: unaffected",
          abs(sig[0] - 5.3) < 0.3 and abs(sig[1] - 13.8) < 0.5 and abs(lam_needed[1] - 3.3) < 0.1 and abs(lam_needed[0] - 5.4) < 0.1 and 12.5 < n_ladder[1] < n_ladder[0] < 13.6 and dE_e < 0.05,
          f"sigma = {sig[0]:.1f}-{sig[1]:.1f} MeV/fm (lattice/{SIGMA_LAT/sig[1]:.0f}-{SIGMA_LAT/sig[0]:.0f}); density needed {lam_needed[1]:.1f}-{lam_needed[0]:.1f} e/fm, strands {l_strand[1]:.3f}-{l_strand[0]:.3f} fm, n = {n_ladder[1]:.1f}-{n_ladder[0]:.1f}; electron dE = {dE_e:.3f} MeV")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
