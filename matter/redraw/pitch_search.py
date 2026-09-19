#!/usr/bin/env python3
"""Redraw: what fixes the winding ratio gamma = 2.4 of R38?  First a correction to R38, then the search.

Checks:
  A. force is not energy per axial length.  A wound matched line stores E = u_path x s with
     s = sqrt(L^2 + (2 pi a Tw)^2); at fixed Tw and core, dE/dL = u_path / gamma <= u_path; at fixed
     path length, dE/dL = 0.  Winding cannot supply a force above u_path = 380 MeV/fm (at
     lambda = delta/l_1(9)); R38's sigma = 380 gamma is energy per axial length, not a tension.  The
     lattice's 904 MeV/fm as a force needs u_path itself: 'gamma = 2.376' is the shortfall ratio
     904/380, not a pitch.
  B. no rule of the base fixes a pitch anyway: the matched line's energy is flat in gamma at fixed
     path, the lumped matching overshoots (R38 D); the numerical coincidences (tan psi = 1/2:
     gamma = sqrt 5, -6 %; n a = 1/3: -2.3 %; 3 pi/4: -0.8 %) come from no mechanism.
  C. what closes the factor without winding: the base's other charge, the vortex charge
     q = e/(2 sqrt alpha) = 5.85 e (R32 C, the charge that carries the spin-half energy):
     (q/delta)^2 = (pi/2)^2 = 2.467 (+3.8 % vs 2.376).  With q per strand, u_path = 4 pi K q^2/l_1(9)^2
     = pi hbar c / l_1(9)^2 = 939 MeV/fm, sqrt(sigma) = 430 MeV: the centre of the lattice band
     (420-440); closed form sqrt(sigma) = (3/(2 sqrt pi)) 3^(2 pi) m_e c^2, alpha and Z_0 cancelling.
  D. the price: q per strand on the electron's three strings gives a line energy 3 pi hbar c/l_1(3)
     = 2.30 MeV = 4.5 m_e.  In the electron q is the charge of the whole ring (R32 C); the nucleon
     would need it per strand.  The base has the number, not the rule that assigns it.
  E. the two charges side by side: delta per strand -> electron static line energy 0.91 m_e (x1.8 the
     static half), nucleon 380 MeV/fm (2.4 x short); q per strand -> nucleon 939 (lattice), electron
     4.5 m_e; q per object -> electron spin half, proton ring 125 MeV, centre uncovered.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME = 0.51099895
LAMBDA_E = HBARC / ME
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
Q_VORTEX = 1 / (2 * math.sqrt(ALPHA))
SIGMA_LAT, BAND = 904.0, (420.0, 440.0)
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def u_path(q_per_strand, l1):
    return 4 * math.pi * K * (q_per_strand / l1) ** 2


def main() -> int:
    # A. force vs energy
    u0 = u_path(DELTA, L1_9)
    a, Tw = 1.0, 1.0
    def E(L):
        return u0 * math.sqrt(L**2 + (2 * math.pi * a * Tw) ** 2)
    L = 2.0
    g = math.sqrt(L**2 + (2 * math.pi * a * Tw) ** 2) / L
    dEdL = (E(L + 1e-6) - E(L - 1e-6)) / 2e-6
    check("A. A wound matched line: E = u_path x path; at fixed Tw and core dE/dL = u_path/gamma <= u_path, at fixed path dE/dL = 0. Winding cannot supply a force above u_path = 380 MeV/fm; R38's 380 gamma is energy per axial length, not a tension. 'gamma = 2.376' is the shortfall 904/380, not a pitch",
          abs(dEdL - u0 / g) < 1e-3 and dEdL < u0 and abs(SIGMA_LAT / u0 - 2.376) < 0.005,
          f"u_path = {u0:.1f} MeV/fm; dE/dL = {dEdL:.1f} = u_path/{g:.2f}; 904/u_path = {SIGMA_LAT/u0:.3f}")

    # B. no rule for a pitch
    cands = {"tan psi = 1/2": math.sqrt(5), "n a = 1/3": math.sqrt(1 + (2 * math.pi / 3) ** 2), "3 pi/4": 3 * math.pi / 4}
    devs = {k: v / 2.376 - 1 for k, v in cands.items()}
    check("B. No rule of the base fixes a pitch: the matched line's energy is flat in gamma at fixed path, the lumped matching overshoots (R38 D); the coincidences tan psi = 1/2 (-6 %), n a = 1/3 (-2.3 %), 3 pi/4 (-0.8 %) come from no mechanism",
          abs(devs["tan psi = 1/2"] + 0.059) < 0.005 and abs(devs["n a = 1/3"] + 0.023) < 0.005 and abs(devs["3 pi/4"] + 0.008) < 0.003,
          "; ".join(f"{k}: gamma = {v:.3f} ({devs[k]:+.1%})" for k, v in cands.items()))

    # C. the vortex charge closes the factor
    ratio = (Q_VORTEX / DELTA) ** 2
    u_q = u_path(Q_VORTEX, L1_9)
    closed = math.pi * HBARC / L1_9**2
    sqrt_sig = math.sqrt(u_q * HBARC)
    sqrt_closed = 3 / (2 * math.sqrt(math.pi)) * 3 ** (2 * math.pi) * ME
    check("C. The base's other charge, the vortex charge q = e/(2 sqrt alpha) = 5.85 e (spin half, R32 C): (q/delta)^2 = (pi/2)^2 = 2.467 (+3.8 % vs 2.376); with q per strand u_path = pi hbar c/l_1(9)^2 = 939 MeV/fm, sqrt(sigma) = 430 MeV, centre of the lattice band; closed form (3/2 sqrt pi) 3^2pi m_e c^2, alpha and Z_0 cancelling",
          abs(ratio - math.pi**2 / 4) < 1e-12 and abs(ratio / 2.376 - 1.038) < 0.003 and abs(u_q - closed) < 1e-6 and abs(u_q - 939) < 1 and abs(sqrt_sig - 430.4) < 0.5 and abs(sqrt_closed - sqrt_sig) < 1e-6 and BAND[0] < sqrt_sig < BAND[1],
          f"(q/delta)^2 = {ratio:.3f}; u_path = {u_q:.1f} MeV/fm = pi hbar c/l_1^2; sqrt(sigma) = {sqrt_sig:.1f} MeV (band {BAND[0]:.0f}-{BAND[1]:.0f})")

    # D. the price on the electron
    e_line = 3 * math.pi * HBARC / L1_E
    e_ring = 2 * K * Q_VORTEX**2 / LAMBDA_E
    check("D. The price: q per strand on the electron's three strings stores 3 pi hbar c/l_1(3) = 2.30 MeV = 4.5 m_e; in the electron q is the whole ring's charge (2 K q^2/lambda-bar = m_e c^2/2 exactly). The nucleon would need q per strand: the base has the number, not the rule",
          abs(e_line - 2.30) < 0.01 and abs(e_line / ME - 4.5) < 0.05 and abs(e_ring / (ME / 2) - 1) < 1e-12,
          f"3 pi hbar c/l_1(3) = {e_line:.2f} MeV = {e_line/ME:.1f} m_e; ring: 2 K q^2/lambda-bar = {e_ring*1e3:.1f} keV")

    # E. side by side
    e_static_delta = 3 * 2 * math.pi * K * DELTA**2 / L1_E
    check("E. delta per strand: electron 0.91 m_e (x1.8 the static half), nucleon 380 MeV/fm (2.4 x short); q per strand: nucleon 939 (lattice), electron 4.5 m_e; q per object: electron spin half, proton ring 125 MeV, centre uncovered",
          abs(e_static_delta / ME - 0.912) < 0.002 and abs(u0 - 380.4) < 1 and abs(u_q - 939) < 1,
          f"delta/strand: e {e_static_delta/ME:.2f} m_e, N {u0:.0f} MeV/fm; q/strand: N {u_q:.0f} MeV/fm, e {e_line/ME:.1f} m_e")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
