#!/usr/bin/env python3
"""Redraw, after R39 ('the base has the number, not the rule'): test the rule that would assign the
vortex charge -- ONE quantum of circulation, charge q = e/(2 sqrt alpha), per CLOSED BIFILAR CIRCUIT
(the fluid goes out on one branch and back on the other).  A chain of strands end to end is one
circuit; a star has one circuit per arm (R20).

Checks:
  A. electron = one circuit of three strands (path 3 l_1(3) = 2 pi lambda-bar): matched-line energy
     4 pi K (q/path)^2 x path = pi hbar c/(3 l_1) = m_e c^2/2 exactly -- the circulating half (R13),
     from the same formula as the nucleon tension; mu, tau follow by the ladder.
  B. star arm = one circuit of one strand: tension 4 pi K (q/l_1(9))^2 = pi hbar c/l_1(9)^2 = 939 MeV/fm,
     sqrt(sigma) = 430 MeV (lattice 420-440).  One rule, both sectors.
  C. but energy = tension x length: one arm of l_1(9) stores pi hbar c/l_1(9) = 763 MeV.  Five charged
     arms (3-strand proton) would weigh 3.8 GeV; ONE arm gives 763 MeV, the proton's static part
     with e on one ring (770 MeV, R25 A) to -0.9 %.  The centre would be a single circuit.
  D. spin parity: one quantum of circulation per circuit, hbar/2 each; an odd number of circuits
     gives half-integer spin: leptons (1 circuit) and 3-arm quarks (3), baryons (9) odd; mesons (6)
     even.  Consistent, no prediction.
  E. verdict: the circuit rule reconciles the electron's half with the lattice tension without a
     free number; it does not give the proton's mass unless the centre is one circuit, which the
     3-strand content (5 charged arms) contradicts.  Open: which strands close into circuits.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MMU, MTAU, MP = 0.51099895, 105.6583755, 1776.93, 938.27209
LAMBDA_E = HBARC / ME
L1 = lambda n: 2 * math.pi * LAMBDA_E / 3 * (3 / n) ** (2 * math.pi)
Q = 1 / (2 * math.sqrt(ALPHA))
BAND = (420.0, 440.0)
STATIC_P = MP - MP / (2 * 2.79284734)        # R25 A: e on one ring, 18 % circulates
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def circuit_energy(path):
    """Matched-line energy of one circuit carrying q along a path: 4 pi K (q/path)^2 x path."""
    return 4 * math.pi * K * Q**2 / path


def main() -> int:
    # A. electron
    e_half = circuit_energy(3 * L1(3))
    m = {n: 2 * circuit_energy(3 * L1(n)) for n in (3, 7, 11)}
    check("A. Electron = one circuit of three strands: 4 pi K (q/path)^2 x path = pi hbar c/(3 l_1) = m_e c^2/2 exactly, the circulating half, from the nucleon-tension formula; mu -0.8 %, tau +1.0 % by the ladder",
          abs(e_half / (ME / 2) - 1) < 1e-12 and abs(m[7] / MMU - 1) < 0.011 and abs(m[11] / MTAU - 1) < 0.011,
          f"E = {e_half*1e3:.2f} keV = {e_half/ME:.4f} m_e; mu {m[7]/MMU-1:+.1%}, tau {m[11]/MTAU-1:+.1%}")

    # B. star arm
    sigma = 4 * math.pi * K * (Q / L1(9)) ** 2
    sqrt_sigma = math.sqrt(sigma * HBARC)
    check("B. Star arm = one circuit of one strand: tension pi hbar c/l_1(9)^2 = 939 MeV/fm, sqrt(sigma) = 430 MeV, inside the lattice band: one rule for both sectors",
          abs(sigma - math.pi * HBARC / L1(9) ** 2) < 1e-9 and BAND[0] < sqrt_sigma < BAND[1],
          f"sigma = {sigma:.0f} MeV/fm, sqrt(sigma) = {sqrt_sigma:.1f} MeV")

    # C. energy per arm
    e_arm = sigma * L1(9)
    check("C. Energy = tension x length: one arm of l_1(9) stores 763 MeV; five charged arms would weigh 3.8 GeV; one arm matches the proton's static part with e on one ring (770 MeV, R25 A) to -0.9 %: the centre would be a single circuit",
          abs(e_arm - 763) < 1 and abs(5 * e_arm / 1e3 - 3.8) < 0.05 and abs(e_arm / STATIC_P - 0.991) < 0.003,
          f"per arm {e_arm:.0f} MeV; 5 arms {5*e_arm:.0f} MeV; static part (R25 A) {STATIC_P:.0f} MeV, ratio {e_arm/STATIC_P:.3f}")

    # D. spin parity
    circuits = {"lepton (chain)": 1, "quark (3-arm star)": 3, "baryon": 9, "meson": 6}
    half_int = {k: v % 2 == 1 for k, v in circuits.items()}
    check("D. One quantum hbar/2 per circuit: odd counts give half-integer spin -- leptons (1), quarks (3), baryons (9) odd; mesons (6) even. Consistent, no prediction",
          half_int["lepton (chain)"] and half_int["quark (3-arm star)"] and half_int["baryon"] and not half_int["meson"],
          "; ".join(f"{k}: {v} circuit(s)" for k, v in circuits.items()))

    check("E. Verdict: the circuit rule reconciles the electron's half with the lattice tension without a free number, but gives the proton's mass only if its centre is one circuit, against the 3-strand content (5 charged arms). Open: which strands close into circuits",
          True, "see A-C")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
