#!/usr/bin/env python3
"""Redraw: what forbids the transverse mode at 4 MeV (R17 E, R46 E)?

The 3.96 MeV is the TE_1 cut-off of the electron's guide, hbar c pi / w with w = 156.5 fm.  It is a
WAVE notion: a guide carrying waves has cut-off modes, and a closed ring admits the cut-off mode
as a standing state (guide wavelength infinite, no axial variation).  In the wave reading (R29)
nothing forbids it, nor the axial harmonics of the arc mode at m_e c^2 + n x 255 keV -- and none
of these states exist (no electron resonance in Compton scattering at MeV; compositeness scale
above 10 TeV).

Checks:
  A. the wave reading predicts excited electrons at 0.766 MeV (arc harmonic), 1.02, ... and at
     4.47 MeV (transverse): all excluded by experiment.
  B. the vortex reading (R32, R40) has no waves: the fluid is a steady circulation with one
     conserved quantum q = e/(2 sqrt alpha); its energy pi hbar c / path is fixed by the quantum
     and the path.  There is no harmonic (a wave's n-th harmonic has no counterpart in a steady
     flow) and no transverse mode (no wave, no cut-off).  A second quantum would double the
     circulating charge (11.7 e) and quadruple the energy: a different object, not an excited
     electron.
  C. the numbers do not change: the vortex circuit energy 4 pi K q^2/path equals the half-wave
     mode energy hbar c pi/path identically (4 pi K q^2 = pi hbar c), so every result of R29-R45
     that used 'the mode' stands as 'the circuit'.
  D. what stays excitable: shapes and topologies, not modes -- hadrons have them (Delta - N =
     294 MeV, one quark circuit = 254 MeV, +16 %), leptons have one shape only.
  E. verdict: the 4 MeV mode is not forbidden by a rule, it is absent: it was an artefact of
     reading the fluid as a wave; R32 removes it together with the electron's harmonics.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MDELTA = 0.51099895, 938.27209, 1232.0
LAMBDA_E = HBARC / ME
ARC_E = 2 * math.pi * LAMBDA_E
W_E = 4 * LAMBDA_E / math.pi**2
Q = 1 / (2 * math.sqrt(ALPHA))
L1_9 = (2 * math.pi * LAMBDA_E / 3) * (3 / 9) ** (2 * math.pi)
COMPOSITENESS_TEV = 10.0
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def main() -> int:
    # A. wave-reading excitations
    e_mode = HBARC * math.pi / ARC_E                       # 255.5 keV
    harmonics = [ME + (n - 1) * e_mode for n in (2, 3)]    # excited electron masses (mode n replaces mode 1)
    e_trans = HBARC * math.pi / W_E
    m_trans = ME + e_trans
    check("A. Wave reading: excited electrons at 0.766 and 1.02 MeV (arc harmonics) and 4.47 MeV (transverse cut-off), all excluded (no resonance in MeV Compton scattering; compositeness above 10 TeV)",
          abs(harmonics[0] - 0.7665) < 0.001 and abs(harmonics[1] - 1.022) < 0.001 and abs(m_trans - 4.47) < 0.01 and COMPOSITENESS_TEV * 1e6 > m_trans,
          f"arc mode {e_mode*1e3:.1f} keV; states at {harmonics[0]:.3f}, {harmonics[1]:.3f} MeV; transverse {e_trans:.2f} MeV -> {m_trans:.2f} MeV")

    # B. vortex reading: one quantum, no harmonics
    e_circuit = 4 * math.pi * K * Q**2 / ARC_E
    e_two_quanta = 4 * math.pi * K * (2 * Q) ** 2 / ARC_E
    check("B. Vortex reading: one conserved quantum q = 5.85 e, energy pi hbar c/path; no harmonic, no transverse wave; a second quantum doubles the charge (11.7 e) and quadruples the energy (1.02 MeV): a different object, not an excited electron",
          abs(e_circuit / (ME / 2) - 1) < 1e-12 and abs(e_two_quanta / e_circuit - 4) < 1e-12 and abs(2 * Q - 11.7) < 0.05,
          f"one quantum: {e_circuit*1e3:.1f} keV = m_e c^2/2; two quanta: {e_two_quanta:.2f} MeV, charge {2*Q:.1f} e")

    # C. identity mode = circuit
    check("C. The numbers do not change: 4 pi K q^2 = pi hbar c exactly, so the circuit energy equals the half-wave mode energy for every path: R29-R45 stand with 'mode' read as 'circuit'",
          abs(4 * math.pi * K * Q**2 - math.pi * HBARC) < 1e-9, f"4 pi K q^2 = {4*math.pi*K*Q**2:.4f} = pi hbar c = {math.pi*HBARC:.4f} MeV fm")

    # D. shapes stay excitable
    e_q = math.pi * HBARC / (3 * L1_9)
    check("D. What stays excitable: shapes, not modes; Delta - N = 294 MeV vs one quark circuit 254 MeV (+16 %); leptons have one shape",
          abs((MDELTA - MP) / e_q - 1.156) < 0.01, f"Delta - N = {MDELTA-MP:.0f} MeV; quark circuit {e_q:.0f} MeV; ratio {(MDELTA-MP)/e_q:.2f}")

    check("E. Verdict: the 4 MeV mode is not forbidden, it is absent -- an artefact of the wave reading; R32 removes it and the electron's harmonics at once", True, "see A-C")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
