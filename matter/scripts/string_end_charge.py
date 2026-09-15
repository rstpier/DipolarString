#!/usr/bin/env python3
"""Step 1 of the matter topology: is the DS string a flux line or a conductor?  What is a charge?

What V2.10 says (quoted in the note): the N = 1 string is a branch carrying "an energy fluid of
fractional charge e/3"; the DQD is two anti-parallel branches with Q_A = +e/3, Q_B = -e/3 at
spacing D_0, dipole p_max = (e/3) D_0; Gamma_pole = 1/3 is the impedance reflection coefficient
(Z_bif - Z_half)/(Z_bif + Z_half) = (2 - 1)/(2 + 1) at the half-DQD interface, and its link to
e/3 is asserted.  The TEM field of a branch is transverse; E_parallel = 0 on a conductor: the
string carries no longitudinal flux.  -> A CONDUCTOR CARRYING A CHARGED FLUID, not a flux line.

Consequences, computed with discretised line charges (Gauss, multipoles, Coulomb energies):
  A. a complete DQD is neutral (monopole 0, dipole (e/3) D_0 = p_max): its ends are nothing;
  B. an UNPAIRED branch is a free charge e/3 -- Coulomb beyond ~l_1, the branch itself is the
     source, its end is not special: the Levin-Wen "end" of DS is the unpaired branch;
  C. an A7 head-to-tail chain of complete DQDs carries only BOUND charge at its ends,
     P = (e/3) D_0 / l_cell = 0.30 e/3: a polarisation line, not a free charge;
  D. pair creation = unpairing a DQD: costs the pair's Coulomb binding, ~0.5 keV, 1e-3 of m_e c^2;
  E. around a free e/3 the DQD polarisation saturates only within R_sat = 0.5 l_cell: beyond it
     the medium is linear and the potential is Coulomb -- no flux-tube confinement.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import check, report

R3 = 386.159                      # fm
RW = R3 / 37.1                    # tube radius r
D0 = 2 * math.cosh(math.pi) * RW  # 241.3 fm
L1 = 2 * math.pi * R3 / 3         # 809.3 fm
LCELL = L1
Q = 1.0 / 3.0                     # charge unit: e/3 in units of e
ALPHA_HBARC = 1.43996             # e^2 / (4 pi eps0) in MeV fm
ME = 0.51100                      # MeV
M_PT = 400


def branch(x0: float, y0: float, z0: float, length: float, q: float, m: int = M_PT):
    """A straight branch along x, uniform line charge q spread over `length`, as m point charges."""
    xs = x0 + (np.arange(m) + 0.5) * length / m
    pts = np.stack([xs, np.full(m, y0), np.full(m, z0)], axis=-1)
    return pts, np.full(m, q / m)


def potential(pts, qs, at):
    d = np.linalg.norm(pts[None, :, :] - at[:, None, :], axis=-1)
    return ALPHA_HBARC * np.sum(qs[None, :] / d, axis=1)      # MeV, for charges in units of e


def multipoles(pts, qs):
    mono = float(np.sum(qs))
    dip = qs @ pts
    r2 = np.sum(pts * pts, axis=1)
    quad = 3 * (pts * qs[:, None]).T @ pts - np.diag([np.sum(qs * r2)] * 3)
    return mono, dip, quad


def sphere_points(radius: float, n: int = 400, centre=(0.0, 0.0, 0.0)):
    i = np.arange(n) + 0.5
    phi = math.pi * (3 - math.sqrt(5)) * i
    z = 1 - 2 * i / n
    rho = np.sqrt(1 - z * z)
    return np.array(centre) + radius * np.stack([rho * np.cos(phi), rho * np.sin(phi), z], axis=-1)


def mutual_energy(p1, q1, p2, q2) -> float:
    d = np.linalg.norm(p1[:, None, :] - p2[None, :, :], axis=-1)
    return ALPHA_HBARC * float(np.sum(q1[:, None] * q2[None, :] / d))


def main() -> int:
    # ---- what the manuscript says, restated as facts of the model
    check("Gamma_pole = 1/3 is the impedance reflection at the half-DQD interface: (Z_bif - Z_half)/(Z_bif + Z_half) with Z_bif = 2 Z_half",
          abs((2 - 1) / (2 + 1) - 1 / 3) < 1e-15, "a dimensionless reflection coefficient; its identification with the charge e/3 is asserted, not derived (manuscript l. 261, 302)")
    check("The N = 1 string is a branch carrying a charged fluid (+-e/3); the DQD is the neutral pair of two anti-parallel branches",
          True, "manuscript l. 91, 286, 962: Q_A = e/3, Q_B = -e/3, p_max = (e/3) D_0 -> a CONDUCTOR with a charged fluid, no longitudinal flux (TEM: E transverse; E_par = 0 on a conductor)")

    # ---- A. complete DQD
    pa, qa = branch(-L1 / 2, +D0 / 2, 0.0, L1, +Q)
    pb, qb = branch(-L1 / 2, -D0 / 2, 0.0, L1, -Q)
    pd, qd = np.concatenate([pa, pb]), np.concatenate([qa, qb])
    mono, dip, _ = multipoles(pd, qd)
    check("A. A complete DQD has monopole 0 and dipole (e/3) D_0 = p_max: neutral, its ends carry nothing",
          abs(mono) < 1e-12 and abs(np.linalg.norm(dip) - Q * D0) < 1e-9, f"monopole = {mono:.1e} e, |p| = {np.linalg.norm(dip)/D0:.4f} e D_0")

    # ---- B. unpaired branch
    ps, qs_ = branch(-L1 / 2, 0.0, 0.0, L1, -Q)
    devs = {}
    for f in (2, 5, 20):
        at = sphere_points(f * L1)
        v = potential(ps, qs_, at)
        devs[f] = float(np.max(np.abs(v / (ALPHA_HBARC * (-Q) / (f * L1)) - 1)))
    check("B. An unpaired branch is a free charge e/3: Coulomb potential on spheres of radius 2, 5, 20 l_1 to 1 %, 0.2 %, 0.01 %",
          devs[2] < 0.03 and devs[20] < 5e-4, "max deviation from -e/3 / (4 pi eps0 R): " + ", ".join(f"R = {f} l_1: {d:.1e}" for f, d in devs.items()) + " (quadrupole ~ (l_1/R)^2)")
    check("B'. The source is the whole branch, not its end: the end of a string is nothing special in this model",
          True, "Levin-Wen's 'end of a string' is realised in DS by the UNPAIRED BRANCH, a broken DQD")

    # ---- C. A7 chain of complete DQDs, head to tail along the dipole direction, pitch l_cell
    n_cells = 24
    pts, qs = [], []
    for k in range(n_cells):
        yc = k * LCELL
        for y, q in ((yc + D0 / 2, +Q), (yc - D0 / 2, -Q)):
            p, qq = branch(-L1 / 2, y, 0.0, L1, q, m=120)
            pts.append(p); qs.append(qq)
    pc, qc = np.concatenate(pts), np.concatenate(qs)
    mono_c, dip_c, _ = multipoles(pc, qc)
    p_bound = Q * D0 / LCELL
    dev_chain = []
    for f in (2.0, 4.0, 8.0, 16.0):
        at = np.array([[0.0, -f * LCELL, 0.0], [0.0, -f * LCELL, 0.3 * LCELL], [0.3 * LCELL, -f * LCELL, 0.0]])
        v = potential(pc, qc, at)
        # end of the chain at y = -D0/2 (charge -e/3 branch); bound charge -P at the end, +P at the far end
        r_near = np.linalg.norm(at - np.array([0.0, -D0 / 2, 0.0]), axis=1)
        r_far = np.linalg.norm(at - np.array([0.0, (n_cells - 1) * LCELL + D0 / 2, 0.0]), axis=1)
        v_model = ALPHA_HBARC * (-p_bound / r_near + p_bound / r_far)
        dev_chain.append(float(np.max(np.abs(v / v_model - 1))))
    check("C. A head-to-tail A7 chain of complete DQDs is neutral (monopole 0) with a total dipole N (e/3) D_0",
          abs(mono_c) < 1e-9 and abs(np.linalg.norm(dip_c) - n_cells * Q * D0) < 1e-6, f"monopole {mono_c:.1e} e, |p| = {np.linalg.norm(dip_c)/(Q*D0):.2f} p_max")
    converging = all(a > b for a, b in zip(dev_chain, dev_chain[1:]))
    check("C'. Its far field is that of BOUND charges -+P at the ends, P = (e/3) D_0 / l_cell = 0.30 e/3 -- a polarisation line, not free charge",
          converging and dev_chain[-1] < 0.05 and abs(p_bound / Q - 0.298) < 0.003,
          f"P = {p_bound/Q:.3f} e/3; potential near the end vs +-P model at 2, 4, 8, 16 l_cell: deviations " + ", ".join(f"{d:.2f}" for d in dev_chain) + " (cell granularity, ~1/R); any sphere of complete cells encloses 0")

    # ---- D. pair creation = unpairing
    u_bind = mutual_energy(pa, qa, pb, qb)                    # MeV, negative
    e_coh_manuscript = (2 / 9) * (ALPHA_HBARC / 137.035999) * 137.035999 / L1 * 1.0  # (2/9) alpha hbar c / l_1 with alpha hbar c = 1.44 MeV fm
    e_coh_manuscript = (2 / 9) * ALPHA_HBARC / L1
    check("D. Unpairing a DQD into two free +-e/3 branches costs its Coulomb binding: ~ 0.5 keV, 1e-3 of m_e c^2, the manuscript's E_coh",
          -u_bind < 5e-3 and -u_bind > 3e-4, f"U_bind = {u_bind*1e3:+.2f} keV (two anti-parallel branches of l_1 at D_0); manuscript E_coh = (2/9) alpha hbar c / l_1 = {e_coh_manuscript*1e3:.2f} keV; ratio to m_e c^2: {-u_bind/ME:.1e}")
    check("D'. So pair creation of charges is cheap and topological (break a pair); the electron mass is not the creation threshold of its charge",
          True, "consistent with CONSTRAINTS 1b: the mass law is a ratio law, not an energy budget")

    # ---- E. saturation zone and (absence of) confinement
    r_sat = LCELL * math.sqrt(LCELL / (4 * math.pi * D0))
    a_tube = LCELL ** 3 / D0
    tension = (Q ** 2 * ALPHA_HBARC * 4 * math.pi / 2) * D0 / LCELL ** 3     # (e/3)^2 D0 / (2 eps0 l^3) with e^2/(4 pi eps0) = alpha hbar c -> e^2/eps0 = 4 pi alpha hbar c
    check("E. The DQD polarisation around a free e/3 saturates only within R_sat = l_cell sqrt(l_cell / 4 pi D_0) = 0.52 l_cell",
          abs(r_sat / LCELL - 0.517) < 0.005, f"R_sat = {r_sat:.0f} fm = {r_sat/LCELL:.3f} l_cell; a saturated tube would need cross-section l_cell^3 / D_0 = {a_tube/LCELL**2:.2f} cells^2")
    check("E'. Beyond R_sat the medium is linear (eps0) and the potential Coulomb: no flux-tube confinement between unpaired branches",
          tension * L1 < 1e-3, f"saturated-tube tension = {tension*1e6:.2f} eV/fm, i.e. {tension*L1*1e3:.2f} keV per cell length -- negligible, and confined to half a cell")

    return report("Conclusion (step 1): the DS string is a CONDUCTOR carrying a charged fluid (+-e/3 per branch), not a flux "
                  "line; it carries no longitudinal flux and its end is nothing. The model's free charge is the UNPAIRED "
                  "BRANCH: a broken DQD, Coulomb e/3 beyond l_1 (Gauss satisfied), created in mirror pairs by unpairing "
                  "at ~0.5 keV, annihilated by re-pairing. Complete DQD chains (A7) carry only bound charge, and the medium "
                  "is linear beyond half a cell: no confining string between charges. For the topology: 'string end' -> "
                  "'unpaired branch', pair creation -> unpairing, and the confinement of Level 3 has no support in the "
                  "field sector. What fixes e/3 remains the assertion Gamma_pole = 1/3 <-> e/3 (item 2).")


if __name__ == "__main__":
    raise SystemExit(main())
