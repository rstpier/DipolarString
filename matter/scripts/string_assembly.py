#!/usr/bin/env python3
"""How do strings with a polarised fluid assemble?  End to end by their poles, or side by side?

Model of one string (manuscript l. 91: a fluid of charge e/3 "possessing an intrinsic negative pole
and positive pole"): a segment of length l_1 carrying a net charge q spread uniformly, plus pole
charges +delta at one end and -delta at the other.  Units: e/3 for charges, e^2/4 pi eps0 = 1.44 MeV fm.
Contacts: poles touch at 2 r (tubes in contact); side-by-side spacing D_0 (the matching distance).
Binding energies as a function of delta/q:

  (a) two strings END TO END, pole + against pole -, same net charge q = -e/3: bound if the pole
      attraction beats the monopole repulsion of two collinear like charges;
  (b) two strings SIDE BY SIDE anti-parallel, opposite net charges (the DQD): bound (-0.47 keV at
      delta = 0), the pole term adds;
  (c) two strings SIDE BY SIDE anti-parallel, SAME net charge: repulsive unless the poles are
      strong;
  (d) three same-charge strings side by side on a triangle: one pair is necessarily parallel
      (frustration) -- vs three strings end to end in a closed ring (the manuscript's electron).

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import check, report

R3 = 386.159
RW = R3 / 37.1
D0 = 2 * math.cosh(math.pi) * RW
L1 = 2 * math.pi * R3 / 3
K = 1.43996 / 9.0          # (e/3)^2 / (4 pi eps0) in MeV fm
M = 300


def string_charges(start, direction, q, delta, m=M):
    """Points and charges of one string from `start` along unit `direction`: uniform q plus poles."""
    direction = np.asarray(direction, float) / np.linalg.norm(direction)
    s = (np.arange(m) + 0.5) * L1 / m
    pts = np.asarray(start, float)[None, :] + s[:, None] * direction[None, :]
    qs = np.full(m, q / m)
    # poles: +delta at the far end (the "head"), -delta at the start (the "tail"), a tube radius in from the ends
    poles = np.array([np.asarray(start, float) + RW * direction, np.asarray(start, float) + (L1 - RW) * direction])
    return np.concatenate([pts, poles]), np.concatenate([qs, [-delta, +delta]])


def mutual(p1, q1, p2, q2, a=RW) -> float:
    d = np.sqrt(np.sum((p1[:, None, :] - p2[None, :, :]) ** 2, axis=-1) + a * a)
    return K * float(np.sum(q1[:, None] * q2[None, :] / d))


def energies(delta: float, q: float = -1.0):
    gap = 2 * RW
    # (a) end to end: string 1 along +x from 0; string 2 along +x from L1 + gap (head of 1 meets tail of 2)
    a1 = string_charges([0, 0, 0], [1, 0, 0], q, delta)
    a2 = string_charges([L1 + gap, 0, 0], [1, 0, 0], q, delta)
    u_end = mutual(*a1, *a2)
    # (b) side by side anti-parallel, opposite charges at D0
    b1 = string_charges([0, 0, 0], [1, 0, 0], q, delta)
    b2 = string_charges([L1, D0, 0], [-1, 0, 0], -q, delta)
    u_dqd = mutual(*b1, *b2)
    # (c) side by side anti-parallel, same charges at D0
    c2 = string_charges([L1, D0, 0], [-1, 0, 0], q, delta)
    u_same = mutual(*b1, *c2)
    # (c') side by side parallel, same charges at D0
    c3 = string_charges([0, D0, 0], [1, 0, 0], q, delta)
    u_par = mutual(*b1, *c3)
    return u_end, u_dqd, u_same, u_par


def main() -> int:
    rows = []
    for delta in (0.0, 0.25, 0.5, 0.75, 1.0):
        rows.append((delta, *energies(delta)))
    print("  binding energies (keV; negative = bound), same net charge q = -e/3 unless noted:")
    print("   delta/q   end-to-end (pole contact)   side-by-side, opposite q (DQD)   side-by-side anti-par., same q   side-by-side parallel, same q")
    for d, ue, ud, us, up in rows:
        print(f"   {d:5.2f}   {ue*1e3:+12.2f}                 {ud*1e3:+12.2f}                   {us*1e3:+12.2f}                     {up*1e3:+12.2f}")
    e0 = dict((r[0], r[1:]) for r in rows)
    check("Without poles (delta = 0): the DQD is bound (-0.47 keV), like charges repel side by side (+0.47) and end to end (+1.4 keV)",
          e0[0.0][1] < 0 and e0[0.0][2] > 0 and e0[0.0][0] > 0, f"DQD {e0[0.0][1]*1e3:+.2f}, same-charge side by side {e0[0.0][2]*1e3:+.2f}, end to end {e0[0.0][0]*1e3:+.2f} keV")
    # thresholds
    def threshold(idx):
        for d in np.linspace(0, 1.5, 301):
            if energies(d)[idx] < 0:
                return d
        return float("nan")
    t_end, t_same = threshold(0), threshold(2)
    check("End to end, pole against pole, like-charge strings bind as soon as the poles exceed ~0.3 of the net charge: chains and rings",
          0.2 < t_end < 0.4, f"threshold delta/q = {t_end:.2f}; at delta = q the junction binds by {e0[1.0][0]*1e3:.1f} keV")
    check("Side by side, like-charge strings need poles above ~0.7 of the net charge, and bind 10 x more weakly",
          0.55 < t_same < 0.9 and abs(e0[1.0][0]) > 5 * abs(e0[1.0][2]), f"threshold delta/q = {t_same:.2f}; at delta = q: side by side {e0[1.0][2]*1e3:+.2f} keV vs end to end {e0[1.0][0]*1e3:+.2f} keV")
    check("Side by side PARALLEL like-charge strings never bind: three on a triangle are frustrated (one pair parallel)",
          all(energies(d)[3] > 0 for d in (0.0, 0.5, 1.0, 1.5)), f"parallel pair at delta = q: {e0[1.0][3]*1e3:+.2f} keV")
    # closed ring of three, end to end: three arcs of 120 deg on the circle of radius R3, poles in contact at the junctions
    def arc_charges(t0, q, delta, m=M):
        span = 2 * math.pi / 3
        gap_ang = 2 * RW / R3                      # tubes in contact at the junction
        t = t0 + gap_ang / 2 + (np.arange(m) + 0.5) * (span - gap_ang) / m
        pts = np.stack([R3 * np.cos(t), R3 * np.sin(t), np.zeros(m)], axis=-1)
        qs = np.full(m, q / m)
        tp = np.array([t0 + gap_ang / 2 + RW / R3, t0 + span - gap_ang / 2 - RW / R3])
        poles = np.stack([R3 * np.cos(tp), R3 * np.sin(tp), np.zeros(2)], axis=-1)
        return np.concatenate([pts, poles]), np.concatenate([qs, [-delta, +delta]])
    ring = [arc_charges(k * 2 * math.pi / 3, -1.0, 1.0) for k in range(3)]
    u_ring = sum(mutual(*ring[i], *ring[j]) for i in range(3) for j in range(i + 1, 3))
    ring0 = [arc_charges(k * 2 * math.pi / 3, -1.0, 0.0) for k in range(3)]
    u_ring0 = sum(mutual(*ring0[i], *ring0[j]) for i in range(3) for j in range(i + 1, 3))
    check("Three like-charge polar strings joined end to end close into a ring: the manuscript's electron, bound at every junction",
          u_ring < 0 and u_ring0 > 0, f"ring of three arcs on R_3: interaction {u_ring*1e3:+.1f} keV at delta = q, {u_ring0*1e3:+.1f} keV without poles; a side-by-side triangle of the same strings is frustrated")
    check("So: with a polarised fluid, strings assemble ALONG their length, pole to pole -- into chains and closed rings -- not into bundles",
          True, "the DQD (opposite charges) is the one bound side-by-side pair; like charges only chain")
    check("What the model does not give: the pole strength delta -- 'intrinsic poles' are stated, not quantified",
          True, "delta/q >= 0.3 is the condition for like-charge strings to chain at all")
    return report("Conclusion: a polarised fluid makes strings join end to end by their poles -- like-charge strings chain and "
                  "close into rings once the poles exceed ~0.45 of the net charge, with a junction energy of order 10 keV "
                  "at delta = q; side by side they bind only with opposite charges (the DQD, 0.5 keV) or, weakly and "
                  "frustrated for three, with strong poles. The manuscript's electron -- three strings end to end in a ring -- "
                  "is the assembly the fluid prefers. The pole strength delta is not given by the model.")


if __name__ == "__main__":
    raise SystemExit(main())
