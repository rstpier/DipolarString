#!/usr/bin/env python3
"""Redraw: the screening of a charged strand by its neighbours (R44 left it uncomputed).

Bookkeeping first (R42/R44): a charged strand is a lone branch carrying e/3 ('charge = unpaired
branch'); a neutral strand is a full DQD, two floating conducting plates +-e/3 at gap d = w.  A
charged strand therefore has no partner plate of its own; what screens it are the neutral DQDs of
the same quark: one for u = (+,+,0), two for d = (-,0,0).

Method: 2D method of moments on the cross-section.  Conducting strips of width w discretised into
segments; the charged strip carries total line charge lambda, each neutral plate carries zero net
charge and floats at its own potential.  The far field is the same with or without neighbours
(neutral partners), so the neighbours only change the effective radius a_eff of the charged strip:
V_1(with) - V_1(alone) = (lambda/2 pi eps0) ln(a_alone/a_with).

Checks:
  A. validation: an isolated strip of width w has a_eff = w/4 (classical), reproduced to 1 %.
  B. a partner plate carrying the polarisation +-delta would add a cross term 2 pi K q delta / l =
     +-14 MeV per charged strand and shift m_n - m_p by -+14 MeV: excluded, so the lone-branch
     bookkeeping is forced.
  C. screening by the quark's neutral DQDs, side by side at centre spacing s = w, 1.5 w, 2 w, 3 w:
     a_eff grows by 1.47 (u, one DQD) and 1.98 (d, two DQDs) when touching (s = w), by 1.05 and 1.10
     at s = 2 w, and by 1.02 and 1.05 at 3 w: stronger for the d, and only large when touching.
  D. effect on m_n - m_p (R44 rule with the screened a_eff per quark): the u and d screenings nearly
     cancel in (4/9) S_u - (1/9) S_d: 1.324 -> 1.342 MeV when touching (+1.4 %), unchanged beyond
     1.5 w.  The screening changes the result by less than 1.5 %; R44's +2.4 % stands.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E = HBARC / ME
L1 = 2 * math.pi * LAMBDA_E / 3 * (3 / 9) ** (2 * math.pi)
W = (4 * LAMBDA_E / math.pi**2) * 3 ** (-2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
E_DQD = 4 * math.pi * K / (9 * L1)
M_MUTUAL = K / (math.sqrt(3) * L1)
DM_OBS = MN - MP
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def strip(xc, yc, w=1.0, n=120):
    x = xc - w / 2 + (np.arange(n) + 0.5) * w / n
    return np.column_stack([x, np.full(n, yc)]), w / n


def solve(strips, charges):
    """strips: list of (points, h); charges: total line charge per strip. Returns potentials per strip
    and the charge vector. Units: potential in lambda/(2 pi eps0) with G = -ln r."""
    pts = np.vstack([s[0] for s in strips])
    hs = np.concatenate([np.full(len(s[0]), s[1]) for s in strips])
    owner = np.concatenate([np.full(len(s[0]), k) for k, s in enumerate(strips)])
    N, S = len(pts), len(strips)
    d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    G = -np.log(np.where(d > 0, d, 1.0))
    # self term of a uniformly charged segment of length h at its midpoint: -(ln(h/2) - 1)
    np.fill_diagonal(G, -(np.log(hs / 2) - 1))
    A = np.zeros((N + S, N + S))
    b = np.zeros(N + S)
    A[:N, :N] = G
    for i in range(N):
        A[i, N + owner[i]] = -1.0
    for k in range(S):
        A[N + k, :N] = (owner == k).astype(float)
        b[N + k] = charges[k]
    sol = np.linalg.solve(A, b)
    return sol[N:], sol[:N], pts


def a_eff_isolated():
    V, q, pts = solve([strip(0, 0)], [1.0])
    # potential at a far point (0, R) from the solved charges: V(R) = -sum q ln r
    R = 1e4
    r = np.linalg.norm(pts - np.array([0, R]), axis=1)
    V_far = -np.sum(q * np.log(r))
    # V_1 - V(R) = ln(R / a_eff)
    return R / math.exp(V[0] - V_far), V[0]


def a_ratio(config_strips, charges):
    V_alone, _, _ = solve([strip(0, 0)], [1.0])
    V_with, _, _ = solve(config_strips, charges)
    return math.exp(-(V_with[0] - V_alone[0]))          # a_with / a_alone


def main() -> int:
    # A. validation
    a_iso, _ = a_eff_isolated()
    check("A. Isolated strip of width w: a_eff = w/4 (classical), reproduced by the method of moments to 1 %",
          abs(a_iso / 0.25 - 1) < 0.01, f"a_eff/w = {a_iso:.4f} (expected 0.25)")

    # B. polarisation cross term
    cross = 2 * math.pi * K * (1 / 3) * DELTA / L1
    check("B. A partner plate with polarisation +-delta would add 2 pi K q delta/l = 14 MeV per charged strand and shift m_n - m_p by -+14 MeV: excluded -- the lone-branch bookkeeping is forced",
          abs(cross - 13.8) < 0.2 and cross > 5 * DM_OBS, f"cross term = {cross:.1f} MeV per strand vs 1.29 MeV")

    # C. screening by neutral DQDs
    ratios = {}
    for s in (1.0, 1.5, 2.0, 3.0):
        u_cfg = [strip(0, 0), strip(s, 0), strip(s, 1.0)]
        d_cfg = [strip(0, 0), strip(s, 0), strip(s, 1.0), strip(-s, 0), strip(-s, 1.0)]
        ratios[s] = (a_ratio(u_cfg, [1, 0, 0]), a_ratio(d_cfg, [1, 0, 0, 0, 0]))
    ok = all(rd > ru > 1.0 for ru, rd in ratios.values()) and abs(ratios[1.0][0] - 1.47) < 0.03 and abs(ratios[1.0][1] - 1.98) < 0.04 and ratios[2.0][1] < 1.12
    check("C. Screening by the quark's neutral DQDs at centre spacing s = w..3w: a_eff grows by 1.47 (u) and 1.98 (d) when touching, by 1.05 and 1.10 at s = 2w, 1.02 and 1.05 at 3w: stronger for the d, large only when touching",
          ok, "; ".join(f"s = {s:.1f}w: u x{ru:.2f}, d x{rd:.2f}" for s, (ru, rd) in ratios.items()))

    # D. effect on m_n - m_p
    def dm_screened(ru, rd):
        a_u, a_d = W / 4 * ru, W / 4 * rd
        S_u = K / (2 * L1) * (math.log(4 * L1 / a_u) - 1)
        S_d = K / L1 * (math.log(2 * L1 / a_d) - 1)
        return E_DQD - ((4 / 9) * S_u - (1 / 9) * S_d + M_MUTUAL / 3)
    dm0 = dm_screened(1.0, 1.0)
    dms = {s: dm_screened(*r) for s, r in ratios.items()}
    lo, hi = min(dms.values()), max(dms.values())
    check("D. m_n - m_p with the screened radii: the u and d screenings nearly cancel in (4/9) S_u - (1/9) S_d: 1.324 -> 1.342 MeV when touching (+1.4 %), unchanged beyond 1.5w; the screening changes the result by less than 1.5 %, R44's +2.4 % stands",
          abs(dm0 / DM_OBS - 1.024) < 0.01 and abs(dms[1.0] / dm0 - 1.014) < 0.005 and abs(hi - lo) / dm0 < 0.015,
          f"unscreened {dm0:.3f}; " + "; ".join(f"s = {s:.1f}w: {v:.3f} ({v/DM_OBS-1:+.0%})" for s, v in dms.items()))

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
