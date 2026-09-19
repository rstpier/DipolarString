#!/usr/bin/env python3
"""Redraw: what fixes the fluid's charge per unit path length lambda (R34 left it free: the lattice
tension with the Z_0-matched pitch needs 3.3-5.4 e/fm at the nucleon)?  Every density the base can
produce is tried.

Candidates:
  A. ring rule + ladder scale: lambda = (e/3)/l_1(n) = 0.41 e/fm at n = 9: 11 x short (needs n = 13).
  B. the loop's own density e/(2 pi R+) = 0.20 e/fm: 20 x short.
  C. the matched pole charge as the fluid's charge per strand (R17: delta = e/(pi sqrt alpha) = 3.73 e):
     lambda = delta / l_1(9) = 4.6 e/fm -- INSIDE the window.  The wound-strand tension is then
     4 pi K lambda^2 x gain, 660-1720 MeV/fm for the matched pitch with cut-off b/a = 37-4, i.e.
     0.7-1.9 x the lattice, and exactly the lattice at b/a = 14.  No number put in by hand: delta
     from Z_0, l_1 from the ladder, the pitch from Z_0.
  D. the same rule on the electron: three strings of 809 fm with delta per branch store a static line
     energy 3 x 2 pi K delta^2 / l_1 = (9/pi^2) m_e c^2 = 0.91 m_e c^2 (alpha cancels).  But the base
     already pays m_e with mode + junctions (R13): the rule overcounts by 466 keV unless the static
     half (255 keV) IS this line energy, which it misses by 1.8.  CONDITIONAL.
  E. the same line energy per charged strand at the nucleon, 2 hbar c/(pi l_1(9)) = 155 MeV: five
     strands 773 MeV (p), four 618 MeV (n); with R29's mode parts (153, 191) p = 926 MeV (-1.3 %),
     n = 809 MeV (-14 %), wrong order.  Recorded, not claimed.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP, MN = 0.51099895, 938.27209, 939.56542
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
RP = 0.8409
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
R_PLUS = 3.749 * LAMBDA_P
SIGMA_LAT = 904.0
WINDOW = (3.3, 5.4)
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def gain(b_over_a):
    na = math.sqrt(2 / math.log(b_over_a))
    return math.pi * na**2


def main() -> int:
    # A. ring rule + ladder
    lam_a = (1 / 3) / L1_9
    check("A. Ring rule + ladder scale: lambda = (e/3)/l_1(9) = 0.41 e/fm, 11 x short of the window (3.3-5.4 e/fm)",
          abs(lam_a - 0.410) < 0.003 and WINDOW[0] / lam_a > 8, f"lambda = {lam_a:.3f} e/fm; window/lambda = {WINDOW[0]/lam_a:.0f}-{WINDOW[1]/lam_a:.0f}")

    # B. loop density
    lam_b = 1 / (2 * math.pi * R_PLUS)
    check("B. The loop's own density e/(2 pi R+) = 0.20 e/fm: 20 x short",
          abs(lam_b - 0.202) < 0.003, f"lambda = {lam_b:.3f} e/fm")

    # C. matched pole charge per strand
    lam_c = DELTA / L1_9
    sig = {b: 4 * math.pi * K * lam_c**2 * gain(b) for b in (37.1, 14.0, 4.0)}
    b_exact = None
    lo, hi = 2.0, 60.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if 4 * math.pi * K * lam_c**2 * gain(mid) > SIGMA_LAT:
            lo = mid
        else:
            hi = mid
    b_exact = lo
    check("C. Matched pole charge per strand: lambda = delta/l_1(9) = e/(pi sqrt alpha l_1) = 4.6 e/fm, inside the window; wound-strand tension 660-1720 MeV/fm for b/a = 37-4 (0.7-1.9 x lattice), exactly the lattice at b/a = 14 -- delta from Z_0, l_1 from the ladder, pitch from Z_0",
          WINDOW[0] < lam_c < WINDOW[1] and abs(sig[37.1] - 660) < 15 and abs(sig[4.0] - 1720) < 40 and abs(b_exact - 14) < 1,
          f"lambda = {lam_c:.2f} e/fm; sigma(b/a = 37.1) = {sig[37.1]:.0f}, (14) = {sig[14.0]:.0f}, (4) = {sig[4.0]:.0f} MeV/fm; lattice at b/a = {b_exact:.1f}")

    # D. electron consistency
    e_line = 3 * 2 * math.pi * K * DELTA**2 / L1_E
    identity = 9 / math.pi**2
    check("D. Same rule on the electron: three strings with delta per branch store 3 x 2 pi K delta^2/l_1 = (9/pi^2) m_e c^2 = 0.91 m_e c^2 (alpha cancels); the base already pays m_e with mode + junctions, so the rule overcounts by 466 keV unless the static half (255 keV) is this energy, missed by 1.8. CONDITIONAL",
          abs(e_line / ME - identity) < 1e-9 and abs(e_line * 1e3 - 466) < 1 and abs(e_line / (ME / 2) - 1.82) < 0.01,
          f"E_line = {e_line*1e3:.0f} keV = {e_line/ME:.3f} m_e; / static half = {e_line/(ME/2):.2f}")

    # E. nucleons
    e_strand = 2 * HBARC / (math.pi * L1_9)
    mode_p, mode_n = 152.5, 190.7            # R29 D
    p_tot, n_tot = 5 * e_strand + mode_p, 4 * e_strand + mode_n
    check("E. Same line energy per charged strand at the nucleon, 2 hbar c/(pi l_1(9)) = 155 MeV: p = 5 x 155 + 153 (mode) = 926 MeV (-1.3 %), n = 4 x 155 + 191 = 809 MeV (-14 %), wrong order. Recorded, not claimed",
          abs(e_strand - 154.6) < 0.5 and abs(p_tot / MP - 0.987) < 0.003 and abs(n_tot / MN - 0.861) < 0.003 and p_tot > n_tot,
          f"per strand {e_strand:.1f} MeV; p = {p_tot:.0f} ({p_tot/MP-1:+.1%}), n = {n_tot:.0f} ({n_tot/MN-1:+.1%})")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
