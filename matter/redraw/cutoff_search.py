#!/usr/bin/env python3
"""Redraw: what fixes the cut-off b/a = 14 of R35 (the outer radius b, where the wound strand's
field is neutralised, over its core radius a, entering the Z_0-matched pitch n a = sqrt(2/ln(b/a)))?

With lambda = delta/l_1(9) (R35), sigma = 4 pi K lambda^2 x pi (n a)^2 = (4 hbar c/(pi l_1^2)) x 2 pi/ln(b/a).
Checks:
  A. the lattice band sqrt(sigma) = 420-440 MeV (894-981 MeV/fm) maps to b/a = 11.4-14.5.
  B. candidates from the base: the partner branch at D, D/r = 2 cosh pi = 23.2, ln = pi exactly ->
     gain 2, sigma = 8 hbar c/(pi l_1^2) = 761 MeV/fm (-16 %); the half-spacing D/2r = cosh pi = 11.6
     -> gain 2.56, sigma = 976 MeV/fm (+8 %, inside the band); R_3/r = 37.1 -> 661; r_p/lambda-bar_p
     = 4 -> 1724; l_1/lambda-bar_p = 3.9 -> 1770; R-/lambda-bar_p = 5.2 -> 1453; (r_e/2)/lambda-bar_p
     = 6.7 -> 1257.  Only the two Z_0-matching ratios come close.
  C. the two matching readings bracket the lattice: 761 (b = D) and 976 (b = D/2); the band lies
     between.  The cut-off is fixed by the same matching that fixes D/r, to within the D-vs-D/2
     identification, i.e. +-15 % on sigma.
  D. closed form with ln(b/a) = pi: sigma = 8 hbar c/(pi l_1(9)^2) = (18 x 3^(4 pi)/pi^3) m_e^2 c^4/(hbar c),
     sqrt(sigma) = 0.76 x 3^(2 pi) m_e c^2 = 388 MeV (lattice 420-440, -8 to -12 %); with b = D/2,
     sqrt(sigma) = 439 MeV.  The strong tension is written in m_e and the ladder exponent 2 pi alone.
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
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
LAM = DELTA / L1_9
BAND = (420.0, 440.0)                      # sqrt(sigma), MeV
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def sigma(b_over_a):
    return 4 * math.pi * K * LAM**2 * 2 * math.pi / math.log(b_over_a)


def b_over_a_for(sig):
    return math.exp(4 * math.pi * K * LAM**2 * 2 * math.pi / sig)


def main() -> int:
    sig_lo, sig_hi = (s**2 / HBARC for s in BAND)
    ba_hi, ba_lo = b_over_a_for(sig_lo), b_over_a_for(sig_hi)
    check("A. Lattice band sqrt(sigma) = 420-440 MeV (894-981 MeV/fm) maps to b/a = 11.4-14.5",
          abs(sig_lo - 894) < 1 and abs(sig_hi - 981) < 1 and abs(ba_lo - 11.4) < 0.1 and abs(ba_hi - 14.5) < 0.1,
          f"sigma = {sig_lo:.0f}-{sig_hi:.0f} MeV/fm -> b/a = {ba_lo:.1f}-{ba_hi:.1f}")

    cands = {"D/r = 2 cosh pi": 2 * math.cosh(math.pi), "D/2r = cosh pi": math.cosh(math.pi), "R_3/r": 37.1,
             "r_p/lambda_p": RP / LAMBDA_P, "l_1/lambda_p": L1_9 / LAMBDA_P, "R-/lambda_p": 5.184,
             "(r_e/2)/lambda_p": (ALPHA * LAMBDA_E / 2) / LAMBDA_P}
    sig = {k: sigma(v) for k, v in cands.items()}
    inside = [k for k, s in sig.items() if sig_lo <= s <= sig_hi]
    check("B. Candidates: partner branch D/r = 2 cosh pi (ln = pi, gain 2) -> 761 MeV/fm (-16 %); half-spacing D/2r = cosh pi -> 976 (+8 %, inside the band); R_3/r -> 661; r_p/lambda_p -> 1724; l_1/lambda_p -> 1770; R-/lambda_p -> 1453; (r_e/2)/lambda_p -> 1257. Only the two Z_0-matching ratios come close",
          abs(sig["D/r = 2 cosh pi"] - 761) < 2 and abs(sig["D/2r = cosh pi"] - 976) < 3 and inside == ["D/2r = cosh pi"],
          "; ".join(f"{k} = {v:.1f} -> {sig[k]:.0f}" for k, v in cands.items()))

    s_D, s_D2 = sig["D/r = 2 cosh pi"], sig["D/2r = cosh pi"]
    check("C. The two matching readings bracket the lattice: 761 (b = D) and 976 (b = D/2), the band 894-981 lies between: the cut-off is fixed by the matching that fixes D/r, to within the D-vs-D/2 identification, +-15 % on sigma",
          s_D < sig_lo < sig_hi < s_D2 * 1.01 and abs(s_D2 / s_D - 1.28) < 0.01,
          f"761 < 894-981 <= 976; ratio D/2 vs D = {s_D2/s_D:.2f}")

    sig_closed = 8 * HBARC / (math.pi * L1_9**2)
    sig_me = 18 * 3 ** (4 * math.pi) / math.pi**3 * ME**2 / HBARC
    sqrt_D, sqrt_D2 = math.sqrt(s_D * HBARC), math.sqrt(s_D2 * HBARC)
    check("D. Closed form with ln(b/a) = pi: sigma = 8 hbar c/(pi l_1(9)^2) = (18 x 3^4pi/pi^3) m_e^2 c^4/hbar c, sqrt(sigma) = 388 MeV (lattice 420-440, -8 to -12 %); with b = D/2, 439 MeV: the strong tension written in m_e and the exponent 2 pi alone",
          abs(sig_closed - s_D) < 1 and abs(sig_me - sig_closed) < 1e-6 and abs(sqrt_D - 387.5) < 1 and abs(sqrt_D2 - 439) < 1,
          f"sigma = {sig_closed:.1f} MeV/fm (ln(2 cosh pi) = {math.log(2*math.cosh(math.pi)):.4f} vs pi, {sig_closed/s_D-1:+.2%}); sqrt(sigma) = {sqrt_D:.0f} MeV (b = D), {sqrt_D2:.0f} MeV (b = D/2); 0.76 x 3^2pi x m_e = {math.sqrt(18/math.pi**3)*3**(2*math.pi)*ME:.0f} MeV")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
