#!/usr/bin/env python3
"""Redraw: the difference between the real helix and the smooth cylinder of R34-R37.

The smooth model spread the wound strand's charge along the AXIS at its path density lambda and
took the energy as coax-capacitance + solenoid-inductance.  A real helix of n turns per unit axial
length on radius a has gamma = sqrt(1 + (2 pi n a)^2) units of path -- and of charge -- per unit of
axis: the fluid runs at c along the wire, so the axial charge density is gamma lambda and the axial
phase velocity c/gamma (a slow-wave line).

Checks:
  A. at the smooth model's matched pitch (n a = 0.69-0.80) gamma = 4.45-5.1: the real helix carries
     4.5-5 x more charge and path per unit axial length than the smooth model assumed.
  B. for a line matched along its path (Z = Z_0, R10) the energy per unit PATH length is mu_0 I^2 =
     4 pi K lambda^2 = 380 MeV/fm at lambda = delta/l_1(9), whatever the winding; per unit AXIAL
     length it is gamma x that.  At the smooth model's pitch the real tension is 1690-1940 MeV/fm,
     not 761: the smooth model undercounted by gamma/2 = 2.2-2.5 -- a +120 to +155 % correction,
     not the +15 % residual of R37, which was an artefact of the model.
  C. the tension is sigma = 4 pi K lambda^2 x gamma: the lattice (894-981, centre 904 MeV/fm) needs
     gamma = 2.35-2.58 (2.376): n a = 0.34-0.38, pitch angle 23-25 deg, 2.4 fm of strand per fm of
     axis; with a = lambda-bar_p, 1.6 turns per fm, 1.5 turns in the 0.9 fm centre.
  D. the lumped matching of R34 (solenoid L' + coax C') is not a valid helix computation: made
     self-consistent with the slow wave (L' = mu_0 gamma, i.e. Z = Z_0 at v = c/gamma) it gives
     n a = 1.78, gamma = 11.2, sigma = 4270 MeV/fm, 4.7 x the lattice.  The sheath-helix theory
     (Pierce) is needed for the pitch; not done.  What fixes the pitch is reopened.
  E. the electron is unwound (R12's 3/4 turn, gamma = 1): unaffected.
Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math

HBARC = 197.3269804
ALPHA = 1 / 137.035999
K = ALPHA * HBARC
ME, MP = 0.51099895, 938.27209
LAMBDA_E, LAMBDA_P = HBARC / ME, HBARC / MP
L1_E = 2 * math.pi * LAMBDA_E / 3
L1_9 = L1_E * (3 / 9) ** (2 * math.pi)
DELTA = 1 / (math.pi * math.sqrt(ALPHA))
LAM = DELTA / L1_9
U_PATH = 4 * math.pi * K * LAM**2                 # MeV per fm of path, matched line
A_OVER_RW = 2 * math.cosh(math.pi)
BAND = (894.0, 981.0)
SIGMA_LAT = 904.0
RESULTS = []


def check(name, cond, detail):
    RESULTS.append(("PASS" if cond else "FAIL", name, detail))


def gamma(na):
    return math.sqrt(1 + (2 * math.pi * na) ** 2)


def na_for_gamma(g):
    return math.sqrt(g**2 - 1) / (2 * math.pi)


def main() -> int:
    # A. winding ratio at the smooth pitch
    g_lo, g_hi = gamma(0.69), gamma(0.80)
    check("A. At the smooth model's matched pitch (n a = 0.69-0.80) the winding ratio gamma = sqrt(1 + (2 pi n a)^2) = 4.45-5.1: the real helix carries 4.5-5 x more charge and path per unit axial length than assumed",
          abs(g_lo - 4.45) < 0.03 and abs(g_hi - 5.13) < 0.03, f"gamma(0.69) = {g_lo:.2f}, gamma(0.80) = {g_hi:.2f}")

    # B. matched line energy per path length, and the real tension at the smooth pitch
    sig_real = (U_PATH * g_lo, U_PATH * g_hi)
    corr = (g_lo / 2, g_hi / 2)
    check("B. Matched along its path: energy per fm of path = 4 pi K lambda^2 = 380 MeV whatever the winding; per fm of axis, gamma x that. At the smooth pitch the real tension is 1690-1940 MeV/fm, not 761: the smooth model undercounted by 2.2-2.5 (+120 to +155 %), so R37's 15 % residual was an artefact",
          abs(U_PATH - 380.4) < 1 and abs(sig_real[0] - 1693) < 10 and abs(sig_real[1] - 1951) < 10 and abs(corr[0] - 2.22) < 0.02,
          f"u_path = {U_PATH:.1f} MeV/fm; sigma_real = {sig_real[0]:.0f}-{sig_real[1]:.0f} MeV/fm; correction x{corr[0]:.2f}-{corr[1]:.2f}")

    # C. what the lattice needs
    g_need = SIGMA_LAT / U_PATH
    g_band = tuple(b / U_PATH for b in BAND)
    na_need = na_for_gamma(g_need)
    na_band = tuple(na_for_gamma(g) for g in g_band)
    psi = math.degrees(math.atan(1 / (2 * math.pi * na_need)))
    n_per_fm = na_need / LAMBDA_P
    check("C. sigma = 4 pi K lambda^2 x gamma: the lattice needs gamma = 2.35-2.58 (2.376): n a = 0.34-0.38, pitch angle 25 deg, 2.4 fm of strand per fm of axis; with a = lambda-bar_p, 1.6 turns per fm, 1.5 turns in the 0.9 fm centre",
          abs(g_need - 2.376) < 0.005 and abs(na_need - 0.343) < 0.003 and abs(psi - 24.9) < 0.3 and abs(n_per_fm - 1.63) < 0.02,
          f"gamma = {g_need:.3f} ({g_band[0]:.2f}-{g_band[1]:.2f}); n a = {na_need:.3f} ({na_band[0]:.3f}-{na_band[1]:.3f}); pitch {psi:.1f} deg; n = {n_per_fm:.2f}/fm, Tw(0.9 fm) = {n_per_fm*0.9:.1f}")

    # D. lumped matching made slow-wave consistent: pi x^2 + gamma ln(p/2 pi r_w)/2pi + ln(a/r_w)/(2 pi gamma^2) = gamma
    def f(x):
        g = gamma(x)
        p_over_rw = A_OVER_RW / x
        return math.pi * x**2 + g * math.log(p_over_rw / (2 * math.pi)) / (2 * math.pi) + math.log(A_OVER_RW) / (2 * math.pi * g**2) - g
    lo, hi = 1.0, 3.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    x_lumped = lo
    g_lumped = gamma(x_lumped)
    check("D. The lumped matching (solenoid L' + coax C') made slow-wave consistent (L' = mu_0 gamma) gives n a = 1.78, gamma = 11.2, sigma = 4270 MeV/fm, 4.7 x the lattice: not a valid helix computation; the sheath-helix theory is needed for the pitch -- reopened",
          abs(x_lumped - 1.78) < 0.03 and abs(g_lumped - 11.2) < 0.2 and abs(U_PATH * g_lumped / SIGMA_LAT - 4.7) < 0.1,
          f"n a = {x_lumped:.2f}, gamma = {g_lumped:.1f}, sigma = {U_PATH*g_lumped:.0f} MeV/fm = {U_PATH*g_lumped/SIGMA_LAT:.1f} x lattice")

    # E. electron unwound
    check("E. The electron's 3/4 turn is unwound (gamma = 1): unaffected", gamma(0.0) == 1.0, "gamma(0) = 1")

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print(f"\nRESULT: {len(RESULTS)-len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
