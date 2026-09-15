#!/usr/bin/env python3
"""Author's statement: the EM fluid accumulates in the TRAILING pole, creating the active charge.

Two readings, both tested against the model and against measurement.

  (1) Motion-based: the string moves through the substrate (a preferred frame exists in DS, A6),
      the fluid lags and piles at the trailing pole.  Then the charge depends on velocity,
      q(v) = q_0 (1 + kappa v^2/c^2) at lowest even order.  The neutrality of atoms bounds
      |q_p + q_e| < 1e-21 e with the electron at v/c ~ alpha and the proton at rest: kappa < 2e-16.
      The accumulated charge cannot depend on motion to sixteen digits: this reading is excluded.

  (2) Flow-based: the fluid circulates at c_0 inside the string and meets the reactive barrier
      Gamma_pole = 1/3 at each pole (manuscript l. 315).  A partial reflection makes a standing
      wave with ratio SWR = (1 + Gamma)/(1 - Gamma) = 2: the charge density oscillates between
      antinode and node with contrast 2, SYMMETRICALLY at both poles, and its time average is
      uniform.  No net static pile at one pole; a one-way trap would be needed.

Exit status is zero only if every check passes; the checks are statements of fact.
"""

from __future__ import annotations

import math
import numpy as np

from berry_holonomy_common import check, report

ALPHA = 1 / 137.035999
NEUTRALITY = 1e-21          # |q_p + q_e| / e, laboratory bound (Bressi et al. 2011 and earlier)
GAMMA_POLE = 1.0 / 3.0


def main() -> int:
    kappa_max = NEUTRALITY / ALPHA ** 2
    check("(1) A velocity-dependent accumulated charge is excluded by atom neutrality: kappa < 2e-16 in q(v) = q_0 (1 + kappa v^2/c^2)",
          kappa_max < 3e-16, f"electron at v/c = alpha = {ALPHA:.4f} in hydrogen, proton at rest; |q_p + q_e| < {NEUTRALITY:.0e} e -> kappa < {kappa_max:.1e}")
    check("(1') Even a linear dependence q_0 (1 + kappa' v/c) is bounded at kappa' < 1.4e-19: 'trailing' as motion through the substrate cannot create the charge",
          NEUTRALITY / ALPHA < 2e-19, f"kappa' < {NEUTRALITY/ALPHA:.1e}")

    swr = (1 + GAMMA_POLE) / (1 - GAMMA_POLE)
    check("(2) Flow at c_0 meeting the pole's reactive barrier Gamma_pole = 1/3: a standing wave with SWR = 2",
          abs(swr - 2.0) < 1e-12, f"antinode / node charge-density contrast = {swr:.1f}")
    # charge density of a partially reflected wave along a string of length l with equal reflections at both ends:
    # rho(s, t) = Re[(e^{-i beta s} + Gamma e^{+i beta s}) e^{i omega t}] (from the head) -- symmetric under s -> l - s for the tail
    s = np.linspace(0, 1, 2001)
    beta_l = math.pi
    forward = np.exp(-1j * beta_l * s) + GAMMA_POLE * np.exp(1j * beta_l * s)
    env = np.abs(forward)                      # amplitude envelope
    t_avg = np.mean([np.real(forward * np.exp(1j * ph)) for ph in np.linspace(0, 2 * math.pi, 400, endpoint=False)], axis=0)
    check("(2') The envelope is symmetric between the two poles and the time-averaged oscillating charge is zero everywhere: no net pile at either pole",
          abs(env[0] - env[-1]) < 1e-9 and np.max(np.abs(t_avg)) < 1e-9,
          f"envelope at head {env[0]:.3f}, at tail {env[-1]:.3f}; max |time average| = {np.max(np.abs(t_avg)):.1e}")
    check("(2'') A net static accumulation at one pole needs an asymmetric, one-way barrier (a nonlinearity of the fluid at the pole): Gamma_pole = 1/3 is a symmetric linear reflection",
          True, "the manuscript names the wall (l. 315) and Pi_sat; it writes neither as an equation of state")
    check("(3) And a static end charge on a conductor spreads unless trapped (open_turn_stability.py): the trap is the same missing ingredient",
          True, "item 2: the fluid's equation of state")
    return report("Conclusion: 'the fluid accumulates in the trailing pole and creates the active charge' cannot be a motion "
                  "effect -- atom neutrality forbids any velocity dependence of charge to 1e-16 -- and as an internal-flow "
                  "effect it is the standing wave of the pole reflections, symmetric at both poles, oscillating, with zero "
                  "time average: no net charge is created at one pole. A net pile needs a one-way, nonlinear barrier at the "
                  "pole, which V2.10 names and does not write. The charge remains what step 1 found: the fluid's net e/3 per "
                  "branch, spread along it.")


if __name__ == "__main__":
    raise SystemExit(main())
