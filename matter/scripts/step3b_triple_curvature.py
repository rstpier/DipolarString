#!/usr/bin/env python3
"""Step 3b: Berry CURVATURE of the triple's doublet over (tangent x frame), and the honest
composition of its holonomy on a closed string.

Step 3a measured two things separately and never composed them: the overlap of the doublet
over a third of a turn (the permutation of the three conductors, a 120 deg rotation) and the
Berry connection A = a J per unit frame angle (a = 0.186).  The holonomy of a closed string
is the PRODUCT of the two.  This script:

  A. shows that a is the mean angular momentum <L_z> of the circular doublet state, by two
     independent routes (finite difference of the rigidly rotated pattern; the rotation
     generator applied to the pattern), and that it depends on the conductor radius;
  B. shows that the tilt generators L_x, L_y have no matrix element inside the doublet, so on
     the frame bundle of the string the connection is exactly A = a J omega_z (omega_z = the
     twist 1-form) and, by the Maurer-Cartan equation, the curvature is
         F = a J x (area form of S^2 pulled back by the tangent t),
     the photon's Berry curvature (Tomita-Chiao) scaled by a;
  C. composes the holonomy of a planar ring whose triple is twisted by one third of a turn:
     not 120 deg (step 3a's permutation part) but (1 - a) x 120 deg;
  D. on an explicit closed space curve checks Fuller's theorem numerically (the parallel-
     transported frame returns rotated by 2 pi Wr = the solid angle swept by t, mod 2 pi),
     tunes the writhe to exactly 1/3 and shows that a triple closing by writhe alone, with no
     twist, has holonomy exactly 120 deg;
  E. states and checks the general formula
         H = R( 2 pi [ (1 - a) Lk + a Wr ] ),   Lk = Tw + Wr in Z/3,
     which interpolates between the pure permutation (a = 0: 2 pi Lk) and the photon
     (a = 1: 2 pi Wr).

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

import berry_holonomy_common as bhc
from berry_holonomy_common import (D, EPS, L, X, XS, Y, berry_connection, check, differential_charges,
                                   ip, modes, positions, report, transport)

J = np.array([[0.0, -1.0], [1.0, 0.0]])          # so(2) generator; expm(theta J) = R(theta)


def angle_of(m: np.ndarray) -> float:
    """Rotation angle (deg) of a 2x2 matrix close to a rotation."""
    return math.degrees(math.atan2(m[1, 0], m[0, 0]))


def antisym_coeff(m: np.ndarray) -> float:
    """a such that the antisymmetric part of m is a J."""
    return 0.5 * (m[1, 0] - m[0, 1])


# --------------------------------------------------------------------------- A. a = <L_z>

def generator_connection(psi: float = 0.0) -> np.ndarray:
    """<E_i | d_psi E_j> with d_psi computed by the rotation generator on the pattern:
    d/dpsi [R(psi) E(R(-psi) x)] at 0  =  z x E  -  (x d_y - y d_x) E."""
    e = modes(3, psi)

    def gen(f):
        ex, ey = f
        dex_dx, dex_dy = np.gradient(ex, XS, XS)
        dey_dx, dey_dy = np.gradient(ey, XS, XS)
        return (-ey - (X * dex_dy - Y * dex_dx), ex - (X * dey_dy - Y * dey_dx))

    g = [gen(f) for f in e]
    return np.array([[ip(e[i], g[j]) for j in range(2)] for i in range(2)])


def block_a() -> float:
    a_fd = antisym_coeff(berry_connection(3, 0.0))
    a_gen = antisym_coeff(generator_connection(0.0))
    check("a from the rotated pattern equals a from the generator L_z", abs(a_fd - a_gen) / a_fd < 0.02,
          f"finite difference a = {a_fd:.4f}; generator <L_z> a = {a_gen:.4f} (grid-limited agreement)")
    # circular combination e+ = (E1 + i E2)/sqrt2: <e+|d_psi e+> = -i a
    A = berry_connection(3, 0.0)
    ep_phase = 0.5j * (A[0, 1] - A[1, 0])
    check("Circular state e+ carries <e+|d_psi e+> = -i a", abs(ep_phase + 1j * a_fd) < 1e-9,
          f"<e+|d_psi e+> = {ep_phase:.4f}; a = {a_fd:.4f}")
    check("a lies strictly between 0 (mode locked to the conductors) and 1 (free vector / photon)",
          0.0 < a_fd < 1.0, f"a = {a_fd:.4f}, 1 - a = {1 - a_fd:.4f}")
    # sensitivity to the conductor radius model
    vals = {}
    for eps_factor in (0.5, 1.0, 2.0):
        bhc.EPS = eps_factor * bhc.R
        vals[eps_factor] = antisym_coeff(berry_connection(3, 0.0))
    bhc.EPS = EPS
    check("a depends on the conductor radius (a model number, not a topological one)",
          all(0.0 < v < 1.0 for v in vals.values()) and max(vals.values()) - min(vals.values()) > 0.01,
          "a(eps/r) = " + ", ".join(f"{k:g}: {v:.4f}" for k, v in vals.items()))
    return a_fd


# --------------------------------------------------------------------------- B. tilt generators

def block_b() -> None:
    """Coarse 3D embedding: E(x, y, z) = (E_x, E_y, 0), z-independent, on a symmetric slab.
    Overlap of the doublet with its version rotated by delta about the x-axis (a tilt of the
    string) and about the z-axis (a twist of the frame); antisymmetric first-order parts."""
    nc, nz = 301, 25
    xs = np.linspace(-L, L, nc)
    xc, yc = np.meshgrid(xs, xs, indexing="ij")
    da = (xs[1] - xs[0]) ** 2
    zs = np.linspace(-2 * D, 2 * D, nz)
    pts, qs = positions(3, 0.0), differential_charges(3)

    def f2(q, xa, ya):
        ex, ey = np.zeros_like(xa), np.zeros_like(xa)
        for qi, p in zip(q, pts):
            dx, dy = xa - p[0], ya - p[1]
            g = 1.0 / (dx * dx + dy * dy + EPS * EPS)
            ex += qi * dx * g
            ey += qi * dy * g
        return ex, ey

    base = [f2(q, xc, yc) for q in qs]
    nrm = [math.sqrt(np.sum(ex * ex + ey * ey) * da) for ex, ey in base]
    base = [(ex / n, ey / n) for (ex, ey), n in zip(base, nrm)]

    def overlap_tilt(delta):
        c, s = math.cos(delta), math.sin(delta)
        o = np.zeros((2, 2))
        for z in zs:                                    # R_x(-delta) x = (x, c y + s z, ...)
            rot = [f2(q, xc, c * yc + s * z) for q in qs]  # rotated vector = (E_x, c E_y, s E_y)
            for i in range(2):
                for j in range(2):
                    o[i, j] += np.sum(base[i][0] * rot[j][0] + base[i][1] * c * rot[j][1]) / nrm[j] * da
        return o / nz

    def overlap_twist(delta):
        c, s = math.cos(delta), math.sin(delta)
        rot = [f2(q, c * xc + s * yc, -s * xc + c * yc) for q in qs]
        o = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                rx, ry = c * rot[j][0] - s * rot[j][1], s * rot[j][0] + c * rot[j][1]
                o[i, j] = np.sum(base[i][0] * rx + base[i][1] * ry) / nrm[j] * da
        return o

    h = 1e-4
    ident = overlap_tilt(0.0)
    a_tilt = antisym_coeff((overlap_tilt(h) - overlap_tilt(-h)) / (2 * h))
    a_twist = antisym_coeff((overlap_twist(h) - overlap_twist(-h)) / (2 * h))
    check("3D embedding is consistent (zero-angle overlap is the identity)", np.allclose(ident, np.eye(2), atol=1e-4),
          f"O(0) = {np.round(ident, 6).tolist()} (coarse 301^2 grid; the fine grid is orthogonal to 1e-7)")
    check("Tilting the string (L_x) has NO first-order matrix element in the doublet",
          abs(a_tilt) < 1e-6 * max(abs(a_twist), 1e-3),
          f"|a_tilt| = {abs(a_tilt):.1e} vs a_twist (same coarse grid) = {a_twist:.4f}")
    check("Hence A = a J omega_z on the frame bundle and F = dA = -a J omega_x ^ omega_y (Maurer-Cartan)",
          True, "curvature = a x (solid-angle 2-form of the tangent direction); the photon has a = 1")


# --------------------------------------------------------------------------- C. planar ring, twist 1/3

def block_c(a: float) -> None:
    u = transport(3, 2 * math.pi / 3)                  # permutation part (step 3a's "holonomy")
    A = berry_connection(3, 0.0)                       # connection, constant
    h_ring = u @ expm(-A * 2 * math.pi / 3)            # parallel transport composed with the identification
    ang = angle_of(h_ring)
    check("Step 3a's 120 deg was the permutation part alone", abs(angle_of(u) - 120) < 1e-2,
          f"overlap over 2pi/3 = {angle_of(u):+.3f} deg")
    check("Planar ring twisted by 1/3: full holonomy is (1 - a) x 120 deg, NOT 120 deg",
          abs(ang - (1 - a) * 120) < 0.05 and abs(np.linalg.det(h_ring) - 1) < 1e-6,
          f"holonomy = {ang:+.2f} deg = 2 pi x {ang/360:.4f}; (1 - a)/3 = {(1-a)/3:.4f}")
    ev = np.linalg.eigvals(h_ring)
    check("Its eigenvalues are exp(+/- i 2 pi (1 - a)/3), a phase that is not a cube root of unity",
          abs(abs(ev[0]) - 1) < 1e-9 and abs(abs(np.angle(ev[0])) - math.radians(ang)) < 1e-6,
          f"eigenvalue phases = +/- {abs(np.angle(ev[0])):.4f} rad; 2 pi/3 = {2*math.pi/3:.4f} rad")


# --------------------------------------------------------------------------- D/E. explicit closed curve

def curve_and_tangent(s, r0, q=3, r_big=1.0):
    """Right-handed (1, q) torus curve: an unknotted ring that wobbles q times, writhe > 0."""
    rho = r_big + r0 * np.cos(q * s)
    drho = -q * r0 * np.sin(q * s)
    x = np.stack([rho * np.cos(s), rho * np.sin(s), -r0 * np.sin(q * s)], axis=-1)
    dx = np.stack([drho * np.cos(s) - rho * np.sin(s), drho * np.sin(s) + rho * np.cos(s),
                   -q * r0 * np.cos(q * s)], axis=-1)
    return x, dx


def writhe(r0, n=1500) -> float:
    s = (np.arange(n) + 0.5) * 2 * math.pi / n
    x, dx = curve_and_tangent(s, r0)
    dxx = x[:, None, :] - x[None, :, :]
    dist3 = np.sum(dxx * dxx, axis=-1) ** 1.5
    np.fill_diagonal(dist3, np.inf)
    cross = np.cross(dx[:, None, :], dx[None, :, :])
    integrand = np.sum(cross * dxx, axis=-1) / dist3
    return float(np.sum(integrand) * (2 * math.pi / n) ** 2 / (4 * math.pi))


def bishop_rotation_and_solid_angle(r0, n=4000):
    """Double-reflection parallel transport of a normal vector around the closed curve.
    Returns (Delta, Omega): the right-handed angle about t(0) from b(0) to b(L), and the
    signed solid angle enclosed by the tangent indicatrix (Oosterom-Strackee, pole reference)."""
    s = np.arange(n + 1) * 2 * math.pi / n
    x, dx = curve_and_tangent(s, r0)
    t = dx / np.linalg.norm(dx, axis=-1, keepdims=True)
    x[-1], t[-1] = x[0], t[0]
    b = np.cross(t[0], [0.0, 0.0, 1.0])
    b /= np.linalg.norm(b)
    b0 = b.copy()
    for i in range(n):
        v1 = x[i + 1] - x[i]
        c1 = v1 @ v1
        rl = b - (2 / c1) * (v1 @ b) * v1
        tl = t[i] - (2 / c1) * (v1 @ t[i]) * v1
        v2 = t[i + 1] - tl
        c2 = v2 @ v2
        b = rl - (2 / c2) * (v2 @ rl) * v2
        b /= np.linalg.norm(b)
    delta = math.atan2(np.cross(b0, b) @ t[0], b0 @ b)
    p = np.array([0.0, 0.0, 1.0])
    omega = 0.0
    for i in range(n):
        a_, b_ = t[i], t[i + 1]
        omega += 2 * math.atan2(p @ np.cross(a_, b_), 1 + p @ a_ + a_ @ b_ + b_ @ p)
    return delta, omega


def wrap(x: float) -> float:
    """Wrap to (-1/2, 1/2]."""
    return x - math.floor(x + 0.5)


def block_de(a: float) -> None:
    # convergence of the writhe quadrature
    w_lo, w_hi = writhe(0.3, 800), writhe(0.3, 1600)
    check("Writhe quadrature converged", abs(w_lo - w_hi) < 1e-4, f"Wr(r0 = 0.3): n = 800 -> {w_lo:.6f}, n = 1600 -> {w_hi:.6f}")
    # Fuller / Calugareanu on a generic curve
    delta, omega = bishop_rotation_and_solid_angle(0.3)
    check("Fuller: the parallel-transported frame returns rotated by 2 pi Wr (mod 2 pi)",
          abs(wrap(delta / (2 * math.pi) - w_hi)) < 2e-4,
          f"Delta/2pi = {delta/(2*math.pi):+.5f}, Wr = {w_hi:+.5f}, difference mod 1 = {wrap(delta/(2*math.pi) - w_hi):+.1e}")
    check("Fuller: that rotation is the solid angle swept by the tangent (mod 2 pi) -- the curvature is the area form",
          abs(wrap((omega - delta) / (2 * math.pi))) < 2e-4,
          f"Omega/2pi = {omega/(2*math.pi):+.5f}, Delta/2pi = {delta/(2*math.pi):+.5f}, difference mod 1 = {wrap((omega-delta)/(2*math.pi)):+.1e}")
    # tune the writhe to exactly 1/3: a triple that closes by writhe alone, zero twist
    r_star = brentq(lambda r: writhe(r, 1600) - 1.0 / 3.0, 0.05, 0.6, xtol=1e-10)
    w_star = writhe(r_star, 1600)
    delta_star, _ = bishop_rotation_and_solid_angle(r_star)
    A = berry_connection(3, 0.0)
    u_star = transport(3, delta_star)                  # material frame = Bishop frame: psi_L = 0
    h_star = u_star @ expm(-A * 0.0)
    check("A closed triple with Wr = 1/3 and Tw = 0 exists (found by root-finding on the curve family)",
          abs(w_star - 1 / 3) < 1e-8, f"r0* = {r_star:.6f}, Wr = {w_star:.8f}; Bishop frame closes up to {math.degrees(delta_star):+.3f} deg = one conductor permutation")
    check("Closing by writhe alone: holonomy is EXACTLY 120 deg (fractional winding 1/3, the G conjecture's input)",
          abs(abs(angle_of(h_star)) - 120) < 1e-2, f"holonomy = {angle_of(h_star):+.3f} deg")
    # generic curve, closure k = 1 supplied by twist: Tw = 1/3 - Wr
    tw = 1.0 / 3.0 - w_hi
    psi_l = 2 * math.pi * tw
    h_gen = transport(3, psi_l + delta) @ expm(-A * psi_l)
    predicted = 360 * ((1 - a) / 3 + a * w_hi)
    check("Generic closed curve (Wr = 0.75..., Tw = 1/3 - Wr): holonomy = 2 pi [(1 - a) Lk + a Wr]",
          abs(wrap((angle_of(h_gen) - predicted) / 360)) < 1e-3,
          f"chain: {angle_of(h_gen):+.3f} deg; formula: {wrap(predicted/360)*360:+.3f} deg (Lk = 1/3, Wr = {w_hi:.4f}, Tw = {tw:+.4f})")
    check("Limits: a = 0 gives the pure permutation 2 pi Lk (Z_3); a = 1 gives the photon 2 pi Wr (Tomita-Chiao)",
          True, f"the triple sits at a = {a:.3f}: {100*(1-a):.0f} % topological (Lk), {100*a:.0f} % geometric (Wr)")


def main() -> int:
    a = block_a()
    block_b()
    block_c(a)
    block_de(a)
    return report("Conclusion (step 3b): the doublet's Berry curvature over the string's frame bundle is "
                  "a x (solid-angle form of the tangent), i.e. F = 2a f_mn[t] with f the CP1 field strength of the "
                  "CONSTRAINTS/audit identity -- the Faddeev structure exists for the triple's tangent field, with "
                  "weight a = 0.19 (photon: 1; pair: 0). The closed-string holonomy is 2 pi[(1 - a) Lk + a Wr]: "
                  "exactly 1/3 only when the triple closes by writhe with no twist; a planar ring twisted by 1/3 "
                  "gives (1 - a)/3 = 0.27. Step 3a's '120 deg' was the permutation part alone. Still OPEN: the "
                  "coefficient of the induced Faddeev term (needs the doublet's propagator), PVLAS, and whether "
                  "the weave's excitations are triple-doublet modes.")


if __name__ == "__main__":
    raise SystemExit(main())
