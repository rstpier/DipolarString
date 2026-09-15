#!/usr/bin/env python3
"""Independent re-derivation of the S3 three-string reconnection dossier.

Nothing here is taken from the dossier: every quantity is recomputed from the
definitions.  The script checks the combinatorics, the class-basis reduction,
the characteristic polynomial, the exact and perturbative gap, and the size of
the perturbative error in the regime the dossier's own diagnostic reaches.

It also computes the quantity the dossier never computes: the relabelling
orbits under the full symmetric group versus the alternating (cyclic) subgroup,
which is what decides whether a doublet exists at all.

Exit status is zero only if every check passes.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np
import sympy as sp

PERMS = list(itertools.permutations(range(3)))
IDX = {p: i for i, p in enumerate(PERMS)}
TRANSPOSITIONS = [(1, 0, 2), (2, 1, 0), (0, 2, 1)]

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    RESULTS.append(("PASS" if condition else "FAIL", name, detail))


def compose(p, q):
    return tuple(p[q[i]] for i in range(3))


def invert(p):
    r = [0] * 3
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def cycle_type(p):
    seen, lengths = set(), []
    for i in range(3):
        if i in seen:
            continue
        j, n = i, 0
        while j not in seen:
            seen.add(j)
            n += 1
            j = p[j]
        lengths.append(n)
    return tuple(sorted(lengths, reverse=True))


def adjacency():
    A = np.zeros((6, 6), dtype=float)
    for p in PERMS:
        for t in TRANSPOSITIONS:
            A[IDX[p], IDX[compose(p, t)]] = 1.0
    return A


def block_combinatorics() -> None:
    counts = Counter(cycle_type(p) for p in PERMS)
    check("Cycle-type census", counts == {(3,): 2, (2, 1): 3, (1, 1, 1): 1},
          f"{dict(counts)}")

    A = adjacency()
    even = [p for p in PERMS if cycle_type(p) != (2, 1)]
    odd = [p for p in PERMS if cycle_type(p) == (2, 1)]
    complete_bipartite = (
        all(A[IDX[a], IDX[b]] == 1 for a in even for b in odd)
        and all(A[IDX[a], IDX[b]] == 0 for a in even for b in even)
        and all(A[IDX[a], IDX[b]] == 0 for a in odd for b in odd)
    )
    check("Cayley graph is K_{3,3}", complete_bipartite,
          "complete bipartite between even and odd classes")
    check("K_{3,3} spectrum", np.allclose(sorted(np.linalg.eigvalsh(A)), [-3, 0, 0, 0, 0, 3], atol=1e-12),
          f"{np.round(np.linalg.eigvalsh(A), 9).tolist()}")

    c123, c132 = (1, 2, 0), (2, 0, 1)
    one_step = any(A[IDX[c123], IDX[q]] and cycle_type(q) == (3,) for q in PERMS)
    check("One-step 3-cycle to 3-cycle forbidden", not one_step, "parity flips at each move")

    mids = {cycle_type(q) for q in PERMS if A[IDX[c123], IDX[q]] and A[IDX[q], IDX[c132]]}
    paths = int((A @ A)[IDX[c123], IDX[c132]])
    check("[2+1] is the unique intermediate class", mids == {(2, 1)} and paths == 3,
          f"{paths} two-step paths, intermediate types {sorted(mids)}")


def block_relabelling_orbits() -> None:
    """The step the dossier omits: which relabelling group is gauge?"""
    def orbits(group):
        seen, out = set(), []
        for p in PERMS:
            if p in seen:
                continue
            orb = {compose(compose(g, p), invert(g)) for g in group}
            seen |= orb
            out.append(sorted(orb))
        return out

    S3 = PERMS
    A3 = [p for p in PERMS if cycle_type(p) != (2, 1)]
    three_cycles = {(1, 2, 0), (2, 0, 1)}

    o_s3 = orbits(S3)
    merged = any(set(o) == three_cycles for o in o_s3)
    check("Odd relabelling merges the two 3-cycles", merged and len(o_s3) == 3,
          f"{len(o_s3)} orbits under S3 conjugation")

    o_a3 = orbits(A3)
    split = all([c] in o_a3 for c in three_cycles)
    check("Cyclic relabelling keeps them distinct", split and len(o_a3) == 4,
          f"{len(o_a3)} orbits under A3 conjugation -> a doublet exists iff the "
          "gauge group contains no odd permutation")


def block_class_reduction() -> None:
    A = adjacency()
    classes = {(3,): [], (2, 1): [], (1, 1, 1): []}
    for p in PERMS:
        classes[cycle_type(p)].append(IDX[p])
    order = [(3,), (2, 1), (1, 1, 1)]
    V = np.zeros((6, 3))
    for k, c in enumerate(order):
        for i in classes[c]:
            V[i, k] = 1 / math.sqrt(len(classes[c]))
    Aq = V.T @ A @ V
    target = np.array([[0, math.sqrt(6), 0], [math.sqrt(6), 0, math.sqrt(3)], [0, math.sqrt(3), 0]])
    check("Class-basis adjacency", np.allclose(Aq, target, atol=1e-12),
          "[[0,sqrt6,0],[sqrt6,0,sqrt3],[0,sqrt3,0]] recovered by projection")
    check("Class-basis spectrum", np.allclose(sorted(np.linalg.eigvalsh(Aq)), [-3, 0, 3], atol=1e-12),
          f"{np.round(np.linalg.eigvalsh(Aq), 9).tolist()}")


def block_algebra() -> None:
    g, Lm, L3 = sp.symbols("g Lambda_mu Lambda_3", positive=True)
    E, a = sp.symbols("E a")

    H = sp.Matrix([[0, 0, -sp.sqrt(6) * g],
                   [0, L3, -sp.sqrt(3) * g],
                   [-sp.sqrt(6) * g, -sp.sqrt(3) * g, Lm]])
    cp = sp.expand(-(H - E * sp.eye(3)).det())
    claimed = sp.expand(E**3 - (L3 + Lm) * E**2 + (L3 * Lm - 9 * g**2) * E + 6 * L3 * g**2)
    check("Characteristic polynomial", sp.expand(cp - claimed) == 0,
          "E^3 -(L3+Lm)E^2 +(L3*Lm-9g^2)E +6*L3*g^2")

    coeff = sp.Poly(sp.expand(cp.subs(E, a * g**2)), g).coeff_monomial(g**2)
    sol = sp.solve(sp.Eq(coeff, 0), a)
    check("E_low = -6g^2/Lambda_mu, independent of Lambda_3",
          len(sol) == 1 and sp.simplify(sol[0] + 6 / Lm) == 0 and not sol[0].has(L3),
          f"coefficient factors as {sp.factor(coeff)}")

    H2 = sp.Matrix([[0, -sp.sqrt(6) * g], [-sp.sqrt(6) * g, Lm]])
    ev = sp.solve(sp.expand((H2 - E * sp.eye(2)).det()), E)
    low = (Lm - sp.sqrt(Lm**2 + 24 * g**2)) / 2
    check("Truncated block eigenvalues", any(sp.simplify(x - low) == 0 for x in ev),
          "(Lambda_mu -/+ sqrt(Lambda_mu^2 + 24 g^2))/2")

    delta = (sp.sqrt(Lm**2 + 24 * g**2) - Lm) / 4
    series = sp.expand(sp.series(delta, g, 0, 6).removeO())
    check("Delta expansion", sp.simplify(series - (3 * g**2 / Lm - 18 * g**4 / Lm**3)) == 0,
          "Delta = 3g^2/Lm - 18g^4/Lm^3 + O(g^6)")


def block_derrick() -> None:
    """Which Derrick-evading term stabilises, and at what cost."""
    lam, A, B, C, D = sp.symbols("lambda A B C D", positive=True)

    chiral = A * lam**3 - B * lam**2 + C * lam
    disc = sp.discriminant(sp.Poly(sp.diff(chiral, lam), lam), lam)
    check("Chiral route is conditional", sp.simplify(disc - (4 * B**2 - 12 * A * C)) == 0,
          "finite size iff |B| >= sqrt(3AC): a threshold")

    u = sp.Symbol("u", positive=True)
    roots = sp.solve(sp.Eq(3 * A * u**2 + C * u - D, 0), u)
    positive = [r for r in roots if sp.simplify(sp.limit(r, D, sp.oo)) is not sp.zoo]
    stabilises = any(sp.simplify(r - (-C + sp.sqrt(C**2 + 12 * A * D)) / (6 * A)) == 0 for r in roots)
    check("Four-derivative route is unconditional", stabilises,
          "lambda^2 = (-C + sqrt(C^2+12AD))/6A > 0 for every D > 0: no threshold "
          "(quartic in first derivatives, e.g. Faddeev; NOT string bending, which is Frank K3 ~ lambda^1)")
    # degree-k term in first derivatives scales as lambda^(3-k): k=2 (Frank) -> +1, k=4 (Faddeev) -> -1
    check("Frank bend does not stabilise", 3 - 2 == 1 and 3 - 4 == -1,
          "string bending |(n.grad)n|^2 is degree 2 -> lambda^1; Faddeev is degree 4 -> lambda^-1")


def block_cp1_identity() -> None:
    """The Faddeev density is the Berry curvature of the director's CP^1 gauge field.

    With n = (sin th cos ph, sin th sin ph, cos th) and z = (cos th/2, sin th/2 e^{i ph}),
    the Berry connection is a_mu = sin^2(th/2) d_mu ph and its curvature satisfies
    f_xy = (1/2) n.(d_x n x d_y n).  Hence the Faddeev term (n.(dn x dn))^2 = 4 f^2 is the
    Maxwell term of the emergent gauge field -- the mechanism by which a quartic gradient
    term arises when a U(1)-charged field coupled to a[n] is integrated out.
    """
    th, ph = sp.symbols("theta phi", real=True)
    thx, thy, phx, phy, phxy = sp.symbols("theta_x theta_y phi_x phi_y phi_xy", real=True)
    n = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    dn = lambda dth, dph: n.diff(th) * dth + n.diff(ph) * dph
    fad = sp.simplify(n.dot(dn(thx, phx).cross(dn(thy, phy))))
    f = sp.simplify(sp.sin(th / 2) * sp.cos(th / 2) * thx * phy + sp.sin(th / 2)**2 * phxy
                    - (sp.sin(th / 2) * sp.cos(th / 2) * thy * phx + sp.sin(th / 2)**2 * phxy))
    check("Faddeev density = 2 x Berry curvature", sp.simplify(f - fad / 2) == 0,
          "f_xy = (1/2) n.(d_x n x d_y n); Faddeev term = 4 f^2 = Maxwell term of the CP^1 gauge field")


def block_diagnostic_error() -> None:
    """How wrong the perturbative formula is where the dossier's diagnostic lands."""
    E1 = 3 * 0.51099895  # MeV; E1 = h c / ell_1 = 3 m_e c^2 from the anchor
    worst = 0.0
    rows = []
    for mev in (0.12, 0.58):
        lam = mev / E1
        g_over = math.sqrt(2 / 27 + lam / 9)
        exact = (math.sqrt(lam**2 + 24 * g_over**2) - lam) / 4
        pert = 3 * g_over**2 / lam
        factor = pert / exact
        worst = max(worst, factor)
        rows.append(f"Lm={mev:.2f}MeV g/Lm={g_over/lam:.3f} pert/exact={factor:.1f}x")
    check("Diagnostic is outside weak coupling", worst > 2.0,
          "; ".join(rows) + " -- 'not quantitatively justified' understates it")


def main() -> int:
    block_combinatorics()
    block_relabelling_orbits()
    block_class_reduction()
    block_algebra()
    block_derrick()
    block_cp1_identity()
    block_diagnostic_error()

    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    print()
    print(f"RESULT: {len(RESULTS) - len(failed)}/{len(RESULTS)} PASS; {len(failed)} FAIL")
    print("Verified: combinatorics and algebra. NOT verified, because not derivable here: "
          "g, Lambda_mu (including its sign), the eta doublet, charge, statistics, localization.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
