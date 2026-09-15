#!/usr/bin/env python3
"""Draws matter/figures/matter_topology.png: the proposed matter topology of the DS medium,
after Levin-Wen (string ends as fermions) and Bilson-Thompson (braids of three twisted ribbons)."""

from __future__ import annotations

import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

C = {"bg": "#fbfaf7", "ink": "#1d1d1b", "s1": "#3b6ea5", "s2": "#c1553b", "s3": "#3f8f5a", "grey": "#9a9a94", "gold": "#b8860b"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": C["grey"]})

fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.6), facecolor=C["bg"])
fig.subplots_adjust(left=0.03, right=0.97, top=0.875, bottom=0.03, wspace=0.1, hspace=0.34)
fig.suptitle("Topologie matière proposée pour le milieu à cordes dipolaires\n"
             "d'après Levin–Wen (bouts de cordes = fermions) et Bilson-Thompson (tresses de trois rubans tordus)",
             fontsize=12, color=C["ink"], y=0.985)


def braid(ax, x0, x1, y0, dy, twists, colors, n_cross=2, lw=2.2):
    """Three strands braided between x0 and x1 (schematic), with chirality marks."""
    n = 3
    xs = np.linspace(x0, x1, 400)
    ys = np.zeros((n, len(xs)))
    order = list(range(n))
    seg = (x1 - x0) / (n_cross + 1)
    for k in range(n):
        ys[k] = y0 + (k - 1) * dy
    # build crossing schedule: swap adjacent strands sequentially
    pos = np.tile(np.arange(n)[:, None], (1, len(xs))).astype(float)
    for c in range(n_cross):
        i = c % 2                       # swap strands at positions i, i+1
        xa, xb = x0 + seg * (c + 0.5), x0 + seg * (c + 1.5)
        mask = (xs >= xa) & (xs <= xb)
        t = (xs[mask] - xa) / (xb - xa)
        smooth = 0.5 - 0.5 * np.cos(math.pi * t)
        # which strand is at position i / i+1 at xa
        cur = pos[:, mask][:, 0]
        a = int(np.where(cur == i)[0][0]); b = int(np.where(cur == i + 1)[0][0])
        pos[a, mask] = i + smooth
        pos[b, mask] = i + 1 - smooth
        pos[a, xs > xb] = i + 1
        pos[b, xs > xb] = i
    for k in range(n):
        yk = y0 + (pos[k] - 1) * dy
        ax.plot(xs, yk, color=colors[k], lw=lw, solid_capstyle="round", zorder=3 + (k % 2))
        # chirality marks (twist) along the strand
        for xm in np.linspace(x0 + 0.08 * (x1 - x0), x1 - 0.08 * (x1 - x0), 3):
            j = int(np.argmin(np.abs(xs - xm)))
            tw = twists[k]
            if tw != 0:
                ax.text(xs[j], yk[j] + 0.16 * dy * (1 if k != 1 else -1.0), "+" if tw > 0 else "−",
                        ha="center", va="center", fontsize=9, color=colors[k], fontweight="bold", zorder=6)
    return pos


# ---------------------------------------------------------------- (1) vacuum and pair creation
ax = axes[0, 0]; ax.set_facecolor(C["bg"]); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 7); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("1 · Le vide : un réseau de cordes fermées.\nCouper une corde crée deux bouts, de charges opposées",
             fontsize=9.5, loc="left", color=C["ink"])
for y in (1.5, 3.0, 4.5, 6.0):
    ax.plot([0.5, 9.5], [y, y], color=C["grey"], lw=1.4, zorder=1)
for x in (2.0, 4.0, 6.0, 8.0):
    ax.plot([x, x], [0.7, 6.6], color=C["grey"], lw=1.4, zorder=1)
# a cut on the line y = 3 between x = 4 and 6
ax.plot([4.55, 5.45], [3.0, 3.0], color=C["bg"], lw=4, zorder=2)
ax.plot(4.5, 3.0, "o", color=C["s2"], ms=9, zorder=4); ax.plot(5.5, 3.0, "o", color=C["s1"], ms=9, zorder=4)
ax.text(4.5, 3.45, "+e/3", ha="center", color=C["s2"], fontsize=9, fontweight="bold")
ax.text(5.5, 3.45, "−e/3", ha="center", color=C["s1"], fontsize=9, fontweight="bold")
ax.annotate("", xy=(4.35, 3.0), xytext=(3.6, 3.0), arrowprops=dict(arrowstyle="->", color=C["s2"], lw=1.5))
ax.annotate("", xy=(5.65, 3.0), xytext=(6.4, 3.0), arrowprops=dict(arrowstyle="<-", color=C["s1"], lw=1.5))
ax.text(5.0, 0.75, "Pas de bout ⇒ pas de charge : le vide est neutre par construction.\n"
                  "Un bout de corde = une source de divergence, de signe fixé par la chiralité.\n"
                  "Création par paires miroir · annihilation = recollement · P agit comme C.",
        ha="center", va="center", fontsize=7.9, color=C["ink"],
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["grey"], lw=0.8))

# ---------------------------------------------------------------- (2) the electron
ax = axes[0, 1]; ax.set_facecolor(C["bg"]); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 7); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("2 · L'électron : une tresse de trois cordes de même chiralité,\nson bout replié en spire partiellement ouverte",
             fontsize=9.5, loc="left", color=C["ink"])
ax.text(0.6, 3.5, "vers\nla trame", ha="center", va="center", fontsize=8, color=C["grey"])
braid(ax, 1.3, 6.0, 3.5, 0.55, twists=(-1, -1, -1), colors=(C["s1"], C["s1"], C["s1"]), n_cross=3)
# the folded end: an open turn
th = np.linspace(-0.15 * math.pi, 1.75 * math.pi, 300)
Rr = 1.35; cx, cy = 7.6, 3.5
for off in (-0.13, 0.0, 0.13):
    ax.plot(cx + (Rr + off) * np.cos(th), cy + (Rr + off) * np.sin(th), color=C["s1"], lw=1.8, zorder=3)
ax.plot([6.0, cx + Rr * math.cos(th[0])], [3.5, cy + Rr * math.sin(th[0])], color=C["s1"], lw=1.8, alpha=0)
ax.plot(cx + Rr * math.cos(th[-1]), cy + Rr * math.sin(th[-1]), "o", color=C["s2"], ms=8, zorder=5)
ax.text(cx + 0.2, cy - Rr - 0.55, "bout libre = charge −e", fontsize=8, color=C["s2"], ha="center")
ax.annotate("", xy=(cx + Rr, cy), xytext=(cx, cy), arrowprops=dict(arrowstyle="<->", color=C["gold"], lw=1.2))
ax.text(cx + 0.55, cy + 0.2, "R", color=C["gold"], fontsize=9)
ax.text(5.0, 0.75, "Charge : 3 × (−e/3), une par corde, portée par le bout. Chiralité : sens de la tresse.\n"
                  "Masse : énergie du mode du bout replié, ∝ 1/ℓ (spire ouverte : βℓ = π, dérivé).\n"
                  "Spin : demi-torsion ⇒ mode antipériodique, L_z = n + ½, signe −1 sous 2π (un axe).",
        ha="center", va="center", fontsize=7.9, color=C["ink"],
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["grey"], lw=0.8))

# ---------------------------------------------------------------- (3) catalogue
ax = axes[1, 0]; ax.set_facecolor(C["bg"]); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 7); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("3 · Catalogue de première génération (Bilson-Thompson) :\ncharge = somme des chiralités / 3",
             fontsize=9.5, loc="left", color=C["ink"])
entries = [("ν", (0, 0, 0), "0"), ("e⁻", (-1, -1, -1), "−1"), ("u", (1, 1, 0), "+2/3"), ("d", (-1, 0, 0), "−1/3")]
for i, (name, tw, q) in enumerate(entries):
    y = 6.1 - i * 1.4
    cols = [C["s2"] if t > 0 else (C["s1"] if t < 0 else C["grey"]) for t in tw]
    braid(ax, 1.2, 4.6, y, 0.32, twists=tw, colors=cols, n_cross=2, lw=1.8)
    ax.text(0.5, y, name, ha="center", va="center", fontsize=11, color=C["ink"], fontweight="bold")
    ax.text(5.1, y, f"q = {q}", va="center", fontsize=9, color=C["ink"])
    if name in ("u", "d"):
        ax.text(6.4, y, "3 couleurs = position\nde la corde impaire", va="center", fontsize=7.6, color=C["s3"])
ax.text(5.0, 0.7, "Antiparticule = tresse miroir (chiralités et sens de tresse inversés).\nPhoton = onde sur une corde sans bout (natif CD).",
        ha="center", va="center", fontsize=7.9, color=C["ink"],
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["grey"], lw=0.8))

# ---------------------------------------------------------------- (4) proton
ax = axes[1, 1]; ax.set_facecolor(C["bg"]); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 7); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("4 · Le proton : trois tresses ouvertes (u, u, d) jointes à leurs bouts,\nun graphe Θ : le confinement est topologique",
             fontsize=9.5, loc="left", color=C["ink"])
jx1, jx2, jy = 2.2, 7.8, 3.6
for k, (lab, tw, yoff) in enumerate((("u", (1, 1, 0), 1.6), ("u", (1, 1, 0), 0.0), ("d", (-1, 0, 0), -1.6))):
    cols = [C["s2"] if t > 0 else (C["s1"] if t < 0 else C["grey"]) for t in tw]
    # curved paths from junction 1 to junction 2
    xs = np.linspace(jx1, jx2, 400)
    bump = yoff * np.sin(math.pi * (xs - jx1) / (jx2 - jx1))
    for j, t in enumerate(tw):
        ax.plot(xs, jy + bump + (j - 1) * 0.16, color=cols[j], lw=1.7, zorder=3)
    ax.text(5.0, jy + yoff * 1.0 + (0.45 if yoff >= 0 else -0.5), lab, ha="center", fontsize=11, fontweight="bold", color=C["ink"])
for jx in (jx1, jx2):
    ax.plot(jx, jy, "o", color=C["ink"], ms=10, zorder=6)
ax.text(jx1 - 0.2, jy - 0.55, "jonction", ha="center", fontsize=8, color=C["grey"])
ax.text(5.0, 0.75, "Un bout de tresse ouverte doit finir sur une autre tresse : les quarks n'existent pas seuls.\n"
                  "Charge totale +2/3 +2/3 −1/3 = +1.\n"
                  "Masse : hors du modèle (pas de secteur fort ; CONSTRAINTS §1a, §1c).",
        ha="center", va="center", fontsize=7.9, color=C["ink"],
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["grey"], lw=0.8))

out = "/home/user/DipolarString/matter/figures/matter_topology.png"
fig.savefig(out, dpi=170, facecolor=C["bg"])
print("wrote", out)
