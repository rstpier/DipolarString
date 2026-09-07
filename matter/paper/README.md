# Article — matter sector

`DS_matter_sector.tex` is a skeleton, not a draft. It compiles clean (2 pages,
zero undefined references) so that the first content pass is content and not
plumbing.

## What is already set up

- The **preamble is byte-identical** to the V2.10 manuscript's, up to the
  hyperref block: same document class, packages, theorem environments, TikZ
  libraries and the `\dqd` / `\Ncd` macros. It builds with the same toolchain
  (`latexmk -pdf`), so nothing has to be rediscovered.
- A `\status{...}` macro. **Every claim carries its maturity level** — one of
  the six: *Derived*, *Conditionally derived*, *Calibrated*, *Postulated*,
  *Reparametrized*, *Closed*. A ratio obtained by fitting a measured mass is
  *Reparametrized*.
- A section skeleton in the order the argument has to run, each section carrying
  as comments the facts and traps it must respect, taken from
  [`../CONSTRAINTS.md`](../CONSTRAINTS.md).
- A bibliography with the two entries any such paper needs.

## Two things to settle before writing

**The scope section comes first, and is not optional.** The main manuscript
withdrew its muon, tau, W, Z and H table entries as *reverse-engineered mass fits
without predictive content*. A matter-sector paper that does not open by
delimiting itself repeats that failure in a new place.

**The provenance paragraph decides what may be claimed.** The main manuscript
states its own limit plainly: *both numerical implementations originate from this
same collaboration; independent third-party verification has not yet been
performed.* A second AI system changes that sentence only under one of two
readings, and they must not be conflated:

- the second system was given the DS manuscript and asked to extend it — the
  prior is shared, and this is **not** independent verification;
- the second system reproduced a result from the structure equations alone,
  without the manuscript's conclusions — closer to independence, and the protocol
  must be described precisely enough for a reader to judge.

Which one applies is a fact about how the work was actually produced. It belongs
in the Acknowledgments, stated, not implied.

## Build

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error DS_matter_sector.tex
```

Build products are gitignored: this is a draft. When it is deposited it moves to
a release folder and ships its build products like the other packages.
