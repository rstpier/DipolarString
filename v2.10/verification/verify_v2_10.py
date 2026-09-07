#!/usr/bin/env python3
"""Verification for DS V2.10.

V2.10 is an incremental package.  It changes three sentences of the manuscript
and nothing else: no equation, no numerical value, no figure, no script.  The
numerical suite is therefore unchanged and is the V2.9.2 one, run here from the
archived V2.9 package; this script adds the four manuscript-hygiene checks,
re-run against the V2.10 source, plus three checks specific to the V2.10
corrections.

Exit status is zero only if the delegated suite passes and every local check
passes.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript" / "DS_model_V2_10_Zenodo.tex"
V29 = ROOT.parent / "v2.9"
DELEGATED = V29 / "verification" / "verify_v291.py"

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    RESULTS.append(("PASS" if condition else "FAIL", name, detail))


def run_delegated() -> tuple[bool, str]:
    """Run the unchanged V2.9.2 numerical suite from the archived package."""
    if not DELEGATED.is_file():
        return False, f"not found: {DELEGATED}"
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run(
        [sys.executable, str(DELEGATED)],
        cwd=str(V29), env=env, capture_output=True, text=True,
    )
    tail = proc.stdout.strip().splitlines()
    summary = next((ln for ln in reversed(tail) if ln.startswith("RESULT:")), "no RESULT line")
    return proc.returncode == 0, summary


def block_manuscript_hygiene(source: str) -> None:
    """The four V2.9.2 hygiene checks, re-run against the V2.10 source."""
    check("Low-impedance equation withdrawn", "eq:vsync" not in source,
          "no synchronization-speed equation label")
    check("A8 source printed", "u - u_0" in source, "corrected source present")
    check("Conditional phase equation printed",
          "eq:phase_action" in source and "eq:z3_coefficient" in source,
          "phase normalization and Z3 coefficient present")
    check("Instanton-to-G gap disclosed", "The two missing bridges to $G$" in source,
          "missing gravitational map stated")


def block_v2_10_corrections(source: str) -> None:
    """The three corrections this revision makes, and the version stamp."""
    check("Version stamped V2.10", "V2.10 (Zenodo edition)" in source,
          "title-page date line")
    check("Primary-target sentence corrected",
          "remains the model's primary target" not in source
          and "absence of PBH evaporation bursts is now the model's primary target" in source,
          "shadow-diameter no longer named as primary target")
    check("Maturity-level count consistent",
          "six levels of theoretical maturity, and records separately" in source,
          "seventh heading declared as a withdrawal record, not a level")
    check("Weave-chirality entry no longer titled resolved",
          "Status of the weave chirality (chiral closure dissolved; the converse question open)" in source
          and "(resolved, with a residual question)" not in source,
          "converse question stated as open")
    check("Document History records V2.10",
          "\\paragraph{V2.10 (September 2026).}" in source, "revision traced")


def main() -> int:
    if not MANUSCRIPT.is_file():
        print(f"FAIL  manuscript not found: {MANUSCRIPT}")
        return 1
    source = MANUSCRIPT.read_text(encoding="utf-8")

    print("=" * 72)
    print("Delegated numerical suite (unchanged from V2.9.2)")
    print("=" * 72)
    ok, summary = run_delegated()
    print(f"  {'PASS' if ok else 'FAIL'}  {DELEGATED.relative_to(ROOT.parent)}: {summary}")

    block_manuscript_hygiene(source)
    block_v2_10_corrections(source)

    print()
    print("=" * 72)
    print("V2.10 source checks")
    print("=" * 72)
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")

    failed = [r for r in RESULTS if r[0] == "FAIL"]
    total = len(RESULTS)
    print()
    print("=" * 72)
    print(f"RESULT: {total - len(failed)}/{total} local PASS; delegated suite "
          f"{'PASS' if ok else 'FAIL'} ({summary})")
    return 0 if (ok and not failed) else 1


if __name__ == "__main__":
    raise SystemExit(main())
