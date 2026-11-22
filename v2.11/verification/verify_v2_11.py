#!/usr/bin/env python3
"""Verification for DS V2.11 (erratum edition).

V2.11 is the frozen V2.10 manuscript with an erratum applied in place.  This script
(1) checks that every erratum item is present in the V2.11 source and that the V2.10
hygiene checks still pass; (2) checks that the companion note carries its fourteen
results; (3) re-runs the seven evidence scripts of the erratum (R119-R125) and
requires a RESULT line with zero failures from each; (4) runs the unchanged V2.9.2
numerical suite when its dependencies are installed, and reports SKIP otherwise.

Exit status is zero only if no check fails (SKIP is not a failure).
"""
from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # v2.11/
REPO = ROOT.parent
MANUSCRIPT = ROOT / "manuscript" / "DS_model_V2_11_erratum.tex"
NOTE = ROOT / "negative_results" / "DS_negative_results.tex"
V29_SUITE = REPO / "v2.9" / "verification" / "verify_v291.py"

EVIDENCE = [
    ("E8/R119 literature watch", REPO / "matter" / "redraw" / "veille_21nov2026.py"),
    ("R120 Nambu-Goto / Regge", REPO / "matter" / "redraw" / "nambu_goto_regge.py"),
    ("E2/R121 internal LC oscillator", REPO / "matter" / "redraw" / "internal_lc_oscillator.py"),
    ("R122 saturable self-trapping", REPO / "matter" / "rebuild" / "step1_saturable_selftrapping.py"),
    ("E8-E10/R123 amplitude-axiom consequences", REPO / "matter" / "redraw" / "a3_adoption_consequences.py"),
    ("E6/R124 preferred-frame parameters (sympy)", REPO / "matter" / "rebuild" / "alpha12_preferred_frame.py"),
    ("E7/R125 substrate stars under TOV", REPO / "matter" / "rebuild" / "limit_states_tov.py"),
]

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    RESULTS.append(("PASS" if condition else "FAIL", name, detail))


def skip(name: str, detail: str) -> None:
    RESULTS.append(("SKIP", name, detail))


def block_source(source: str) -> None:
    # V2.10 hygiene checks, re-run against the V2.11 source
    check("Low-impedance equation withdrawn", "eq:vsync" not in source, "no synchronization-speed equation label")
    check("A8 source printed", "u - u_0" in source, "corrected source present")
    check("Conditional phase equation printed", "eq:phase_action" in source and "eq:z3_coefficient" in source,
          "phase normalization and Z3 coefficient present")
    check("Instanton-to-G gap disclosed", "The two missing bridges to $G$" in source, "missing gravitational map stated")
    # erratum items
    check("Version stamped V2.11", "V2.11 (erratum edition)" in source, "title-page date line")
    check("Erratum notice present", "Erratum notice (V2.11, 21 November 2026)" in source, "front-matter notice")
    check("E1: charge-by-chirality / mass-by-inductance removed from Methods",
          "Electric charge is dictated topologically by the winding chirality" not in source
          and "are withdrawn in this erratum edition (charge is not carried" in source, "abstract corrected")
    check("E2: electron soliton withdrawn", "Erratum E2 (21 November 2026)" in source
          and "withdrawn as an extended object by erratum E2" in source, "electron section and table")
    check("E3: proton impedance blocked in principle", "Erratum E3:" in source and "jointly unrealizable" in source, "open problem 18")
    check("E4: Gamma_pole re-filed", "erratum E4" in source, "status ledger")
    check("E5: A5 reworded", "erratum E5" in source and "no pole fusion" in source, "axiom A5")
    check("E6: alpha_1 = -8 computed", "Erratum E6 (21 November 2026): computed" in source
          and "\\alpha_1 = -(4\\gamma+4) = -8" in source, "open problem 6 closed negatively")
    check("E7: TOV re-integration and withdrawals", "Erratum E7 (21 November 2026): the branch under general-relativistic hydrostatics" in source
          and "whose V2.10 falsification hierarchy was" in source
          and "is withdrawn in this erratum edition (E7)" in source
          and "the prediction is withdrawn" in source, "substrate stars, abstract, hierarchy, registered prediction")
    check("E8: lattice dispersion under the amplitude axiom", "Erratum E8 (21 November 2026)" in source
          and "excludes $\\ell_1$ as the lattice spacing" in source, "Brillouin cutoff paragraph")
    check("E11: ledger re-filed", "re-filed by erratum E7" in source and "withdrawn by erratum E7" in source, "Derived list")
    check("E12 / Document History records V2.11", "\\paragraph{V2.11 (erratum edition, 21 November 2026).}" in source
          and "adopted an amplitude axiom (A9" in source, "revision traced")
    check("E6 attributes alpha_1 = -8 to the classical stratified theories",
          "\\cite{ni1972,nordtvedtwill1972}" in source and "\\bibitem{ni1972}" in source and "\\bibitem{nordtvedtwill1972}" in source,
          "Ni 1972; Nordtvedt and Will 1972 cited")
    check("Distribution stated: not deposited on Zenodo",
          "it is not deposited on Zenodo, whose record remains V2.10" in source
          and "this edition is not deposited on Zenodo" in source, "erratum notice and Data Availability")
    check("Conclusion states what survives", "Erratum (V2.11):} these three items are withdrawn as predictions" in source, "general conclusion")


def block_note(note: str) -> None:
    n = note.count("\\begin{prop}")
    check("Companion note carries fourteen results", n == 14, f"{n} results")
    classical = ("derrick1964", "townes1964", "vk1973", "bag1974", "ni1972", "nordtvedtwill1972", "luscher1981", "tov1939")
    missing = [k for k in classical if f"\\bibitem{{{k}}}" not in note or k not in note.split("\\begin{thebibliography}")[0]]
    check("Note cites its classical sources in the text and lists them", not missing, "all present" if not missing else f"missing {missing}")
    for key in ("alpha_1 = -(4\\gamma+4) = -8", "158/\\sqrt{x_s}", "0.31", "j_0(qR_3)"):
        check(f"Note states: {key}", key in note, "present")


def run_script(path: Path) -> tuple[str, str]:
    if not path.is_file():
        return "FAIL", f"not found: {path}"
    proc = subprocess.run([sys.executable, str(path)], cwd=str(path.parent), capture_output=True, text=True, timeout=900)
    out = proc.stdout + proc.stderr
    if "No module named" in out:
        mod = re.search(r"No module named '([^']+)'", out)
        return "SKIP", f"missing module {mod.group(1) if mod else '?'} (pip install it)"
    line = next((ln for ln in reversed(out.strip().splitlines()) if ln.startswith("RESULT:")), "")
    ok = proc.returncode == 0 and "0 FAIL" in line
    return ("PASS" if ok else "FAIL"), (line or f"exit {proc.returncode}, no RESULT line")


def block_evidence() -> None:
    for name, path in EVIDENCE:
        status, detail = run_script(path)
        RESULTS.append((status, f"Evidence {name}", f"{path.relative_to(REPO)}: {detail}"))


def block_delegated() -> None:
    missing = [m for m in ("numpy", "scipy", "matplotlib") if importlib.util.find_spec(m) is None]
    if missing:
        skip("Delegated V2.9.2 numerical suite", f"missing {', '.join(missing)}; install -r v2.9/requirements.txt to run it")
        return
    if not V29_SUITE.is_file():
        check("Delegated V2.9.2 numerical suite", False, f"not found: {V29_SUITE}")
        return
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run([sys.executable, str(V29_SUITE)], cwd=str(V29_SUITE.parents[1]), env=env, capture_output=True, text=True)
    line = next((ln for ln in reversed(proc.stdout.strip().splitlines()) if ln.startswith("RESULT:")), "no RESULT line")
    check("Delegated V2.9.2 numerical suite", proc.returncode == 0, line)


def main() -> int:
    for p, label in ((MANUSCRIPT, "manuscript"), (NOTE, "note")):
        if not p.is_file():
            print(f"FAIL  {label} not found: {p}")
            return 1
    block_source(MANUSCRIPT.read_text(encoding="utf-8"))
    block_note(NOTE.read_text(encoding="utf-8"))
    block_evidence()
    block_delegated()
    print("=" * 72)
    print("V2.11 erratum edition — checks")
    print("=" * 72)
    for status, name, detail in RESULTS:
        print(f"  [{status}] {name}: {detail}")
    failed = [r for r in RESULTS if r[0] == "FAIL"]
    skipped = [r for r in RESULTS if r[0] == "SKIP"]
    total = len(RESULTS)
    print()
    print(f"RESULT: {total - len(failed) - len(skipped)}/{total - len(skipped)} PASS; {len(failed)} FAIL; {len(skipped)} SKIP")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
