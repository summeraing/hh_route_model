from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "smoke_tests"


def run(cmd: list[str]) -> None:
    print("RUN", " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def assert_close(value: float, expected: float, tol: float = 1e-6) -> None:
    if not math.isclose(value, expected, rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"Expected {expected}, observed {value}")


def check_sox_route_graph() -> None:
    out = OUT / "sox_route_graph"
    run([
        sys.executable,
        "code/route_graph/run_route_graph_case.py",
        "--input", "data/sox_transfer/sox_evidence_units_expanded_v2.csv",
        "--case-id", "sox_hgt_transfer",
        "--donors", "lineage_core,hgt_pathway,mobile_context",
        "--roles", "conserved_metabolic_backbone,variable_energy_module,mobility_boundary",
        "--prespec", "lineage_core=conserved_metabolic_backbone,hgt_pathway=variable_energy_module,mobile_context=mobility_boundary",
        "--out", str(out),
        "--iterations", "200",
        "--enforce-gates",
    ])
    summary = pd.read_csv(out / "summary.csv").iloc[0]
    if int(summary["n_route_eligible_units"]) != 65:
        raise AssertionError(summary.to_dict())
    if int(summary["prespecified_rank"]) != 1:
        raise AssertionError(summary.to_dict())
    assert_close(float(summary["prespecified_source_equal"]), 1.0)
    assert_close(float(summary["margin_vs_best_alternative"]), 0.3087121212121212, tol=1e-9)


def check_permutation_equivalence() -> None:
    run([sys.executable, "code/route_graph/test_permutation_equivalence.py"])


def check_sox_af3_gates() -> None:
    out_md = OUT / "sox_af3_gate_check.md"
    out_csv = OUT / "sox_af3_gate_check.csv"
    run([
        sys.executable,
        "code/af3_postprocess/evaluate_sox_af3_gates.py",
        "--summary-csv", "data/sox_af3/sox_af3_combined_primary_boundary.csv",
        "--output-md", str(out_md),
        "--output-csv", str(out_csv),
    ])
    text = out_md.read_text(encoding="utf-8")
    if "STRUCTURAL_GATE_PASS" not in text:
        raise AssertionError(text)


def check_sox_genome_atlas_release() -> None:
    run([
        sys.executable,
        "code/sox_genome_atlas/validate_released_results.py",
        "--data-root", "data/sox_genome_atlas",
    ])


def check_eukaryogenesis_core_metrics() -> None:
    run([sys.executable, "code/core/reproduce_core_metrics.py"])
    run([sys.executable, "code/core/reproduce_continuous_route_diagnostic.py"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    check_eukaryogenesis_core_metrics()
    check_sox_route_graph()
    check_permutation_equivalence()
    check_sox_af3_gates()
    check_sox_genome_atlas_release()
    print("SMOKE_TESTS_PASS")


if __name__ == "__main__":
    main()
