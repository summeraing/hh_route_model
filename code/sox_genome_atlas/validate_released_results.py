from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.data_root

    qc = read_tsv(root / "04_candidate_qc" / "candidate_qc_summary.tsv")
    metrics = dict(zip(qc["metric"], qc["value"], strict=True))
    require(int(metrics["candidate_genomes"]) == 427, "candidate genome count changed")
    require(int(metrics["gtdb_mapped"]) == 398, "GTDB-mapped genome count changed")

    route = pd.read_csv(root / "15_atlas_route_fixed" / "summary.csv").iloc[0]
    require(int(route["prespecified_rank"]) == 5, "frozen symmetric route rank changed")
    require(bool(route["source_null_non_degenerate"]), "source null became degenerate")
    require(bool(route["source_layer_null_non_degenerate"]), "source-layer null became degenerate")

    fixed = read_tsv(root / "39_generax_calibrated_summary" / "generax_family_summary.tsv")
    selected = read_tsv(root / "44_generax_modelselected_summary" / "generax_family_summary.tsv")
    require(set(fixed["family"]) == set(selected["family"]), "tree-model family sets differ")
    for name, table in (("fixed", fixed), ("ModelFinder", selected)):
        rates = table["transfer_events_per_100_tips"]
        require(len(rates) == 6, f"{name} family count changed")
        require(rates.between(30, 50).all(), f"{name} transfer-rate range changed")

    mobile = read_tsv(root / "27_matched_mobile_boundary" / "matched_mobile_boundary_tests.tsv")
    require((mobile["observed"] < 0).all(), "matched mobile-control direction changed")
    composition = read_tsv(root / "30_matched_composition_shift" / "matched_composition_tests.tsv")
    require((composition["observed"] < 0).all(), "matched composition-control direction changed")

    cross = read_tsv(root / "46_cross_layer_modelselected" / "sox_cross_layer_association_tests.tsv")
    completeness = cross.loc[cross["outcome"].eq("module_completeness")].iloc[0]
    require(float(completeness["observed_spearman_rho"]) > 0, "module-completeness direction changed")

    print("SOX_GENOME_ATLAS_RELEASE_VALIDATION_PASS")


if __name__ == "__main__":
    main()
