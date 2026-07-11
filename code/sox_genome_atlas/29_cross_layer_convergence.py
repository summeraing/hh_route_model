from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def within_group_permutation(
    frame: pd.DataFrame,
    x: str,
    y: str,
    group: str,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    subset = frame[[x, y, group]].dropna().copy()
    observed = float(subset[x].corr(subset[y], method="spearman"))
    null = np.empty(permutations)
    x_values = subset[x].to_numpy(float)
    group_indices = list(subset.groupby(group, sort=False).indices.values())
    for iteration in range(permutations):
        shuffled = x_values.copy()
        for indices in group_indices:
            if len(indices) > 1:
                shuffled[indices] = rng.permutation(shuffled[indices])
        null[iteration] = pd.Series(shuffled).corr(subset[y].reset_index(drop=True), method="spearman")
    bootstrap = np.empty(permutations)
    for iteration in range(permutations):
        sampled = subset.iloc[rng.integers(0, len(subset), size=len(subset))]
        bootstrap[iteration] = sampled[x].corr(sampled[y], method="spearman")
    return {
        "genomes": len(subset),
        "observed_spearman_rho": observed,
        "permutation_p_two_sided": float((1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)),
        "null_q025": float(np.nanquantile(null, 0.025)),
        "null_q975": float(np.nanquantile(null, 0.975)),
        "bootstrap_q025": float(np.nanquantile(bootstrap, 0.025)),
        "bootstrap_q975": float(np.nanquantile(bootstrap, 0.975)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generax-branches", type=Path, required=True)
    parser.add_argument("--operon-catalog", type=Path, required=True)
    parser.add_argument("--composition", type=Path, required=True)
    parser.add_argument("--matched-mobile", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    branches = pd.read_csv(args.generax_branches, sep="\t")
    tips = branches.loc[branches["is_terminal_genome"].astype(bool)].copy()
    hgt = (
        tips.groupby("node", as_index=False)
        .agg(
            hgt_recipient_fraction_mean=("transfer_recipient_sample_fraction", "mean"),
            hgt_recipient_fraction_max=("transfer_recipient_sample_fraction", "max"),
            mean_transfer_events_per_sample=("mean_T_events_per_sample", "mean"),
            sox_families_in_generax=("family", "nunique"),
        )
        .rename(columns={"node": "assembly_accession"})
    )

    operons = pd.read_csv(args.operon_catalog, sep="\t", dtype=str).fillna("")
    operons["strand_coherence"] = pd.to_numeric(
        operons["strand_coherence"], errors="coerce"
    )
    signature_frequency = operons["architecture_signature"].value_counts()
    operons["architecture_frequency"] = operons["architecture_signature"].map(signature_frequency)
    operons["architecture_rarity"] = -np.log10(operons["architecture_frequency"] / len(operons))
    architecture = (
        operons.groupby("assembly_accession", as_index=False)
        .agg(
            architecture_rarity_mean=("architecture_rarity", "mean"),
            architecture_rarity_max=("architecture_rarity", "max"),
            operon_architectures=("architecture_signature", "nunique"),
            sox_clusters=("cluster_id", "nunique"),
            mean_strand_coherence=("strand_coherence", "mean"),
        )
    )

    composition = pd.read_csv(args.composition, sep="\t", low_memory=False)
    composition = composition.loc[composition["organizational_role"].eq("sox_energy_module")].copy()
    for column in ["gene_gc_deviation", "local_gc_deviation"]:
        composition[column] = pd.to_numeric(composition[column], errors="coerce")
        composition[f"absolute_{column}"] = composition[column].abs()
    composition_summary = (
        composition.groupby("assembly_accession", as_index=False)
        .agg(
            mean_abs_sox_gene_gc_deviation=("absolute_gene_gc_deviation", "mean"),
            max_abs_sox_gene_gc_deviation=("absolute_gene_gc_deviation", "max"),
            mean_abs_sox_local_gc_deviation=("absolute_local_gc_deviation", "mean"),
            max_abs_sox_local_gc_deviation=("absolute_local_gc_deviation", "max"),
        )
    )

    mobile = pd.read_csv(args.matched_mobile, sep="\t")
    mobile = mobile.loc[mobile["window_bp"].eq(10000)].copy()
    mobile_summary = (
        mobile.groupby("assembly_accession", as_index=False)
        .agg(
            matched_mobile_count_difference=("count_difference", "mean"),
            matched_mobile_presence_difference=("presence_difference", "mean"),
        )
    )

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.drop_duplicates("assembly_accession")
    metadata["module_completeness"] = pd.to_numeric(metadata["module_completeness"], errors="coerce")
    keep = [
        "assembly_accession",
        "gtdb_class",
        "gtdb_order",
        "module_class",
        "module_completeness",
    ]
    integrated = metadata[keep].merge(hgt, on="assembly_accession", how="left")
    integrated = integrated.merge(architecture, on="assembly_accession", how="left")
    integrated = integrated.merge(composition_summary, on="assembly_accession", how="left")
    integrated = integrated.merge(mobile_summary, on="assembly_accession", how="left")

    outcomes = [
        "architecture_rarity_mean",
        "architecture_rarity_max",
        "mean_abs_sox_gene_gc_deviation",
        "mean_abs_sox_local_gc_deviation",
        "module_completeness",
        "matched_mobile_count_difference",
        "matched_mobile_presence_difference",
    ]
    predictor = "hgt_recipient_fraction_mean"
    test_rows: list[dict[str, object]] = []
    for outcome in outcomes:
        stats = within_group_permutation(
            integrated,
            predictor,
            outcome,
            "gtdb_order",
            args.permutations,
            rng,
        )
        test_rows.append({"predictor": predictor, "outcome": outcome, **stats})
    tests = pd.DataFrame(test_rows)

    class_summary = (
        integrated.groupby("gtdb_class", as_index=False)
        .agg(
            genomes=("assembly_accession", "nunique"),
            mean_hgt_recipient_fraction=("hgt_recipient_fraction_mean", "mean"),
            mean_architecture_rarity=("architecture_rarity_mean", "mean"),
            mean_abs_local_gc_deviation=("mean_abs_sox_local_gc_deviation", "mean"),
            mean_mobile_difference=("matched_mobile_count_difference", "mean"),
        )
    )

    integrated.to_csv(args.outdir / "sox_cross_layer_genome_table.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "sox_cross_layer_association_tests.tsv", sep="\t", index=False)
    class_summary.to_csv(args.outdir / "sox_cross_layer_class_summary.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "status": "post-boundary cross-layer sensitivity",
                "predictor": predictor,
                "taxonomic_control": "shuffle predictor within GTDB order",
                "mobile_window_bp": 10000,
                "permutations": args.permutations,
                "seed": args.seed,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"genomes={len(integrated)} tests={len(tests)}")


if __name__ == "__main__":
    main()
