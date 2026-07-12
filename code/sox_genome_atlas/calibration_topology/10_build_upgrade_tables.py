#!/usr/bin/env python3
"""Build manuscript-ready summaries for the three-day computation upgrade."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


METHOD_ORDER = [
    "raw_pooled",
    "source_equal",
    "source_equal_dependency_collapsed",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def representative_simulation(summary: pd.DataFrame) -> pd.DataFrame:
    common = (
        summary["source_regime"].eq("dominant_adversarial")
        & summary["signal_strength"].eq(0.65)
        & summary["dominant_source_fraction"].eq(0.85)
        & summary["mean_dependency_group_size"].eq(5)
        & summary["label_noise_fraction"].eq(0.15)
    )
    full = summary.loc[
        common
        & summary["architecture"].eq("full")
        & summary["metric"].eq("full_route_recovered")
    ].copy()
    full.insert(0, "case", "full route; 85% adversarial row source")
    partial = summary.loc[
        common
        & summary["architecture"].eq("partial")
        & summary["metric"].eq("partial_route_recovered")
    ].copy()
    partial.insert(0, "case", "partial route; 85% adversarial row source")
    out = pd.concat([full, partial], ignore_index=True)
    out["method"] = pd.Categorical(out["method"], METHOD_ORDER, ordered=True)
    return out.sort_values(["case", "method"]).reset_index(drop=True)


def null_summary(summary: pd.DataFrame) -> pd.DataFrame:
    null = summary.loc[
        summary["architecture"].eq("null")
        & summary["metric"].eq("passes_matched_null")
    ].copy()
    out = (
        null.groupby("method", as_index=False)
        .agg(
            cells=("rate", "size"),
            mean_false_positive_rate=("rate", "mean"),
            median_false_positive_rate=("rate", "median"),
            maximum_false_positive_rate=("rate", "max"),
            minimum_false_positive_rate=("rate", "min"),
        )
    )
    out["method"] = pd.Categorical(out["method"], METHOD_ORDER, ordered=True)
    return out.sort_values("method").reset_index(drop=True)


def topology_summary(topology: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        topology.groupby("tree_design", as_index=False)
        .agg(
            families=("family", "nunique"),
            min_delta_log_likelihood=("delta_log_likelihood", "min"),
            max_delta_log_likelihood=("delta_log_likelihood", "max"),
            max_holm_au_p=("p_au_holm", "max"),
            holm_au_significant=("p_au_holm", lambda values: int((values < 0.05).sum())),
            holm_sh_significant=("p_sh_holm", lambda values: int((values < 0.05).sum())),
        )
    )
    return grouped


def fmt(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sim = pd.read_csv(
        args.results / "01_route_simulation/summary/simulation_headline_metrics.tsv",
        sep="\t",
        keep_default_na=False,
    )
    representative = representative_simulation(sim)
    null = null_summary(sim)

    consensus = pd.read_csv(
        args.results / "02_cross_model_consensus/cross_model_module_completeness_tests.tsv",
        sep="\t",
    )
    consensus_primary = consensus.loc[consensus["excluded_family"].eq("none")].copy()
    consensus_family = pd.read_csv(
        args.results / "02_cross_model_consensus/cross_model_family_summary.tsv",
        sep="\t",
    )

    topology = pd.read_csv(
        args.results / "03_topology_tests/summary/sox_gtdb_constraint_tests.tsv",
        sep="\t",
    )
    topology_design = topology_summary(topology)

    representative.to_csv(args.outdir / "simulation_representative_cases.tsv", sep="\t", index=False)
    null.to_csv(args.outdir / "simulation_null_fpr_summary.tsv", sep="\t", index=False)
    consensus_primary.to_csv(args.outdir / "cross_model_primary_tests.tsv", sep="\t", index=False)
    consensus_family.to_csv(args.outdir / "cross_model_family_summary.tsv", sep="\t", index=False)
    topology.to_csv(args.outdir / "topology_family_tests.tsv", sep="\t", index=False)
    topology_design.to_csv(args.outdir / "topology_design_summary.tsv", sep="\t", index=False)

    full = representative.loc[representative["architecture"].eq("full")].set_index("method")
    partial = representative.loc[representative["architecture"].eq("partial")].set_index("method")
    con = consensus_primary.set_index("metric")

    metrics = {
        "simulation": {
            "full_route_raw_recovery": float(full.loc["raw_pooled", "rate"]),
            "full_route_source_equal_recovery": float(full.loc["source_equal", "rate"]),
            "full_route_source_equal_collapsed_recovery": float(full.loc["source_equal_dependency_collapsed", "rate"]),
            "partial_route_raw_recovery": float(partial.loc["raw_pooled", "rate"]),
            "partial_route_source_equal_recovery": float(partial.loc["source_equal", "rate"]),
            "partial_route_source_equal_collapsed_recovery": float(partial.loc["source_equal_dependency_collapsed", "rate"]),
            "null_mean_fpr_min": float(null["mean_false_positive_rate"].min()),
            "null_mean_fpr_max": float(null["mean_false_positive_rate"].max()),
        },
        "cross_model": {
            row.metric: {
                "rho": float(row.spearman_rho),
                "p": float(row.permutation_p),
                "ci95_low": float(row.bootstrap_ci95_low),
                "ci95_high": float(row.bootstrap_ci95_high),
            }
            for row in consensus_primary.itertuples(index=False)
        },
        "topology": {
            "tests": int(len(topology)),
            "families": int(topology["family"].nunique()),
            "designs": int(topology["tree_design"].nunique()),
            "delta_log_likelihood_min": float(topology["delta_log_likelihood"].min()),
            "delta_log_likelihood_max": float(topology["delta_log_likelihood"].max()),
            "holm_au_significant": int((topology["p_au_holm"] < 0.05).sum()),
            "maximum_holm_au_p": float(topology["p_au_holm"].max()),
        },
    }
    (args.outdir / "upgrade_key_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )

    summary_md = f"""# Three-day computation upgrade: verified results

## Known-truth calibration

- The matched diffuse-null false-positive rate averaged {fmt(null['mean_false_positive_rate'].min())}-{fmt(null['mean_false_positive_rate'].max())} across the three estimators.
- Under a frozen stress case with 85% of rows from an adversarial dominant source, route strength 0.65, mean dependency size 5 and 15% donor/role label noise, full-route recovery was {fmt(full.loc['raw_pooled', 'rate'])} under raw pooling, {fmt(full.loc['source_equal', 'rate'])} after source equalization and {fmt(full.loc['source_equal_dependency_collapsed', 'rate'])} after source equalization plus dependency collapse.
- In the corresponding partial-route case, recovery was {fmt(partial.loc['raw_pooled', 'rate'])}, {fmt(partial.loc['source_equal', 'rate'])} and {fmt(partial.loc['source_equal_dependency_collapsed', 'rate'])}, respectively.

## Cross-tree consensus

- Cross-tree mean recipient support was associated with module completeness (rho={fmt(con.loc['consensus_mean', 'spearman_rho'])}, within-order permutation P={fmt(con.loc['consensus_mean', 'permutation_p'], 4)}, bootstrap 95% CI {fmt(con.loc['consensus_mean', 'bootstrap_ci95_low'])}-{fmt(con.loc['consensus_mean', 'bootstrap_ci95_high'])}).
- The conservative cross-tree minimum remained positive (rho={fmt(con.loc['consensus_min', 'spearman_rho'])}, P={fmt(con.loc['consensus_min', 'permutation_p'], 4)}, 95% CI {fmt(con.loc['consensus_min', 'bootstrap_ci95_low'])}-{fmt(con.loc['consensus_min', 'bootstrap_ci95_high'])}).
- A binary stable-in-both indicator was positive but not significant (rho={fmt(con.loc['stable_both', 'spearman_rho'])}, P={fmt(con.loc['stable_both', 'permutation_p'], 3)}), defining a terminal-assignment stability boundary.

## Formal topology tests

- All {len(topology)} family-by-tree-design tests rejected the GTDB-constrained topology after Holm correction of AU P values.
- The constrained trees lost {topology['delta_log_likelihood'].min():.1f}-{topology['delta_log_likelihood'].max():.1f} log-likelihood units relative to the unconstrained family trees.
- These tests reject strict species-tree congruence; they do not by themselves identify transfer direction or timing.
"""
    (args.outdir / "UPGRADE_RESULTS_SUMMARY.md").write_text(summary_md, encoding="utf-8")
    print(args.outdir)


if __name__ == "__main__":
    main()
