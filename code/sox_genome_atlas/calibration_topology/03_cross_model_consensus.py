#!/usr/bin/env python3
"""Conservative cross-tree consensus analysis for SOX transfer recipients."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, rankdata, spearmanr


FAMILIES = ("soxa", "soxb", "soxc", "soxx", "soxy", "soxz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--fixed-support", required=True)
    parser.add_argument("--mfp-support", required=True)
    parser.add_argument("--genome-table", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def holm_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * p[idx])
        running = max(running, value)
        adjusted[idx] = running
    return pd.Series(adjusted, index=values.index)


def within_order_permutation(
    frame: pd.DataFrame,
    outcome: str,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    x_rank = rankdata(frame["module_completeness"].to_numpy(dtype=float))
    y_rank = rankdata(frame[outcome].to_numpy(dtype=float))
    x_centered = x_rank - x_rank.mean()
    y_centered = y_rank - y_rank.mean()
    denominator = float(np.sqrt(np.dot(x_centered, x_centered) * np.dot(y_centered, y_centered)))
    observed = float(np.dot(x_centered, y_centered) / denominator)
    null = np.empty(permutations, dtype=float)
    groups = [idx.to_numpy() for _, idx in frame.groupby("gtdb_order").groups.items()]
    for i in range(permutations):
        shuffled = x_rank.copy()
        for idx in groups:
            shuffled[idx] = rng.permutation(shuffled[idx])
        null[i] = float(np.dot(shuffled - shuffled.mean(), y_centered) / denominator)
    p = (1.0 + float(np.sum(np.abs(null) >= abs(observed)))) / (permutations + 1.0)
    return observed, p


def within_order_bootstrap(
    frame: pd.DataFrame,
    outcome: str,
    bootstraps: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    chunks = [chunk for _, chunk in frame.groupby("gtdb_order", sort=False)]
    values = np.empty(bootstraps, dtype=float)
    for i in range(bootstraps):
        sampled = []
        for chunk in chunks:
            take = rng.integers(0, len(chunk), size=len(chunk))
            sampled.append(chunk.iloc[take])
        b = pd.concat(sampled, ignore_index=True)
        values[i] = float(spearmanr(b["module_completeness"], b[outcome]).statistic)
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def load_support(path: str, label: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame = frame.loc[frame["is_terminal_genome"].astype(str).str.lower().isin(["true", "1"])]
    frame = frame.loc[frame["family"].isin(FAMILIES), [
        "family", "node", "transfer_recipient_sample_fraction"
    ]].copy()
    return frame.rename(columns={
        "node": "assembly_accession",
        "transfer_recipient_sample_fraction": label,
    })


def genome_summary(long: pd.DataFrame, excluded_family: str | None = None) -> pd.DataFrame:
    frame = long if excluded_family is None else long.loc[long["family"] != excluded_family]
    return (
        frame.groupby(["assembly_accession", "gtdb_order", "module_completeness"], as_index=False)
        .agg(
            consensus_min=("consensus_min", "mean"),
            consensus_mean=("consensus_mean", "mean"),
            stable_both=("stable_both", "mean"),
            families=("family", "nunique"),
        )
    )


def main() -> None:
    args = parse_args()
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    cfg = spec["cross_model_consensus"]
    fixed = load_support(args.fixed_support, "support_fixed")
    mfp = load_support(args.mfp_support, "support_mfp")
    long = fixed.merge(mfp, on=["family", "assembly_accession"], how="inner", validate="one_to_one")
    metadata = pd.read_csv(args.genome_table, sep="\t")[[
        "assembly_accession", "gtdb_order", "module_completeness"
    ]].drop_duplicates("assembly_accession")
    long = long.merge(metadata, on="assembly_accession", how="inner", validate="many_to_one")
    long["consensus_min"] = long[["support_fixed", "support_mfp"]].min(axis=1)
    long["consensus_mean"] = long[["support_fixed", "support_mfp"]].mean(axis=1)
    long["stable_both"] = ((long["support_fixed"] >= 0.5) & (long["support_mfp"] >= 0.5)).astype(float)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    long.to_csv(outdir / "cross_model_family_genome_support.tsv", sep="\t", index=False)

    rng = np.random.default_rng(int(cfg["seed"]))
    metrics = ("consensus_min", "consensus_mean", "stable_both")
    result_rows = []
    for excluded in [None, *FAMILIES]:
        summary = genome_summary(long, excluded)
        if excluded is None:
            summary.to_csv(outdir / "cross_model_genome_consensus.tsv", sep="\t", index=False)
        for metric in metrics:
            rho, p = within_order_permutation(summary, metric, int(cfg["permutations"]), rng)
            if excluded is None:
                ci_lo, ci_hi = within_order_bootstrap(summary, metric, int(cfg["bootstraps"]), rng)
            else:
                ci_lo, ci_hi = float("nan"), float("nan")
            result_rows.append({
                "excluded_family": "none" if excluded is None else excluded,
                "metric": metric,
                "n_genomes": len(summary),
                "n_orders": summary["gtdb_order"].nunique(),
                "spearman_rho": rho,
                "permutation_p": p,
                "bootstrap_ci95_low": ci_lo,
                "bootstrap_ci95_high": ci_hi,
            })

    results = pd.DataFrame(result_rows)
    results["holm_p_within_metric"] = np.nan
    for metric in metrics:
        idx = (results["metric"] == metric) & (results["excluded_family"] != "none")
        results.loc[idx, "holm_p_within_metric"] = holm_adjust(results.loc[idx, "permutation_p"])
    results.to_csv(outdir / "cross_model_module_completeness_tests.tsv", sep="\t", index=False)

    family_summary = (
        long.groupby("family", as_index=False)
        .agg(
            n_genome_pairs=("assembly_accession", "nunique"),
            mean_fixed=("support_fixed", "mean"),
            mean_mfp=("support_mfp", "mean"),
            mean_consensus_min=("consensus_min", "mean"),
            stable_both_count=("stable_both", "sum"),
        )
    )
    family_summary.to_csv(outdir / "cross_model_family_summary.tsv", sep="\t", index=False)
    params = {
        "spec_id": spec["spec_id"],
        "permutations": int(cfg["permutations"]),
        "bootstraps": int(cfg["bootstraps"]),
        "family_genome_rows": len(long),
        "genomes": long["assembly_accession"].nunique(),
    }
    (outdir / "analysis_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    print(json.dumps(params, indent=2))


if __name__ == "__main__":
    main()
