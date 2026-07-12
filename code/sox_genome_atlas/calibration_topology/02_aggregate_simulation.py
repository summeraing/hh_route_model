#!/usr/bin/env python3
"""Aggregate known-truth simulation outputs and calibrate matched null gates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def binomial_interval(successes: int, total: int, alpha: float = 0.05) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    lower = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, total - successes + 1))
    upper = 1.0 if successes == total else float(beta.ppf(1 - alpha / 2, successes + 1, total - successes))
    return lower, upper


def main() -> None:
    args = parse_args()
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    files = sorted(Path(args.input_dir).glob("scenario_*.csv"))
    if not files:
        raise SystemExit("No scenario CSV files found")
    # Preserve the literal architecture label "null"; pandas otherwise treats
    # it as a missing-value token and silently removes the calibration stratum.
    frame = pd.concat(
        [pd.read_csv(path, keep_default_na=False) for path in files],
        ignore_index=True,
    )

    structural = [
        "method",
        "dominant_source_fraction",
        "mean_dependency_group_size",
        "label_noise_fraction",
    ]
    null = frame.loc[frame["architecture"] == "null"].copy()
    calibration_fraction = float(spec["simulation"]["null_calibration_fraction"])
    cutoff = int(spec["simulation"]["replicates_per_scenario"] * calibration_fraction)
    null_cal = null.loc[null["replicate"] < cutoff]
    null_eval = null.loc[null["replicate"] >= cutoff]
    quantile = float(spec["simulation"]["null_quantile"])
    gates = (
        null_cal.groupby(structural, as_index=False)["margin"]
        .quantile(quantile)
        .rename(columns={"margin": "matched_null_margin_q975"})
    )

    evaluated = frame.merge(gates, on=structural, how="left", validate="many_to_one")
    evaluated["rank1"] = evaluated["prespecified_rank"].eq(1)
    evaluated["passes_matched_null"] = evaluated["margin"] > evaluated["matched_null_margin_q975"]
    evaluated["full_route_recovered"] = (
        evaluated["rank1"]
        & evaluated["passes_matched_null"]
        & evaluated["full_components_correct"].eq(1)
    )
    evaluated["partial_route_recovered"] = (
        evaluated["rank1"]
        & evaluated["passes_matched_null"]
        & evaluated["partial_components_correct"].eq(1)
    )

    eval_frame = pd.concat([
        evaluated.loc[evaluated["architecture"] != "null"],
        evaluated.loc[(evaluated["architecture"] == "null") & (evaluated["replicate"] >= cutoff)],
    ], ignore_index=True)
    group_cols = [
        "architecture",
        "signal_strength",
        "source_regime",
        "method",
        "dominant_source_fraction",
        "mean_dependency_group_size",
        "label_noise_fraction",
    ]
    summary_rows = []
    for key, chunk in eval_frame.groupby(group_cols, dropna=False, sort=True):
        base = dict(zip(group_cols, key))
        for metric in (
            "rank1",
            "passes_matched_null",
            "full_route_recovered",
            "partial_route_recovered",
        ):
            successes = int(chunk[metric].sum())
            total = int(len(chunk))
            lo, hi = binomial_interval(successes, total)
            summary_rows.append({
                **base,
                "metric": metric,
                "successes": successes,
                "replicates": total,
                "rate": successes / total if total else float("nan"),
                "ci95_low": lo,
                "ci95_high": hi,
                "mean_margin": float(chunk["margin"].mean()),
                "median_margin": float(chunk["margin"].median()),
            })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summary_rows)
    gates.to_csv(outdir / "matched_null_margin_gates.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "simulation_recovery_summary.tsv", sep="\t", index=False)
    evaluated.to_csv(outdir / "simulation_replicate_results.tsv.gz", sep="\t", index=False, compression="gzip")

    headline = summary.loc[
        ((summary["architecture"] == "null") & (summary["metric"] == "passes_matched_null"))
        | ((summary["architecture"] == "full") & (summary["metric"] == "full_route_recovered"))
        | ((summary["architecture"] == "partial") & (summary["metric"] == "partial_route_recovered"))
    ].copy()
    headline.to_csv(outdir / "simulation_headline_metrics.tsv", sep="\t", index=False)
    metadata = {
        "spec_id": spec["spec_id"],
        "scenario_files": len(files),
        "replicate_rows": len(evaluated),
        "null_calibration_replicates_per_cell": cutoff,
        "null_evaluation_replicates_per_cell": int(spec["simulation"]["replicates_per_scenario"]) - cutoff,
    }
    (outdir / "analysis_parameters.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
