#!/usr/bin/env python3
"""Leave-one-SOX-family sensitivity for cross-layer module completeness.

This is a fast post-boundary robustness test.  It does not alter or rescue the
frozen symmetric route model; it asks whether the localized HGT-to-module
association depends on one reconciled gene family.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def within_order_test(frame: pd.DataFrame, permutations: int, rng: np.random.Generator) -> dict:
    cols = ["recipient_mean", "module_completeness", "gtdb_order"]
    x = frame[cols].dropna().copy()
    observed = float(x["recipient_mean"].corr(x["module_completeness"], method="spearman"))
    values = x["recipient_mean"].to_numpy(float)
    groups = list(x.groupby("gtdb_order", sort=False).indices.values())
    null = np.empty(permutations)
    for i in range(permutations):
        shuffled = values.copy()
        for idx in groups:
            if len(idx) > 1:
                shuffled[idx] = rng.permutation(shuffled[idx])
        null[i] = pd.Series(shuffled).corr(
            x["module_completeness"].reset_index(drop=True), method="spearman"
        )
    return {
        "genomes": len(x),
        "orders": x["gtdb_order"].nunique(),
        "observed_spearman_rho": observed,
        "permutation_p_two_sided": float(
            (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
        ),
        "null_q025": float(np.nanquantile(null, 0.025)),
        "null_q975": float(np.nanquantile(null, 0.975)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--branches", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tree-model", default="fixed")
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260711)
    args = ap.parse_args()

    branches = pd.read_csv(args.branches, sep="\t")
    tips = branches.loc[branches["is_terminal_genome"].astype(bool)].copy()
    tips["family"] = tips["family"].str.lower()
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.drop_duplicates("assembly_accession")
    metadata["module_completeness"] = pd.to_numeric(
        metadata["module_completeness"], errors="coerce"
    )
    metadata = metadata[["assembly_accession", "gtdb_order", "module_completeness"]]
    families = sorted(tips["family"].dropna().unique())

    rows = []
    rng = np.random.default_rng(args.seed)
    omission_sets = [("none", set())] + [(family, {family}) for family in families]
    for label, omitted in omission_sets:
        subset = tips.loc[~tips["family"].isin(omitted)].copy()
        aggregate = (
            subset.groupby("node", as_index=False)
            .agg(
                recipient_mean=("transfer_recipient_sample_fraction", "mean"),
                families_retained=("family", "nunique"),
            )
            .rename(columns={"node": "assembly_accession"})
        )
        integrated = metadata.merge(aggregate, on="assembly_accession", how="inner")
        stats = within_order_test(integrated, args.permutations, rng)
        rows.append(
            {
                "tree_model": args.tree_model,
                "omitted_family": label,
                "families_available": len(families) - len(omitted),
                **stats,
            }
        )

    # Family-specific estimates diagnose heterogeneity without treating them as independent tests.
    for family in families:
        subset = tips.loc[tips["family"].eq(family), [
            "node", "transfer_recipient_sample_fraction"
        ]].rename(
            columns={
                "node": "assembly_accession",
                "transfer_recipient_sample_fraction": "recipient_mean",
            }
        )
        integrated = metadata.merge(subset, on="assembly_accession", how="inner")
        stats = within_order_test(integrated, args.permutations, rng)
        rows.append(
            {
                "tree_model": args.tree_model,
                "omitted_family": f"family_only:{family}",
                "families_available": 1,
                **stats,
            }
        )

    result = pd.DataFrame(rows)
    result["holm_p_leave_one_family"] = np.nan
    mask = result["omitted_family"].isin(families)
    p = result.loc[mask, "permutation_p_two_sided"].astype(float)
    ordered = p.sort_values()
    adjusted = []
    running = 0.0
    m = len(ordered)
    for rank, (idx, value) in enumerate(ordered.items()):
        running = max(running, (m - rank) * value)
        adjusted.append((idx, min(1.0, running)))
    for idx, value in adjusted:
        result.loc[idx, "holm_p_leave_one_family"] = value

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    print(args.out)


if __name__ == "__main__":
    main()
