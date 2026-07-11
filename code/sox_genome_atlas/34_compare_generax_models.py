from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", type=Path, required=True)
    parser.add_argument("--model-selected", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fixed_family = pd.read_csv(args.fixed / "generax_family_summary.tsv", sep="\t")
    selected_family = pd.read_csv(args.model_selected / "generax_family_summary.tsv", sep="\t")
    family = fixed_family.merge(
        selected_family,
        on="family",
        suffixes=("_fixed_lg", "_model_selected"),
        validate="one_to_one",
    )
    family["transfer_event_difference"] = family["T_model_selected"] - family["T_fixed_lg"]
    family["transfer_event_relative_difference"] = (
        family["transfer_event_difference"] / family["T_fixed_lg"].replace(0, pd.NA)
    )

    fixed_branch = pd.read_csv(args.fixed / "generax_per_branch_support.tsv", sep="\t")
    selected_branch = pd.read_csv(args.model_selected / "generax_per_branch_support.tsv", sep="\t")
    columns = ["family", "node", "is_terminal_genome", "transfer_recipient_sample_fraction"]
    branches = fixed_branch[columns].merge(
        selected_branch[columns],
        on=["family", "node", "is_terminal_genome"],
        suffixes=("_fixed_lg", "_model_selected"),
        how="inner",
        validate="one_to_one",
    )
    # ``is_terminal_genome`` is part of the merge key and therefore remains
    # unsuffixed in pandas.
    terminal = branches.loc[branches["is_terminal_genome"].astype(bool)].copy()
    comparison_rows: list[dict[str, object]] = []
    for family_name, group in terminal.groupby("family"):
        fixed_positive = set(
            group.loc[
                group["transfer_recipient_sample_fraction_fixed_lg"] >= 0.5, "node"
            ]
        )
        selected_positive = set(
            group.loc[
                group["transfer_recipient_sample_fraction_model_selected"] >= 0.5, "node"
            ]
        )
        union = fixed_positive | selected_positive
        comparison_rows.append(
            {
                "family": family_name,
                "shared_terminal_tips": len(group),
                "recipient_fraction_spearman": group[
                    "transfer_recipient_sample_fraction_fixed_lg"
                ].corr(
                    group["transfer_recipient_sample_fraction_model_selected"],
                    method="spearman",
                ),
                "stable_recipient_tips_fixed_lg": len(fixed_positive),
                "stable_recipient_tips_model_selected": len(selected_positive),
                "stable_recipient_jaccard": len(fixed_positive & selected_positive) / len(union)
                if union
                else 1.0,
            }
        )
    comparison = pd.DataFrame(comparison_rows)

    family.to_csv(args.outdir / "generax_family_model_comparison.tsv", sep="\t", index=False)
    terminal.to_csv(args.outdir / "generax_terminal_model_comparison.tsv", sep="\t", index=False)
    comparison.to_csv(args.outdir / "generax_recipient_stability_comparison.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "fixed_model": "LG+F+R6",
                "alternative": "ModelFinder-selected model per Sox family",
                "stable_recipient_threshold": 0.5,
                "interpretation": "Post-boundary model sensitivity only.",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"families={len(family)} terminal_rows={len(terminal)}")


if __name__ == "__main__":
    main()
