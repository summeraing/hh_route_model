from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def read_stats(path: Path) -> tuple[float, float]:
    first_line = path.read_text().splitlines()[0]
    values = [float(value) for value in first_line.split()]
    if len(values) != 2:
        raise ValueError(f"Expected two likelihoods in {path}, found {values}")
    return values[0], values[1]


def read_rates(path: Path, expected: int) -> list[float]:
    values = [float(value) for value in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", path.read_text())]
    if len(values) < expected:
        raise ValueError(f"Could not parse reconciliation rates from {path}")
    return values[:expected]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dtl-root", type=Path, required=True)
    parser.add_argument("--dl-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for family in ["soxa", "soxb", "soxc", "soxx", "soxy", "soxz"]:
        dtl_family = args.dtl_root / family
        dl_family = args.dl_root / family
        dtl_joint, dtl_rec = read_stats(dtl_family / "results" / family / "stats.txt")
        dl_joint, dl_rec = read_stats(dl_family / "results" / family / "stats.txt")
        dtl_rates = read_rates(dtl_family / "gene_optimization_0" / "dtl_rates.txt", expected=3)
        dl_rates = read_rates(dl_family / "gene_optimization_0" / "dtl_rates.txt", expected=2)
        rows.append(
            {
                "family": family,
                "dtl_joint_likelihood": dtl_joint,
                "dl_joint_likelihood": dl_joint,
                "delta_joint_likelihood_dtl_minus_dl": dtl_joint - dl_joint,
                "dtl_reconciliation_likelihood": dtl_rec,
                "dl_reconciliation_likelihood": dl_rec,
                "delta_reconciliation_likelihood_dtl_minus_dl": dtl_rec - dl_rec,
                "twice_delta_reconciliation_likelihood": 2 * (dtl_rec - dl_rec),
                "dtl_duplication_rate": dtl_rates[0],
                "dtl_loss_rate": dtl_rates[1],
                "dtl_transfer_rate": dtl_rates[2],
                "dl_duplication_rate": dl_rates[0],
                "dl_loss_rate": dl_rates[1],
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "generax_dtl_vs_dl_likelihood.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "software": "GeneRax 2.1.3",
                "strategy": "EVAL",
                "models": ["UndatedDTL", "UndatedDL"],
                "formal_likelihood_ratio_test": False,
                "reason": "The transfer-rate null is on a model boundary; delta likelihoods are reported descriptively.",
                "known_software_issue": "GeneRax 2.1.3 segfaulted while writing per-species event counts after valid likelihood and rate files were written.",
            },
            indent=2,
        )
        + "\n"
    )
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
