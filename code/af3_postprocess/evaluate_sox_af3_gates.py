from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


GATE_SCORE = {"missing": -1, "low": 0, "moderate": 1, "high": 2}
REQUIRED_PRIORITY1_JOBS = {
    "SOX_AF3_POS_001",
    "SOX_AF3_POS_002",
    "SOX_AF3_FUNC_003",
    "SOX_AF3_POS_006",
    "SOX_AF3_NEG_010",
    "SOX_AF3_NEG_011",
    "SOX_AF3_NEG_012",
}


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def gate_score(value: str) -> int:
    return GATE_SCORE.get((value or "").strip().lower(), -1)


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate SOX-AF3 positive/negative gate separation.")
    parser.add_argument("--summary-csv", required=True, type=Path)
    parser.add_argument("--output-md", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    with args.summary_csv.open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))

    for row in rows:
        group = row.get("comparison_group", "").lower()
        route = row.get("route_class", "").lower()
        if "negative" in group or "negative" in route or "decoy" in route:
            row["gate_family"] = "negative"
        elif "functional" in group or "functional" in route:
            row["gate_family"] = "functional"
        else:
            row["gate_family"] = "positive"
        row["gate_score"] = str(gate_score(row.get("observed_gate", "")))

    families: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        families.setdefault(row["gate_family"], []).append(row)

    summary_rows = []
    for family, frecs in sorted(families.items()):
        iptms = [as_float(r.get("median_iptm", "")) for r in frecs]
        iptms = [v for v in iptms if v is not None]
        scores = [int(r["gate_score"]) for r in frecs]
        summary_rows.append({
            "gate_family": family,
            "n_jobs": len(frecs),
            "median_iptm": statistics.median(iptms) if iptms else "",
            "mean_gate_score": statistics.mean(scores) if scores else "",
            "jobs": ";".join(r.get("job_id", "") for r in frecs),
            "observed_gates": ";".join(r.get("observed_gate", "") for r in frecs),
        })

    pos_scores = [int(r["gate_score"]) for r in families.get("positive", [])]
    func_scores = [int(r["gate_score"]) for r in families.get("functional", [])]
    neg_scores = [int(r["gate_score"]) for r in families.get("negative", [])]

    pos_ok = pos_scores and statistics.median(pos_scores) >= 1
    neg_ok = neg_scores and statistics.median(neg_scores) <= 0
    func_ok = (not func_scores) or statistics.median(func_scores) >= 1

    observed_jobs = {row.get("job_id", "") for row in rows}
    missing_priority1 = sorted(REQUIRED_PRIORITY1_JOBS - observed_jobs)
    has_positive = bool(pos_scores)
    has_negative = bool(neg_scores)
    has_functional = bool(func_scores)

    if not rows:
        decision = "Pending: no completed SOX-AF3 jobs were parsed; do not make a structural claim."
    elif missing_priority1 or not (has_positive and has_negative and has_functional):
        decision = (
            "Partial screen: wait for remaining priority-1 SOX-AF3 jobs before deciding "
            "Vita/iMeta structural claims."
        )
    elif pos_ok and func_ok and neg_ok:
        decision = "Vita-first possible: structural gates separate."
    elif pos_scores and neg_scores and statistics.median(pos_scores) > statistics.median(neg_scores):
        decision = "iMeta-first: partial structural separation; keep SOX-AF3 as support or Supplementary."
    else:
        decision = "Do not claim SOX structural validation in the main text."

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = ["gate_family", "n_jobs", "median_iptm", "mean_gate_score", "jobs", "observed_gates"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    md = [
        "# SOX-AF3 gate evaluation",
        "",
        f"- Input summary: `{args.summary_csv}`",
        f"- Decision: **{decision}**",
        f"- Parsed jobs: {len(observed_jobs)}",
        f"- Missing priority-1 jobs: {', '.join(missing_priority1) if missing_priority1 else 'none'}",
        "",
        "## Family summary",
        "",
        "| gate family | n jobs | median ipTM | mean gate score | observed gates |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        md.append(
            f"| {row['gate_family']} | {row['n_jobs']} | {row['median_iptm']} | "
            f"{row['mean_gate_score']} | {row['observed_gates']} |"
        )
    md.extend([
        "",
        "Gate score: low=0, moderate=1, high=2. Missing=-1.",
        "",
        "Interpretation rule: Vita-style claims require cognate positives and functional assemblies to be moderate/high while mobile-boundary negatives remain low.",
    ])
    args.output_md.write_text("\n".join(md) + "\n", encoding="utf-8")

    print(decision)
    print(f"output_md={args.output_md}")
    print(f"output_csv={args.output_csv}")


if __name__ == "__main__":
    main()
