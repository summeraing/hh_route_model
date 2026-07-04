from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any


def find_summary_jsons(root: Path) -> list[Path]:
    patterns = ["*summary_confidences*.json", "*ranking*.json", "*scores*.json"]
    found: list[Path] = []
    for pattern in patterns:
        found.extend(root.rglob(pattern))
    return sorted(p for p in set(found) if p.name != "confidences.json")


def flatten_numeric(d: dict[str, Any]) -> dict[str, float]:
    out = {}
    for key in ["iptm", "ipTM", "ptm", "pTM", "ranking_score", "rankingScore"]:
        val = d.get(key)
        if isinstance(val, (int, float)):
            out[key.lower().replace("rankingscore", "ranking_score")] = float(val)
    return out


def infer_job_id(path: Path, known_job_ids: set[str] | None = None) -> str:
    if known_job_ids:
        known_by_lower = {job_id.lower(): job_id for job_id in known_job_ids}
        parts = list(path.parts)
        for part in reversed(parts):
            hit = known_by_lower.get(part.lower())
            if hit:
                return hit
        text = str(path)
        text_lower = text.lower()
        hits = [job_id for job_id in known_job_ids if job_id.lower() in text_lower]
        if hits:
            return sorted(hits, key=len, reverse=True)[0]

    name = path.name
    for suffix in ["_summary_confidences", "_ranking", "_scores"]:
        if suffix in name:
            return name.split(suffix)[0]
    return path.parent.name


def gate_from_iptm(values: list[float]) -> str:
    if not values:
        return "missing"
    med = statistics.median(values)
    high_frac = sum(v >= 0.60 for v in values) / len(values)
    mod_frac = sum(v >= 0.40 for v in values) / len(values)
    if med >= 0.60 and high_frac >= 0.60:
        return "high"
    if med >= 0.40 and mod_frac >= 0.60:
        return "moderate"
    return "low"


def q1_q3(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], values[0]
    qs = statistics.quantiles(values, n=4, method="inclusive")
    return qs[0], qs[2]


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize AF3 confidence JSONs for SOX validation jobs.")
    parser.add_argument("--af3-output-dir", required=True, type=Path)
    parser.add_argument("--target-matrix", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    args = parser.parse_args()

    target_rows = {}
    with args.target_matrix.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            target_rows[row["job_id"]] = row

    records = []
    for path in find_summary_jsons(args.af3_output_dir):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        nums = flatten_numeric(data if isinstance(data, dict) else {})
        if not nums:
            continue
        job_id = infer_job_id(path, set(target_rows))
        records.append({
            "job_id": job_id,
            "file": str(path),
            "iptm": nums.get("iptm"),
            "ptm": nums.get("ptm"),
            "ranking_score": nums.get("ranking_score"),
        })

    grouped: dict[str, list[dict[str, Any]]] = {}
    for rec in records:
        grouped.setdefault(rec["job_id"], []).append(rec)

    rows = []
    for job_id, recs in sorted(grouped.items()):
        iptms = [r["iptm"] for r in recs if isinstance(r.get("iptm"), float)]
        ptms = [r["ptm"] for r in recs if isinstance(r.get("ptm"), float)]
        ranks = [r["ranking_score"] for r in recs if isinstance(r.get("ranking_score"), float)]
        iqr = q1_q3(iptms)
        base = target_rows.get(job_id, {})
        rows.append({
            "job_id": job_id,
            "comparison_group": base.get("comparison_group", ""),
            "route_class": base.get("route_class", ""),
            "chain_set": base.get("chain_set", ""),
            "expected_gate": base.get("expected_gate", ""),
            "n_models": len(recs),
            "median_iptm": statistics.median(iptms) if iptms else "",
            "q1_iptm": iqr[0] if iqr[0] is not None else "",
            "q3_iptm": iqr[1] if iqr[1] is not None else "",
            "median_ptm": statistics.median(ptms) if ptms else "",
            "median_ranking_score": statistics.median(ranks) if ranks else "",
            "observed_gate": gate_from_iptm(iptms),
        })

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "job_id", "comparison_group", "route_class", "chain_set", "expected_gate",
        "n_models", "median_iptm", "q1_iptm", "q3_iptm", "median_ptm",
        "median_ranking_score", "observed_gate",
    ]
    with args.output_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"summary_jsons={len(records)}")
    print(f"jobs={len(rows)}")
    print(f"output_csv={args.output_csv}")


if __name__ == "__main__":
    main()
