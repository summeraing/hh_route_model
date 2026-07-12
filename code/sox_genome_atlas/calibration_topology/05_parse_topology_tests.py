#!/usr/bin/env python3
"""Parse IQ-TREE user-tree topology tests and apply family-wise Holm correction."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def holm(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    out = np.empty_like(values, dtype=float)
    running = 0.0
    m = len(values)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * float(values[idx]))
        running = max(running, value)
        out[idx] = running
    return out


def parse_user_trees(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip().startswith("USER TREES")), None)
    if start is None:
        raise ValueError(f"USER TREES section missing: {path}")
    header_idx = next((i for i in range(start, min(start + 30, len(lines))) if "Tree" in lines[i] and "logL" in lines[i]), None)
    if header_idx is None:
        raise ValueError(f"Topology-test header missing: {path}")
    header = re.split(r"\s+", lines[header_idx].strip())
    rows = []
    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("-"):
            if rows:
                break
            continue
        fields = re.split(r"\s+", stripped)
        if not fields[0].isdigit():
            if rows:
                break
            continue
        rows.append(fields)
    if len(rows) < 2:
        raise ValueError(f"Expected at least two tested trees: {path}")

    # IQ-TREE appends +/- significance markers as separate columns. Parse the
    # stable numeric order documented for topology tests.
    parsed = []
    for fields in rows:
        numeric = []
        for token in fields:
            try:
                numeric.append(float(token))
            except ValueError:
                continue
        if len(numeric) < 8:
            raise ValueError(f"Could not parse topology row: {fields}")
        parsed.append({
            "tree_index": int(numeric[0]),
            "log_likelihood": numeric[1],
            "delta_log_likelihood": numeric[2],
            "bp_rell": numeric[3],
            "p_kh": numeric[4],
            "p_sh": numeric[5],
            "p_wkh": numeric[6] if len(numeric) >= 10 else np.nan,
            "p_wsh": numeric[7] if len(numeric) >= 10 else np.nan,
            "c_elw": numeric[-2],
            "p_au": numeric[-1],
        })
    return pd.DataFrame(parsed)


def main() -> None:
    args = parse_args()
    root = Path(args.root)
    rows = []
    for design_dir in sorted(root.iterdir()):
        if not design_dir.is_dir():
            continue
        for family_dir in sorted(design_dir.iterdir()):
            report = family_dir / "topology_test.iqtree"
            if not report.exists():
                continue
            parsed = parse_user_trees(report)
            parsed["tree_design"] = design_dir.name
            parsed["family"] = family_dir.name
            parsed["tree_label"] = ["unconstrained_ML", "GTDB_constrained_ML"] + [
                f"tree_{i}" for i in range(3, len(parsed) + 1)
            ]
            rows.append(parsed)
    if not rows:
        raise SystemExit("No topology-test reports found")
    frame = pd.concat(rows, ignore_index=True)
    constrained = frame.loc[frame["tree_label"] == "GTDB_constrained_ML"].copy()
    for design in constrained["tree_design"].unique():
        idx = constrained["tree_design"] == design
        for col in ("p_sh", "p_wsh", "p_au"):
            constrained.loc[idx, f"{col}_holm"] = holm(constrained.loc[idx, col].to_numpy(dtype=float))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(outdir / "sox_topology_test_all_trees.tsv", sep="\t", index=False)
    constrained.to_csv(outdir / "sox_gtdb_constraint_tests.tsv", sep="\t", index=False)
    summary = {
        "reports": int(frame.groupby(["tree_design", "family"]).ngroups),
        "families": sorted(frame["family"].unique().tolist()),
        "tree_designs": sorted(frame["tree_design"].unique().tolist()),
        "constrained_rejected_au_005": int((constrained["p_au"] < 0.05).sum()),
        "constrained_rejected_au_holm_005": int((constrained["p_au_holm"] < 0.05).sum()),
    }
    (outdir / "analysis_parameters.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
