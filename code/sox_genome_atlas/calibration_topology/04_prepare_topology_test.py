#!/usr/bin/env python3
"""Relabel a pruned GTDB topology with gene-sequence tip names."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--species-tree", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mapping = {}
    for line in Path(args.mapping).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        gene, genome = line.split()[:2]
        if genome in mapping:
            raise SystemExit(f"Non-single-copy mapping for {genome}")
        mapping[genome] = gene

    tree = Path(args.species_tree).read_text(encoding="utf-8").strip()
    observed = set(re.findall(r"(?<=[(,])([^():;,]+)(?=:)", tree))
    missing = sorted(observed - set(mapping))
    if missing:
        raise SystemExit(f"Species-tree tips missing from mapping: {missing[:10]}")
    for genome in sorted(observed, key=len, reverse=True):
        tree = re.sub(
            rf"(?<=[(,]){re.escape(genome)}(?=:)",
            mapping[genome],
            tree,
        )
    gene_tips = set(re.findall(r"(?<=[(,])([^():;,]+)(?=:)", tree))
    expected = {mapping[x] for x in observed}
    if gene_tips != expected:
        raise SystemExit("Relabelled constraint tree failed tip-set audit")

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(tree + "\n", encoding="utf-8")
    output.with_suffix(".json").write_text(json.dumps({
        "species_tips": len(observed),
        "gene_tips": len(gene_tips),
        "mapping_rows": len(mapping),
    }, indent=2), encoding="utf-8")
    print(f"constraint_tips={len(gene_tips)}")


if __name__ == "__main__":
    main()
