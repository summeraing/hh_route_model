from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


EVENT_RE = re.compile(r"^(S|SL|D|T|TL|L|Leaf):(\d+)$")
VALUE_RE = re.compile(r"(S|SL|D|T|TL|L)=(\d+)")


def read_event_counts(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    for line in path.read_text().splitlines():
        match = EVENT_RE.match(line.strip())
        if match:
            result[match.group(1)] = int(match.group(2))
    return result


def read_per_species(path: Path, family: str, samples: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        node = fields[0]
        values = {key: int(value) for key, value in VALUE_RE.findall(line)}
        rows.append(
            {
                "family": family,
                "node": node,
                "is_terminal_genome": node.startswith(("GCF_", "GCA_")),
                **{f"{key}_count": values.get(key, 0) for key in ["S", "SL", "D", "T", "TL", "L"]},
                **{
                    f"mean_{key}_events_per_sample": values.get(key, 0) / samples
                    for key in ["S", "SL", "D", "T", "TL", "L"]
                },
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generax-root", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--discordance", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=100)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.drop_duplicates("assembly_accession").set_index("assembly_accession")
    discordance = pd.read_csv(args.discordance, sep="\t")
    discordance["family"] = discordance["family"].astype(str).str.lower()

    family_rows: list[dict[str, object]] = []
    per_species_frames: list[pd.DataFrame] = []
    edge_rows: list[dict[str, object]] = []

    for family_dir in sorted(path for path in args.generax_root.iterdir() if path.is_dir()):
        family = family_dir.name.lower()
        rec = family_dir / "reconciliations"
        event_file = rec / f"{family}_eventCounts.txt"
        species_file = family_dir / "per_species_event_counts.txt"
        if not event_file.exists() or not species_file.exists():
            continue
        counts = read_event_counts(event_file)
        per_species = read_per_species(species_file, family, args.samples)
        tip = per_species.loc[per_species["is_terminal_genome"]].copy()

        sample_files = sorted(
            path for path in rec.glob(f"{family}_*_transfers.txt")
            if re.search(rf"{re.escape(family)}_(\d+)_transfers\.txt$", path.name)
        )
        edge_counter: Counter[tuple[str, str]] = Counter()
        recipient_sample_counter: Counter[str] = Counter()
        for path in sample_files:
            seen: set[tuple[str, str]] = set()
            for line in path.read_text().splitlines():
                fields = line.split()
                if len(fields) >= 2:
                    seen.add((fields[0], fields[1]))
            edge_counter.update(seen)
            recipient_sample_counter.update({recipient for _, recipient in seen})

        denominator = max(1, len(sample_files))
        per_species["transfer_recipient_sample_fraction"] = (
            per_species["node"].map(recipient_sample_counter).fillna(0).astype(int) / denominator
        )
        per_species_frames.append(per_species)
        tip = per_species.loc[per_species["is_terminal_genome"]].copy()

        for (donor, recipient), support_count in edge_counter.items():
            donor_tip = donor.startswith(("GCF_", "GCA_"))
            recipient_tip = recipient.startswith(("GCF_", "GCA_"))
            row: dict[str, object] = {
                "family": family,
                "donor_node": donor,
                "recipient_node": recipient,
                "donor_is_terminal": donor_tip,
                "recipient_is_terminal": recipient_tip,
                "sample_support_count": support_count,
                "sample_support": support_count / max(1, len(sample_files)),
            }
            if donor_tip and donor in metadata.index:
                row["donor_class"] = metadata.at[donor, "gtdb_class"]
                row["donor_order"] = metadata.at[donor, "gtdb_order"]
            if recipient_tip and recipient in metadata.index:
                row["recipient_class"] = metadata.at[recipient, "gtdb_class"]
                row["recipient_order"] = metadata.at[recipient, "gtdb_order"]
            edge_rows.append(row)

        leaves = counts.get("Leaf", len(tip))
        family_rows.append(
            {
                "family": family,
                **counts,
                "transfer_events_per_100_tips": 100 * counts.get("T", 0) / max(1, leaves),
                "terminal_branches": len(tip),
                "terminal_recipient_fraction_ge_0.5": int(
                    (tip["transfer_recipient_sample_fraction"] >= 0.5).sum()
                ),
                "terminal_recipient_fraction_ge_0.9": int(
                    (tip["transfer_recipient_sample_fraction"] >= 0.9).sum()
                ),
                "median_terminal_recipient_sample_fraction": float(
                    tip["transfer_recipient_sample_fraction"].median()
                ),
                "max_terminal_recipient_sample_fraction": float(
                    tip["transfer_recipient_sample_fraction"].max()
                ),
                "max_terminal_mean_transfer_events_per_sample": float(
                    tip["mean_T_events_per_sample"].max()
                ),
                "reconciliation_samples": len(sample_files),
            }
        )

    family_summary = pd.DataFrame(family_rows).merge(discordance, on="family", how="left")
    per_species_all = pd.concat(per_species_frames, ignore_index=True)
    edges = pd.DataFrame(edge_rows).fillna("")
    if not edges.empty:
        terminal_edges = edges.loc[edges["donor_is_terminal"] & edges["recipient_is_terminal"]].copy()
        taxonomic_edges = (
            terminal_edges.groupby(
                ["family", "donor_class", "recipient_class"], dropna=False
            )
            .agg(
                terminal_edges=("sample_support", "size"),
                mean_sample_support=("sample_support", "mean"),
                max_sample_support=("sample_support", "max"),
            )
            .reset_index()
        )
    else:
        terminal_edges = pd.DataFrame()
        taxonomic_edges = pd.DataFrame()

    family_summary.to_csv(args.outdir / "generax_family_summary.tsv", sep="\t", index=False)
    per_species_all.to_csv(args.outdir / "generax_per_branch_support.tsv", sep="\t", index=False)
    edges.to_csv(args.outdir / "generax_sampled_transfer_edges.tsv", sep="\t", index=False)
    terminal_edges.to_csv(args.outdir / "generax_terminal_transfer_edges.tsv", sep="\t", index=False)
    taxonomic_edges.to_csv(args.outdir / "generax_taxonomic_transfer_summary.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "software": "GeneRax 2.0.4",
                "model": "UndatedDTL",
                "strategy": "EVAL",
                "unrooted_gene_tree": True,
                "reconciliation_samples": args.samples,
                "interpretation": "Post-boundary sensitivity; does not alter frozen donor labels or route ranking.",
            },
            indent=2,
        )
        + "\n"
    )
    try:
        with pd.ExcelWriter(args.outdir / "generax_dtl_source_data.xlsx") as writer:
            family_summary.to_excel(writer, sheet_name="family_summary", index=False)
            per_species_all.to_excel(writer, sheet_name="branch_support", index=False)
            edges.to_excel(writer, sheet_name="sampled_edges", index=False)
            terminal_edges.to_excel(writer, sheet_name="terminal_edges", index=False)
            taxonomic_edges.to_excel(writer, sheet_name="taxonomic_edges", index=False)
    except ImportError:
        pass
    print(f"families={len(family_summary)} branches={len(per_species_all)} sampled_edges={len(edges)}")


if __name__ == "__main__":
    main()
