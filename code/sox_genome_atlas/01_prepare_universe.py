from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


LEVELS = {"Complete Genome", "Chromosome"}
CATEGORIES = {"reference genome", "representative genome"}
GROUPS = {"bacteria", "archaea"}


def read_summary(path: Path) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    rows: list[list[str]] = []
    malformed: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        first = handle.readline()
        header = next(csv.reader([handle.readline()], delimiter="\t"))
        header = [column.lstrip("#") for column in header]
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), start=3):
            if len(fields) != len(header):
                malformed.append(
                    {
                        "line_number": line_number,
                        "observed_fields": len(fields),
                        "expected_fields": len(header),
                        "assembly_accession": fields[0] if fields else "",
                    }
                )
                continue
            rows.append(fields)
    return pd.DataFrame(rows, columns=header).fillna("na"), malformed


def usable(frame: pd.DataFrame) -> pd.Series:
    excluded = frame["excluded_from_refseq"].str.lower().isin({"na", "", "none"})
    return (
        frame["group"].str.lower().isin(GROUPS)
        & frame["version_status"].eq("latest")
        & frame["assembly_level"].isin(LEVELS)
        & frame["ftp_path"].str.startswith("https://")
        & excluded
    )


def add_urls(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    stem = out["ftp_path"].str.rsplit("/", n=1).str[-1]
    out["gff_url"] = out["ftp_path"] + "/" + stem + "_genomic.gff.gz"
    out["protein_url"] = out["ftp_path"] + "/" + stem + "_protein.faa.gz"
    out["genome_url"] = out["ftp_path"] + "/" + stem + "_genomic.fna.gz"
    return out


def summarize(name: str, frame: pd.DataFrame) -> pd.DataFrame:
    summary = (
        frame.groupby(["group", "refseq_category", "assembly_level"], dropna=False)
        .size()
        .rename("assemblies")
        .reset_index()
    )
    summary.insert(0, "universe", name)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    frame, malformed = read_summary(args.assembly_summary)
    broad = add_urls(frame.loc[usable(frame)].copy())
    core = broad.loc[broad["refseq_category"].isin(CATEGORIES)].copy()

    columns = [
        "assembly_accession",
        "refseq_category",
        "taxid",
        "species_taxid",
        "organism_name",
        "assembly_level",
        "group",
        "genome_size",
        "gc_percent",
        "scaffold_count",
        "contig_count",
        "annotation_provider",
        "annotation_date",
        "ftp_path",
        "gff_url",
        "protein_url",
        "genome_url",
    ]

    core[columns].to_csv(args.outdir / "refseq_core_universe.tsv", sep="\t", index=False)
    broad[columns].to_csv(args.outdir / "refseq_broad_universe.tsv", sep="\t", index=False)
    core["assembly_accession"].to_csv(
        args.outdir / "refseq_core_accessions.txt", index=False, header=False
    )

    overview = pd.concat(
        [summarize("core_reference_representative", core), summarize("broad_complete_chromosome", broad)],
        ignore_index=True,
    )
    overview.to_csv(args.outdir / "universe_summary.tsv", sep="\t", index=False)
    pd.DataFrame(malformed).to_csv(
        args.outdir / "malformed_assembly_summary_rows.tsv", sep="\t", index=False
    )

    print(f"all_refseq_rows={len(frame)}")
    print(f"malformed_rows_skipped={len(malformed)}")
    print(f"broad_complete_or_chromosome={len(broad)}")
    print(f"core_reference_or_representative={len(core)}")
    print(overview.to_string(index=False))


if __name__ == "__main__":
    main()
