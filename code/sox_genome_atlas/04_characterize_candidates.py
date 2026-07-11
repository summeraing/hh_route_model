from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_taxonomy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", names=["gtdb_accession", "gtdb_lineage"], dtype=str)
    frame["assembly_accession"] = frame["gtdb_accession"].str.replace(
        r"^(RS_|GB_)", "", regex=True
    )
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    split = frame["gtdb_lineage"].str.split(";", expand=True)
    for index, rank in enumerate(ranks):
        frame[f"gtdb_{rank}"] = split[index] if index in split.columns else ""
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--bac-taxonomy", type=Path, required=True)
    parser.add_argument("--ar-taxonomy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str).fillna("")
    hits = pd.read_csv(args.hits, sep="\t", dtype=str).fillna("")
    universe = pd.read_csv(args.universe, sep="\t", dtype=str).fillna("")
    taxonomy = pd.concat(
        [read_taxonomy(args.bac_taxonomy), read_taxonomy(args.ar_taxonomy)], ignore_index=True
    ).drop_duplicates("assembly_accession")

    candidates = candidates.merge(universe, on="assembly_accession", how="left")
    candidates = candidates.merge(taxonomy, on="assembly_accession", how="left")

    gene_sets = candidates["sox_genes"].map(lambda value: set(filter(None, value.split(","))))
    candidates["has_soxB"] = gene_sets.map(lambda genes: "soxb" in genes)
    candidates["has_carrier_YZ"] = gene_sets.map(lambda genes: bool(genes & {"soxy", "soxz"}))
    candidates["has_cytochrome_AX"] = gene_sets.map(lambda genes: bool(genes & {"soxa", "soxx"}))
    candidates["has_dehydrogenase_CD"] = gene_sets.map(lambda genes: bool(genes & {"soxc", "soxd"}))
    candidates["module_completeness"] = (
        candidates[["has_soxB", "has_carrier_YZ", "has_cytochrome_AX", "has_dehydrogenase_CD"]]
        .sum(axis=1)
        .astype(int)
    )
    candidates["module_class"] = pd.cut(
        candidates["module_completeness"],
        bins=[-1, 2, 3, 4],
        labels=["minimal", "near_complete", "complete"],
    ).astype(str)
    candidates["gtdb_mapped"] = candidates["gtdb_lineage"].fillna("").ne("")

    core_hits = hits.loc[hits["feature_class"].eq("sox_core")].copy()
    core_hits["product_lower"] = core_hits["product"].str.lower()
    suspect = core_hits.loc[
        core_hits["product_lower"].str.contains(
            "sarcosine oxidase|superoxide|hypothetical protein", regex=True, na=False
        )
    ].copy()

    taxon_summary = (
        candidates.groupby(["gtdb_phylum", "gtdb_class", "module_class"], dropna=False)
        .agg(
            genomes=("assembly_accession", "nunique"),
            orders=("gtdb_order", "nunique"),
            genera=("gtdb_genus", "nunique"),
            median_sox_genes=("n_distinct_sox_genes", lambda x: pd.to_numeric(x).median()),
        )
        .reset_index()
        .sort_values("genomes", ascending=False)
    )

    summary = pd.DataFrame(
        [
            {"metric": "candidate_genomes", "value": len(candidates)},
            {"metric": "gtdb_mapped", "value": int(candidates["gtdb_mapped"].sum())},
            {"metric": "complete_modules", "value": int((candidates["module_class"] == "complete").sum())},
            {"metric": "near_complete_modules", "value": int((candidates["module_class"] == "near_complete").sum())},
            {"metric": "minimal_modules", "value": int((candidates["module_class"] == "minimal").sum())},
            {"metric": "distinct_gtdb_phyla", "value": int(candidates["gtdb_phylum"].nunique())},
            {"metric": "distinct_gtdb_orders", "value": int(candidates["gtdb_order"].nunique())},
            {"metric": "suspect_sox_products", "value": len(suspect)},
        ]
    )

    candidates.to_csv(args.outdir / "sox_candidate_genomes_annotated.tsv", sep="\t", index=False)
    candidates["assembly_accession"].to_csv(
        args.outdir / "sox_candidate_accessions.txt", index=False, header=False
    )
    taxon_summary.to_csv(args.outdir / "sox_taxonomic_distribution.tsv", sep="\t", index=False)
    suspect.to_csv(args.outdir / "suspect_sox_products.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "candidate_qc_summary.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nTop taxonomic strata")
    print(taxon_summary.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
