from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


SOX_GENES = {"soxa", "soxx", "soxy", "soxz", "soxb", "soxc", "soxd"}


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def canonical_gene(gene: str, product: str) -> str:
    token = re.sub(r"[^a-z0-9]", "", gene.lower())
    if "sarcosine oxidase" in product.lower():
        return ""
    if token in SOX_GENES:
        return token
    match = re.search(r"\bSox([ABCDXYZ])\b", product, flags=re.IGNORECASE)
    if match and re.search(r"sulfur|thiosulf|sulfane", product, flags=re.IGNORECASE):
        return f"sox{match.group(1).lower()}"
    return ""


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks)
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current:
            sequences[current] = "".join(chunks)
    return sequences


def protein_key(attrs: dict[str, str]) -> str:
    protein_id = attrs.get("protein_id", "")
    if protein_id:
        return protein_id
    cds_id = attrs.get("ID", "")
    return re.sub(r"^cds-", "", cds_id)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--candidate-accessions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_dir = args.outdir / "fasta"
    fasta_dir.mkdir(exist_ok=True)
    single_copy_dir = args.outdir / "fasta_single_copy"
    single_copy_dir.mkdir(exist_ok=True)

    candidates = set(args.candidate_accessions.read_text().split())
    data_root = args.package / "ncbi_dataset" / "data"
    family_records: dict[str, list[tuple[str, str]]] = defaultdict(list)
    manifest_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, str]] = []

    for accession in sorted(candidates):
        folder = data_root / accession
        gff_files = list(folder.glob("*.gff")) + list(folder.glob("*.gff.gz"))
        faa_files = list(folder.glob("*.faa")) + list(folder.glob("*.faa.gz"))
        if not gff_files or not faa_files:
            missing_rows.append({"assembly_accession": accession, "reason": "missing_gff_or_faa"})
            continue
        sequences = read_fasta(faa_files[0])
        with open_text(gff_files[0]) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                attrs = attributes(fields[8])
                product = attrs.get("product", "")
                family = canonical_gene(attrs.get("gene", ""), product)
                if not family:
                    continue
                key = protein_key(attrs)
                if attrs.get("pseudo", "").lower() == "true" or attrs.get("partial", "").lower() == "true":
                    missing_rows.append(
                        {
                            "assembly_accession": accession,
                            "reason": "pseudogene_or_partial_cds_excluded",
                            "protein_id": key,
                            "family": family,
                        }
                    )
                    continue
                sequence = sequences.get(key, "")
                if not sequence:
                    missing_rows.append(
                        {
                            "assembly_accession": accession,
                            "reason": "protein_id_not_found",
                            "protein_id": key,
                            "family": family,
                        }
                    )
                    continue
                locus = attrs.get("locus_tag", key)
                sequence_id = f"{accession}__{family}__{key}__{locus}"
                family_records[family].append((sequence_id, sequence))
                manifest_rows.append(
                    {
                        "assembly_accession": accession,
                        "family": family,
                        "protein_id": key,
                        "locus_tag": locus,
                        "seqid": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "product": product,
                        "sequence_length": len(sequence),
                        "sequence_id": sequence_id,
                    }
                )

    family_summary = []
    for family, records in sorted(family_records.items()):
        out = fasta_dir / f"{family}.faa"
        with out.open("w", encoding="ascii") as handle:
            for identifier, sequence in records:
                handle.write(f">{identifier}\n")
                for index in range(0, len(sequence), 80):
                    handle.write(sequence[index : index + 80] + "\n")
        by_genome: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for identifier, sequence in records:
            by_genome[identifier.split("__", 1)[0]].append((identifier, sequence))
        single_copy_records = [items[0] for items in by_genome.values() if len(items) == 1]
        single_out = single_copy_dir / f"{family}.faa"
        with single_out.open("w", encoding="ascii") as handle:
            for identifier, sequence in single_copy_records:
                handle.write(f">{identifier}\n")
                for index in range(0, len(sequence), 80):
                    handle.write(sequence[index : index + 80] + "\n")
        family_summary.append(
            {
                "family": family,
                "sequences": len(records),
                "genomes": len(by_genome),
                "single_copy_genomes": len(single_copy_records),
                "multicopy_genomes": sum(len(items) > 1 for items in by_genome.values()),
                "fasta": str(out),
                "single_copy_fasta": str(single_out),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(args.outdir / "sox_sequence_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(missing_rows).to_csv(args.outdir / "sequence_extraction_failures.tsv", sep="\t", index=False)
    summary = pd.DataFrame(family_summary)
    summary.to_csv(args.outdir / "sox_family_sequence_counts.tsv", sep="\t", index=False)
    summary.loc[summary["single_copy_genomes"] >= 20, "family"].to_csv(
        args.outdir / "tree_families.txt", index=False, header=False
    )
    print(summary.to_string(index=False))
    print(f"manifest_rows={len(manifest_rows)} failures={len(missing_rows)}")


if __name__ == "__main__":
    main()
