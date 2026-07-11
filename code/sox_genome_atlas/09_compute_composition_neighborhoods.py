from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


SOX_GENES = {"soxa", "soxx", "soxy", "soxz", "soxb", "soxc", "soxd"}
SUPPORT_GENES = {"fcca", "fccb", "sqr", "dsra", "dsrb", "apra", "aprb", "sat", "hdra", "hdrb", "hdrc"}
MOBILE_TERMS = re.compile(
    r"transposase|integrase|relaxase|insertion sequence|mobile element|conjugative",
    flags=re.IGNORECASE,
)


def attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def canonical_gene(gene: str, product: str) -> str:
    token = re.sub(r"[^a-z0-9]", "", gene.lower())
    product_lower = product.lower()
    if "sarcosine oxidase" in product_lower:
        return ""
    if token in SOX_GENES | SUPPORT_GENES:
        if token == "apra" and "metalloproteinase" in product_lower:
            return ""
        return token
    match = re.search(r"\bSox([ABCDXYZ])\b", product, flags=re.IGNORECASE)
    if match and re.search(r"sulfur|thiosulf|sulfane", product, flags=re.IGNORECASE):
        return f"sox{match.group(1).lower()}"
    return ""


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name:
            sequences[name] = "".join(chunks).upper()
    return sequences


def gc(sequence: str) -> float:
    valid = sum(base in "ACGT" for base in sequence)
    if not valid:
        return float("nan")
    return sum(base in "GC" for base in sequence) / valid


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if end_a < start_b:
        return start_b - end_a
    if end_b < start_a:
        return start_a - end_b
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--candidate-accessions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--window", type=int, default=10000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    candidates = set(args.candidate_accessions.read_text().split())
    data_root = args.package / "ncbi_dataset" / "data"
    rows: list[dict[str, object]] = []
    genome_rows: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []

    for index, accession in enumerate(sorted(candidates), start=1):
        folder = data_root / accession
        gffs = list(folder.glob("*.gff")) + list(folder.glob("*.gff.gz"))
        fnas = list(folder.glob("*.fna")) + list(folder.glob("*.fna.gz"))
        if not gffs or not fnas:
            failures.append({"assembly_accession": accession, "reason": "missing_gff_or_genome"})
            continue
        sequences = read_fasta(fnas[0])
        assembly_sequence = "".join(sequences.values())
        genome_gc = gc(assembly_sequence)
        all_cds: list[dict[str, object]] = []
        targets: list[dict[str, object]] = []
        excluded_pseudogene_or_partial_sulfur_features = 0
        with open_text(gffs[0]) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                annotation = attrs(fields[8])
                gene_raw = annotation.get("gene", "")
                product = annotation.get("product", "")
                token = canonical_gene(gene_raw, product)
                gene_simple = re.sub(r"[^a-z0-9]", "", gene_raw.lower())
                excluded = (
                    annotation.get("pseudo", "").lower() == "true"
                    or annotation.get("partial", "").lower() == "true"
                )
                if excluded:
                    if token in SOX_GENES | SUPPORT_GENES:
                        excluded_pseudogene_or_partial_sulfur_features += 1
                    continue
                feature = {
                    "assembly_accession": accession,
                    "seqid": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene": token or gene_simple,
                    "locus_tag": annotation.get("locus_tag", annotation.get("ID", "")),
                    "product": product,
                    "is_mobile_element": bool(
                        gene_simple not in {"reca", "xerc", "xerd"}
                        and MOBILE_TERMS.search(f"{gene_raw} {product}")
                    ),
                }
                if token in SOX_GENES:
                    feature["organizational_role"] = "sox_energy_module"
                    feature["role_basis"] = "sulfur_sox_family_identity"
                    targets.append(feature)
                elif token in SUPPORT_GENES:
                    feature["organizational_role"] = "recipient_sulfur_backbone"
                    feature["role_basis"] = "sulfur_backbone_family_identity"
                    targets.append(feature)
                all_cds.append(feature)

        mobile = [feature for feature in all_cds if feature["is_mobile_element"]]

        sox_by_contig: dict[str, list[dict[str, object]]] = {}
        for feature in targets:
            if feature["organizational_role"] == "sox_energy_module":
                sox_by_contig.setdefault(str(feature["seqid"]), []).append(feature)
        cluster_lookup: dict[tuple[str, str], str] = {}
        cluster_intervals: list[dict[str, object]] = []
        for seqid, features in sox_by_contig.items():
            ordered = sorted(features, key=lambda item: int(item["start"]))
            clusters: list[list[dict[str, object]]] = []
            for feature in ordered:
                if not clusters or int(feature["start"]) - int(clusters[-1][-1]["end"]) > 20000:
                    clusters.append([])
                clusters[-1].append(feature)
            for cluster_index, cluster in enumerate(clusters, start=1):
                cluster_id = f"{accession}:{seqid}:sox_cluster_{cluster_index}"
                cluster_start = min(int(feature["start"]) for feature in cluster)
                cluster_end = max(int(feature["end"]) for feature in cluster)
                cluster_intervals.append(
                    {
                        "seqid": seqid,
                        "cluster_id": cluster_id,
                        "start": cluster_start,
                        "end": cluster_end,
                    }
                )
                for feature in cluster:
                    cluster_lookup[(str(feature["seqid"]), str(feature["locus_tag"]))] = cluster_id

        # A boundary role is defined only by position outside a SOX cluster. Mobile
        # annotation remains a separate donor-evidence feature and is never used to
        # assign the organizational role.
        boundary_by_locus: dict[tuple[str, str], dict[str, object]] = {}
        target_loci = {(str(feature["seqid"]), str(feature["locus_tag"])) for feature in targets}
        for feature in all_cds:
            key = (str(feature["seqid"]), str(feature["locus_tag"]))
            if key in target_loci:
                continue
            compatible = [cluster for cluster in cluster_intervals if cluster["seqid"] == feature["seqid"]]
            if not compatible:
                continue
            nearest_cluster = min(
                compatible,
                key=lambda cluster: interval_distance(
                    int(feature["start"]), int(feature["end"]), int(cluster["start"]), int(cluster["end"])
                ),
            )
            distance = interval_distance(
                int(feature["start"]),
                int(feature["end"]),
                int(nearest_cluster["start"]),
                int(nearest_cluster["end"]),
            )
            if distance == 0 or distance > args.window:
                continue
            boundary = dict(feature)
            boundary["organizational_role"] = "mobility_boundary"
            boundary["role_basis"] = "positional_flank_outside_sox_cluster"
            boundary["sox_cluster_id"] = nearest_cluster["cluster_id"]
            boundary["distance_to_sox_cluster_bp"] = distance
            boundary["boundary_side"] = (
                "left" if int(feature["end"]) < int(nearest_cluster["start"]) else "right"
            )
            previous = boundary_by_locus.get(key)
            if previous is None or distance < int(previous["distance_to_sox_cluster_bp"]):
                boundary_by_locus[key] = boundary

        selected_features = targets + list(boundary_by_locus.values())

        for feature in targets:
            key = (str(feature["seqid"]), str(feature["locus_tag"]))
            if key in cluster_lookup:
                feature["sox_cluster_id"] = cluster_lookup[key]
                feature["distance_to_sox_cluster_bp"] = 0
            else:
                compatible = [cluster for cluster in cluster_intervals if cluster["seqid"] == feature["seqid"]]
                if compatible:
                    nearest_cluster = min(
                        compatible,
                        key=lambda cluster: interval_distance(
                            int(feature["start"]), int(feature["end"]), int(cluster["start"]), int(cluster["end"])
                        ),
                    )
                    feature["sox_cluster_id"] = nearest_cluster["cluster_id"]
                    feature["distance_to_sox_cluster_bp"] = interval_distance(
                        int(feature["start"]),
                        int(feature["end"]),
                        int(nearest_cluster["start"]),
                        int(nearest_cluster["end"]),
                    )
                else:
                    feature["sox_cluster_id"] = ""
                    feature["distance_to_sox_cluster_bp"] = None
            feature["boundary_side"] = ""

        for feature in selected_features:
            seqid = str(feature["seqid"])
            sequence = sequences.get(seqid, "")
            start = int(feature["start"])
            end = int(feature["end"])
            gene_sequence = sequence[max(0, start - 1) : end]
            local_sequence = sequence[max(0, start - 1 - args.window) : min(len(sequence), end + args.window)]
            same_contig_mobile = [item for item in mobile if item["seqid"] == seqid]
            nearest = None
            distance = None
            if same_contig_mobile:
                nearest = min(
                    same_contig_mobile,
                    key=lambda item: min(abs(int(item["start"]) - end), abs(start - int(item["end"]))),
                )
                distance = min(abs(int(nearest["start"]) - end), abs(start - int(nearest["end"])))
            gene_gc = gc(gene_sequence)
            local_gc = gc(local_sequence)
            rows.append(
                {
                    **feature,
                    "gene_gc": gene_gc,
                    "local_20kb_gc": local_gc,
                    "genome_gc": genome_gc,
                    "gene_gc_deviation": gene_gc - genome_gc,
                    "local_gc_deviation": local_gc - genome_gc,
                    "nearest_mobile_distance_bp": distance,
                    "nearest_mobile_locus_tag": "" if nearest is None else nearest["locus_tag"],
                    "mobile_within_20kb": distance is not None and distance <= 20000,
                }
            )
        genome_rows.append(
            {
                "assembly_accession": accession,
                "genome_gc": genome_gc,
                "genome_length": len(assembly_sequence),
                "target_features": len(targets),
                "sox_features": sum(item["organizational_role"] == "sox_energy_module" for item in targets),
                "backbone_features": sum(item["organizational_role"] == "recipient_sulfur_backbone" for item in targets),
                "mobile_features": len(mobile),
                "sox_clusters": len(set(cluster_lookup.values())),
                "boundary_context_features": len(boundary_by_locus),
                "mobile_boundary_context_features": sum(
                    bool(feature["is_mobile_element"]) for feature in boundary_by_locus.values()
                ),
                "excluded_pseudogene_or_partial_sulfur_features": excluded_pseudogene_or_partial_sulfur_features,
            }
        )
        if index % 50 == 0:
            print(f"processed={index}/{len(candidates)}", flush=True)

    pd.DataFrame(rows).to_csv(args.outdir / "sox_composition_neighborhood_features.tsv", sep="\t", index=False)
    pd.DataFrame(genome_rows).to_csv(args.outdir / "sox_genome_architecture_summary.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(args.outdir / "composition_failures.tsv", sep="\t", index=False)
    print(f"genomes={len(genome_rows)} feature_rows={len(rows)} failures={len(failures)}")


if __name__ == "__main__":
    main()
