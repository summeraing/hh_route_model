from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


SOX_GENES = {"soxa", "soxx", "soxy", "soxz", "soxb", "soxc", "soxd"}
SUPPORT_GENES = {
    "fcca",
    "fccb",
    "sqr",
    "dsra",
    "dsrb",
    "apra",
    "aprb",
    "sat",
    "hdra",
    "hdrb",
    "hdrc",
}
MOBILE_TERMS = re.compile(
    r"transposase|integrase|recombinase|resolvase|relaxase|insertion sequence|"
    r"mobile element|conjugative|site-specific recombination",
    flags=re.IGNORECASE,
)


def canonical_gene(gene: str, product: str) -> str:
    token = re.sub(r"[^a-z0-9]", "", gene.lower())
    product_lower = product.lower()
    if "sarcosine oxidase" in product_lower:
        return token
    if token in SOX_GENES | SUPPORT_GENES:
        return token
    match = re.search(r"\bSox([ABCDXYZ])\b", product, flags=re.IGNORECASE)
    if match and re.search(r"sulfur|thiosulf|sulfane", product, flags=re.IGNORECASE):
        return f"sox{match.group(1).lower()}"
    return token


@dataclass
class Feature:
    accession: str
    seqid: str
    start: int
    end: int
    strand: str
    gene: str
    locus_tag: str
    product: str
    feature_class: str
    is_pseudogene: bool
    is_partial: bool


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        result[key] = unquote(value)
    return result


def classify(gene: str, product: str) -> str:
    token = canonical_gene(gene, product)
    product_lower = product.lower()
    if "sarcosine oxidase" in product_lower:
        return "other_cds"
    if token in SOX_GENES:
        return "sox_core"
    if token in SUPPORT_GENES and not (
        token == "apra" and "metalloproteinase" in product_lower
    ):
        return "sulfur_support"
    if token in {"reca", "xerc", "xerd"}:
        return "other_cds"
    if MOBILE_TERMS.search(f"{gene} {product}"):
        return "mobile_context"
    return "other_cds"


def assembly_from_path(path: Path) -> str:
    for part in reversed(path.parts):
        if part.startswith(("GCF_", "GCA_")):
            return part
    return path.parent.name


def read_gff(path: Path) -> list[Feature]:
    accession = assembly_from_path(path)
    features: list[Feature] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attributes(fields[8])
            gene_raw = attrs.get("gene", "")
            product = attrs.get("product", "")
            feature_class = classify(gene_raw, product)
            gene = canonical_gene(gene_raw, product) or gene_raw
            is_pseudogene = attrs.get("pseudo", "").lower() == "true"
            is_partial = attrs.get("partial", "").lower() == "true"
            features.append(
                Feature(
                    accession=accession,
                    seqid=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    gene=gene,
                    locus_tag=attrs.get("locus_tag", attrs.get("ID", "")),
                    product=product,
                    feature_class=feature_class,
                    is_pseudogene=is_pseudogene,
                    is_partial=is_partial,
                )
            )
    return features


def nearest_mobile(feature: Feature, mobile: list[Feature]) -> tuple[int | None, str]:
    candidates = [item for item in mobile if item.seqid == feature.seqid]
    if not candidates:
        return None, ""
    nearest = min(
        candidates,
        key=lambda item: min(abs(item.start - feature.end), abs(feature.start - item.end)),
    )
    distance = min(abs(nearest.start - feature.end), abs(feature.start - nearest.end))
    return distance, nearest.locus_tag


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-list", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    paths = [Path(line.strip()) for line in args.gff_list.read_text().splitlines() if line.strip()]
    hit_rows: list[dict[str, object]] = []
    all_sulfur_rows: list[dict[str, object]] = []
    genome_rows: list[dict[str, object]] = []
    parse_failures: list[dict[str, str]] = []

    for index, path in enumerate(paths, start=1):
        try:
            features = read_gff(path)
        except Exception as exc:
            parse_failures.append({"path": str(path), "error": repr(exc)})
            continue
        targets = [feature for feature in features if feature.feature_class != "other_cds"]
        active_targets = [
            feature for feature in targets if not feature.is_pseudogene and not feature.is_partial
        ]
        all_sox = [feature for feature in targets if feature.feature_class == "sox_core"]
        sox = [feature for feature in active_targets if feature.feature_class == "sox_core"]
        mobile = [feature for feature in active_targets if feature.feature_class == "mobile_context"]
        sox_names = {re.sub(r"[^a-z0-9]", "", feature.gene.lower()) for feature in sox}
        has_carrier_or_cytochrome = bool(sox_names & {"soxx", "soxy", "soxz"})
        is_candidate = "soxb" in sox_names and has_carrier_or_cytochrome and len(sox_names) >= 3

        for feature in targets:
            if feature.feature_class not in {"sox_core", "sulfur_support"}:
                continue
            distance, nearest_tag = nearest_mobile(feature, mobile)
            all_sulfur_rows.append(
                {
                    **feature.__dict__,
                    "is_strict_sox_candidate_genome": is_candidate,
                    "nearest_mobile_distance_bp": distance,
                    "nearest_mobile_locus_tag": nearest_tag,
                    "within_20kb_mobile": distance is not None and distance <= 20000,
                }
            )

        if is_candidate:
            for feature in active_targets:
                distance, nearest_tag = nearest_mobile(feature, mobile)
                hit_rows.append(
                    {
                        **feature.__dict__,
                        "nearest_mobile_distance_bp": distance,
                        "nearest_mobile_locus_tag": nearest_tag,
                        "within_20kb_mobile": distance is not None and distance <= 20000,
                    }
                )
            genome_rows.append(
                {
                    "assembly_accession": sox[0].accession,
                    "n_cds": len(features),
                    "n_sox_hits": len(sox),
                    "n_excluded_pseudogene_or_partial_sox_hits": len(all_sox) - len(sox),
                    "n_distinct_sox_genes": len(sox_names),
                    "sox_genes": ",".join(sorted(sox_names)),
                    "n_mobile_features": len(mobile),
                    "sox_with_mobile_within_20kb": sum(
                        1
                        for feature in sox
                        if (nearest_mobile(feature, mobile)[0] is not None)
                        and nearest_mobile(feature, mobile)[0] <= 20000
                    ),
                }
            )
        if index % 500 == 0:
            print(f"processed={index}/{len(paths)} candidates={len(genome_rows)}", flush=True)

    hits = pd.DataFrame(hit_rows)
    all_sulfur = pd.DataFrame(all_sulfur_rows)
    genomes = pd.DataFrame(genome_rows)
    failures = pd.DataFrame(parse_failures)
    hits.to_csv(args.outdir / "sox_candidate_target_hits.tsv", sep="\t", index=False)
    all_sulfur.to_csv(args.outdir / "all_core_sulfur_hits.tsv", sep="\t", index=False)
    genomes.to_csv(args.outdir / "sox_candidate_genomes.tsv", sep="\t", index=False)
    failures.to_csv(args.outdir / "gff_parse_failures.tsv", sep="\t", index=False)

    print(f"gff_files={len(paths)}")
    print(f"candidate_genomes={len(genomes)}")
    print(f"target_hits_in_candidates={len(hits)}")
    print(f"sulfur_marker_hits_in_full_universe={len(all_sulfur)}")
    if not all_sulfur.empty:
        excluded = all_sulfur["is_pseudogene"].astype(bool) | all_sulfur["is_partial"].astype(bool)
        print(f"excluded_pseudogene_or_partial_sulfur_hits={int(excluded.sum())}")
    print(f"parse_failures={len(failures)}")


if __name__ == "__main__":
    main()
