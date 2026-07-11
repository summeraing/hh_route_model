from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd


SOX_GENES = {"soxa", "soxx", "soxy", "soxz", "soxb", "soxc", "soxd"}
MOBILE_TERMS = re.compile(
    r"transposase|integrase|relaxase|insertion sequence|mobile element|conjugative",
    flags=re.IGNORECASE,
)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def canonical_sox(gene: str, product: str) -> str:
    token = re.sub(r"[^a-z0-9]", "", gene.lower())
    if "sarcosine oxidase" in product.lower():
        return ""
    if token in SOX_GENES:
        return token
    match = re.search(r"\bSox([ABCDXYZ])\b", product, flags=re.IGNORECASE)
    if match and re.search(r"sulfur|thiosulf|sulfane", product, flags=re.IGNORECASE):
        return f"sox{match.group(1).lower()}"
    return ""


def interval_distance(a0: int, a1: int, b0: int, b1: int) -> int:
    if a1 < b0:
        return b0 - a1
    if b1 < a0:
        return a0 - b1
    return 0


def mobile_count(mobile: list[tuple[int, int]], start: int, end: int, window: int) -> int:
    return sum(interval_distance(left, right, start, end) <= window for left, right in mobile)


def sign_flip_test(values: np.ndarray, permutations: int, rng: np.random.Generator) -> dict[str, float]:
    values = values[np.isfinite(values)]
    observed = float(values.mean()) if len(values) else float("nan")
    if not len(values):
        return {
            "observed": observed,
            "p_greater": float("nan"),
            "p_less": float("nan"),
            "p_two_sided": float("nan"),
            "null_q025": float("nan"),
            "null_q975": float("nan"),
        }
    null = np.empty(permutations)
    for index in range(permutations):
        null[index] = np.mean(values * rng.choice([-1.0, 1.0], size=len(values)))
    return {
        "observed": observed,
        "p_greater": float((1 + np.sum(null >= observed)) / (permutations + 1)),
        "p_less": float((1 + np.sum(null <= observed)) / (permutations + 1)),
        "p_two_sided": float((1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--candidate-accessions", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--controls-per-cluster", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260714)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    windows = [5000, 10000, 20000]
    max_window = max(windows)
    data_root = args.package / "ncbi_dataset" / "data"
    accessions = args.candidate_accessions.read_text().split()
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.drop_duplicates("assembly_accession").set_index("assembly_accession")

    rows: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    for accession in accessions:
        folder = data_root / accession
        gffs = list(folder.glob("*.gff")) + list(folder.glob("*.gff.gz"))
        if not gffs:
            failures.append({"assembly_accession": accession, "reason": "missing_gff"})
            continue
        contig_max: dict[str, int] = {}
        sox: dict[str, list[tuple[int, int, str]]] = {}
        mobile: dict[str, list[tuple[int, int]]] = {}
        with open_text(gffs[0]) as handle:
            for line in handle:
                if line.startswith("##sequence-region"):
                    fields = line.split()
                    if len(fields) >= 4:
                        contig_max[fields[1]] = int(fields[3])
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                ann = parse_attrs(fields[8])
                if ann.get("pseudo", "").lower() == "true" or ann.get("partial", "").lower() == "true":
                    continue
                seqid = fields[0]
                start, end = int(fields[3]), int(fields[4])
                contig_max[seqid] = max(contig_max.get(seqid, 0), end)
                gene = ann.get("gene", "")
                product = ann.get("product", "")
                token = canonical_sox(gene, product)
                if token:
                    sox.setdefault(seqid, []).append((start, end, token))
                gene_simple = re.sub(r"[^a-z0-9]", "", gene.lower())
                if gene_simple not in {"reca", "xerc", "xerd"} and MOBILE_TERMS.search(f"{gene} {product}"):
                    mobile.setdefault(seqid, []).append((start, end))

        cluster_index = 0
        for seqid, features in sox.items():
            ordered = sorted(features)
            clusters: list[list[tuple[int, int, str]]] = []
            for feature in ordered:
                if not clusters or feature[0] - clusters[-1][-1][1] > 20000:
                    clusters.append([])
                clusters[-1].append(feature)
            intervals = [(min(x[0] for x in cluster), max(x[1] for x in cluster)) for cluster in clusters]
            contig_length = contig_max.get(seqid, 0)
            for cluster, (start, end) in zip(clusters, intervals):
                cluster_index += 1
                cluster_id = f"{accession}:{seqid}:sox_cluster_{cluster_index}"
                span = end - start + 1
                controls: list[tuple[int, int]] = []
                attempts = 0
                while len(controls) < args.controls_per_cluster and attempts < 5000:
                    attempts += 1
                    if contig_length <= span + 2:
                        break
                    control_start = int(rng.integers(1, contig_length - span + 1))
                    control_end = control_start + span - 1
                    if any(interval_distance(control_start, control_end, a, b) <= 2 * max_window for a, b in intervals):
                        continue
                    controls.append((control_start, control_end))
                if len(controls) < 5:
                    failures.append({"assembly_accession": accession, "reason": f"insufficient_controls:{cluster_id}"})
                    continue
                for window in windows:
                    observed_count = mobile_count(mobile.get(seqid, []), start, end, window)
                    for control_id, (control_start, control_end) in enumerate(controls, start=1):
                        control_count = mobile_count(mobile.get(seqid, []), control_start, control_end, window)
                        rows.append(
                            {
                                "assembly_accession": accession,
                                "gtdb_class": metadata.at[accession, "gtdb_class"] if accession in metadata.index else "",
                                "gtdb_order": metadata.at[accession, "gtdb_order"] if accession in metadata.index else "",
                                "seqid": seqid,
                                "cluster_id": cluster_id,
                                "sox_genes": ";".join(x[2] for x in cluster),
                                "window_bp": window,
                                "control_id": control_id,
                                "observed_mobile_count": observed_count,
                                "control_mobile_count": control_count,
                                "observed_mobile_present": observed_count > 0,
                                "control_mobile_present": control_count > 0,
                                "count_difference": observed_count - control_count,
                                "presence_difference": int(observed_count > 0) - int(control_count > 0),
                            }
                        )

    pairs = pd.DataFrame(rows)
    if pairs.empty:
        raise RuntimeError("No matched control rows were generated")
    cluster = (
        pairs.groupby(["assembly_accession", "cluster_id", "window_bp"], as_index=False)
        .agg(
            observed_mobile_count=("observed_mobile_count", "first"),
            control_mobile_count_mean=("control_mobile_count", "mean"),
            observed_mobile_present=("observed_mobile_present", "first"),
            control_mobile_present_mean=("control_mobile_present", "mean"),
            controls=("control_id", "nunique"),
        )
    )
    cluster["count_difference"] = cluster["observed_mobile_count"] - cluster["control_mobile_count_mean"]
    cluster["presence_difference"] = cluster["observed_mobile_present"].astype(float) - cluster["control_mobile_present_mean"]
    genome = (
        cluster.groupby(["assembly_accession", "window_bp"], as_index=False)
        .agg(
            count_difference=("count_difference", "mean"),
            presence_difference=("presence_difference", "mean"),
            clusters=("cluster_id", "nunique"),
        )
    )

    summary_rows: list[dict[str, object]] = []
    for window, group in genome.groupby("window_bp"):
        for metric in ["count_difference", "presence_difference"]:
            stats = sign_flip_test(group[metric].to_numpy(float), args.permutations, rng)
            summary_rows.append(
                {
                    "window_bp": window,
                    "metric": metric,
                    "genomes": len(group),
                    "clusters": int(cluster.loc[cluster["window_bp"].eq(window), "cluster_id"].nunique()),
                    **stats,
                }
            )
    summary = pd.DataFrame(summary_rows)
    class_summary = (
        pairs.groupby(["gtdb_class", "window_bp"], as_index=False)
        .agg(
            genomes=("assembly_accession", "nunique"),
            clusters=("cluster_id", "nunique"),
            observed_mobile_prevalence=("observed_mobile_present", "mean"),
            matched_control_prevalence=("control_mobile_present", "mean"),
            mean_count_difference=("count_difference", "mean"),
        )
    )

    pairs.to_csv(args.outdir / "matched_mobile_boundary_pairs.tsv.gz", sep="\t", index=False, compression="gzip")
    cluster.to_csv(args.outdir / "matched_mobile_boundary_cluster_summary.tsv", sep="\t", index=False)
    genome.to_csv(args.outdir / "matched_mobile_boundary_genome_summary.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "matched_mobile_boundary_tests.tsv", sep="\t", index=False)
    class_summary.to_csv(args.outdir / "matched_mobile_boundary_class_summary.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(args.outdir / "matched_mobile_boundary_failures.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "controls_per_cluster": args.controls_per_cluster,
                "windows_bp": windows,
                "control_exclusion_distance_bp": 2 * max_window,
                "permutations": args.permutations,
                "seed": args.seed,
                "dependency_unit": "genome mean across SOX clusters",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"pairs={len(pairs)} clusters={cluster['cluster_id'].nunique()} genomes={genome['assembly_accession'].nunique()} failures={len(failures)}")


if __name__ == "__main__":
    main()
