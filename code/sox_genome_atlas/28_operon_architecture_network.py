from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SOX = {"soxa", "soxx", "soxy", "soxz", "soxb", "soxc", "soxd"}


def entropy(labels: pd.Series) -> float:
    probabilities = labels.value_counts(normalize=True).to_numpy(float)
    return float(-(probabilities * np.log(probabilities)).sum()) if len(probabilities) else 0.0


def normalized_mutual_information(left: pd.Series, right: pd.Series) -> float:
    frame = pd.DataFrame({"left": left.astype(str), "right": right.astype(str)})
    joint = frame.value_counts(normalize=True)
    p_left = frame["left"].value_counts(normalize=True)
    p_right = frame["right"].value_counts(normalize=True)
    mi = 0.0
    for (a, b), p_ab in joint.items():
        mi += p_ab * np.log(p_ab / (p_left[a] * p_right[b]))
    denominator = np.sqrt(entropy(frame["left"]) * entropy(frame["right"]))
    return float(mi / denominator) if denominator > 0 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260715)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    features = pd.read_csv(args.features, sep="\t", dtype=str).fillna("")
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.drop_duplicates("assembly_accession")
    sox = features.loc[
        features["organizational_role"].eq("sox_energy_module")
        & features["gene"].str.lower().isin(SOX)
        & features["sox_cluster_id"].ne("")
    ].copy()
    sox["start"] = pd.to_numeric(sox["start"])

    cluster_rows: list[dict[str, object]] = []
    cluster_gene_lists: dict[str, list[str]] = {}
    for cluster_id, group in sox.groupby("sox_cluster_id"):
        ordered = group.sort_values("start")
        genes = ordered["gene"].str.lower().tolist()
        minus_fraction = ordered["strand"].eq("-").mean()
        if minus_fraction > 0.5:
            genes = list(reversed(genes))
        signature = "-".join(genes)
        cluster_gene_lists[cluster_id] = genes
        cluster_rows.append(
            {
                "cluster_id": cluster_id,
                "assembly_accession": ordered["assembly_accession"].iloc[0],
                "seqid": ordered["seqid"].iloc[0],
                "genes": ";".join(genes),
                "architecture_signature": signature,
                "genes_in_cluster": len(genes),
                "distinct_sox_families": len(set(genes)),
                "majority_strand": "-" if minus_fraction > 0.5 else "+",
                "strand_coherence": max(minus_fraction, 1 - minus_fraction),
            }
        )
    clusters = pd.DataFrame(cluster_rows).merge(metadata, on="assembly_accession", how="left")

    architecture = (
        clusters.groupby("architecture_signature", as_index=False)
        .agg(
            clusters=("cluster_id", "nunique"),
            genomes=("assembly_accession", "nunique"),
            classes=("gtdb_class", "nunique"),
            orders=("gtdb_order", "nunique"),
            median_gene_count=("genes_in_cluster", "median"),
            mean_strand_coherence=("strand_coherence", "mean"),
        )
        .sort_values(["clusters", "architecture_signature"], ascending=[False, True])
    )

    observed_edges: Counter[tuple[str, str]] = Counter()
    clusters_with_pair: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    genomes_with_pair: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    cluster_to_genome = clusters.set_index("cluster_id")["assembly_accession"].to_dict()
    for cluster_id, genes in cluster_gene_lists.items():
        for left, right in zip(genes, genes[1:]):
            pair = tuple(sorted((left, right)))
            observed_edges[pair] += 1
            clusters_with_pair[pair].add(cluster_id)
            genomes_with_pair[pair].add(cluster_to_genome[cluster_id])

    active_genes = sorted(set(sox["gene"].str.lower()))
    pairs = sorted(
        {tuple(sorted((a, b))) for a in active_genes for b in active_genes if a != b}
    )
    null_counts = {pair: np.zeros(args.permutations, dtype=int) for pair in pairs}
    null_genome_counts = {pair: np.zeros(args.permutations, dtype=int) for pair in pairs}
    for permutation in range(args.permutations):
        counts: Counter[tuple[str, str]] = Counter()
        genome_pairs: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
        for cluster_id, genes in cluster_gene_lists.items():
            shuffled = list(rng.permutation(genes))
            for left, right in zip(shuffled, shuffled[1:]):
                pair = tuple(sorted((left, right)))
                counts[pair] += 1
                genome_pairs[pair].add(cluster_to_genome[cluster_id])
        for pair in pairs:
            null_counts[pair][permutation] = counts[pair]
            null_genome_counts[pair][permutation] = len(genome_pairs[pair])

    edge_rows: list[dict[str, object]] = []
    for pair in pairs:
        observed = observed_edges[pair]
        null = null_counts[pair]
        observed_genomes = len(genomes_with_pair[pair])
        genome_null = null_genome_counts[pair]
        edge_rows.append(
            {
                "gene_a": pair[0],
                "gene_b": pair[1],
                "observed_adjacent_clusters": observed,
                "null_mean": float(null.mean()),
                "null_q025": float(np.quantile(null, 0.025)),
                "null_q975": float(np.quantile(null, 0.975)),
                "enrichment_ratio": float(observed / null.mean()) if null.mean() > 0 else float("nan"),
                "permutation_p_greater": float((1 + np.sum(null >= observed)) / (args.permutations + 1)),
                "observed_adjacent_genomes": observed_genomes,
                "genome_null_mean": float(genome_null.mean()),
                "genome_null_q025": float(np.quantile(genome_null, 0.025)),
                "genome_null_q975": float(np.quantile(genome_null, 0.975)),
                "genome_enrichment_ratio": float(observed_genomes / genome_null.mean())
                if genome_null.mean() > 0
                else float("nan"),
                "genome_permutation_p_greater": float(
                    (1 + np.sum(genome_null >= observed_genomes)) / (args.permutations + 1)
                ),
            }
        )
    edges = pd.DataFrame(edge_rows)

    largest_per_genome = (
        clusters.sort_values(
            ["assembly_accession", "distinct_sox_families", "genes_in_cluster", "cluster_id"],
            ascending=[True, False, False, True],
        )
        .drop_duplicates("assembly_accession")
        .copy()
    )
    subsets = {
        "all_clusters": clusters,
        "multigene_clusters": clusters.loc[clusters["genes_in_cluster"] >= 2],
        "near_complete_or_complete_clusters": clusters.loc[clusters["distinct_sox_families"] >= 5],
        "largest_cluster_per_genome": largest_per_genome,
    }
    association_rows: list[dict[str, object]] = []
    for subset_name, subset in subsets.items():
        for rank in ["gtdb_class", "gtdb_order"]:
            valid = subset.loc[subset[rank].ne("")].copy()
            observed = normalized_mutual_information(valid["architecture_signature"], valid[rank])
            null = np.empty(args.permutations)
            labels = valid[rank].to_numpy(copy=True)
            for index in range(args.permutations):
                null[index] = normalized_mutual_information(
                    valid["architecture_signature"],
                    pd.Series(rng.permutation(labels), index=valid.index),
                )
            association_rows.append(
                {
                    "subset": subset_name,
                    "taxonomic_rank": rank,
                    "clusters": len(valid),
                    "genomes": valid["assembly_accession"].nunique(),
                    "categories": valid[rank].nunique(),
                    "observed_nmi": observed,
                    "null_mean": float(null.mean()),
                    "null_q975": float(np.quantile(null, 0.975)),
                    "permutation_p_greater": float(
                        (1 + np.sum(null >= observed)) / (args.permutations + 1)
                    ),
                }
            )
    association = pd.DataFrame(association_rows)

    nodes = (
        sox.groupby("gene", as_index=False)
        .agg(features=("locus_tag", "size"), genomes=("assembly_accession", "nunique"), clusters=("sox_cluster_id", "nunique"))
        .rename(columns={"gene": "node"})
    )
    clusters.to_csv(args.outdir / "sox_operon_architecture_catalog.tsv", sep="\t", index=False)
    architecture.to_csv(args.outdir / "sox_operon_architecture_frequency.tsv", sep="\t", index=False)
    nodes.to_csv(args.outdir / "sox_operon_network_nodes.tsv", sep="\t", index=False)
    edges.to_csv(args.outdir / "sox_operon_network_edges.tsv", sep="\t", index=False)
    association.to_csv(args.outdir / "sox_architecture_taxonomy_association.tsv", sep="\t", index=False)
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(
            {
                "cluster_definition": "SOX genes separated by no more than 20 kb on the same contig",
                "orientation_rule": "reverse order when the majority strand is negative",
                "adjacency_null": "within-cluster gene-order permutation preserving gene content and cluster size",
                "permutations": args.permutations,
                "seed": args.seed,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"clusters={len(clusters)} architectures={len(architecture)} edges={len(edges)}")


if __name__ == "__main__":
    main()
