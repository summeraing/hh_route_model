from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

import dendropy
from dendropy.calculate import treecompare


def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8")


def normalize_species_label(label: str) -> str:
    if label.startswith(("RS_", "GB_")):
        return label[3:]
    return label


def normalize_gene_label(label: str) -> str:
    return label.split("__", 1)[0]


def normalize_tree(tree: dendropy.Tree, mode: str) -> None:
    seen: set[str] = set()
    for leaf in tree.leaf_node_iter():
        original = leaf.taxon.label
        label = normalize_species_label(original) if mode == "species" else normalize_gene_label(original)
        if label in seen:
            raise ValueError(f"duplicate normalized tip: {label}")
        seen.add(label)
        leaf.taxon.label = label


def coherence(tree: dendropy.Tree, rank_by_accession: dict[str, str], k: int = 10) -> dict[str, float]:
    matrix = tree.phylogenetic_distance_matrix()
    taxa = [leaf.taxon for leaf in tree.leaf_node_iter()]
    output: dict[str, float] = {}
    for taxon in taxa:
        rank = rank_by_accession.get(taxon.label, "")
        others = sorted(
            ((matrix.distance(taxon, other), other) for other in taxa if other is not taxon),
            key=lambda item: item[0],
        )[: min(k, len(taxa) - 1)]
        if not rank or not others:
            output[taxon.label] = float("nan")
            continue
        output[taxon.label] = sum(rank_by_accession.get(other.label, "") == rank for _, other in others) / len(others)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-tree", type=Path, required=True)
    parser.add_argument("--candidate-metadata", type=Path, required=True)
    parser.add_argument("--gene-tree-root", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.candidate_metadata, sep="\t", dtype=str).fillna("")
    order_by_accession = dict(zip(metadata["assembly_accession"], metadata["gtdb_order"]))
    class_by_accession = dict(zip(metadata["assembly_accession"], metadata["gtdb_class"]))
    species_text = read_text(args.species_tree)
    families = [line.strip() for line in args.families.read_text().splitlines() if line.strip()]
    global_rows: list[dict[str, object]] = []
    tip_rows: list[dict[str, object]] = []

    for family in families:
        treefile = args.gene_tree_root / family / f"{family}.treefile"
        if not treefile.exists():
            global_rows.append({"family": family, "status": "missing_gene_tree"})
            continue
        species = dendropy.Tree.get(data=species_text, schema="newick", preserve_underscores=True)
        gene = dendropy.Tree.get(path=treefile, schema="newick", preserve_underscores=True)
        normalize_tree(species, "species")
        normalize_tree(gene, "gene")
        species_labels = {leaf.taxon.label for leaf in species.leaf_node_iter()}
        gene_labels = {leaf.taxon.label for leaf in gene.leaf_node_iter()}
        common = species_labels & gene_labels
        if len(common) < 10:
            global_rows.append({"family": family, "status": "fewer_than_10_shared_tips", "shared_tips": len(common)})
            continue
        species.retain_taxa_with_labels(common)
        gene.retain_taxa_with_labels(common)
        namespace = dendropy.TaxonNamespace(sorted(common))
        species.migrate_taxon_namespace(namespace)
        gene.migrate_taxon_namespace(namespace)
        species.is_rooted = False
        gene.is_rooted = False
        species.encode_bipartitions()
        gene.encode_bipartitions()
        rf = treecompare.symmetric_difference(species, gene)
        max_rf = max(1, 2 * (len(common) - 3))
        weighted_rf = treecompare.weighted_robinson_foulds_distance(species, gene)

        gene_order = coherence(gene, order_by_accession)
        species_order = coherence(species, order_by_accession)
        gene_class = coherence(gene, class_by_accession)
        species_class = coherence(species, class_by_accession)
        for accession in sorted(common):
            tip_rows.append(
                {
                    "family": family,
                    "assembly_accession": accession,
                    "gtdb_order": order_by_accession.get(accession, ""),
                    "gtdb_class": class_by_accession.get(accession, ""),
                    "gene_tree_order_coherence_k10": gene_order.get(accession),
                    "species_tree_order_coherence_k10": species_order.get(accession),
                    "order_discordance": species_order.get(accession, float("nan")) - gene_order.get(accession, float("nan")),
                    "gene_tree_class_coherence_k10": gene_class.get(accession),
                    "species_tree_class_coherence_k10": species_class.get(accession),
                    "class_discordance": species_class.get(accession, float("nan")) - gene_class.get(accession, float("nan")),
                }
            )
        global_rows.append(
            {
                "family": family,
                "status": "ok",
                "shared_tips": len(common),
                "symmetric_difference": rf,
                "normalized_rf": rf / max_rf,
                "weighted_rf": weighted_rf,
                "mean_order_discordance": pd.Series(
                    [row["order_discordance"] for row in tip_rows if row["family"] == family]
                ).mean(),
                "mean_class_discordance": pd.Series(
                    [row["class_discordance"] for row in tip_rows if row["family"] == family]
                ).mean(),
            }
        )

    global_frame = pd.DataFrame(global_rows)
    tip_frame = pd.DataFrame(tip_rows)
    global_frame.to_csv(args.outdir / "gene_species_global_discordance.tsv", sep="\t", index=False)
    tip_frame.to_csv(args.outdir / "gene_species_tip_discordance.tsv", sep="\t", index=False)
    print(global_frame.to_string(index=False))
    print(f"tip_rows={len(tip_frame)}")


if __name__ == "__main__":
    main()
