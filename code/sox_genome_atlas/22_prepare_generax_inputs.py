from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
from pathlib import Path

import dendropy
import pandas as pd


def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8")


def normalize_species(label: str) -> str:
    return label[3:] if label.startswith(("RS_", "GB_")) else label


def gene_species(label: str) -> str:
    return label.split("__", 1)[0]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    sequence: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(sequence)))
                name = line[1:].strip().split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
        if name:
            records.append((name, "".join(sequence)))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-tree", type=Path, required=True)
    parser.add_argument("--gene-tree-root", type=Path, required=True)
    parser.add_argument("--alignment-root", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    spec_bytes = args.spec.read_bytes()
    spec = json.loads(spec_bytes)
    spec_hash = hashlib.sha256(spec_bytes).hexdigest()

    species_tree = dendropy.Tree.get(
        data=read_text(args.species_tree), schema="newick", preserve_underscores=True
    )
    seen: set[str] = set()
    for leaf in species_tree.leaf_node_iter():
        label = normalize_species(leaf.taxon.label)
        if label in seen:
            raise ValueError(f"Duplicate normalized species-tree label: {label}")
        seen.add(label)
        leaf.taxon.label = label

    families = [line.strip() for line in args.families.read_text().splitlines() if line.strip()]
    requested = set(spec["families"])
    families = [family for family in families if family in requested]
    if set(families) != requested:
        raise ValueError(f"Missing requested families: {sorted(requested - set(families))}")

    family_records: list[dict[str, object]] = []
    union_species: set[str] = set()
    blocks = ["[FAMILIES]"]
    for family in families:
        source_alignment = args.alignment_root / family / f"{family}.trimmed.faa"
        source_tree = args.gene_tree_root / family / f"{family}.treefile"
        if not source_alignment.exists() or not source_tree.exists():
            raise FileNotFoundError(
                f"Missing GeneRax input for {family}: {source_alignment} or {source_tree}"
            )
        records = read_fasta(source_alignment)
        retained = [(name, sequence) for name, sequence in records if gene_species(name) in seen]
        retained_labels = [name for name, _ in retained]
        retained_species = [gene_species(name) for name in retained_labels]
        if len(retained) < 3:
            raise ValueError(f"{family} has fewer than three GTDB-tree-mapped sequences")
        if len(retained_species) != len(set(retained_species)):
            raise ValueError(f"{family} is not single-copy after GTDB-tree filtering")

        family_dir = args.outdir / family
        family_dir.mkdir(exist_ok=True)
        alignment = family_dir / f"{family}.gtdb_shared.trimmed.faa"
        write_fasta(retained, alignment)

        gene_tree = dendropy.Tree.get(
            path=source_tree, schema="newick", preserve_underscores=True
        )
        gene_tree.retain_taxa_with_labels(set(retained_labels))
        gene_tree.is_rooted = False
        tree = family_dir / f"{family}.gtdb_shared.tree"
        gene_tree.write(
            path=tree,
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True,
        )

        mapping = family_dir / f"{family}.mapping.link"
        mapping.write_text(
            "\n".join(f"{gene} {gene_species(gene)}" for gene in retained_labels) + "\n",
            encoding="utf-8",
        )
        family_species_tree = species_tree.clone(depth=2)
        family_species_tree.retain_taxa_with_labels(set(retained_species))
        family_species_tree.is_rooted = True
        family_species_tree.resolve_polytomies(
            update_bipartitions=False,
            rng=random.Random(int(spec["random_seed"]) + len(family_records)),
        )
        family_species_path = family_dir / f"{family}.gtdb_rooted_binary.tree"
        family_species_tree.write(
            path=family_species_path,
            schema="newick",
            suppress_rooting=True,
            suppress_internal_node_labels=True,
            unquoted_underscores=True,
        )
        family_file = family_dir / f"{family}.generax_family.txt"
        family_file.write_text(
            "\n".join(
                [
                    "[FAMILIES]",
                    f"- {family}",
                    f"starting_gene_tree = {tree}",
                    f"alignment = {alignment}",
                    f"mapping = {mapping}",
                    f"subst_model = {spec['substitution_model_for_evaluation']}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        union_species.update(retained_species)
        blocks.extend(
            [
                f"- {family}",
                f"starting_gene_tree = {tree}",
                f"alignment = {alignment}",
                f"mapping = {mapping}",
                f"subst_model = {spec['substitution_model_for_evaluation']}",
            ]
        )
        family_records.append(
            {
                "family": family,
                "input_sequences": len(records),
                "retained_gtdb_tree_sequences": len(retained),
                "excluded_nonrepresentative_or_unmapped": len(records) - len(retained),
                "species": len(set(retained_species)),
                "alignment": str(alignment),
                "starting_gene_tree": str(tree),
                "mapping": str(mapping),
                "family_species_tree": str(family_species_path),
                "family_file": str(family_file),
            }
        )

    species_tree.retain_taxa_with_labels(union_species)
    species_tree.is_rooted = True
    species_tree.resolve_polytomies(
        update_bipartitions=False, rng=random.Random(int(spec["random_seed"]))
    )
    pruned_tree = args.outdir / "gtdb_r232_sox_union_rooted_binary.tree"
    species_tree.write(
        path=pruned_tree,
        schema="newick",
        suppress_rooting=True,
        suppress_internal_node_labels=True,
        unquoted_underscores=True,
    )
    (args.outdir / "generax_families.txt").write_text(
        "\n".join(blocks) + "\n", encoding="utf-8"
    )
    pd.DataFrame(family_records).to_csv(
        args.outdir / "generax_input_manifest.tsv", sep="\t", index=False
    )
    (args.outdir / "analysis_spec_used.json").write_bytes(spec_bytes)
    (args.outdir / "analysis_spec_sha256.txt").write_text(
        spec_hash + "\n", encoding="utf-8"
    )
    print(
        f"families={len(families)} union_species={len(union_species)} "
        f"species_tree_tips={len(list(species_tree.leaf_node_iter()))}"
    )
    print(pd.DataFrame(family_records).to_string(index=False))
    print(f"spec_sha256={spec_hash}")


if __name__ == "__main__":
    main()
