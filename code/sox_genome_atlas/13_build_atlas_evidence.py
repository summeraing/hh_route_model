from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def read_taxonomy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", names=["gtdb_accession", "gtdb_lineage"], dtype=str)
    frame["assembly_accession"] = frame["gtdb_accession"].str.replace(r"^(RS_|GB_)", "", regex=True)
    split = frame["gtdb_lineage"].str.split(";", expand=True)
    frame["gtdb_class"] = split[2]
    frame["gtdb_order"] = split[3]
    return frame[["assembly_accession", "gtdb_class", "gtdb_order"]]


def donor_from_phylogeny(value: float | None, spec: dict) -> tuple[str, str]:
    if value is None or pd.isna(value):
        return "uncertain", "missing"
    threshold = spec["donor_thresholds"]["phylogeny"]
    if value >= threshold["hgt_like"]:
        return "hgt_like", "high" if value >= threshold["hgt_like_high"] else "medium"
    if value <= threshold["vertical_like"]:
        return "vertical_like", "high" if value <= threshold["vertical_like_high"] else "medium"
    return "uncertain", "low"


def donor_from_distribution(prevalence: float | None, spec: dict) -> tuple[str, str]:
    if prevalence is None or pd.isna(prevalence):
        return "uncertain", "missing"
    threshold = spec["donor_thresholds"]["order_distribution"]
    if prevalence <= threshold["hgt_like_max_prevalence"]:
        return "hgt_like", "high" if prevalence <= threshold["hgt_like_high_max_prevalence"] else "medium"
    if prevalence >= threshold["vertical_like_min_prevalence"]:
        return "vertical_like", "high" if prevalence >= threshold["vertical_like_high_min_prevalence"] else "medium"
    return "uncertain", "low"


def donor_from_composition(value: float | None, spec: dict) -> tuple[str, str]:
    if value is None or pd.isna(value):
        return "uncertain", "missing"
    absolute = abs(value)
    threshold = spec["donor_thresholds"]["composition"]
    if absolute >= threshold["hgt_like_min_absolute_deviation"]:
        return "hgt_like", "high" if absolute >= threshold["hgt_like_high_min_absolute_deviation"] else "medium"
    if absolute <= threshold["vertical_like_max_absolute_deviation"]:
        return "vertical_like", "high" if absolute <= threshold["vertical_like_high_max_absolute_deviation"] else "medium"
    return "uncertain", "low"


def donor_from_mobility(is_mobile_element: bool, distance: float | None, spec: dict) -> tuple[str, str]:
    if is_mobile_element:
        return "mobile_context", "high"
    if distance is not None and not pd.isna(distance) and distance <= spec["donor_thresholds"]["mobile_proximity_bp"]:
        return "hgt_like", "medium"
    return "uncertain", "low"


def aggregate_donor(group: pd.DataFrame, spec: dict) -> tuple[str, str, int, int]:
    layer = str(group["evidence_layer"].iloc[0])
    if layer == "mobile_proximity":
        if group["is_mobile_element"].astype(bool).any():
            return "mobile_context", "high", int(group["is_mobile_element"].astype(bool).sum()), len(group)
        if group["donor_class"].eq("hgt_like").any():
            return "hgt_like", "medium", int(group["donor_class"].eq("hgt_like").sum()), len(group)
        return "uncertain", "low", 0, len(group)

    eligible = group.loc[group["donor_class"].ne("uncertain"), "donor_class"]
    if eligible.empty:
        return "uncertain", "missing", 0, len(group)
    counts = eligible.value_counts()
    winners = counts.index[counts.eq(counts.iloc[0])].tolist()
    if len(winners) != 1:
        return "uncertain", "tie", int(counts.iloc[0]), len(group)
    winner = str(winners[0])
    vote_fraction = float(counts.iloc[0] / len(group))
    aggregation = spec["aggregation"]
    confidence = (
        "high"
        if vote_fraction >= aggregation["high_vote_fraction"]
        else "medium"
        if vote_fraction >= aggregation["medium_vote_fraction"]
        else "low"
    )
    return winner, confidence, int(counts.iloc[0]), len(group)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--all-sulfur-hits", type=Path, required=True)
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--candidate-metadata", type=Path, required=True)
    parser.add_argument("--tip-discordance", type=Path, required=True)
    parser.add_argument("--bac-taxonomy", type=Path, required=True)
    parser.add_argument("--ar-taxonomy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    spec_bytes = args.spec.read_bytes()
    spec = json.loads(spec_bytes)
    spec_sha256 = hashlib.sha256(spec_bytes).hexdigest()

    features = pd.read_csv(args.features, sep="\t", dtype=str).fillna("")
    all_hits = pd.read_csv(args.all_sulfur_hits, sep="\t", dtype=str).fillna("")
    if "assembly_accession" not in all_hits.columns and "accession" in all_hits.columns:
        all_hits = all_hits.rename(columns={"accession": "assembly_accession"})
    required_hit_columns = {"assembly_accession", "gene"}
    missing_hit_columns = required_hit_columns - set(all_hits.columns)
    if missing_hit_columns:
        raise ValueError(f"all-sulfur-hits missing columns: {sorted(missing_hit_columns)}")
    for flag in ["is_pseudogene", "is_partial"]:
        if flag in all_hits.columns:
            all_hits = all_hits.loc[
                ~all_hits[flag].str.strip().str.lower().isin({"1", "true", "yes", "y"})
            ].copy()
    universe = pd.read_csv(args.universe, sep="\t", dtype=str).fillna("")
    candidates = pd.read_csv(args.candidate_metadata, sep="\t", dtype=str).fillna("")
    discordance = pd.read_csv(args.tip_discordance, sep="\t", dtype=str).fillna("")
    taxonomy = pd.concat(
        [read_taxonomy(args.bac_taxonomy), read_taxonomy(args.ar_taxonomy)], ignore_index=True
    ).drop_duplicates("assembly_accession")

    universe_tax = universe[["assembly_accession"]].merge(taxonomy, on="assembly_accession", how="left")
    order_denominator = universe_tax.groupby("gtdb_order")["assembly_accession"].nunique().to_dict()
    hit_tax = all_hits.merge(taxonomy, on="assembly_accession", how="left")
    hit_tax["family"] = hit_tax["gene"].str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
    order_presence = (
        hit_tax.groupby(["family", "gtdb_order"])["assembly_accession"].nunique().to_dict()
    )

    candidates = candidates[["assembly_accession", "gtdb_class", "gtdb_order"]].drop_duplicates()
    features = features.merge(candidates, on="assembly_accession", how="left")
    discordance["order_discordance"] = pd.to_numeric(discordance["order_discordance"], errors="coerce")
    discordance_map = discordance.set_index(["family", "assembly_accession"])["order_discordance"].to_dict()

    output: list[dict[str, object]] = []
    for _, feature in features.iterrows():
        accession = feature["assembly_accession"]
        family = str(feature["gene"]).lower()
        role = feature["organizational_role"]
        order = feature.get("gtdb_order", "")
        denominator = order_denominator.get(order, 0)
        presence_key = (family, order)
        prevalence = (
            order_presence[presence_key] / denominator
            if denominator and presence_key in order_presence
            else float("nan")
        )
        phylo = discordance_map.get((family, accession), float("nan"))
        gene_composition = pd.to_numeric(feature.get("gene_gc_deviation", ""), errors="coerce")
        local_composition = pd.to_numeric(feature.get("local_gc_deviation", ""), errors="coerce")
        mobility_distance = pd.to_numeric(feature.get("nearest_mobile_distance_bp", ""), errors="coerce")
        cluster_id = str(feature.get("sox_cluster_id", ""))
        if role == "recipient_sulfur_backbone":
            operon_dependency = f"{accession}:recipient_sulfur_backbone"
        elif cluster_id:
            operon_dependency = f"{cluster_id}:{role}"
        else:
            operon_dependency = f"{accession}:{role}:unclustered"
        locus_dependency = f"{accession}:{feature.get('locus_tag','')}"
        homolog_dependency = f"{role}:{family}"
        is_mobile_element = str(feature.get("is_mobile_element", "")).strip().lower() in {
            "1", "true", "yes", "y"
        }

        layers = [
            ("GTDB_R232", "gene_tree_discordance", phylo, donor_from_phylogeny(phylo, spec)),
            ("GTDB_R232", "order_distribution", prevalence, donor_from_distribution(prevalence, spec)),
            (
                "RefSeq_PGAP",
                "gene_composition",
                gene_composition,
                donor_from_composition(gene_composition, spec),
            ),
            (
                "RefSeq_PGAP",
                "local_composition",
                local_composition,
                donor_from_composition(local_composition, spec),
            ),
            (
                "RefSeq_PGAP",
                "mobile_proximity",
                mobility_distance,
                donor_from_mobility(is_mobile_element, mobility_distance, spec),
            ),
        ]
        for source, layer, value, (donor, confidence) in layers:
            output.append(
                {
                    "case_id": "sox_genome_atlas",
                    "unit_id": f"ATLAS-{len(output)+1:07d}",
                    "source_id": source,
                    "source_type": "public_genome_database",
                    "evidence_layer": layer,
                    "evidence_unit": "genome_feature_metric",
                    "module_id": feature.get("sox_cluster_id", ""),
                    "gene_symbol": family,
                    "taxon_or_genome": accession,
                    "donor_class": donor,
                    "functional_class": role,
                    "route_eligible": donor != "uncertain",
                    "dependency_group": operon_dependency,
                    "genome_dependency_group": f"{accession}:{role}",
                    "operon_dependency_group": operon_dependency,
                    "locus_dependency_group": locus_dependency,
                    "homolog_dependency_group": homolog_dependency,
                    "confidence": confidence,
                    "metric_value": value,
                    "gtdb_class": feature.get("gtdb_class", ""),
                    "gtdb_order": order,
                    "provenance": f"{source}:{layer}:{accession}:{feature.get('locus_tag','')}",
                    "role_basis": feature.get("role_basis", ""),
                    "is_mobile_element": is_mobile_element,
                    "notes": "Donor evidence and organizational role are assigned from separate feature sets.",
                    "analysis_spec_sha256": spec_sha256,
                }
            )

    feature_metrics = pd.DataFrame(output)
    unit_rows: list[dict[str, object]] = []
    group_columns = [
        "source_id",
        "evidence_layer",
        "operon_dependency_group",
        "functional_class",
    ]
    for _, group in feature_metrics.groupby(group_columns, dropna=False, sort=True):
        donor, confidence, supporting_features, total_features = aggregate_donor(group, spec)
        metric_values = pd.to_numeric(group["metric_value"], errors="coerce").dropna()
        families = sorted(set(filter(None, group["gene_symbol"].astype(str))))
        loci = sorted(set(filter(None, group["locus_dependency_group"].astype(str))))
        unit_rows.append(
            {
                "case_id": "sox_genome_atlas",
                "unit_id": f"ATLAS-UNIT-{len(unit_rows)+1:07d}",
                "source_id": group["source_id"].iloc[0],
                "source_type": "public_genome_database",
                "evidence_layer": group["evidence_layer"].iloc[0],
                "evidence_unit": "operon_role_summary",
                "module_id": group["module_id"].replace("", pd.NA).dropna().iloc[0]
                if group["module_id"].replace("", pd.NA).notna().any()
                else "",
                "gene_symbol": ";".join(families),
                "taxon_or_genome": group["taxon_or_genome"].iloc[0],
                "donor_class": donor,
                "functional_class": group["functional_class"].iloc[0],
                "route_eligible": donor != "uncertain",
                "dependency_group": group["operon_dependency_group"].iloc[0],
                "genome_dependency_group": group["genome_dependency_group"].iloc[0],
                "operon_dependency_group": group["operon_dependency_group"].iloc[0],
                "homolog_dependency_group": ";".join(
                    sorted(set(group["homolog_dependency_group"].astype(str)))
                ),
                "confidence": confidence,
                "metric_value": float(metric_values.median()) if not metric_values.empty else float("nan"),
                "supporting_feature_metrics": supporting_features,
                "total_feature_metrics": total_features,
                "n_loci": len(loci),
                "gtdb_class": group["gtdb_class"].iloc[0],
                "gtdb_order": group["gtdb_order"].iloc[0],
                "provenance": ";".join(sorted(set(group["provenance"].astype(str)))),
                "role_basis": ";".join(sorted(set(group["role_basis"].astype(str)))),
                "contains_mobile_element": bool(group["is_mobile_element"].astype(bool).any()),
                "notes": "Operon-level donor call aggregated independently within a frozen evidence modality.",
                "analysis_spec_sha256": spec_sha256,
            }
        )

    evidence = pd.DataFrame(unit_rows)
    strict = evidence.loc[evidence["route_eligible"]].copy()
    matrix = pd.crosstab(strict["donor_class"], strict["functional_class"])
    layer_matrix = (
        strict.groupby(["source_id", "evidence_layer", "donor_class", "functional_class"])
        .size()
        .rename("rows")
        .reset_index()
    )
    summary = pd.DataFrame(
        [
            {"metric": "all_feature_metric_rows", "value": len(feature_metrics)},
            {"metric": "all_operon_role_units", "value": len(evidence)},
            {"metric": "route_eligible_operon_role_units", "value": len(strict)},
            {"metric": "dependency_groups", "value": strict["dependency_group"].nunique()},
            {"metric": "genomes", "value": strict["taxon_or_genome"].nunique()},
            {"metric": "sources", "value": strict["source_id"].nunique()},
            {"metric": "layers", "value": strict["evidence_layer"].nunique()},
            {"metric": "analysis_spec_sha256", "value": spec_sha256},
        ]
    )

    feature_metrics.to_csv(args.outdir / "sox_atlas_feature_metrics.tsv", sep="\t", index=False)
    feature_metrics.to_csv(args.outdir / "sox_atlas_feature_metrics.csv", index=False)
    evidence.to_csv(args.outdir / "sox_atlas_evidence_units.tsv", sep="\t", index=False)
    evidence.to_csv(args.outdir / "sox_atlas_evidence_units.csv", index=False)
    matrix.to_csv(args.outdir / "sox_atlas_donor_role_matrix.tsv", sep="\t")
    layer_matrix.to_csv(args.outdir / "sox_atlas_source_layer_matrix.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "sox_atlas_evidence_summary.tsv", sep="\t", index=False)
    (args.outdir / "analysis_spec_used.json").write_bytes(spec_bytes)
    print(summary.to_string(index=False))
    print("\nStrict donor-role matrix")
    print(matrix.to_string())


if __name__ == "__main__":
    main()
