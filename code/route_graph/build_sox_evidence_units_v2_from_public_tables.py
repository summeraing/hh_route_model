from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
CASE = ROOT / "15_iMETA_ROUTE_B_DUAL_CASE_TRANSFER_20260701" / "01_HGT_SOX_CASE"
TMP = ROOT / "sox_tmp_tables"
OUT = CASE / "curated_seed" / "sox_evidence_units_expanded_v2.csv"


COLUMNS = [
    "case_id",
    "unit_id",
    "source_id",
    "source_type",
    "evidence_layer",
    "evidence_unit",
    "module_id",
    "gene_symbol",
    "taxon_or_genome",
    "donor_class",
    "functional_class",
    "route_eligible",
    "dependency_group",
    "confidence",
    "evidence_text",
    "provenance",
    "notes",
]


def norm_presence(x: object) -> bool:
    s = str(x).strip()
    return s.startswith("+")


def safe_id(text: object) -> str:
    s = str(text).strip()
    for ch in [" ", ".", "(", ")", "/", "\\", ":", ";", ",", "'", '"', "−", "+"]:
        s = s.replace(ch, "_")
    while "__" in s:
        s = s.replace("__", "_")
    return s.strip("_")[:80]


def row(
    unit_id: str,
    source_id: str,
    source_type: str,
    evidence_layer: str,
    evidence_unit: str,
    module_id: str,
    gene_symbol: str,
    taxon_or_genome: str,
    donor_class: str,
    functional_class: str,
    route_eligible: bool,
    dependency_group: str,
    confidence: str,
    evidence_text: str,
    provenance: str,
    notes: str,
) -> dict[str, object]:
    return {
        "case_id": "sox_hgt_transfer",
        "unit_id": unit_id,
        "source_id": source_id,
        "source_type": source_type,
        "evidence_layer": evidence_layer,
        "evidence_unit": evidence_unit,
        "module_id": module_id,
        "gene_symbol": gene_symbol,
        "taxon_or_genome": taxon_or_genome,
        "donor_class": donor_class,
        "functional_class": functional_class,
        "route_eligible": bool(route_eligible),
        "dependency_group": dependency_group,
        "confidence": confidence,
        "evidence_text": evidence_text,
        "provenance": provenance,
        "notes": notes,
    }


def build_gregersen_rows() -> list[dict[str, object]]:
    table = pd.read_csv(TMP / "gregersen2011_table_0.csv")
    gene_names = table.iloc[0, 5:].tolist()
    data = table.iloc[1:].copy()

    # Gregersen et al. state that thiosulfate utilization by SOX in GSB was
    # acquired horizontally from Proteobacteria and that the conserved
    # soxJXYZAKBW cluster was also exchanged within GSB lineages. Other sulfur
    # systems are retained as support unless the text explicitly supports a
    # route assignment.
    strict_hgt_genes = {"sox"}
    support_genes = set(str(g) for g in gene_names) - strict_hgt_genes

    rows: list[dict[str, object]] = []
    n = 1
    for _, r in data.iterrows():
        strain = str(r["Strain"]).strip()
        acc = str(r["Accession Number"]).strip()
        taxon = strain if acc in {"nan", "n.a.", "n.a"} else f"{strain} [{acc}]"
        strain_id = safe_id(strain)
        for col, gene in zip(table.columns[5:], gene_names):
            gene = str(gene).strip()
            if not gene or gene == "nan":
                continue
            if not norm_presence(r[col]):
                continue
            if gene in strict_hgt_genes:
                rows.append(
                    row(
                        unit_id=f"SOX-V2-GREG-{n:04d}",
                        source_id="GREGERSEN2011",
                        source_type="genomic_synthesis",
                        evidence_layer="table1_gene_presence",
                        evidence_unit="strain_gene_presence",
                        module_id="GSB_sox_cluster",
                        gene_symbol=gene,
                        taxon_or_genome=taxon,
                        donor_class="hgt_pathway",
                        functional_class="variable_energy_module",
                        route_eligible=True,
                        dependency_group=f"GREGERSEN2011:Table1:{strain_id}:{gene}",
                        confidence="medium",
                        evidence_text="Table 1 reports SOX presence in this strain; text states SOX thiosulfate utilization in GSB was horizontally acquired from Proteobacteria.",
                        provenance="Gregersen2011 Table 1; doi:10.3389/fmicb.2011.00116",
                        notes="Auto-expanded from PMC table; verify against article table before final submission.",
                    )
                )
                n += 1
            elif gene in support_genes:
                rows.append(
                    row(
                        unit_id=f"SOX-V2-GREG-{n:04d}",
                        source_id="GREGERSEN2011",
                        source_type="genomic_synthesis",
                        evidence_layer="table1_gene_presence_support",
                        evidence_unit="strain_gene_presence",
                        module_id="GSB_sulfur_metabolism_support",
                        gene_symbol=gene,
                        taxon_or_genome=taxon,
                        donor_class="uncertain",
                        functional_class="support",
                        route_eligible=False,
                        dependency_group=f"GREGERSEN2011:Table1:{strain_id}:{gene}",
                        confidence="medium",
                        evidence_text="Table 1 reports sulfur-metabolism gene presence; donor-route label not assigned in the automated expansion.",
                        provenance="Gregersen2011 Table 1; doi:10.3389/fmicb.2011.00116",
                        notes="Support row retained for transparency; not counted in strict route support.",
                    )
                )
                n += 1

    # Explicit mobile-element context from the source text around Chlorobium
    # phaeovibrioides acquisition. Kept as source-level rows because the HTML
    # table does not expose the exact mobility genes.
    mobile_items = ["mobile_genetic_element", "sox_cluster_boundary", "within_GSB_exchange"]
    for item in mobile_items:
        rows.append(
            row(
                unit_id=f"SOX-V2-GREG-{n:04d}",
                source_id="GREGERSEN2011",
                source_type="genomic_synthesis",
                evidence_layer="text_mobile_context",
                evidence_unit="source_text_event",
                module_id="GSB_sox_cluster_mobility",
                gene_symbol=item,
                taxon_or_genome="Chlorobium_phaeovibrioides_DSM_265",
                donor_class="mobile_context",
                functional_class="mobility_boundary",
                route_eligible=True,
                dependency_group=f"GREGERSEN2011:text_mobile_context:{item}",
                confidence="medium",
                evidence_text="Source text links conserved soxJXYZAKBW acquisition to a mobile genetic element and within-GSB exchange.",
                provenance="Gregersen2011 text; doi:10.3389/fmicb.2011.00116",
                notes="Text-derived event row; final curation should confirm exact wording.",
            )
        )
        n += 1
    return rows


def build_berben_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    n = 1
    # Direct statements from Berben et al. text: fccAB, soxAX and soxB are
    # found in all 75 genomes; truncated sox and soe are present in all strains;
    # hdr is detected in 73 genomes; dsr occurs in only six strains.
    core_genes = [
        ("fccAB", "fccA", "all_75_genomes"),
        ("fccAB", "fccB", "all_75_genomes"),
        ("soxAX", "soxA", "all_75_genomes"),
        ("soxAX", "soxX", "all_75_genomes"),
        ("soxB", "soxB", "all_75_genomes"),
        ("truncated_sox", "soxA", "all_75_genomes"),
        ("truncated_sox", "soxB", "all_75_genomes"),
        ("truncated_sox", "soxX", "all_75_genomes"),
        ("truncated_sox", "soxY", "all_75_genomes"),
        ("truncated_sox", "soxZ", "all_75_genomes"),
        ("soeABC", "soeA", "all_75_genomes"),
        ("soeABC", "soeB", "all_75_genomes"),
        ("soeABC", "soeC", "all_75_genomes"),
        ("hdr", "hdrC1", "73_genomes"),
        ("hdr", "hdrB1", "73_genomes"),
        ("hdr", "hdrA", "73_genomes"),
        ("hdr", "hyp", "73_genomes"),
        ("hdr", "hdrC2", "73_genomes"),
        ("hdr", "hdrB2", "73_genomes"),
        ("hdr", "lbpA1", "73_genomes"),
        ("hdr", "dsrE", "73_genomes"),
        ("hdr", "lbpA2", "73_genomes"),
    ]
    for module, gene, scope in core_genes:
        rows.append(
            row(
                unit_id=f"SOX-V2-BERB-{n:04d}",
                source_id="BERBEN2019",
                source_type="comparative_genomics",
                evidence_layer="genome_distribution_text",
                evidence_unit="gene_presence_summary",
                module_id=module,
                gene_symbol=gene,
                taxon_or_genome=f"Thioalkalivibrio_{scope}",
                donor_class="lineage_core",
                functional_class="conserved_metabolic_backbone",
                route_eligible=True,
                dependency_group=f"BERBEN2019:{module}:{gene}:{scope}",
                confidence="medium",
                evidence_text="Berben et al. report conserved sulfur-oxidation genes across 75 Thioalkalivibrio genomes or 73-genome hdr distribution.",
                provenance="Berben2019 text/table; doi:10.3389/fmicb.2019.00160",
                notes="Auto-expanded from text statements and Table 1 module names; verify exact gene list before final submission.",
            )
        )
        n += 1

    variable_events = [
        ("second_soxB_copy_species_1_15", "soxB_copy2", "species_1_to_15_clade"),
        ("second_soxB_copy_ALJD", "soxB_copy2", "ALJD_separate_event"),
        ("dsrABC_rare_presence", "dsrABC", "six_genomes"),
    ]
    for module, gene, scope in variable_events:
        rows.append(
            row(
                unit_id=f"SOX-V2-BERB-{n:04d}",
                source_id="BERBEN2019",
                source_type="comparative_genomics",
                evidence_layer="copy_number_or_variable_distribution",
                evidence_unit="event_summary",
                module_id=module,
                gene_symbol=gene,
                taxon_or_genome=f"Thioalkalivibrio_{scope}",
                donor_class="hgt_pathway",
                functional_class="variable_energy_module",
                route_eligible=True,
                dependency_group=f"BERBEN2019:{module}:{gene}:{scope}",
                confidence="low" if "dsr" in gene else "medium",
                evidence_text="Berben et al. discuss second soxB copies as possible duplication/HGT events and report rare dsr distribution.",
                provenance="Berben2019 text; doi:10.3389/fmicb.2019.00160",
                notes="Variable event row; dsr row should be downgraded to support if final curation does not justify HGT/pathway donor label.",
            )
        )
        n += 1
    return rows


def main() -> None:
    seed = pd.read_csv(CASE / "curated_seed" / "sox_evidence_units_seed_v1.csv")
    rows = seed.to_dict("records")
    rows.extend(build_gregersen_rows())
    rows.extend(build_berben_rows())
    df = pd.DataFrame(rows, columns=COLUMNS)
    df = df.drop_duplicates(subset=["unit_id"], keep="first")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, encoding="utf-8-sig")
    route = df[df["route_eligible"].astype(str).str.lower().isin(["true", "1", "yes"])]
    print(f"Wrote {OUT}")
    print(f"rows={len(df)} route_eligible={len(route)} sources={route['source_id'].nunique()} dependency_groups={route['dependency_group'].nunique()}")
    print(route["source_id"].value_counts().to_string())


if __name__ == "__main__":
    main()

