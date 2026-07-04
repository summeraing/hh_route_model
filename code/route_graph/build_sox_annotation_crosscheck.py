from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
CASE = ROOT / "15_iMETA_ROUTE_B_DUAL_CASE_TRANSFER_20260701" / "01_HGT_SOX_CASE"
V4 = ROOT / "16_iMETA_ROUTE_B_DRAFT_v4_DUAL_CASE_SOX_TRANSFER"
SOURCE_OUT = V4 / "04_Source_Data_DUAL_CASE"

INPUT = CASE / "curated_seed" / "sox_evidence_units_expanded_v2.csv"
# Keep generated files under a short path. Pandas/openpyxl can still hit
# Windows MAX_PATH limits under the full manuscript directory.
OUTDIR = ROOT / "sox_annotation_crosscheck_v1"
STRESS = OUTDIR / "sox_evidence_units_expanded_v2_annotation_stress.csv"
WORKBOOK = OUTDIR / "SOX_annotation_crosscheck_v1.xlsx"


METABOLIC_PREFIXES = (
    "sox",
    "dsr",
    "hdr",
    "soe",
    "fcc",
    "sqr",
    "apr",
    "sat",
)

MOBILE_KEYWORDS = (
    "mobile",
    "transposase",
    "integrase",
    "recombinase",
    "resolvase",
    "plasmid",
    "insertion",
    "island",
    "boundary",
    "border",
    "repeat",
    "reverse_transcriptase",
)

HGT_KEYWORDS = (
    "hgt",
    "horizontal",
    "horizontally",
    "acquired",
    "exchanged",
    "gain",
    "loss",
    "gain/loss",
    "variable",
    "copy",
    "discordant",
    "phylogeny",
    "phylogen",
    "distribution implications",
    "rare",
)

CORE_KEYWORDS = (
    "all_75",
    "73_genomes",
    "conserved",
    "present across",
    "across analysed",
    "across analyzed",
    "most analysed",
    "most analyzed",
    "stable repertoire",
    "genome_distribution",
    "genome cluster",
    "reported across",
)


def clean_text(*parts: object) -> str:
    return " ".join("" if pd.isna(p) else str(p) for p in parts).lower()


def starts_with_metabolic_gene(token: str) -> bool:
    token = re.sub(r"[^a-z0-9_]+", "", token.lower())
    return any(token.startswith(prefix) for prefix in METABOLIC_PREFIXES)


def has_any(text: str, keywords: tuple[str, ...]) -> bool:
    return any(k in text for k in keywords)


def annotation_axis(row: pd.Series) -> str:
    text = clean_text(
        row.get("gene_symbol"),
        row.get("module_id"),
        row.get("evidence_layer"),
        row.get("evidence_text"),
        row.get("notes"),
    )
    gene = clean_text(row.get("gene_symbol"))
    module = clean_text(row.get("module_id"))
    if has_any(text, MOBILE_KEYWORDS):
        return "mobility_context_signal"
    if starts_with_metabolic_gene(gene) or starts_with_metabolic_gene(module):
        return "sulfur_metabolism_signal"
    if "support" in clean_text(row.get("functional_class")):
        return "support_only"
    return "unresolved_annotation"


def function_supported(row: pd.Series, axis: str) -> bool:
    functional_class = clean_text(row.get("functional_class"))
    if functional_class in {"support", "uncertain"}:
        return True
    if functional_class == "mobility_boundary":
        return axis == "mobility_context_signal"
    if functional_class in {"conserved_metabolic_backbone", "variable_energy_module"}:
        return axis == "sulfur_metabolism_signal"
    return False


def donor_supported(row: pd.Series, axis: str) -> bool:
    donor_class = clean_text(row.get("donor_class"))
    text = clean_text(
        row.get("evidence_layer"),
        row.get("evidence_text"),
        row.get("taxon_or_genome"),
        row.get("module_id"),
        row.get("notes"),
    )
    if donor_class in {"uncertain", "support"}:
        return True
    if donor_class == "mobile_context":
        return axis == "mobility_context_signal" or has_any(text, MOBILE_KEYWORDS)
    if donor_class == "hgt_pathway":
        return has_any(text, HGT_KEYWORDS) or row.get("evidence_layer") in {
            "phylogeny",
            "gain_loss",
            "table1_gene_presence",
            "copy_number_or_variable_distribution",
            "gene_distribution",
        }
    if donor_class == "lineage_core":
        return has_any(text, CORE_KEYWORDS) or row.get("evidence_layer") in {
            "genome_distribution_text",
            "genome_cluster",
        }
    return False


def confidence_rank(value: object) -> int:
    s = clean_text(value)
    return {"low": 0, "medium": 1, "high": 2}.get(s, 1)


def add_crosscheck(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    axes = []
    function_ok = []
    donor_ok = []
    retained = []
    flags = []

    for _, row in out.iterrows():
        axis = annotation_axis(row)
        f_ok = function_supported(row, axis)
        d_ok = donor_supported(row, axis)
        is_route = str(row.get("route_eligible")).lower() in {"true", "1", "yes"}
        low_conf = confidence_rank(row.get("confidence")) == 0
        retain = bool(is_route and f_ok and d_ok and not low_conf)

        row_flags: list[str] = []
        if is_route and not f_ok:
            row_flags.append("functional_label_not_supported_by_gene_module_annotation")
        if is_route and not d_ok:
            row_flags.append("donor_label_not_supported_by_text_or_layer_signal")
        if is_route and low_conf:
            row_flags.append("low_confidence_route_row")
        if not row_flags:
            row_flags.append("none")

        axes.append(axis)
        function_ok.append(f_ok)
        donor_ok.append(d_ok)
        retained.append(retain)
        flags.append(";".join(row_flags))

    out["annotation_axis"] = axes
    out["functional_label_supported_by_annotation"] = function_ok
    out["donor_label_supported_by_layer_text"] = donor_ok
    out["annotation_stress_route_eligible"] = retained
    out["annotation_crosscheck_flag"] = flags
    return out


def build_stress_table(checked: pd.DataFrame) -> pd.DataFrame:
    stress = checked.drop(
        columns=[
            "annotation_axis",
            "functional_label_supported_by_annotation",
            "donor_label_supported_by_layer_text",
            "annotation_stress_route_eligible",
            "annotation_crosscheck_flag",
        ]
    ).copy()
    stress["route_eligible"] = checked["annotation_stress_route_eligible"].astype(bool)
    stress.loc[~stress["route_eligible"], "notes"] = (
        stress.loc[~stress["route_eligible"], "notes"].fillna("")
        + " | annotation-stress rerun excludes this row from strict route scoring."
    ).str.strip()
    return stress


def summarize(checked: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    is_route = checked["route_eligible"].astype(str).str.lower().isin(["true", "1", "yes"])
    summary = pd.DataFrame(
        [
            {
                "metric": "total_rows",
                "value": len(checked),
            },
            {
                "metric": "strict_route_rows_original",
                "value": int(is_route.sum()),
            },
            {
                "metric": "route_rows_with_supported_function_label",
                "value": int((is_route & checked["functional_label_supported_by_annotation"]).sum()),
            },
            {
                "metric": "route_rows_with_supported_donor_label",
                "value": int((is_route & checked["donor_label_supported_by_layer_text"]).sum()),
            },
            {
                "metric": "route_rows_retained_after_annotation_stress",
                "value": int(checked["annotation_stress_route_eligible"].sum()),
            },
            {
                "metric": "route_rows_flagged_for_review",
                "value": int((is_route & (checked["annotation_crosscheck_flag"] != "none")).sum()),
            },
        ]
    )

    by_class = (
        checked[is_route]
        .groupby(["donor_class", "functional_class", "annotation_axis"], dropna=False)
        .agg(
            rows=("unit_id", "count"),
            retained_after_stress=("annotation_stress_route_eligible", "sum"),
            function_supported=("functional_label_supported_by_annotation", "sum"),
            donor_supported=("donor_label_supported_by_layer_text", "sum"),
        )
        .reset_index()
    )

    by_source = (
        checked[is_route]
        .groupby(["source_id"], dropna=False)
        .agg(
            route_rows=("unit_id", "count"),
            retained_after_stress=("annotation_stress_route_eligible", "sum"),
            flagged=("annotation_crosscheck_flag", lambda x: int((x != "none").sum())),
        )
        .reset_index()
    )
    return summary, by_class, by_source


def main() -> None:
    df = pd.read_csv(INPUT)
    checked = add_crosscheck(df)
    stress = build_stress_table(checked)
    summary, by_class, by_source = summarize(checked)
    flagged = checked[checked["annotation_crosscheck_flag"] != "none"].copy()

    rules = pd.DataFrame(
        [
            {"rule_group": "metabolic_prefixes", "values": ", ".join(METABOLIC_PREFIXES)},
            {"rule_group": "mobile_keywords", "values": ", ".join(MOBILE_KEYWORDS)},
            {"rule_group": "hgt_keywords", "values": ", ".join(HGT_KEYWORDS)},
            {"rule_group": "core_keywords", "values": ", ".join(CORE_KEYWORDS)},
            {
                "rule_group": "interpretation",
                "values": "Rule-based cross-check is a transparent annotation stress test, not a replacement for manual source curation.",
            },
        ]
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    checked.to_csv(OUTDIR / "sox_annotation_crosscheck_rows_v1.csv", index=False, encoding="utf-8-sig")
    stress.to_csv(STRESS, index=False, encoding="utf-8-sig")
    summary.to_csv(OUTDIR / "sox_annotation_crosscheck_summary_v1.csv", index=False, encoding="utf-8-sig")
    flagged.to_csv(OUTDIR / "sox_annotation_crosscheck_flagged_rows_v1.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(WORKBOOK, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="summary", index=False)
        by_class.to_excel(xw, sheet_name="by_route_class", index=False)
        by_source.to_excel(xw, sheet_name="by_source", index=False)
        flagged.to_excel(xw, sheet_name="flagged_rows", index=False)
        rules.to_excel(xw, sheet_name="rules", index=False)
        checked.to_excel(xw, sheet_name="row_crosscheck", index=False)

    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    copy_targets = {
        WORKBOOK: "SOX_annotation_crosscheck.xlsx",
        STRESS: "SOX_annotation_stress.csv",
        OUTDIR / "sox_annotation_crosscheck_summary_v1.csv": "SOX_annotation_summary.csv",
        OUTDIR / "sox_annotation_crosscheck_flagged_rows_v1.csv": "SOX_annotation_flagged_rows.csv",
    }
    for src, dst_name in copy_targets.items():
        dst = SOURCE_OUT / dst_name
        dst.write_bytes(src.read_bytes())

    print(f"Wrote {WORKBOOK}")
    print(f"Wrote {STRESS}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
