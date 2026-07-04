from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
V4 = ROOT / "16_iMETA_ROUTE_B_DRAFT_v4_DUAL_CASE_SOX_TRANSFER"
SRC = V4 / "04_Source_Data_DUAL_CASE"
BASE = SRC / "Source_Data_iMeta_RouteB_v4_base.xlsx"
SOX = ROOT / "sox_v2_results"
SOX_STRESS = ROOT / "sox_annotation_stress_results"
SOX_ANN = ROOT / "sox_annotation_crosscheck_v1"
OUT = SRC / "Source_Data_iMeta_RouteB_v4_consolidated.xlsx"


def safe_sheet(name: str) -> str:
    bad = set(r"[]:*?/\\")
    cleaned = "".join("_" if c in bad else c for c in name)
    return cleaned[:31]


def write_df(writer: pd.ExcelWriter, name: str, df: pd.DataFrame) -> None:
    df.to_excel(writer, sheet_name=safe_sheet(name), index=False)


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def main() -> None:
    if not BASE.exists():
        raise FileNotFoundError(BASE)

    wb = load_workbook(BASE)
    wb.save(OUT)

    index_rows = [
        {
            "sheet": "Index_v4",
            "description": "Top-level index for the consolidated v4 Source Data workbook.",
        },
        {
            "sheet": "SOX_evidence_v2",
            "description": "Machine-readable SOX transfer-case evidence-unit table.",
        },
        {
            "sheet": "SOX_summary",
            "description": "SOX route-graph summary for the original v2 table.",
        },
        {
            "sheet": "SOX_routes",
            "description": "All donor-role route scores for the original SOX v2 table.",
        },
        {
            "sheet": "SOX_leave_source",
            "description": "Leave-one-source sensitivity for the original SOX v2 table.",
        },
        {
            "sheet": "SOX_perm_nulls",
            "description": "Source-only and source-by-layer permutation nulls for the original SOX v2 table.",
        },
        {
            "sheet": "SOX_ann_summary",
            "description": "Annotation cross-check summary.",
        },
        {
            "sheet": "SOX_ann_flagged",
            "description": "Rows flagged by the SOX annotation cross-check.",
        },
        {
            "sheet": "SOX_stress_summary",
            "description": "Route-graph summary after excluding low-confidence or unsupported SOX rows from strict scoring.",
        },
        {
            "sheet": "SOX_stress_leave",
            "description": "Leave-one-source sensitivity after annotation-stress downgrade.",
        },
        {
            "sheet": "SOX_stress_perm",
            "description": "Permutation nulls after annotation-stress downgrade.",
        },
    ]

    with pd.ExcelWriter(OUT, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        write_df(writer, "Index_v4", pd.DataFrame(index_rows))
        write_df(writer, "SOX_evidence_v2", read_csv(SRC / "SOX_evidence_units_expanded_v2.csv"))
        write_df(writer, "SOX_summary", read_csv(SOX / "summary.csv"))
        write_df(writer, "SOX_routes", read_csv(SOX / "route_scores.csv"))
        write_df(writer, "SOX_leave_source", read_csv(SOX / "leave_one_source.csv"))
        write_df(writer, "SOX_perm_nulls", read_csv(SOX / "permutation_nulls.csv"))
        write_df(writer, "SOX_ann_summary", read_csv(SOX_ANN / "sox_annotation_crosscheck_summary_v1.csv"))
        write_df(writer, "SOX_ann_flagged", read_csv(SOX_ANN / "sox_annotation_crosscheck_flagged_rows_v1.csv"))
        write_df(writer, "SOX_stress_summary", read_csv(SOX_STRESS / "summary.csv"))
        write_df(writer, "SOX_stress_leave", read_csv(SOX_STRESS / "leave_one_source.csv"))
        write_df(writer, "SOX_stress_perm", read_csv(SOX_STRESS / "permutation_nulls.csv"))

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
