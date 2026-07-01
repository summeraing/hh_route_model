from __future__ import annotations

from itertools import permutations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
FULL_TABLE = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "v77_clean" / "01_MAIN_FIGURES" / "Fig1" / "source_package" / "Fig1B_selected_original_network" / "source_combined_v4_flat_export.csv"
SOURCE_XLSX = ROOT / "10_WORKING_COPY_FROM_NC12" / "04_Source_Data" / "Source_Data_complete.xlsx"
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_row_level_adversarial_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
ROLE_LABELS = {"scaffold": "scaffold", "energetic_incorporation": "energy", "transition": "transition"}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}

COL = {
    "host": "#356B8A",
    "symbiont": "#E0942A",
    "other": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
}

SOURCE_LABELS = {
    "KU2015": "Ku 2015",
    "SANTANA2025": "Santana 2025",
    "ZHANG2025": "Zhang 2025",
    "TOB2026": "Tobiasson 2026",
    "EME2023": "Eme 2023",
}


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def norm(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = s.replace("injection", "energetic_incorporation")
    s = s.replace("energy", "energetic_incorporation")
    return s


def load_full_with_audit_overlay():
    df = pd.read_csv(FULL_TABLE)
    df["row_index"] = np.arange(len(df))
    df["donor_original"] = df["donor_class"].map(norm)
    df["functional_original"] = df["functional_class"].map(norm)
    df["donor_final"] = df["donor_original"]
    df["functional_final"] = df["functional_original"]
    df["audit_overlay_applied"] = False

    audit = pd.read_excel(SOURCE_XLSX, sheet_name="audit_overlay")
    for _, row in audit.iterrows():
        idx = int(row["combined_index"])
        if 0 <= idx < len(df):
            df.loc[idx, "donor_final"] = norm(row["adjudicated_donor_class"])
            df.loc[idx, "functional_final"] = norm(row["adjudicated_functional_class"])
            df.loc[idx, "audit_overlay_applied"] = True
    return df


def source_equal_support(df, mapping):
    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    per_source = []
    for _, g in work.groupby("study_id", sort=False):
        ok = [mapping.get(d) == f for d, f in zip(g["donor_final"], g["functional_final"])]
        if ok:
            per_source.append(float(np.mean(ok)))
    return float(np.mean(per_source)), work


def route_table(df):
    rows = []
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        support, work = source_equal_support(df, mapping)
        raw = float(np.mean([mapping.get(d) == f for d, f in zip(work["donor_final"], work["functional_final"])]))
        rows.append({
            "route_map": "; ".join(f"{d}->{ROLE_LABELS[mapping[d]]}" for d in DONORS),
            "host_route": mapping["host"],
            "symbiont_route": mapping["symbiont"],
            "other_route": mapping["other"],
            "is_prespecified": mapping == PRESPEC,
            "n_route_eligible_rows": len(work),
            "n_sources": work["study_id"].nunique(),
            "source_equal_support": support,
            "raw_support_descriptive": raw,
        })
    out = pd.DataFrame(rows).sort_values("source_equal_support", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def row_level_adversary(df, best_alt):
    support_pres, work = source_equal_support(df, PRESPEC)
    support_alt, _ = source_equal_support(df, best_alt)
    initial_margin = support_pres - support_alt
    n_sources = work["study_id"].nunique()
    source_sizes = work.groupby("study_id").size().to_dict()

    candidates = []
    for _, row in work.iterrows():
        d = row["donor_final"]
        f = row["functional_final"]
        n_s = source_sizes[row["study_id"]]
        old = int(PRESPEC.get(d) == f) - int(best_alt.get(d) == f)
        flipped_role = best_alt.get(d)
        new = int(PRESPEC.get(d) == flipped_role) - int(best_alt.get(d) == flipped_role)
        reduction = (old - new) / (n_sources * n_s)
        if reduction > 0:
            candidates.append({
                "row_index": row["row_index"],
                "study_id": row["study_id"],
                "evidence_layer": row["evidence_layer"],
                "evidence_unit": row["evidence_unit"],
                "unit_id": row["unit_id"],
                "unit_label": row["unit_label"],
                "donor_final": d,
                "functional_final": f,
                "counterfactual_functional_label": flipped_role,
                "old_margin_contribution": old / (n_sources * n_s),
                "new_margin_contribution": new / (n_sources * n_s),
                "source_equal_margin_reduction": reduction,
                "source_route_eligible_n": n_s,
                "audit_overlay_applied": bool(row["audit_overlay_applied"]),
            })
    cand = pd.DataFrame(candidates).sort_values(
        ["source_equal_margin_reduction", "study_id", "row_index"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    cand["cumulative_margin_reduction"] = cand["source_equal_margin_reduction"].cumsum()
    cand["margin_after_flip"] = initial_margin - cand["cumulative_margin_reduction"]
    needed = int((cand["margin_after_flip"] <= 0).idxmax() + 1) if (cand["margin_after_flip"] <= 0).any() else np.nan
    selected = cand.iloc[:needed].copy() if not pd.isna(needed) else cand.iloc[0:0].copy()

    by_source = selected.groupby("study_id").agg(
        targeted_flips=("row_index", "count"),
        cumulative_reduction=("source_equal_margin_reduction", "sum"),
    ).reset_index() if len(selected) else pd.DataFrame(columns=["study_id", "targeted_flips", "cumulative_reduction"])

    summary = pd.DataFrame([{
        "test": "row_level_source_equal_targeted_flip_after_audit_overlay",
        "full_table_rows": len(df),
        "route_eligible_rows_after_audit_overlay": len(work),
        "audit_overlay_rows_applied": int(df["audit_overlay_applied"].sum()),
        "n_sources": n_sources,
        "prespecified_source_equal_support": support_pres,
        "nearest_alternative_source_equal_support": support_alt,
        "initial_source_equal_margin": initial_margin,
        "nearest_alternative_map": "; ".join(f"{d}->{ROLE_LABELS[best_alt[d]]}" for d in DONORS),
        "minimum_targeted_row_flips_to_erase_margin": needed,
        "fraction_of_route_eligible_rows": needed / len(work) if not pd.isna(needed) else np.nan,
        "fraction_of_full_table_rows": needed / len(df) if not pd.isna(needed) else np.nan,
        "selected_flip_reduction": selected["source_equal_margin_reduction"].sum() if len(selected) else 0,
        "final_margin_after_selected_flips": selected["margin_after_flip"].iloc[-1] if len(selected) else initial_margin,
        "interpretation": "Rows are traceable but dependent. This is a targeted row-level adversarial influence calculation under source-equal weighting after overlaying adjudicated audit labels where available.",
    }])
    return summary, cand, selected, by_source, work


def make_figure(routes, summary, cand, by_source, work):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.8), constrained_layout=True)

    ax = axs[0, 0]
    plot = routes.sort_values("source_equal_support", ascending=True)
    colors = [COL["host"] if x else COL["light"] for x in plot["is_prespecified"]]
    ax.barh(np.arange(len(plot)), plot["source_equal_support"], color=colors, edgecolor="white")
    ax.set_yticks(np.arange(len(plot)), [f"R{int(r)}" for r in plot["rank"]])
    ax.set_xlabel("source-equal route support")
    ax.set_title("audit-overlaid row table", fontsize=9, weight="bold")
    for y, (_, r) in enumerate(plot.iterrows()):
        ax.text(r["source_equal_support"] + 0.01, y, f"{r['source_equal_support']:.2f}", va="center", fontsize=6.5)
    ax.text(-0.16, 1.06, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[0, 1]
    sel_n = int(summary["minimum_targeted_row_flips_to_erase_margin"].iloc[0])
    curve = cand.iloc[:max(sel_n + 30, min(len(cand), 200))].copy()
    ax.plot(np.arange(1, len(curve) + 1), curve["margin_after_flip"], color=COL["red"], lw=1.7)
    ax.axhline(0, color="#111827", lw=0.8)
    ax.axvline(sel_n, color=COL["red"], ls="--", lw=1.0)
    ax.set_xlabel("targeted row flips")
    ax.set_ylabel("margin after flips")
    ax.set_title("row-level adversarial erasure", fontsize=9, weight="bold")
    ax.text(sel_n, curve["margin_after_flip"].max() * 0.55, f"threshold\n{sel_n} rows",
            ha="left", va="center", fontsize=6.5, color=COL["red"])
    ax.text(-0.16, 1.06, "b", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 0]
    if len(by_source):
        order = by_source.sort_values("targeted_flips", ascending=True)
        labels = [SOURCE_LABELS.get(x, x) for x in order["study_id"]]
        ax.barh(labels, order["targeted_flips"], color=COL["grey"], edgecolor="white")
        for y, (_, r) in enumerate(order.iterrows()):
            ax.text(r["targeted_flips"] + 0.3, y, str(int(r["targeted_flips"])), va="center", fontsize=6.5)
    ax.set_xlabel("selected flips")
    ax.set_title("where the adversary acts", fontsize=9, weight="bold")
    ax.text(-0.16, 1.06, "c", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 1]
    source_counts = work.groupby("study_id").size().sort_values()
    labels = [SOURCE_LABELS.get(x, x) for x in source_counts.index]
    ax.barh(labels, source_counts.values, color="#D8DEE6", edgecolor="white")
    ax.set_xscale("log")
    ax.set_xlabel("route-eligible rows, log scale")
    ax.set_title("source-size context", fontsize=9, weight="bold")
    ax.text(0.04, 0.06,
            f"threshold: {sel_n} targeted rows\n"
            f"{summary['fraction_of_route_eligible_rows'].iloc[0]:.1%} of route-space rows",
            transform=ax.transAxes, fontsize=6.6, color=COL["grey"],
            bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#D6DEE6"})
    ax.text(-0.16, 1.06, "d", transform=ax.transAxes, fontsize=12, weight="bold")

    fig.savefig(OUT / "row_level_adversarial_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "row_level_adversarial_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "row_level_adversarial_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    routes = route_table(df)
    best_alt_row = routes[~routes["is_prespecified"]].iloc[0]
    best_alt = {
        "host": best_alt_row["host_route"],
        "symbiont": best_alt_row["symbiont_route"],
        "other": best_alt_row["other_route"],
    }
    summary, candidates, selected, by_source, work = row_level_adversary(df, best_alt)

    with pd.ExcelWriter(OUT / "row_level_adversarial_v1.xlsx", engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="summary", index=False)
        routes.to_excel(xw, sheet_name="route_weights", index=False)
        candidates.to_excel(xw, sheet_name="candidate_flips", index=False)
        selected.to_excel(xw, sheet_name="selected_flips", index=False)
        by_source.to_excel(xw, sheet_name="selected_by_source", index=False)
        work.to_excel(xw, sheet_name="route_eligible_full_table", index=False)

    summary.to_csv(OUT / "row_level_adversarial_summary.csv", index=False)
    routes.to_csv(OUT / "row_level_route_weights.csv", index=False)
    candidates.to_csv(OUT / "row_level_candidate_flips.csv", index=False)
    selected.to_csv(OUT / "row_level_selected_flips.csv", index=False)
    by_source.to_csv(OUT / "row_level_selected_by_source.csv", index=False)

    make_figure(routes, summary, candidates, by_source, work)
    (OUT / "README_row_level_adversarial_v1.md").write_text(
        "# Row-level adversarial route erasure\n\n"
        "This analysis uses the full 1,522-row evidence table and overlays adjudicated labels "
        "for the 616 independently audited rows using `combined_index` as the row pointer.\n\n"
        f"- Route-eligible rows after audit overlay: {int(summary['route_eligible_rows_after_audit_overlay'].iloc[0])}\n"
        f"- Audit-overlaid rows: {int(summary['audit_overlay_rows_applied'].iloc[0])}\n"
        f"- Prespecified source-equal support: {summary['prespecified_source_equal_support'].iloc[0]:.3f}\n"
        f"- Nearest alternative source-equal support: {summary['nearest_alternative_source_equal_support'].iloc[0]:.3f}\n"
        f"- Minimum targeted row flips to erase margin: {int(summary['minimum_targeted_row_flips_to_erase_margin'].iloc[0])}\n"
        f"- Fraction of route-eligible rows: {summary['fraction_of_route_eligible_rows'].iloc[0]:.1%}\n\n"
        "Rows are traceable but dependent; this is an adversarial influence calculation under "
        "source-equal weighting, not an independent-observation p value.\n",
        encoding="utf-8",
    )
    print(OUT)


if __name__ == "__main__":
    main()
