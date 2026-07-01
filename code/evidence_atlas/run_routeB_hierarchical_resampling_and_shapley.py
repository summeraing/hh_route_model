from __future__ import annotations

import math
from itertools import combinations, permutations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
FULL_TABLE = (
    Path.cwd()
    / "TOP_JOURNAL_REBUILD_ANALYSES_v1"
    / "v77_clean"
    / "01_MAIN_FIGURES"
    / "Fig1"
    / "source_package"
    / "Fig1B_selected_original_network"
    / "source_combined_v4_flat_export.csv"
)
SOURCE_XLSX = ROOT / "10_WORKING_COPY_FROM_NC12" / "04_Source_Data" / "Source_Data_complete.xlsx"
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_hierarchical_resampling_shapley_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
DONOR_CODE = {d: i for i, d in enumerate(DONORS)}
ROLE_CODE = {r: i for i, r in enumerate(ROLES)}
ROLE_LABEL = {"scaffold": "scaffold", "energetic_incorporation": "energy", "transition": "transition"}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}
SOURCE_LABEL = {
    "KU2015": "Ku 2015",
    "SANTANA2025": "Santana 2025",
    "ZHANG2025": "Zhang 2025",
    "TOB2026": "Tobiasson 2026",
    "EME2023": "Eme 2023",
}
LAYER_LABEL = {
    "orthogroup_donor_profile": "orthogroup donor profile",
    "module_summary": "module summary",
    "taxon_occurrence_enrichment": "taxon occurrence",
    "prefix_group_summary": "prefix groups",
    "group_sign_test": "group sign",
    "esp_distribution": "ESP distribution",
    "extant_genome_size_summary": "genome size",
    "phylogenomic_topology_support": "phylogenomic topology",
    "ogt_prediction": "OGT prediction",
    "ancestral_metabolic_reconstruction": "ancestral metabolism",
    "compartment_summary": "compartment summary",
    "raw_variant_rerun": "raw variant rerun",
    "ancestor_pathway_frequency": "ancestor pathway",
    "clade_pathway_prevalence": "clade pathway",
    "esp_prevalence": "ESP prevalence",
    "genome_dynamics_summary": "genome dynamics",
    "proteome_size_estimate": "proteome size",
}

COL = {
    "blue": "#356B8A",
    "orange": "#E0942A",
    "purple": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
}

N_BOOT = 1000
N_SHAPLEY_LAYER = 3000
RNG_SEED = 20260701

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 7,
    }
)


def norm(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = s.replace("injection", "energetic_incorporation")
    s = s.replace("energy", "energetic_incorporation")
    return s


def route_maps():
    out = []
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        arr = np.array([ROLE_CODE[mapping[d]] for d in DONORS], dtype=int)
        label = "; ".join(f"{d}->{ROLE_LABEL[mapping[d]]}" for d in DONORS)
        short = "_".join(ROLE_LABEL[mapping[d]] for d in DONORS)
        out.append((short, label, arr, mapping == PRESPEC))
    return out


ROUTE_MAPS = route_maps()
PRESPEC_ARR = np.array([ROLE_CODE[PRESPEC[d]] for d in DONORS], dtype=int)


def load_full_with_audit_overlay() -> tuple[pd.DataFrame, pd.DataFrame]:
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

    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    work["donor_code"] = work["donor_final"].map(DONOR_CODE).astype(int)
    work["role_code"] = work["functional_final"].map(ROLE_CODE).astype(int)
    work["source_label"] = work["study_id"].map(lambda x: SOURCE_LABEL.get(str(x), str(x)))
    return work.reset_index(drop=True), audit


def route_scores_for_subset(df: pd.DataFrame, role_codes: np.ndarray | None = None) -> pd.DataFrame:
    if role_codes is None:
        role_codes = df["role_code"].to_numpy()
    donor_codes = df["donor_code"].to_numpy()
    sources = df["study_id"].to_numpy()
    rows = []
    for short, label, arr, is_prespec in ROUTE_MAPS:
        per_source = []
        for s in np.unique(sources):
            idx = sources == s
            if idx.sum():
                per_source.append(float(np.mean(arr[donor_codes[idx]] == role_codes[idx])))
        support = float(np.mean(per_source)) if per_source else np.nan
        raw = float(np.mean(arr[donor_codes] == role_codes)) if len(df) else np.nan
        rows.append(
            {
                "route_short": short,
                "route_map": label,
                "is_prespecified": is_prespec,
                "source_equal_support": support,
                "raw_support_descriptive": raw,
            }
        )
    out = pd.DataFrame(rows).sort_values("source_equal_support", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def margin_for_subset(df: pd.DataFrame, role_codes: np.ndarray | None = None) -> tuple[float, float, float, int]:
    if df.empty or df["study_id"].nunique() == 0:
        return np.nan, np.nan, np.nan, 99
    routes = route_scores_for_subset(df, role_codes=role_codes)
    pres = routes[routes["is_prespecified"]].iloc[0]
    best_alt = routes[~routes["is_prespecified"]].iloc[0]
    return (
        float(pres["source_equal_support"]),
        float(best_alt["source_equal_support"]),
        float(pres["source_equal_support"] - best_alt["source_equal_support"]),
        int(pres["rank"]),
    )


def build_label_noise_lookup(audit: pd.DataFrame) -> dict[str, dict[int, np.ndarray]]:
    lookup: dict[str, dict[int, list[int]]] = {"donor": {i: [] for i in range(3)}, "role": {i: [] for i in range(3)}}
    for _, row in audit.iterrows():
        adj_d = norm(row.get("adjudicated_donor_class"))
        adj_f = norm(row.get("adjudicated_functional_class"))
        if adj_d in DONOR_CODE:
            for col in ["coder1_donor_class", "coder2_donor_class", "adjudicated_donor_class"]:
                val = norm(row.get(col))
                if val in DONOR_CODE:
                    lookup["donor"][DONOR_CODE[adj_d]].append(DONOR_CODE[val])
        if adj_f in ROLE_CODE:
            for col in ["coder1_functional_class", "coder2_functional_class", "adjudicated_functional_class"]:
                val = norm(row.get(col))
                if val in ROLE_CODE:
                    lookup["role"][ROLE_CODE[adj_f]].append(ROLE_CODE[val])
    out: dict[str, dict[int, np.ndarray]] = {"donor": {}, "role": {}}
    for kind in ["donor", "role"]:
        for k, vals in lookup[kind].items():
            out[kind][k] = np.array(vals if vals else [k], dtype=int)
    return out


def draw_noisy_codes(
    donor_codes: np.ndarray,
    role_codes: np.ndarray,
    lookup: dict[str, dict[int, np.ndarray]],
    rng: np.random.Generator,
    mode: str,
    audit_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    donor = donor_codes.copy()
    role = role_codes.copy()
    if mode == "adjudicated_labels":
        return donor, role

    target = np.where(audit_mask)[0] if mode == "audit_rows_coder_draw" else np.arange(len(donor_codes))
    for code in range(len(DONORS)):
        idx = target[donor_codes[target] == code]
        if len(idx):
            donor[idx] = rng.choice(lookup["donor"][code], size=len(idx), replace=True)
    for code in range(len(ROLES)):
        idx = target[role_codes[target] == code]
        if len(idx):
            role[idx] = rng.choice(lookup["role"][code], size=len(idx), replace=True)
    return donor, role


def hierarchical_bootstrap(df: pd.DataFrame, audit: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RNG_SEED)
    sources = np.array(sorted(df["study_id"].unique()))
    src_layers = {
        s: np.array(sorted(df.loc[df["study_id"] == s, "evidence_layer"].dropna().unique()))
        for s in sources
    }
    idx_by_source = {s: df.index[df["study_id"] == s].to_numpy() for s in sources}
    idx_by_source_layer = {
        (s, layer): df.index[(df["study_id"] == s) & (df["evidence_layer"] == layer)].to_numpy()
        for s in sources
        for layer in src_layers[s]
    }
    donor_base = df["donor_code"].to_numpy()
    role_base = df["role_code"].to_numpy()
    audit_mask = df["audit_overlay_applied"].to_numpy(dtype=bool)
    lookup = build_label_noise_lookup(audit)

    designs = [
        ("source_cluster", "adjudicated_labels"),
        ("source_layer_cluster", "adjudicated_labels"),
        ("source_layer_row", "adjudicated_labels"),
        ("source_layer_row", "audit_rows_coder_draw"),
        ("source_layer_row", "projected_coder_noise"),
    ]
    draws = []
    for design, noise_mode in designs:
        for i in range(N_BOOT):
            donor_noisy, role_noisy = draw_noisy_codes(donor_base, role_base, lookup, rng, noise_mode, audit_mask)
            block_scores = []
            sampled_sources = rng.choice(sources, size=len(sources), replace=True)
            for s in sampled_sources:
                if design == "source_cluster":
                    idx = idx_by_source[s]
                    d = donor_noisy[idx]
                    r = role_noisy[idx]
                    block_scores.append([float(np.mean(arr[d] == r)) for _, _, arr, _ in ROUTE_MAPS])
                else:
                    layers = src_layers[s]
                    sampled_layers = rng.choice(layers, size=len(layers), replace=True)
                    layer_idx = []
                    for layer in sampled_layers:
                        base_idx = idx_by_source_layer[(s, layer)]
                        if design == "source_layer_row":
                            layer_idx.append(rng.choice(base_idx, size=len(base_idx), replace=True))
                        else:
                            layer_idx.append(base_idx)
                    idx = np.concatenate(layer_idx) if layer_idx else np.array([], dtype=int)
                    if len(idx):
                        d = donor_noisy[idx]
                        r = role_noisy[idx]
                        block_scores.append([float(np.mean(arr[d] == r)) for _, _, arr, _ in ROUTE_MAPS])
            scores = np.nanmean(np.asarray(block_scores), axis=0)
            order = np.argsort(scores)[::-1]
            pres_idx = [j for j, (_, _, _, isp) in enumerate(ROUTE_MAPS) if isp][0]
            best_alt = [j for j in order if j != pres_idx][0]
            draws.append(
                {
                    "iteration": i + 1,
                    "resampling_design": design,
                    "label_noise_mode": noise_mode,
                    "prespecified_support": scores[pres_idx],
                    "best_alternative_support": scores[best_alt],
                    "margin_vs_best_alternative": scores[pres_idx] - scores[best_alt],
                    "prespecified_rank": int(np.where(order == pres_idx)[0][0] + 1),
                    "best_alternative_route": ROUTE_MAPS[best_alt][1],
                }
            )

    draws_df = pd.DataFrame(draws)
    summary = (
        draws_df.groupby(["resampling_design", "label_noise_mode"])
        .agg(
            iterations=("iteration", "count"),
            rank1_probability=("prespecified_rank", lambda x: float(np.mean(np.asarray(x) == 1))),
            margin_mean=("margin_vs_best_alternative", "mean"),
            margin_q025=("margin_vs_best_alternative", lambda x: float(np.quantile(x, 0.025))),
            margin_q50=("margin_vs_best_alternative", lambda x: float(np.quantile(x, 0.5))),
            margin_q975=("margin_vs_best_alternative", lambda x: float(np.quantile(x, 0.975))),
            support_mean=("prespecified_support", "mean"),
            support_q025=("prespecified_support", lambda x: float(np.quantile(x, 0.025))),
            support_q975=("prespecified_support", lambda x: float(np.quantile(x, 0.975))),
        )
        .reset_index()
    )
    return summary, draws_df


def source_subset_and_shapley(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sources = sorted(df["study_id"].unique())
    subset_rows = []
    margin_cache: dict[frozenset[str], float] = {}
    for r in range(1, len(sources) + 1):
        for combo in combinations(sources, r):
            sub = df[df["study_id"].isin(combo)]
            pres, alt, margin, rank = margin_for_subset(sub)
            key = frozenset(combo)
            margin_cache[key] = margin
            subset_rows.append(
                {
                    "source_subset": " + ".join(combo),
                    "n_sources": len(combo),
                    "n_rows": len(sub),
                    "prespecified_support": pres,
                    "best_alternative_support": alt,
                    "margin_vs_best_alternative": margin,
                    "prespecified_rank": rank,
                    "retains_rank1": rank == 1,
                }
            )
    empty = frozenset()
    margin_cache[empty] = 0.0
    shapley = []
    n = len(sources)
    for s in sources:
        contribs = []
        others = [x for x in sources if x != s]
        for r in range(0, len(others) + 1):
            for combo in combinations(others, r):
                S = frozenset(combo)
                weight = 1.0 / (n * math.comb(n - 1, r))
                contribs.append(weight * (margin_cache[S | {s}] - margin_cache[S]))
        shapley.append({"study_id": s, "study_label": SOURCE_LABEL.get(s, s), "shapley_margin_contribution": float(np.sum(contribs))})
    return pd.DataFrame(subset_rows), pd.DataFrame(shapley).sort_values("shapley_margin_contribution", ascending=False)


def layer_shapley(df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED + 1)
    layers = np.array(sorted(df["evidence_layer"].unique()), dtype=object)
    sources = np.array(sorted(df["study_id"].unique()), dtype=object)
    layer_idx = {layer: i for i, layer in enumerate(layers)}
    source_idx = {source: i for i, source in enumerate(sources)}
    counts = np.zeros((len(layers), len(sources), len(DONORS), len(ROLES)), dtype=float)
    for _, row in df.iterrows():
        counts[
            layer_idx[row["evidence_layer"]],
            source_idx[row["study_id"]],
            int(row["donor_code"]),
            int(row["role_code"]),
        ] += 1.0

    route_arrays = [arr for _, _, arr, _ in ROUTE_MAPS]
    pres_idx = [i for i, (_, _, _, isp) in enumerate(ROUTE_MAPS) if isp][0]
    contrib = {str(layer): [] for layer in layers}

    def margin_from_layer_mask(active_mask: np.ndarray) -> float:
        if not active_mask.any():
            return 0.0
        summed = counts[active_mask].sum(axis=0)
        total_by_source = summed.sum(axis=(1, 2))
        active_sources = total_by_source > 0
        if not active_sources.any():
            return 0.0
        supports = []
        for arr in route_arrays:
            matched_by_source = np.zeros(len(sources), dtype=float)
            for d in range(len(DONORS)):
                matched_by_source += summed[:, d, arr[d]]
            per_source = matched_by_source[active_sources] / total_by_source[active_sources]
            supports.append(float(np.mean(per_source)))
        supports = np.asarray(supports)
        best_alt = np.max(np.delete(supports, pres_idx))
        return float(supports[pres_idx] - best_alt)

    for _ in range(N_SHAPLEY_LAYER):
        order = rng.permutation(np.arange(len(layers)))
        active_mask = np.zeros(len(layers), dtype=bool)
        prev = 0.0
        for idx in order:
            active_mask[idx] = True
            cur = margin_from_layer_mask(active_mask)
            contrib[str(layers[idx])].append(cur - prev)
            prev = cur
    rows = []
    for layer, vals in contrib.items():
        arr = np.asarray(vals)
        sub = df[df["evidence_layer"] == layer]
        rows.append(
            {
                "evidence_layer": layer,
                "n_rows": len(sub),
                "n_sources": sub["study_id"].nunique(),
                "mean_shapley_margin_contribution": float(np.mean(arr)),
                "q025": float(np.quantile(arr, 0.025)),
                "q50": float(np.quantile(arr, 0.5)),
                "q975": float(np.quantile(arr, 0.975)),
            }
        )
    return pd.DataFrame(rows).sort_values("mean_shapley_margin_contribution", ascending=False)


def cross_scale_ledger(df: pd.DataFrame, hier_summary: pd.DataFrame, source_subsets: pd.DataFrame) -> pd.DataFrame:
    src_rank = source_subsets[source_subsets["n_sources"] >= 2]["retains_rank1"].mean()
    source_layer = hier_summary[
        (hier_summary["resampling_design"] == "source_layer_row")
        & (hier_summary["label_noise_mode"] == "adjudicated_labels")
    ].iloc[0]
    projected = hier_summary[
        (hier_summary["resampling_design"] == "source_layer_row")
        & (hier_summary["label_noise_mode"] == "projected_coder_noise")
    ].iloc[0]
    rows = [
            {
                "evidence_axis": "row-level source-equal route map",
                "primary_quantity": "full-table margin over nearest alternative",
                "value": margin_for_subset(df)[2],
                "direction": "supports route-resolved architecture",
                "boundary": "traceable rows are dependent; source-equal weighting is required",
            },
            {
                "evidence_axis": "hierarchical source-layer-row resampling",
                "primary_quantity": "rank-1 probability",
                "value": source_layer["rank1_probability"],
                "direction": "tests source/layer dependence",
                "boundary": "cluster resampling, not independent-observation inference",
            },
            {
                "evidence_axis": "projected coder-noise resampling",
                "primary_quantity": "rank-1 probability under measured label noise",
                "value": projected["rank1_probability"],
                "direction": "tests coding uncertainty",
                "boundary": "projects audit-subset noise to all route rows",
            },
            {
                "evidence_axis": "source-coalition analysis",
                "primary_quantity": "rank-1 fraction among multi-source coalitions",
                "value": src_rank,
                "direction": "tests dependence on a single source",
                "boundary": "single-source subsets remain diagnostic anchors rather than complete tests",
            },
    ]

    try:
        fes = pd.read_excel(SOURCE_XLSX, sheet_name="fes_compartment")
        fixed = fes[fes["run_label"].astype(str).str.contains("fixed", case=False, na=False)].copy()
        if fixed.empty:
            fixed = fes.copy()
        mito = fixed.loc[fixed["compartment_class"] == "mitochondrial_integration", "mean_alpha_fraction"]
        cia = fixed.loc[fixed["compartment_class"] == "cytosolic_CIA", "mean_alpha_fraction"]
        if len(mito) and len(cia):
            rows.append(
                {
                    "evidence_axis": "Fe-S compartment phylogeny",
                    "primary_quantity": "mitochondrial-integration alpha fraction minus cytosolic CIA alpha fraction",
                    "value": float(mito.mean() - cia.mean()),
                    "direction": "supports localized energetic incorporation rather than wholesale Fe-S replacement",
                    "boundary": "marker trees are support-bearing sequence screens, not complete ancestral reconstructions",
                }
            )
    except Exception:
        pass

    try:
        af3 = pd.read_excel(SOURCE_XLSX, sheet_name="AF3_integrated_metrics")
        primary = af3[af3["use_in_claim"].astype(str).str.contains("primary", case=False, na=False)].copy()
        controls = af3[
            af3["route_or_control"].astype(str).str.contains("decoy|control|unrelated", case=False, na=False)
            | af3["use_in_claim"].astype(str).str.contains("negative", case=False, na=False)
        ].copy()
        if len(primary) and len(controls):
            rows.append(
                {
                    "evidence_axis": "AF3 structural-interactivity screen",
                    "primary_quantity": "mean primary-interface ipTM minus mean negative-control ipTM",
                    "value": float(primary["ipTM"].mean() - controls["ipTM"].mean()),
                    "direction": "supports confidence-gated local structural compatibility",
                    "boundary": "AF3 interfaces are plausibility screens and gate calibrations, not ancient complex reconstructions",
                }
            )
    except Exception:
        pass

    return pd.DataFrame(rows)


def make_figure(hier_summary: pd.DataFrame, hier_draws: pd.DataFrame, source_shapley: pd.DataFrame, layer_shap: pd.DataFrame) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.8), constrained_layout=True)

    ax = axs[0, 0]
    order = [
        ("source_cluster", "adjudicated_labels"),
        ("source_layer_cluster", "adjudicated_labels"),
        ("source_layer_row", "adjudicated_labels"),
        ("source_layer_row", "audit_rows_coder_draw"),
        ("source_layer_row", "projected_coder_noise"),
    ]
    labels = ["source\ncluster", "source+\nlayer", "source+\nlayer+\nrow", "audit\nnoise", "projected\nnoise"]
    data = [
        hier_draws[
            (hier_draws["resampling_design"] == d) & (hier_draws["label_noise_mode"] == n)
        ]["margin_vs_best_alternative"].to_numpy()
        for d, n in order
    ]
    parts = ax.violinplot(data, showmedians=False, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("#DCE8F2")
        body.set_edgecolor("#356B8A")
        body.set_alpha(0.9)
    for i, vals in enumerate(data, start=1):
        q = np.quantile(vals, [0.025, 0.5, 0.975])
        ax.plot([i - 0.22, i + 0.22], [q[1], q[1]], color="#111827", lw=1.2)
        ax.plot([i, i], [q[0], q[2]], color="#356B8A", lw=1)
    ax.axhline(0, color="#111827", lw=0.8)
    ax.set_xticks(np.arange(1, len(labels) + 1), labels)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.set_ylabel("margin vs nearest alternative")
    ax.text(-0.16, 1.05, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[0, 1]
    sub = hier_summary.copy()
    sub["label"] = [
        "source cluster",
        "source+layer",
        "source+layer+row",
        "audit noise",
        "projected noise",
    ]
    ax.barh(np.arange(len(sub)), sub["rank1_probability"], color=[COL["blue"], COL["blue"], COL["blue"], COL["purple"], COL["purple"]], alpha=0.9)
    ax.set_yticks(np.arange(len(sub)), sub["label"])
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("rank-1 probability")
    for y, v in enumerate(sub["rank1_probability"]):
        ax.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=6.7)
    ax.text(-0.16, 1.05, "b", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 0]
    ss = source_shapley.sort_values("shapley_margin_contribution")
    ax.barh(ss["study_label"], ss["shapley_margin_contribution"], color=COL["grey"], edgecolor="white")
    ax.axvline(0, color="#111827", lw=0.7)
    ax.set_xlabel("source Shapley contribution to margin")
    ax.text(-0.16, 1.05, "c", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 1]
    ls = layer_shap.sort_values("mean_shapley_margin_contribution", ascending=False).head(10)
    ls = ls.sort_values("mean_shapley_margin_contribution")
    ax.barh(
        [LAYER_LABEL.get(x, x) for x in ls["evidence_layer"]],
        ls["mean_shapley_margin_contribution"],
        color=COL["green"],
        edgecolor="white",
        alpha=0.9,
    )
    ax.axvline(0, color="#111827", lw=0.7)
    ax.set_xlabel("layer Shapley contribution to margin")
    ax.text(-0.16, 1.05, "d", transform=ax.transAxes, fontsize=12, weight="bold")

    fig.savefig(OUT / "hierarchical_resampling_shapley_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "hierarchical_resampling_shapley_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "hierarchical_resampling_shapley_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df, audit = load_full_with_audit_overlay()
    hier_summary, hier_draws = hierarchical_bootstrap(df, audit)
    source_subsets, source_shapley = source_subset_and_shapley(df)
    layer_shap = layer_shapley(df)
    ledger = cross_scale_ledger(df, hier_summary, source_subsets)

    with pd.ExcelWriter(OUT / "hierarchical_resampling_shapley_v1.xlsx", engine="openpyxl") as writer:
        hier_summary.to_excel(writer, sheet_name="hierarchical_summary", index=False)
        hier_draws.sample(min(10000, len(hier_draws)), random_state=RNG_SEED).to_excel(
            writer, sheet_name="hierarchical_draws_sample", index=False
        )
        source_subsets.to_excel(writer, sheet_name="source_coalitions", index=False)
        source_shapley.to_excel(writer, sheet_name="source_shapley", index=False)
        layer_shap.to_excel(writer, sheet_name="layer_shapley", index=False)
        ledger.to_excel(writer, sheet_name="cross_scale_ledger", index=False)

    for name, table in [
        ("hierarchical_summary", hier_summary),
        ("hierarchical_draws", hier_draws),
        ("source_coalitions", source_subsets),
        ("source_shapley", source_shapley),
        ("layer_shapley", layer_shap),
        ("cross_scale_ledger", ledger),
    ]:
        table.to_csv(OUT / f"{name}.csv", index=False, encoding="utf-8-sig")

    (OUT / "README_hierarchical_resampling_shapley_v1.md").write_text(
        "# Hierarchical resampling and contribution analysis\n\n"
        "This analysis adds two computational checks to the route-resolved evidence-atlas rebuild.\n\n"
        "1. Hierarchical resampling treats public evidence as clustered by source, evidence layer "
        "and row, rather than as independent row-level observations. Three cluster designs are "
        "reported: source-only, source+layer and source+layer+row. Two label-noise variants then "
        "draw labels from the independent-coding audit: one only for audited rows and one projected "
        "to all route-eligible rows.\n\n"
        "2. Source and evidence-layer Shapley analyses quantify marginal contributions to the "
        "route margin. Source Shapley is exact over all source coalitions; layer Shapley uses "
        f"{N_SHAPLEY_LAYER:,} random layer orderings.\n\n"
        "Outputs are computational stress tests, not independent-observation p values. They are "
        "designed to expose source dependence, row dependence, label-noise sensitivity and evidence "
        "contribution boundaries.\n",
        encoding="utf-8",
    )

    make_figure(hier_summary, hier_draws, source_shapley, layer_shap)

    print("Wrote", OUT)
    print(hier_summary.to_string(index=False))
    print(source_shapley.to_string(index=False))
    print(layer_shap.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
