from __future__ import annotations

import math
from itertools import permutations
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
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_predictive_permutation_evidence_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
DONOR_CODE = {d: i for i, d in enumerate(DONORS)}
ROLE_CODE = {r: i for i, r in enumerate(ROLES)}
ROLE_LABEL = {
    "scaffold": "scaffold",
    "energetic_incorporation": "energy",
    "transition": "transition",
}
SOURCE_LABELS = {
    "KU2015": "Ku 2015",
    "SANTANA2025": "Santana 2025",
    "ZHANG2025": "Zhang 2025",
    "TOB2026": "Tobiasson 2026",
    "EME2023": "Eme 2023",
}
LAYER_LABELS = {
    "group_sign_test": "group sign",
    "prefix_group_summary": "prefix groups",
    "genome_dynamics_summary": "genome dynamics",
    "proteome_size_estimate": "proteome size",
    "phylogenomic_topology_support": "phylogenomic topology",
    "ogt_prediction": "OGT prediction",
    "extant_genome_size_summary": "extant genome size",
    "ancestral_metabolic_reconstruction": "ancestral metabolism",
    "esp_distribution": "ESP distribution",
    "taxon_occurrence_enrichment": "taxon occurrence",
    "esp_prevalence": "ESP prevalence",
    "orthogroup_donor_profile": "orthogroup donor profile",
    "module_summary": "module summary",
    "clade_pathway_prevalence": "clade pathway",
    "compartment_summary": "compartment summary",
    "ancestor_pathway_frequency": "ancestor pathway",
    "raw_variant_rerun": "raw variant rerun",
}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}

COL = {
    "blue": "#356B8A",
    "orange": "#E0942A",
    "purple": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 7.0,
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
    maps = []
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        label = "; ".join(f"{d}->{ROLE_LABEL[mapping[d]]}" for d in DONORS)
        short = "_".join(ROLE_LABEL[mapping[d]] for d in DONORS)
        maps.append((short, label, mapping, mapping == PRESPEC))
    return maps


def mapping_array(mapping):
    arr = np.zeros(len(DONORS), dtype=int)
    for d in DONORS:
        arr[DONOR_CODE[d]] = ROLE_CODE[mapping[d]]
    return arr


ROUTE_MAPS = route_maps()
PRESPEC_ARR = mapping_array(PRESPEC)


def load_full_with_audit_overlay():
    df = pd.read_csv(FULL_TABLE)
    df["row_index"] = np.arange(len(df))
    df["donor_final"] = df["donor_class"].map(norm)
    df["functional_final"] = df["functional_class"].map(norm)
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
    work["source_label"] = work["study_id"].map(lambda x: SOURCE_LABELS.get(str(x), str(x)))
    return work.reset_index(drop=True)


def source_equal_support_codes(donor_code, role_code, source_id, map_arr):
    vals = []
    for s in np.unique(source_id):
        idx = source_id == s
        vals.append(np.mean(map_arr[donor_code[idx]] == role_code[idx]))
    return float(np.mean(vals))


def route_scores_for_arrays(donor_code, role_code, source_id):
    rows = []
    for short, label, mapping, is_prespec in ROUTE_MAPS:
        arr = mapping_array(mapping)
        rows.append(
            {
                "route_short": short,
                "route_map": label,
                "is_prespecified": is_prespec,
                "source_equal_support": source_equal_support_codes(donor_code, role_code, source_id, arr),
                "raw_support_descriptive": float(np.mean(arr[donor_code] == role_code)),
            }
        )
    out = pd.DataFrame(rows).sort_values("source_equal_support", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def row_weights(df):
    counts = df.groupby("study_id")["row_index"].transform("count").astype(float)
    return 1.0 / counts


def fit_models(train, alpha=1.0):
    donor = train["donor_code"].to_numpy()
    role = train["role_code"].to_numpy()
    weights = row_weights(train).to_numpy()
    weights = weights / weights.sum()

    models = {}
    models["diffuse_uniform"] = np.ones((3, 3), dtype=float) / 3.0

    role_counts = np.zeros(3, dtype=float)
    for r in range(3):
        role_counts[r] = weights[role == r].sum()
    role_prob = (role_counts + alpha) / (role_counts.sum() + 3 * alpha)
    models["donor_independent_global_role"] = np.tile(role_prob, (3, 1))

    counts = np.zeros((3, 3), dtype=float)
    for d in range(3):
        for r in range(3):
            counts[d, r] = weights[(donor == d) & (role == r)].sum()
    models["saturated_donor_specific"] = (counts + alpha) / (counts + alpha).sum(axis=1, keepdims=True)

    for short, label, mapping, is_prespec in ROUTE_MAPS:
        arr = mapping_array(mapping)
        matched = arr[donor] == role
        theta = float(weights[matched].sum() / weights.sum())
        q = np.zeros((3, 3), dtype=float)
        for d in range(3):
            for r in range(3):
                q[d, r] = theta if arr[d] == r else (1 - theta) / 2
        models[f"route_shared_{short}"] = q

        qd = np.zeros((3, 3), dtype=float)
        for d in range(3):
            idx = donor == d
            if idx.sum() == 0:
                theta_d = theta
            else:
                wd = weights[idx] / weights[idx].sum()
                theta_d = float(wd[arr[donor[idx]] == role[idx]].sum())
            for r in range(3):
                qd[d, r] = theta_d if arr[d] == r else (1 - theta_d) / 2
        models[f"route_donor_specific_{short}"] = qd
    return models


def cross_entropy(test, q):
    donor = test["donor_code"].to_numpy()
    role = test["role_code"].to_numpy()
    weights = row_weights(test).to_numpy()
    weights = weights / weights.sum()
    probs = np.clip(q[donor, role], 1e-12, 1.0)
    return float(np.sum(weights * -np.log2(probs)))


def predictive_cv(df, block_col, label):
    rows = []
    blocks = [x for x in sorted(df[block_col].dropna().unique()) if df[block_col].eq(x).sum() > 0]
    for block in blocks:
        train = df[df[block_col] != block].copy()
        test = df[df[block_col] == block].copy()
        if len(train) == 0 or len(test) == 0:
            continue
        models = fit_models(train)
        for model, q in models.items():
            rows.append(
                {
                    "cv_family": label,
                    "heldout_block": block,
                    "heldout_label": SOURCE_LABELS.get(str(block), str(block)) if block_col == "study_id" else str(block),
                    "n_train": len(train),
                    "n_test": len(test),
                    "model": model,
                    "is_prespecified_route_model": model in {
                        "route_shared_scaffold_energy_transition",
                        "route_donor_specific_scaffold_energy_transition",
                    },
                    "is_route_family": model.startswith("route_shared_") or model.startswith("route_donor_specific_"),
                    "cross_entropy_bits": cross_entropy(test, q),
                }
            )
    out = pd.DataFrame(rows)
    out["rank_all_models"] = out.groupby(["cv_family", "heldout_block"])["cross_entropy_bits"].rank(method="min")
    out["rank_route_shared_only"] = np.nan
    shared = out["model"].str.startswith("route_shared_")
    out.loc[shared, "rank_route_shared_only"] = out[shared].groupby(["cv_family", "heldout_block"])[
        "cross_entropy_bits"
    ].rank(method="min")
    return out.sort_values(["cv_family", "heldout_block", "cross_entropy_bits"]).reset_index(drop=True)


def summarize_cv(cv):
    rows = []
    for fam, g in cv.groupby("cv_family"):
        wide = g.pivot_table(index="heldout_block", columns="model", values="cross_entropy_bits", aggfunc="first")
        pres = "route_shared_scaffold_energy_transition"
        pres_ds = "route_donor_specific_scaffold_energy_transition"
        for base in ["diffuse_uniform", "donor_independent_global_role", "saturated_donor_specific"]:
            if base in wide and pres in wide:
                delta = wide[base] - wide[pres]
                rows.append(
                    {
                        "cv_family": fam,
                        "comparison": f"{pres} vs {base}",
                        "n_blocks": int(delta.notna().sum()),
                        "mean_bits_saved_by_prespecified": float(delta.mean()),
                        "median_bits_saved_by_prespecified": float(delta.median()),
                        "blocks_where_prespecified_better": int((delta > 0).sum()),
                        "fraction_blocks_better": float((delta > 0).mean()),
                    }
                )
        shared = g[g["model"].str.startswith("route_shared_")].copy()
        pres_rank = shared[shared["model"].eq(pres)]["rank_route_shared_only"]
        rows.append(
            {
                "cv_family": fam,
                "comparison": "prespecified rank among six shared-theta route maps",
                "n_blocks": int(pres_rank.notna().sum()),
                "mean_bits_saved_by_prespecified": np.nan,
                "median_bits_saved_by_prespecified": np.nan,
                "blocks_where_prespecified_better": int((pres_rank == 1).sum()),
                "fraction_blocks_better": float((pres_rank == 1).mean()),
            }
        )
        if pres_ds in wide and pres in wide:
            delta = wide[pres] - wide[pres_ds]
            rows.append(
                {
                    "cv_family": fam,
                    "comparison": f"{pres_ds} vs {pres}",
                    "n_blocks": int(delta.notna().sum()),
                    "mean_bits_saved_by_prespecified": float(delta.mean()),
                    "median_bits_saved_by_prespecified": float(delta.median()),
                    "blocks_where_prespecified_better": int((delta > 0).sum()),
                    "fraction_blocks_better": float((delta > 0).mean()),
                }
            )
    return pd.DataFrame(rows)


def stratified_permutation(
    df,
    n_iter=10000,
    seed=20260701,
    stratify_cols=("study_id", "evidence_layer"),
    null_label="source_by_layer",
):
    rng = np.random.default_rng(seed)
    donor = df["donor_code"].to_numpy()
    role_obs = df["role_code"].to_numpy()
    source = df["study_id"].astype("category").cat.codes.to_numpy()
    strata_keys = df[list(stratify_cols)].astype(str).agg("||".join, axis=1).astype("category")
    strata = strata_keys.cat.codes.to_numpy()
    strata_indices = [np.where(strata == s)[0] for s in np.unique(strata)]

    observed_routes = route_scores_for_arrays(donor, role_obs, source)
    obs_pres = observed_routes.loc[observed_routes["is_prespecified"], "source_equal_support"].iloc[0]
    obs_best_alt = observed_routes.loc[~observed_routes["is_prespecified"], "source_equal_support"].max()
    obs_margin = obs_pres - obs_best_alt

    draws = []
    for i in range(n_iter):
        role = role_obs.copy()
        for idx in strata_indices:
            if len(idx) > 1:
                role[idx] = rng.permutation(role[idx])
        routes = route_scores_for_arrays(donor, role, source)
        pres = routes.loc[routes["is_prespecified"], "source_equal_support"].iloc[0]
        best_alt = routes.loc[~routes["is_prespecified"], "source_equal_support"].max()
        rank = int(routes.loc[routes["is_prespecified"], "rank"].iloc[0])
        draws.append(
            {
                "iteration": i + 1,
                "null_model": null_label,
                "prespecified_support": pres,
                "best_alternative_support": best_alt,
                "margin_vs_best_alternative": pres - best_alt,
                "prespecified_rank": rank,
            }
        )
    draw_df = pd.DataFrame(draws)
    summary = pd.DataFrame(
        [
            {
                "null_model": "permute_functional_labels_within_source_by_evidence_layer",
                "stratification": " + ".join(stratify_cols),
                "iterations": n_iter,
                "observed_prespecified_support": obs_pres,
                "observed_best_alternative_support": obs_best_alt,
                "observed_margin_vs_best_alternative": obs_margin,
                "null_margin_mean": float(draw_df["margin_vs_best_alternative"].mean()),
                "null_margin_q950": float(draw_df["margin_vs_best_alternative"].quantile(0.95)),
                "null_margin_q975": float(draw_df["margin_vs_best_alternative"].quantile(0.975)),
                "empirical_p_margin_ge_observed": float(
                    (1 + (draw_df["margin_vs_best_alternative"] >= obs_margin).sum()) / (n_iter + 1)
                ),
                "null_rank1_fraction": float((draw_df["prespecified_rank"] == 1).mean()),
                "interpretation": "Functional labels were shuffled within the stated blocks, preserving block composition while breaking donor-role association inside each block.",
            }
        ]
    )
    summary["null_model"] = null_label
    return observed_routes, summary, draw_df


def make_figure(cv_source, cv_layer, cv_summary, perm_summary, perm_draws):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))

    ax = axes[0, 0]
    models = [
        "diffuse_uniform",
        "donor_independent_global_role",
        "route_shared_scaffold_energy_transition",
        "saturated_donor_specific",
    ]
    labels = ["diffuse", "donor-independent", "route", "donor-specific"]
    colors = [COL["light"], "#B8C1CC", COL["blue"], COL["grey"]]
    source_blocks = cv_source["heldout_label"].drop_duplicates().tolist()
    x = np.arange(len(source_blocks))
    width = 0.18
    for i, (model, lab, color) in enumerate(zip(models, labels, colors)):
        vals = []
        for block in source_blocks:
            sub = cv_source[(cv_source["heldout_label"] == block) & (cv_source["model"] == model)]
            vals.append(sub["cross_entropy_bits"].iloc[0] if len(sub) else np.nan)
        ax.bar(x + (i - 1.5) * width, vals, width=width, color=color, edgecolor="white", label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(source_blocks, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("held-out CE (bits)")
    ax.set_title("Leave-one-source prediction", fontsize=8.5, weight="bold")
    ax.legend(frameon=False, fontsize=6, ncol=2)

    ax = axes[0, 1]
    layer = cv_layer.copy()
    pres = layer[layer["model"].eq("route_shared_scaffold_energy_transition")][
        ["heldout_block", "heldout_label", "cross_entropy_bits", "rank_route_shared_only"]
    ].rename(columns={"cross_entropy_bits": "route_ce", "rank_route_shared_only": "route_rank"})
    base = layer[layer["model"].eq("donor_independent_global_role")][["heldout_block", "cross_entropy_bits"]].rename(
        columns={"cross_entropy_bits": "donor_independent_ce"}
    )
    ldat = pres.merge(base, on="heldout_block", how="left")
    ldat["bits_saved_vs_donor_independent"] = ldat["donor_independent_ce"] - ldat["route_ce"]
    ldat = ldat.sort_values("bits_saved_vs_donor_independent")
    bar_colors = [COL["red"] if v < 0 else COL["green"] for v in ldat["bits_saved_vs_donor_independent"]]
    ax.barh(np.arange(len(ldat)), ldat["bits_saved_vs_donor_independent"], color=bar_colors, edgecolor="white")
    ax.axvline(0, color="#111827", lw=0.8)
    ax.set_yticks(np.arange(len(ldat)))
    ax.set_yticklabels([LAYER_LABELS.get(x, x) for x in ldat["heldout_label"]], fontsize=6.2)
    ax.set_xlabel("bits saved by route model")
    ax.set_title("Leave-one-layer prediction", fontsize=8.5, weight="bold")

    ax = axes[1, 0]
    source_vals = perm_draws[perm_draws["null_model"].eq("source_only")]["margin_vs_best_alternative"]
    layer_vals = perm_draws[perm_draws["null_model"].eq("source_by_layer")]["margin_vs_best_alternative"]
    ax.hist(source_vals, bins=50, color="#D8DEE6", edgecolor="white", label="source-only")
    ax.hist(layer_vals, bins=50, color="#7A58A6", alpha=0.24, edgecolor="none", label="source+layer")
    obs = perm_summary["observed_margin_vs_best_alternative"].iloc[0]
    q975_source = perm_summary.loc[perm_summary["null_model"].eq("source_only"), "null_margin_q975"].iloc[0]
    q975_layer = perm_summary.loc[perm_summary["null_model"].eq("source_by_layer"), "null_margin_q975"].iloc[0]
    ax.axvline(q975_source, color=COL["grey"], lw=1.0, ls=":")
    ax.axvline(q975_layer, color=COL["purple"], lw=1.0, ls=":")
    ax.axvline(obs, color=COL["red"], lw=1.5)
    ax.set_xlabel("margin vs best alternative")
    ax.set_ylabel("permutations")
    ax.set_title("Source-layer permutation null", fontsize=8.5, weight="bold")
    ax.text(
        0.98,
        0.95,
        f"obs={obs:.3f}\nq97.5 src={q975_source:.3f}\nq97.5 layer={q975_layer:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.5,
    )
    ax.legend(frameon=False, fontsize=6, loc="upper left")

    ax = axes[1, 1]
    plot = cv_summary[cv_summary["comparison"].str.contains(" vs diffuse_uniform| vs donor_independent_global_role")]
    plot = plot.copy()
    plot["short"] = plot["cv_family"].replace(
        {"leave_one_layer": "leave layer", "leave_one_source": "leave source"}
    ) + " vs " + plot["comparison"].str.extract(r"vs (.*)$")[0].replace(
        {"diffuse_uniform": "diffuse", "donor_independent_global_role": "donor-independent"}
    )
    plot = plot.sort_values("mean_bits_saved_by_prespecified")
    ax.barh(np.arange(len(plot)), plot["mean_bits_saved_by_prespecified"], color=COL["blue"], edgecolor="white")
    ax.axvline(0, color="#111827", lw=0.8)
    ax.set_yticks(np.arange(len(plot)))
    ax.set_yticklabels(plot["short"], fontsize=6.2)
    ax.set_xlabel("mean bits saved")
    ax.set_title("Predictive information gain", fontsize=8.5, weight="bold")

    for label, ax in zip(["a", "b", "c", "d"], axes.ravel()):
        ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=12, weight="bold", va="top")

    fig.tight_layout(w_pad=2.0, h_pad=2.0)
    fig.savefig(OUT / "predictive_permutation_evidence_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "predictive_permutation_evidence_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "predictive_permutation_evidence_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    cv_source = predictive_cv(df, "study_id", "leave_one_source")
    cv_layer = predictive_cv(df, "evidence_layer", "leave_one_layer")
    cv_summary = summarize_cv(pd.concat([cv_source, cv_layer], ignore_index=True))
    observed_routes, perm_source_summary, perm_source_draws = stratified_permutation(
        df,
        n_iter=10000,
        seed=20260701,
        stratify_cols=("study_id",),
        null_label="source_only",
    )
    _, perm_layer_summary, perm_layer_draws = stratified_permutation(
        df,
        n_iter=10000,
        seed=20260702,
        stratify_cols=("study_id", "evidence_layer"),
        null_label="source_by_layer",
    )
    perm_summary = pd.concat([perm_source_summary, perm_layer_summary], ignore_index=True)
    perm_draws = pd.concat([perm_source_draws, perm_layer_draws], ignore_index=True)

    xlsx = OUT / "predictive_permutation_evidence_v1.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        cv_source.to_excel(writer, sheet_name="predictive_cv_source", index=False)
        cv_layer.to_excel(writer, sheet_name="predictive_cv_layer", index=False)
        cv_summary.to_excel(writer, sheet_name="predictive_cv_summary", index=False)
        observed_routes.to_excel(writer, sheet_name="observed_route_scores", index=False)
        perm_summary.to_excel(writer, sheet_name="stratified_perm_summary", index=False)
        perm_draws.to_excel(writer, sheet_name="stratified_perm_draws", index=False)

    for name, frame in [
        ("predictive_cv_source.csv", cv_source),
        ("predictive_cv_layer.csv", cv_layer),
        ("predictive_cv_summary.csv", cv_summary),
        ("observed_route_scores.csv", observed_routes),
        ("stratified_permutation_summary.csv", perm_summary),
        ("stratified_permutation_draws.csv", perm_draws),
    ]:
        frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")

    make_figure(cv_source, cv_layer, cv_summary, perm_summary, perm_draws)

    readme = OUT / "README_predictive_permutation_evidence_v1.md"
    readme.write_text(
        "# Predictive and permutation computational evidence v1\n\n"
        "This analysis adds two tests that are intentionally different from whole-table route scoring.\n\n"
        "1. Leave-one-source and leave-one-layer predictive cross-validation. Route models are fitted "
        "on training blocks and evaluated on held-out sources or evidence layers by cross entropy.\n"
        "2. A source-by-evidence-layer stratified permutation null. Functional labels are shuffled "
        "inside each source/layer block while donor labels and source/layer composition are preserved.\n\n"
        "These tests ask whether the route map has predictive and non-random structure beyond a "
        "descriptive recoding of the pooled evidence table. Row entries remain traceable but dependent; "
        "therefore cross-entropy and permutation results are interpreted as computational stress tests, "
        "not as independent-observation causal tests.\n",
        encoding="utf-8",
    )

    print(xlsx)
    print("\nPredictive CV summary")
    print(cv_summary.to_string(index=False))
    print("\nPermutation summary")
    print(perm_summary.to_string(index=False))


if __name__ == "__main__":
    main()
