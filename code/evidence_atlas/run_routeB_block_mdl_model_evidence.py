from __future__ import annotations

from itertools import permutations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
FULL_TABLE = REPO_ROOT / "data" / "evidence_atlas" / "evidence_units_flat.csv"
CURRENT_SOURCE_XLSX = REPO_ROOT / "data" / "evidence_atlas" / "Source_Data_route_resolved_evidence_atlas_complete.xlsx"
OUT = REPO_ROOT / "outputs" / "routeB_block_mdl_model_evidence_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}
ROLE_SHORT = {"scaffold": "scaffold", "energetic_incorporation": "energy", "transition": "transition"}
DONOR_LABEL = {"host": "host", "symbiont": "alpha", "other": "other"}
COL = {
    "host": "#356B8A",
    "symbiont": "#E0942A",
    "other": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
    "blue": "#356B8A",
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


def load_full_with_audit_overlay() -> pd.DataFrame:
    df = pd.read_csv(FULL_TABLE)
    df["row_index"] = np.arange(len(df))
    df["donor_original"] = df["donor_class"].map(norm)
    df["functional_original"] = df["functional_class"].map(norm)
    df["donor_final"] = df["donor_original"]
    df["functional_final"] = df["functional_original"]
    df["audit_overlay_applied"] = False

    audit = pd.read_excel(CURRENT_SOURCE_XLSX, sheet_name="audit_overlay")
    for _, row in audit.iterrows():
        idx = int(row["combined_index"])
        if 0 <= idx < len(df):
            df.loc[idx, "donor_final"] = norm(row["adjudicated_donor_class"])
            df.loc[idx, "functional_final"] = norm(row["adjudicated_functional_class"])
            df.loc[idx, "audit_overlay_applied"] = True

    for col in ["study_id", "evidence_layer", "evidence_unit", "unit_id", "unit_label", "module_family"]:
        if col not in df.columns:
            df[col] = ""

    df["source_layer"] = df["study_id"].astype(str) + "::" + df["evidence_layer"].astype(str)
    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    return work.reset_index(drop=True)


def route_maps() -> list[dict[str, str]]:
    return [dict(zip(DONORS, perm)) for perm in permutations(ROLES)]


def map_name(mapping: dict[str, str]) -> str:
    return "; ".join(f"{DONOR_LABEL[d]}->{ROLE_SHORT[mapping[d]]}" for d in DONORS)


def row_weights(df: pd.DataFrame, source_equal: bool = True) -> np.ndarray:
    if not source_equal:
        return np.ones(len(df), dtype=float) / max(len(df), 1)
    counts = df.groupby("study_id")["study_id"].transform("size").astype(float)
    n_sources = df["study_id"].nunique()
    return (1.0 / counts / max(n_sources, 1)).to_numpy(dtype=float)


def weighted_count_matrix(df: pd.DataFrame, weights: np.ndarray | None = None) -> np.ndarray:
    if weights is None:
        weights = row_weights(df, source_equal=True)
    mat = np.zeros((len(DONORS), len(ROLES)), dtype=float)
    d_idx = {d: i for i, d in enumerate(DONORS)}
    r_idx = {r: i for i, r in enumerate(ROLES)}
    for w, d, r in zip(weights, df["donor_final"], df["functional_final"]):
        mat[d_idx[d], r_idx[r]] += float(w)
    return mat


def source_equal_joint(df: pd.DataFrame) -> pd.DataFrame:
    mats = []
    for _, g in df.groupby("study_id", sort=False):
        tab = pd.crosstab(g["donor_final"], g["functional_final"]).reindex(index=DONORS, columns=ROLES, fill_value=0)
        if tab.values.sum() > 0:
            mats.append(tab / tab.values.sum())
    joint = sum(mats) / len(mats)
    return joint.reindex(index=DONORS, columns=ROLES)


def fit_probabilities(df: pd.DataFrame, model: str, mapping: dict[str, str] | None = None, alpha: float = 0.5) -> np.ndarray:
    counts = weighted_count_matrix(df)
    q = np.zeros_like(counts, dtype=float)
    if model == "diffuse_uniform":
        q[:] = 1 / len(ROLES)
    elif model == "donor_independent_global":
        role = counts.sum(axis=0) + alpha
        role = role / role.sum()
        q[:] = role[None, :]
    elif model == "saturated_donor_role":
        q = counts + alpha
        q = q / q.sum(axis=1, keepdims=True)
    elif model in {"route_shared_theta", "route_donor_specific_theta"}:
        assert mapping is not None
        for i, d in enumerate(DONORS):
            on_idx = ROLES.index(mapping[d])
            denom = counts[i].sum()
            if model == "route_shared_theta":
                continue
            theta = (counts[i, on_idx] + alpha) / (denom + alpha * len(ROLES))
            q[i, :] = (1 - theta) / (len(ROLES) - 1)
            q[i, on_idx] = theta
        if model == "route_shared_theta":
            on = 0.0
            total = counts.sum()
            for i, d in enumerate(DONORS):
                on += counts[i, ROLES.index(mapping[d])]
            theta = (on + alpha) / (total + alpha * len(ROLES))
            for i, d in enumerate(DONORS):
                on_idx = ROLES.index(mapping[d])
                q[i, :] = (1 - theta) / (len(ROLES) - 1)
                q[i, on_idx] = theta
    else:
        raise ValueError(f"Unknown model: {model}")
    q = np.clip(q, 1e-12, 1)
    return q / q.sum(axis=1, keepdims=True)


def model_specs() -> list[dict]:
    specs = [
        {"model_id": "diffuse_uniform", "family": "diffuse_uniform", "k": 0, "mapping": None},
        {"model_id": "donor_independent_global", "family": "donor_independent_global", "k": 2, "mapping": None},
        {"model_id": "saturated_donor_role", "family": "saturated_donor_role", "k": 6, "mapping": None},
    ]
    for mapping in route_maps():
        route = "prespecified" if mapping == PRESPEC else "alternative"
        specs.append(
            {
                "model_id": "shared_theta__" + map_name(mapping),
                "family": "route_shared_theta",
                "route_type": route,
                "k": 1,
                "mapping": mapping,
            }
        )
        specs.append(
            {
                "model_id": "donor_specific_theta__" + map_name(mapping),
                "family": "route_donor_specific_theta",
                "route_type": route,
                "k": 3,
                "mapping": mapping,
            }
        )
    return specs


def score_df(df: pd.DataFrame, q: np.ndarray, source_equal: bool = True) -> float:
    weights = row_weights(df, source_equal=source_equal)
    d_idx = {d: i for i, d in enumerate(DONORS)}
    r_idx = {r: i for i, r in enumerate(ROLES)}
    losses = []
    for w, d, r in zip(weights, df["donor_final"], df["functional_final"]):
        losses.append(float(w) * (-np.log2(q[d_idx[d], r_idx[r]])))
    return float(np.sum(losses) / np.sum(weights))


def block_cv(df: pd.DataFrame, block_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    specs = model_specs()
    for block, test in df.groupby(block_col, sort=False):
        train = df[df[block_col] != block].copy()
        if train.empty or test.empty:
            continue
        for spec in specs:
            q = fit_probabilities(train, spec["family"], mapping=spec.get("mapping"))
            rows.append(
                {
                    "block_design": "leave_" + block_col,
                    "heldout_block": block,
                    "n_train": len(train),
                    "n_test": len(test),
                    "model_id": spec["model_id"],
                    "model_family": spec["family"],
                    "route_type": spec.get("route_type", "not_route"),
                    "k_parameters": spec["k"],
                    "heldout_log_loss_bits": score_df(test, q, source_equal=False),
                }
            )
    scores = pd.DataFrame(rows)
    summary = (
        scores.groupby(["block_design", "model_id", "model_family", "route_type", "k_parameters"], dropna=False)
        .agg(
            mean_heldout_log_loss_bits=("heldout_log_loss_bits", "mean"),
            median_heldout_log_loss_bits=("heldout_log_loss_bits", "median"),
            sd_heldout_log_loss_bits=("heldout_log_loss_bits", "std"),
            n_blocks=("heldout_block", "nunique"),
            total_test_rows=("n_test", "sum"),
        )
        .reset_index()
        .sort_values(["block_design", "mean_heldout_log_loss_bits"])
    )
    summary["rank_within_design"] = summary.groupby("block_design")["mean_heldout_log_loss_bits"].rank(method="first")
    return scores, summary


def full_mdl(df: pd.DataFrame, block_col: str = "source_layer") -> tuple[pd.DataFrame, pd.DataFrame]:
    specs = model_specs()
    blocks = list(df.groupby(block_col, sort=False))
    rows = []
    prob_rows = []
    for spec in specs:
        q = fit_probabilities(df, spec["family"], mapping=spec.get("mapping"))
        block_loss_sum = 0.0
        block_n = 0
        for _, g in blocks:
            block_loss_sum += score_df(g, q, source_equal=False)
            block_n += 1
        penalty = 0.5 * spec["k"] * np.log2(max(block_n, 2))
        mdl = block_loss_sum + penalty
        rows.append(
            {
                "block_design": block_col,
                "model_id": spec["model_id"],
                "model_family": spec["family"],
                "route_type": spec.get("route_type", "not_route"),
                "k_parameters": spec["k"],
                "n_blocks": block_n,
                "sum_block_log_loss_bits": block_loss_sum,
                "complexity_penalty_bits": penalty,
                "mdl_bits": mdl,
                "mean_block_log_loss_bits": block_loss_sum / block_n,
            }
        )
        for i, d in enumerate(DONORS):
            for j, r in enumerate(ROLES):
                prob_rows.append(
                    {
                        "model_id": spec["model_id"],
                        "model_family": spec["family"],
                        "route_type": spec.get("route_type", "not_route"),
                        "donor_class": d,
                        "functional_role": r,
                        "probability": q[i, j],
                    }
                )
    mdl = pd.DataFrame(rows).sort_values("mdl_bits").reset_index(drop=True)
    mdl["rank_by_mdl"] = np.arange(1, len(mdl) + 1)
    mdl["delta_mdl_bits"] = mdl["mdl_bits"] - mdl["mdl_bits"].min()
    weights = np.exp2(-mdl["delta_mdl_bits"])
    mdl["mdl_weight"] = weights / weights.sum()
    return mdl, pd.DataFrame(prob_rows)


def block_bayesian_bootstrap(df: pd.DataFrame, n_iter: int = 10000, seed: int = 20260701) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    block_mats = []
    block_meta = []
    for block, g in df.groupby("source_layer", sort=False):
        mat = weighted_count_matrix(g, weights=np.ones(len(g)) / len(g))
        if mat.sum() == 0:
            continue
        block_mats.append(mat / mat.sum())
        block_meta.append({"source_layer": block, "n_rows": len(g), "study_id": g["study_id"].iloc[0]})
    mats = np.stack(block_mats, axis=0)
    maps = route_maps()
    pres_idx = [i for i, m in enumerate(maps) if m == PRESPEC][0]
    draws = []
    for i in range(n_iter):
        w = rng.dirichlet(np.ones(len(mats)))
        joint = (mats * w[:, None, None]).sum(axis=0)
        scores = []
        for mapping in maps:
            vals = [joint[DONORS.index(d), ROLES.index(mapping[d])] for d in DONORS]
            scores.append(float(np.mean(vals)))
        order = np.argsort(scores)[::-1]
        rank = int(np.where(order == pres_idx)[0][0] + 1)
        best_alt = max(v for j, v in enumerate(scores) if j != pres_idx)
        draws.append(
            {
                "draw": i,
                "prespecified_support": scores[pres_idx],
                "best_alternative_support": best_alt,
                "margin_vs_best_alternative": scores[pres_idx] - best_alt,
                "prespecified_rank": rank,
            }
        )
    d = pd.DataFrame(draws)
    summary = pd.DataFrame(
        [
            {
                "bootstrap_design": "source_layer_dirichlet_block_weights",
                "n_iter": n_iter,
                "n_blocks": len(mats),
                "prespecified_rank1_probability": float((d["prespecified_rank"] == 1).mean()),
                "margin_mean": float(d["margin_vs_best_alternative"].mean()),
                "margin_q025": float(d["margin_vs_best_alternative"].quantile(0.025)),
                "margin_q50": float(d["margin_vs_best_alternative"].quantile(0.5)),
                "margin_q975": float(d["margin_vs_best_alternative"].quantile(0.975)),
                "support_mean": float(d["prespecified_support"].mean()),
                "support_q025": float(d["prespecified_support"].quantile(0.025)),
                "support_q975": float(d["prespecified_support"].quantile(0.975)),
                "interpretation": "Dirichlet weights are placed on source-layer blocks, not on individual rows. This estimates rank stability under uncertain block influence.",
            }
        ]
    )
    return summary, d


def plot_outputs(mdl: pd.DataFrame, cv_summary: pd.DataFrame, boot_draws: pd.DataFrame, out_png: Path):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    ax = axes[0, 0]
    top = mdl.head(10).iloc[::-1]
    colors = [COL["blue"] if rt == "prespecified" else COL["light"] for rt in top["route_type"]]
    labels = [
        m.replace("donor_specific_theta__", "donor-specific: ")
        .replace("shared_theta__", "shared: ")
        .replace("donor_independent_global", "donor-independent")
        .replace("diffuse_uniform", "diffuse")
        .replace("saturated_donor_role", "saturated")
        for m in top["model_id"]
    ]
    labels = [x[:42] + ("..." if len(x) > 42 else "") for x in labels]
    ax.barh(np.arange(len(top)), top["delta_mdl_bits"], color=colors, edgecolor="white")
    ax.set_yticks(np.arange(len(top)), labels)
    ax.set_xlabel("delta MDL (bits; lower is better)")
    ax.text(0.02, 0.96, "a", transform=ax.transAxes, weight="bold", fontsize=11, va="top")

    ax = axes[0, 1]
    show = cv_summary[cv_summary["block_design"].isin(["leave_study_id", "leave_source_layer"])].copy()
    keep = show["model_id"].isin(
        [
            "diffuse_uniform",
            "donor_independent_global",
            "saturated_donor_role",
            "shared_theta__host->scaffold; alpha->energy; other->transition",
            "donor_specific_theta__host->scaffold; alpha->energy; other->transition",
        ]
    )
    show = show[keep].copy()
    name_map = {
        "diffuse_uniform": "diffuse",
        "donor_independent_global": "donor-independent",
        "saturated_donor_role": "saturated",
        "shared_theta__host->scaffold; alpha->energy; other->transition": "route shared",
        "donor_specific_theta__host->scaffold; alpha->energy; other->transition": "route donor-specific",
    }
    show["short"] = show["model_id"].map(name_map)
    pivot = show.pivot(index="short", columns="block_design", values="mean_heldout_log_loss_bits")
    order = ["diffuse", "donor-independent", "route shared", "route donor-specific", "saturated"]
    pivot = pivot.reindex(order)
    x = np.arange(len(pivot))
    width = 0.36
    ax.bar(x - width / 2, pivot["leave_study_id"], width, color=COL["grey"], label="leave source")
    ax.bar(x + width / 2, pivot["leave_source_layer"], width, color=COL["blue"], label="leave source-layer")
    ax.set_xticks(x, pivot.index, rotation=35, ha="right")
    ax.set_ylabel("held-out loss (bits)")
    ax.legend(frameon=False, fontsize=6)
    ax.text(0.02, 0.96, "b", transform=ax.transAxes, weight="bold", fontsize=11, va="top")

    ax = axes[1, 0]
    ax.hist(boot_draws["margin_vs_best_alternative"], bins=45, color=COL["blue"], alpha=0.85)
    ax.axvline(0, color=COL["red"], lw=1.0)
    q025, q50, q975 = boot_draws["margin_vs_best_alternative"].quantile([0.025, 0.5, 0.975])
    ax.axvline(q50, color="black", lw=1.0)
    ax.axvspan(q025, q975, color=COL["blue"], alpha=0.14, lw=0)
    ax.set_xlabel("block-weight route margin")
    ax.set_ylabel("Dirichlet block draws")
    ax.text(0.02, 0.96, "c", transform=ax.transAxes, weight="bold", fontsize=11, va="top")

    ax = axes[1, 1]
    subset = cv_summary[cv_summary["model_id"].str.contains("host->scaffold; alpha->energy; other->transition", regex=False)].copy()
    subset = subset[subset["block_design"].eq("leave_source_layer")]
    subset["short"] = subset["model_family"].map(
        {
            "route_shared_theta": "shared theta",
            "route_donor_specific_theta": "donor-specific theta",
        }
    )
    vals = subset[["short", "mean_heldout_log_loss_bits"]].dropna()
    ax.scatter(vals["mean_heldout_log_loss_bits"], vals["short"], s=55, color=COL["blue"], zorder=3)
    baselines = cv_summary[cv_summary["block_design"].eq("leave_source_layer")]
    for model, color in [("diffuse_uniform", COL["grey"]), ("donor_independent_global", COL["red"]), ("saturated_donor_role", COL["green"])]:
        y = baselines.loc[baselines["model_id"].eq(model), "mean_heldout_log_loss_bits"]
        if len(y):
            ax.axvline(float(y.iloc[0]), color=color, lw=1.1, alpha=0.85, label=model.replace("_", " "))
    ax.set_xlabel("leave-source-layer loss (bits)")
    ax.set_ylabel("")
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    ax.text(0.02, 0.96, "d", transform=ax.transAxes, weight="bold", fontsize=11, va="top")

    fig.tight_layout(w_pad=1.2, h_pad=1.4)
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    cv_source_scores, cv_source_summary = block_cv(df, "study_id")
    cv_layer_scores, cv_layer_summary = block_cv(df, "source_layer")
    cv_scores = pd.concat([cv_source_scores, cv_layer_scores], ignore_index=True)
    cv_summary = pd.concat([cv_source_summary, cv_layer_summary], ignore_index=True)
    mdl, model_probs = full_mdl(df, block_col="source_layer")
    boot_summary, boot_draws = block_bayesian_bootstrap(df, n_iter=10000)

    pres_shared = "shared_theta__host->scaffold; alpha->energy; other->transition"
    pres_ds = "donor_specific_theta__host->scaffold; alpha->energy; other->transition"
    compact = []
    for design in ["leave_study_id", "leave_source_layer"]:
        sub = cv_summary[cv_summary["block_design"].eq(design)]
        base_diffuse = float(sub.loc[sub["model_id"].eq("diffuse_uniform"), "mean_heldout_log_loss_bits"].iloc[0])
        base_ind = float(sub.loc[sub["model_id"].eq("donor_independent_global"), "mean_heldout_log_loss_bits"].iloc[0])
        for mid in [pres_shared, pres_ds, "saturated_donor_role"]:
            if mid in set(sub["model_id"]):
                row = sub[sub["model_id"].eq(mid)].iloc[0]
                compact.append(
                    {
                        "block_design": design,
                        "model_id": mid,
                        "mean_heldout_log_loss_bits": float(row["mean_heldout_log_loss_bits"]),
                        "rank_within_design": int(row["rank_within_design"]),
                        "bits_saved_vs_diffuse": base_diffuse - float(row["mean_heldout_log_loss_bits"]),
                        "bits_saved_vs_donor_independent": base_ind - float(row["mean_heldout_log_loss_bits"]),
                    }
                )
    compact = pd.DataFrame(compact)

    summary = pd.DataFrame(
        [
            {
                "n_rows": len(df),
                "n_sources": df["study_id"].nunique(),
                "n_source_layers": df["source_layer"].nunique(),
                "best_mdl_model": mdl.iloc[0]["model_id"],
                "best_mdl_delta_bits": float(mdl.iloc[0]["delta_mdl_bits"]),
                "prespecified_shared_theta_mdl_rank": int(mdl.loc[mdl["model_id"].eq(pres_shared), "rank_by_mdl"].iloc[0]),
                "prespecified_donor_specific_mdl_rank": int(mdl.loc[mdl["model_id"].eq(pres_ds), "rank_by_mdl"].iloc[0]),
                "prespecified_shared_theta_delta_mdl_bits": float(mdl.loc[mdl["model_id"].eq(pres_shared), "delta_mdl_bits"].iloc[0]),
                "prespecified_donor_specific_delta_mdl_bits": float(mdl.loc[mdl["model_id"].eq(pres_ds), "delta_mdl_bits"].iloc[0]),
                "block_bootstrap_rank1_probability": float(boot_summary["prespecified_rank1_probability"].iloc[0]),
                "block_bootstrap_margin_q025": float(boot_summary["margin_q025"].iloc[0]),
                "block_bootstrap_margin_q50": float(boot_summary["margin_q50"].iloc[0]),
                "block_bootstrap_margin_q975": float(boot_summary["margin_q975"].iloc[0]),
                "interpretation": "Model evidence is evaluated at source and source-layer block levels, with complexity penalty reported as MDL bits. It complements route-support, permutation and clustering diagnostics.",
            }
        ]
    )

    xlsx = OUT / "block_mdl_model_evidence_v1.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        compact.to_excel(writer, sheet_name="compact_results", index=False)
        mdl.to_excel(writer, sheet_name="mdl_table", index=False)
        cv_summary.to_excel(writer, sheet_name="cv_summary", index=False)
        cv_scores.to_excel(writer, sheet_name="cv_block_scores", index=False)
        model_probs.to_excel(writer, sheet_name="model_probabilities", index=False)
        boot_summary.to_excel(writer, sheet_name="block_bootstrap_summary", index=False)
        boot_draws.head(5000).to_excel(writer, sheet_name="block_bootstrap_draws", index=False)

    for name, frame in {
        "summary.csv": summary,
        "compact_results.csv": compact,
        "mdl_table.csv": mdl,
        "cv_summary.csv": cv_summary,
        "cv_block_scores.csv": cv_scores,
        "model_probabilities.csv": model_probs,
        "block_bootstrap_summary.csv": boot_summary,
        "block_bootstrap_draws_sample.csv": boot_draws.head(5000),
    }.items():
        frame.to_csv(OUT / name, index=False)

    plot_outputs(mdl, cv_summary, boot_draws, OUT / "block_mdl_model_evidence_v1.png")

    (OUT / "README_block_mdl_model_evidence_v1.md").write_text(
        "# Block-level MDL and held-out model evidence v1\n\n"
        "This analysis compares donor-role models using conservative source and source-layer blocks rather than treating "
        "all row-level evidence entries as independent observations. It reports leave-source and leave-source-layer "
        "held-out log loss, an MDL-style complexity penalty at the source-layer level, and a source-layer Dirichlet "
        "block-weight bootstrap for prespecified route rank stability.\n\n"
        "The analysis is a model-evidence diagnostic for the route-resolved atlas. It is not a new primary biological "
        "dataset and should be interpreted together with source-aware permutation, dependency-collapse, coding-audit, "
        "Fe-S phylogeny and AF3 structural-interactivity checks.\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
