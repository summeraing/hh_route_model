from __future__ import annotations

from itertools import permutations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


REPO_ROOT = Path(__file__).resolve().parents[2]
FULL_TABLE = REPO_ROOT / "data" / "evidence_atlas" / "evidence_units_flat.csv"
CURRENT_SOURCE_XLSX = REPO_ROOT / "data" / "evidence_atlas" / "Source_Data_route_resolved_evidence_atlas_complete.xlsx"
OUT = REPO_ROOT / "outputs" / "routeB_unsupervised_route_recovery_v1"
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
    "scaffold": "#356B8A",
    "energetic_incorporation": "#E0942A",
    "transition": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
}
RNG = np.random.default_rng(20260701)


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

    for col in [
        "study_id",
        "evidence_layer",
        "evidence_unit",
        "unit_id",
        "unit_label",
        "module_family",
        "compartment",
        "source_table",
        "notes",
    ]:
        if col not in df.columns:
            df[col] = ""
    for col in ["host_score", "symbiont_score", "other_score", "dominance_margin", "alpha_localization"]:
        if col not in df.columns:
            df[col] = 0.0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    df["source_layer"] = df["study_id"].astype(str) + "::" + df["evidence_layer"].astype(str)
    df["text_features"] = (
        df[
            [
                "evidence_layer",
                "evidence_unit",
                "unit_id",
                "unit_label",
                "module_family",
                "compartment",
                "source_table",
                "notes",
            ]
        ]
        .fillna("")
        .astype(str)
        .agg(" | ".join, axis=1)
    )
    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    return work.reset_index(drop=True)


def source_equal_joint(df: pd.DataFrame, role_col: str = "functional_final") -> pd.DataFrame:
    mats = []
    for _, g in df.groupby("study_id", sort=False):
        tab = pd.crosstab(g["donor_final"], g[role_col]).reindex(index=DONORS, columns=ROLES, fill_value=0)
        if tab.values.sum() == 0:
            continue
        mats.append(tab / tab.values.sum())
    return sum(mats) / len(mats)


def source_equal_joint_from_arrays(donors: np.ndarray, roles: np.ndarray, sources: np.ndarray) -> pd.DataFrame:
    mats = []
    for source in pd.unique(sources):
        idx = np.where(sources == source)[0]
        if len(idx) == 0:
            continue
        tab = pd.crosstab(pd.Series(donors[idx]), pd.Series(roles[idx])).reindex(index=DONORS, columns=ROLES, fill_value=0)
        if tab.values.sum() == 0:
            continue
        mats.append(tab / tab.values.sum())
    return sum(mats) / len(mats)


def standardized_residuals(joint: pd.DataFrame) -> pd.DataFrame:
    p = joint.to_numpy(float)
    expected = p.sum(axis=1, keepdims=True) @ p.sum(axis=0, keepdims=True)
    z = (p - expected) / np.sqrt(np.clip(expected, 1e-12, None))
    return pd.DataFrame(z, index=joint.index, columns=joint.columns)


def route_maps() -> list[dict[str, str]]:
    return [dict(zip(DONORS, perm)) for perm in permutations(ROLES)]


def map_name(mapping: dict[str, str]) -> str:
    return "; ".join(f"{DONOR_LABEL[d]}->{ROLE_SHORT[mapping[d]]}" for d in DONORS)


def residual_route_table(z: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mapping in route_maps():
        vals = [float(z.loc[d, mapping[d]]) for d in DONORS]
        rows.append(
            {
                "route_map": map_name(mapping),
                "is_prespecified": mapping == PRESPEC,
                "mean_standardized_residual": float(np.mean(vals)),
                "minimum_component_residual": float(np.min(vals)),
                "component_residuals": "; ".join(f"{DONOR_LABEL[d]}:{z.loc[d, mapping[d]]:.3f}" for d in DONORS),
            }
        )
    out = pd.DataFrame(rows).sort_values("mean_standardized_residual", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def correspondence_coordinates(joint: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    p = joint.to_numpy(float)
    r = p.sum(axis=1)
    c = p.sum(axis=0)
    expected = np.outer(r, c)
    s = (p - expected) / np.sqrt(np.outer(r, c))
    u, singular, vt = np.linalg.svd(s, full_matrices=False)
    row_coord = (u[:, :2] * singular[:2]) / np.sqrt(r[:, None])
    col_coord = (vt.T[:, :2] * singular[:2]) / np.sqrt(c[:, None])
    rows = pd.DataFrame(row_coord, index=DONORS, columns=["CA1", "CA2"]).reset_index(names="node")
    rows["node_type"] = "donor"
    cols = pd.DataFrame(col_coord, index=ROLES, columns=["CA1", "CA2"]).reset_index(names="node")
    cols["node_type"] = "role"
    eig = singular**2
    return rows, cols, eig / eig.sum()


def modularity_for_partition(joint: pd.DataFrame, mapping: dict[str, str]) -> float:
    nodes = DONORS + ROLES
    idx = {n: i for i, n in enumerate(nodes)}
    a = np.zeros((len(nodes), len(nodes)), dtype=float)
    for d in DONORS:
        for f in ROLES:
            w = float(joint.loc[d, f])
            a[idx[d], idx[f]] = w
            a[idx[f], idx[d]] = w
    m2 = a.sum()
    k = a.sum(axis=1)
    communities = {}
    for i, d in enumerate(DONORS):
        communities[d] = i
        communities[mapping[d]] = i
    q = 0.0
    for i, ni in enumerate(nodes):
        for j, nj in enumerate(nodes):
            if communities[ni] == communities[nj]:
                q += a[i, j] - k[i] * k[j] / m2
    return float(q / m2)


def modularity_route_table(joint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mapping in route_maps():
        rows.append(
            {
                "route_map": map_name(mapping),
                "is_prespecified": mapping == PRESPEC,
                "weighted_modularity": modularity_for_partition(joint, mapping),
            }
        )
    out = pd.DataFrame(rows).sort_values("weighted_modularity", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def block_indices(df: pd.DataFrame, block: str) -> list[np.ndarray]:
    if block == "source":
        groups = df.groupby("study_id").indices
    elif block == "source_layer":
        groups = df.groupby(["study_id", "evidence_layer"]).indices
    else:
        raise ValueError(block)
    return [np.asarray(v, dtype=int) for v in groups.values()]


def permute_array_within_blocks(values: np.ndarray, groups: list[np.ndarray]) -> np.ndarray:
    shuffled = values.copy()
    for idx in groups:
        vals = shuffled[idx].copy()
        RNG.shuffle(vals)
        shuffled[idx] = vals
    return shuffled


def permutation_tests(df: pd.DataFrame, n_perm: int = 1000) -> pd.DataFrame:
    observed_joint = source_equal_joint(df)
    observed_z = standardized_residuals(observed_joint)
    observed_resid = residual_route_table(observed_z)
    observed_mod = modularity_route_table(observed_joint)
    pres_resid = float(observed_resid.loc[observed_resid["is_prespecified"], "mean_standardized_residual"].iloc[0])
    best_alt_resid = float(observed_resid.loc[~observed_resid["is_prespecified"], "mean_standardized_residual"].max())
    resid_margin = pres_resid - best_alt_resid
    pres_mod = float(observed_mod.loc[observed_mod["is_prespecified"], "weighted_modularity"].iloc[0])
    best_alt_mod = float(observed_mod.loc[~observed_mod["is_prespecified"], "weighted_modularity"].max())
    mod_margin = pres_mod - best_alt_mod

    rows = []
    donors = df["donor_final"].to_numpy()
    roles = df["functional_final"].to_numpy()
    sources = df["study_id"].to_numpy()
    group_map = {"source": block_indices(df, "source"), "source_layer": block_indices(df, "source_layer")}
    for design in ["source", "source_layer"]:
        for i in range(n_perm):
            role_perm = permute_array_within_blocks(roles, group_map[design])
            joint = source_equal_joint_from_arrays(donors, role_perm, sources)
            z = standardized_residuals(joint)
            rt = residual_route_table(z)
            mt = modularity_route_table(joint)
            pp = float(rt.loc[rt["is_prespecified"], "mean_standardized_residual"].iloc[0])
            ba = float(rt.loc[~rt["is_prespecified"], "mean_standardized_residual"].max())
            qm = float(mt.loc[mt["is_prespecified"], "weighted_modularity"].iloc[0])
            qba = float(mt.loc[~mt["is_prespecified"], "weighted_modularity"].max())
            rows.append(
                {
                    "design": design,
                    "iteration": i,
                    "prespecified_residual_score": pp,
                    "residual_margin_vs_best_alt": pp - ba,
                    "prespecified_modularity": qm,
                    "modularity_margin_vs_best_alt": qm - qba,
                }
            )
    draws = pd.DataFrame(rows)
    summary_rows = []
    for design, sub in draws.groupby("design"):
        for metric, obs in [
            ("prespecified_residual_score", pres_resid),
            ("residual_margin_vs_best_alt", resid_margin),
            ("prespecified_modularity", pres_mod),
            ("modularity_margin_vs_best_alt", mod_margin),
        ]:
            vals = sub[metric].to_numpy(float)
            summary_rows.append(
                {
                    "design": design,
                    "metric": metric,
                    "observed": obs,
                    "null_mean": float(np.mean(vals)),
                    "null_q025": float(np.quantile(vals, 0.025)),
                    "null_q50": float(np.quantile(vals, 0.5)),
                    "null_q975": float(np.quantile(vals, 0.975)),
                    "empirical_p_greater_equal": float((np.sum(vals >= obs) + 1) / (len(vals) + 1)),
                }
            )
    return draws, pd.DataFrame(summary_rows)


def best_cluster_accuracy(true_labels: np.ndarray, cluster_labels: np.ndarray) -> float:
    true_codes = pd.Categorical(true_labels, categories=ROLES).codes
    clusters = np.unique(cluster_labels)
    cost = np.zeros((len(ROLES), len(clusters)), dtype=int)
    for i in range(len(ROLES)):
        for j, cl in enumerate(clusters):
            cost[i, j] = -np.sum((true_codes == i) & (cluster_labels == cl))
    row_ind, col_ind = linear_sum_assignment(cost)
    correct = -cost[row_ind, col_ind].sum()
    return float(correct / len(true_labels))


def simple_kmeans(x: np.ndarray, n_clusters: int = 3, n_init: int = 60, max_iter: int = 200) -> np.ndarray:
    """Small deterministic k-means to avoid platform-specific sklearn/threadpool issues."""
    x = np.asarray(x, dtype=float)
    x = (x - x.mean(axis=0, keepdims=True)) / np.clip(x.std(axis=0, keepdims=True), 1e-8, None)
    best_labels = None
    best_inertia = np.inf
    for _ in range(n_init):
        centers = x[RNG.choice(len(x), size=n_clusters, replace=False)].copy()
        labels = np.zeros(len(x), dtype=int)
        for _it in range(max_iter):
            dist = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
            new_labels = dist.argmin(axis=1)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for k in range(n_clusters):
                if np.any(labels == k):
                    centers[k] = x[labels == k].mean(axis=0)
                else:
                    centers[k] = x[RNG.integers(0, len(x))]
        inertia = float(((x - centers[labels]) ** 2).sum())
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels.copy()
    return best_labels


def unsupervised_cluster_diagnostics(df: pd.DataFrame, n_perm: int = 500) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_sets = {
        "text_numeric_no_donor": [
            ("text", TfidfVectorizer(lowercase=True, ngram_range=(1, 2), min_df=2, max_features=1200, token_pattern=r"(?u)\b[\w:/.-]+\b"), "text_features"),
            ("numeric", StandardScaler(), ["host_score", "symbiont_score", "other_score", "dominance_margin", "alpha_localization"]),
        ],
        "text_numeric_source_no_donor": [
            ("text", TfidfVectorizer(lowercase=True, ngram_range=(1, 2), min_df=2, max_features=1200, token_pattern=r"(?u)\b[\w:/.-]+\b"), "text_features"),
            ("numeric", StandardScaler(), ["host_score", "symbiont_score", "other_score", "dominance_margin", "alpha_localization"]),
            ("source_layer", OneHotEncoder(handle_unknown="ignore"), ["study_id", "evidence_layer"]),
        ],
        "text_numeric_source_donor": [
            ("text", TfidfVectorizer(lowercase=True, ngram_range=(1, 2), min_df=2, max_features=1200, token_pattern=r"(?u)\b[\w:/.-]+\b"), "text_features"),
            ("numeric", StandardScaler(), ["host_score", "symbiont_score", "other_score", "dominance_margin", "alpha_localization"]),
            ("source_layer", OneHotEncoder(handle_unknown="ignore"), ["study_id", "evidence_layer"]),
            ("donor", OneHotEncoder(handle_unknown="ignore"), ["donor_final"]),
        ],
    }
    obs_rows = []
    assign_tables = []
    perm_rows = []
    y = df["functional_final"].to_numpy()
    groups = block_indices(df, "source_layer")
    for feature_name, transformers in feature_sets.items():
        pre = ColumnTransformer(transformers=transformers, sparse_threshold=0.3)
        pipe = Pipeline([("pre", pre), ("svd", TruncatedSVD(n_components=10, random_state=20260701))])
        emb = pipe.fit_transform(df)
        clusters = simple_kmeans(emb, n_clusters=3, n_init=60)
        nmi = normalized_mutual_info_score(y, clusters)
        ari = adjusted_rand_score(y, clusters)
        acc = best_cluster_accuracy(y, clusters)
        obs_rows.append(
            {
                "feature_set": feature_name,
                "n_rows": len(df),
                "n_clusters": 3,
                "role_nmi": nmi,
                "role_adjusted_rand": ari,
                "best_matched_role_accuracy": acc,
            }
        )
        assign = df[["row_index", "study_id", "evidence_layer", "donor_final", "functional_final", "unit_label"]].copy()
        assign["feature_set"] = feature_name
        assign["cluster"] = clusters
        assign_tables.append(assign)

        for i in range(n_perm):
            yp = permute_array_within_blocks(y, groups)
            perm_rows.append(
                {
                    "feature_set": feature_name,
                    "iteration": i,
                    "role_nmi": normalized_mutual_info_score(yp, clusters),
                    "role_adjusted_rand": adjusted_rand_score(yp, clusters),
                    "best_matched_role_accuracy": best_cluster_accuracy(yp, clusters),
                }
            )
    obs = pd.DataFrame(obs_rows)
    perms = pd.DataFrame(perm_rows)
    summaries = []
    for feature_name, sub in perms.groupby("feature_set"):
        obs_sub = obs[obs["feature_set"] == feature_name].iloc[0]
        for metric in ["role_nmi", "role_adjusted_rand", "best_matched_role_accuracy"]:
            vals = sub[metric].to_numpy(float)
            summaries.append(
                {
                    "feature_set": feature_name,
                    "metric": metric,
                    "observed": float(obs_sub[metric]),
                    "null_mean": float(np.mean(vals)),
                    "null_q025": float(np.quantile(vals, 0.025)),
                    "null_q50": float(np.quantile(vals, 0.5)),
                    "null_q975": float(np.quantile(vals, 0.975)),
                    "empirical_p_greater_equal": float((np.sum(vals >= obs_sub[metric]) + 1) / (len(vals) + 1)),
                }
            )
    return obs, pd.DataFrame(summaries), pd.concat(assign_tables, ignore_index=True)


def plot_results(
    joint: pd.DataFrame,
    coords_d: pd.DataFrame,
    coords_r: pd.DataFrame,
    inertia: np.ndarray,
    residual_routes: pd.DataFrame,
    modularity_routes: pd.DataFrame,
    perm_summary: pd.DataFrame,
    cluster_obs: pd.DataFrame,
    cluster_perm_summary: pd.DataFrame,
):
    fig = plt.figure(figsize=(180 / 25.4, 170 / 25.4), dpi=320)
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, top=0.94, bottom=0.10, wspace=0.42, hspace=0.42)

    ax = fig.add_subplot(gs[0, 0])
    label_offsets = {
        "host": (-0.13, 0.06),
        "symbiont": (0.00, 0.12),
        "other": (0.10, 0.07),
        "scaffold": (-0.13, -0.13),
        "energetic_incorporation": (-0.10, 0.02),
        "transition": (-0.12, 0.03),
    }
    for _, row in coords_d.iterrows():
        ax.scatter(row["CA1"], row["CA2"], s=72, c=COL[row["node"]], edgecolor="white", linewidth=0.8, zorder=3)
        dx, dy = label_offsets[row["node"]]
        ax.text(row["CA1"] + dx, row["CA2"] + dy, DONOR_LABEL[row["node"]], ha="center", va="center", fontweight="bold")
    for _, row in coords_r.iterrows():
        ax.scatter(row["CA1"], row["CA2"], s=78, marker="s", c=COL[row["node"]], edgecolor="white", linewidth=0.8, alpha=0.80, zorder=3)
        dx, dy = label_offsets[row["node"]]
        ax.text(row["CA1"] + dx, row["CA2"] + dy, ROLE_SHORT[row["node"]], ha="center", va="center", fontweight="bold")
    for d, f in PRESPEC.items():
        a = coords_d[coords_d["node"] == d].iloc[0]
        b = coords_r[coords_r["node"] == f].iloc[0]
        ax.plot([a["CA1"], b["CA1"]], [a["CA2"], b["CA2"]], c=COL[d], lw=1.0, alpha=0.65)
    ax.axhline(0, color="#D8DEE6", lw=0.6)
    ax.axvline(0, color="#D8DEE6", lw=0.6)
    ax.margins(x=0.22, y=0.24)
    ax.set_xlabel(f"CA1 ({inertia[0] * 100:.1f}% inertia)")
    ax.set_ylabel(f"CA2 ({inertia[1] * 100:.1f}% inertia)")
    ax.text(-0.10, 1.05, "a", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = fig.add_subplot(gs[0, 1])
    z = standardized_residuals(joint)
    im = ax.imshow(z.loc[DONORS, ROLES].to_numpy(float), cmap="RdBu_r", vmin=-np.max(np.abs(z.values)), vmax=np.max(np.abs(z.values)))
    ax.set_xticks(np.arange(len(ROLES)), [ROLE_SHORT[r] for r in ROLES], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(DONORS)), [DONOR_LABEL[d] for d in DONORS])
    for i, d in enumerate(DONORS):
        for j, f in enumerate(ROLES):
            val = z.loc[d, f]
            weight = "bold" if PRESPEC[d] == f else "normal"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, fontweight=weight)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("standardized residual")
    cb.ax.tick_params(labelsize=6)
    ax.text(-0.10, 1.05, "b", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = fig.add_subplot(gs[1, 0])
    metrics = [
        ("residual_margin_vs_best_alt", "residual margin"),
        ("modularity_margin_vs_best_alt", "modularity margin"),
    ]
    x = np.arange(len(metrics))
    obs = []
    low = []
    high = []
    for metric, _ in metrics:
        row = perm_summary[(perm_summary["design"] == "source_layer") & (perm_summary["metric"] == metric)].iloc[0]
        obs.append(row["observed"])
        low.append(row["null_q025"])
        high.append(row["null_q975"])
    ax.bar(x, obs, width=0.45, color=[COL["host"], COL["green"]], edgecolor="white", linewidth=0.8)
    for xi, lo, hi in zip(x, low, high):
        ax.plot([xi - 0.28, xi + 0.28], [hi, hi], color=COL["red"], lw=1.2)
        ax.fill_between([xi - 0.28, xi + 0.28], [lo, lo], [hi, hi], color=COL["red"], alpha=0.10, linewidth=0)
    ax.axhline(0, color="#2F3742", lw=0.6)
    ax.set_xticks(x, [label for _, label in metrics])
    ax.set_ylabel("observed value; red = source-layer null q2.5-q97.5")
    ax.text(-0.10, 1.05, "c", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = fig.add_subplot(gs[1, 1])
    summary = cluster_perm_summary[cluster_perm_summary["metric"].isin(["role_nmi", "best_matched_role_accuracy"])].copy()
    order = ["text_numeric_no_donor", "text_numeric_source_no_donor", "text_numeric_source_donor"]
    labels = ["text+scores", "text+scores+source", "text+scores+source+donor"]
    width = 0.34
    for k, metric in enumerate(["role_nmi", "best_matched_role_accuracy"]):
        vals = []
        q975 = []
        for feat in order:
            row = summary[(summary["feature_set"] == feat) & (summary["metric"] == metric)].iloc[0]
            vals.append(row["observed"])
            q975.append(row["null_q975"])
        pos = np.arange(len(order)) + (k - 0.5) * width
        ax.bar(pos, vals, width=width, color=[COL["host"], COL["symbiont"]][k], alpha=0.82, edgecolor="white", label="NMI" if metric == "role_nmi" else "matched accuracy")
        for p, q in zip(pos, q975):
            ax.plot([p - width * 0.35, p + width * 0.35], [q, q], color=COL["red"], lw=1.0)
    ax.set_xticks(np.arange(len(order)), labels, rotation=25, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("cluster-role agreement")
    ax.legend(frameon=False, loc="upper left", fontsize=6)
    ax.text(-0.10, 1.05, "d", transform=ax.transAxes, fontsize=12, fontweight="bold")

    fig.savefig(OUT / "unsupervised_route_recovery_v1.png", dpi=450)
    fig.savefig(OUT / "unsupervised_route_recovery_v1.pdf")
    fig.savefig(OUT / "unsupervised_route_recovery_v1.svg")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    joint = source_equal_joint(df)
    z = standardized_residuals(joint)
    residual_routes = residual_route_table(z)
    modularity_routes = modularity_route_table(joint)
    coords_d, coords_r, inertia = correspondence_coordinates(joint)
    perm_draws, perm_summary = permutation_tests(df, n_perm=1000)
    cluster_obs, cluster_perm_summary, cluster_assignments = unsupervised_cluster_diagnostics(df, n_perm=500)

    observed_summary = pd.DataFrame(
        [
            {
                "n_rows": len(df),
                "n_sources": df["study_id"].nunique(),
                "n_layers": df["evidence_layer"].nunique(),
                "prespecified_residual_rank": int(residual_routes.loc[residual_routes["is_prespecified"], "rank"].iloc[0]),
                "prespecified_residual_score": float(residual_routes.loc[residual_routes["is_prespecified"], "mean_standardized_residual"].iloc[0]),
                "prespecified_modularity_rank": int(modularity_routes.loc[modularity_routes["is_prespecified"], "rank"].iloc[0]),
                "prespecified_modularity": float(modularity_routes.loc[modularity_routes["is_prespecified"], "weighted_modularity"].iloc[0]),
            }
        ]
    )

    outputs = {
        "unsupervised_observed_summary.csv": observed_summary,
        "source_equal_joint_matrix.csv": joint.reset_index(names="donor_class"),
        "standardized_residual_matrix.csv": z.reset_index(names="donor_class"),
        "correspondence_donor_coordinates.csv": coords_d,
        "correspondence_role_coordinates.csv": coords_r,
        "correspondence_inertia.csv": pd.DataFrame({"axis": np.arange(1, len(inertia) + 1), "fraction_inertia": inertia}),
        "residual_route_scores.csv": residual_routes,
        "modularity_route_scores.csv": modularity_routes,
        "permutation_summary.csv": perm_summary,
        "permutation_draws_sample.csv": perm_draws.groupby("design", group_keys=False).head(1000),
        "cluster_observed.csv": cluster_obs,
        "cluster_permutation_summary.csv": cluster_perm_summary,
        "cluster_assignments.csv": cluster_assignments,
    }
    for name, table in outputs.items():
        table.to_csv(OUT / name, index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(OUT / "unsupervised_route_recovery_v1.xlsx", engine="openpyxl") as writer:
        for name, table in outputs.items():
            sheet = name.replace(".csv", "")[:31]
            table.to_excel(writer, sheet_name=sheet, index=False)

    plot_results(joint, coords_d, coords_r, inertia, residual_routes, modularity_routes, perm_summary, cluster_obs, cluster_perm_summary)

    readme = (
        "# Route-B unsupervised route-recovery diagnostics\n\n"
        "This analysis asks whether the donor-role architecture is visible without using the prespecified route score as the only diagnostic. "
        "It computes source-equal donor-role correspondence residuals, weighted donor-role graph modularity and unsupervised text/numeric cluster recovery. "
        "Permutation nulls shuffle functional roles within source or source-layer blocks, preserving public-source structure while breaking donor-role pairing.\n\n"
        "Interpretation: these diagnostics are not new biological proof. They are route-independent structure checks showing whether the same organization is recoverable as geometry, graph modularity and feature-space clustering.\n"
    )
    (OUT / "README_unsupervised_route_recovery_v1.md").write_text(readme, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
