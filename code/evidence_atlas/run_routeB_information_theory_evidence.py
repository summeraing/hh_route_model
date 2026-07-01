from __future__ import annotations

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
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_information_theory_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
DONOR_CODE = {d: i for i, d in enumerate(DONORS)}
ROLE_CODE = {r: i for i, r in enumerate(ROLES)}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}
N_PERM = 5000
RNG_SEED = 20260701

COL = {
    "blue": "#356B8A",
    "orange": "#E0942A",
    "purple": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
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

    audit = pd.read_excel(SOURCE_XLSX, sheet_name="audit_overlay")
    for _, row in audit.iterrows():
        idx = int(row["combined_index"])
        if 0 <= idx < len(df):
            df.loc[idx, "donor_final"] = norm(row["adjudicated_donor_class"])
            df.loc[idx, "functional_final"] = norm(row["adjudicated_functional_class"])

    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    work["donor_code"] = work["donor_final"].map(DONOR_CODE).astype(int)
    work["role_code"] = work["functional_final"].map(ROLE_CODE).astype(int)
    work["source_layer"] = work["study_id"].astype(str) + "::" + work["evidence_layer"].astype(str)
    return work.reset_index(drop=True)


def entropy(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.nan
    p = counts[counts > 0] / total
    return float(-(p * np.log2(p)).sum())


def mutual_information(donor: np.ndarray, role: np.ndarray) -> float:
    mat = np.zeros((len(DONORS), len(ROLES)), dtype=float)
    for d, r in zip(donor, role):
        mat[int(d), int(r)] += 1
    total = mat.sum()
    if total <= 0:
        return np.nan
    p = mat / total
    pd_ = p.sum(axis=1, keepdims=True)
    pr = p.sum(axis=0, keepdims=True)
    nz = p > 0
    return float((p[nz] * np.log2(p[nz] / (pd_ @ pr)[nz])).sum())


def normalized_mi(donor: np.ndarray, role: np.ndarray) -> float:
    mi = mutual_information(donor, role)
    h = entropy(np.bincount(role.astype(int), minlength=len(ROLES)))
    if not np.isfinite(mi) or not np.isfinite(h) or h <= 0:
        return np.nan
    return float(mi / h)


def block_equal_metric(df: pd.DataFrame, block_col: str, metric_fn) -> float:
    vals = []
    for _, sub in df.groupby(block_col):
        if len(sub) >= 2 and sub["donor_code"].nunique() >= 2 and sub["role_code"].nunique() >= 2:
            vals.append(metric_fn(sub["donor_code"].to_numpy(), sub["role_code"].to_numpy()))
    return float(np.nanmean(vals)) if vals else np.nan


def route_support(df: pd.DataFrame) -> float:
    mapping = np.array([ROLE_CODE[PRESPEC[d]] for d in DONORS], dtype=int)
    per_source = []
    for _, sub in df.groupby("study_id"):
        per_source.append(float(np.mean(mapping[sub["donor_code"].to_numpy()] == sub["role_code"].to_numpy())))
    return float(np.mean(per_source))


def observed_metrics(df: pd.DataFrame) -> pd.DataFrame:
    donor = df["donor_code"].to_numpy()
    role = df["role_code"].to_numpy()
    rows = [
        {
            "metric": "raw_mutual_information_bits",
            "value": mutual_information(donor, role),
            "interpretation": "Unconditional donor-role mutual information; descriptive because rows are dependent.",
        },
        {
            "metric": "raw_normalized_mi_by_role_entropy",
            "value": normalized_mi(donor, role),
            "interpretation": "Unconditional MI divided by role entropy.",
        },
        {
            "metric": "source_equal_mi_bits",
            "value": block_equal_metric(df, "study_id", mutual_information),
            "interpretation": "Mean within-source donor-role mutual information.",
        },
        {
            "metric": "source_equal_normalized_mi",
            "value": block_equal_metric(df, "study_id", normalized_mi),
            "interpretation": "Mean within-source normalized donor-role MI.",
        },
        {
            "metric": "source_layer_equal_mi_bits",
            "value": block_equal_metric(df, "source_layer", mutual_information),
            "interpretation": "Mean within-source-layer donor-role MI over non-degenerate blocks.",
        },
        {
            "metric": "source_equal_prespecified_route_support",
            "value": route_support(df),
            "interpretation": "Source-equal prespecified route support, included for comparison.",
        },
    ]
    return pd.DataFrame(rows)


def shuffle_within_blocks(df: pd.DataFrame, block_col: str, rng: np.random.Generator) -> pd.DataFrame:
    out = df.copy()
    shuffled = out["role_code"].to_numpy().copy()
    for _, idx in out.groupby(block_col).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        shuffled[idx] = rng.permutation(shuffled[idx])
    out["role_code"] = shuffled
    return out


def permutation_null(df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    designs = [
        ("source_role_shuffle", "study_id"),
        ("source_layer_role_shuffle", "source_layer"),
    ]
    for design, block_col in designs:
        for i in range(N_PERM):
            shuf = shuffle_within_blocks(df, block_col, rng)
            rows.append(
                {
                    "design": design,
                    "iteration": i,
                    "source_equal_mi_bits": block_equal_metric(shuf, "study_id", mutual_information),
                    "source_equal_normalized_mi": block_equal_metric(shuf, "study_id", normalized_mi),
                    "source_layer_equal_mi_bits": block_equal_metric(shuf, "source_layer", mutual_information),
                    "source_equal_route_support": route_support(shuf),
                }
            )
    return pd.DataFrame(rows)


def summarize_null(obs: pd.DataFrame, nulls: pd.DataFrame) -> pd.DataFrame:
    lookup = obs.set_index("metric")["value"].to_dict()
    metric_map = {
        "source_equal_mi_bits": "source_equal_mi_bits",
        "source_equal_normalized_mi": "source_equal_normalized_mi",
        "source_layer_equal_mi_bits": "source_layer_equal_mi_bits",
        "source_equal_prespecified_route_support": "source_equal_route_support",
    }
    rows = []
    for design, sub in nulls.groupby("design"):
        for obs_key, null_col in metric_map.items():
            obs_val = lookup[obs_key]
            vals = sub[null_col].dropna().to_numpy()
            rows.append(
                {
                    "design": design,
                    "metric": obs_key,
                    "observed": obs_val,
                    "null_mean": float(vals.mean()) if len(vals) else np.nan,
                    "null_q025": float(np.quantile(vals, 0.025)) if len(vals) else np.nan,
                    "null_q975": float(np.quantile(vals, 0.975)) if len(vals) else np.nan,
                    "empirical_p_greater_equal": float((np.sum(vals >= obs_val) + 1) / (len(vals) + 1)) if len(vals) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def make_figure(obs: pd.DataFrame, null_summary: pd.DataFrame, nulls: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2), gridspec_kw={"wspace": 0.42})
    obs_lookup = obs.set_index("metric")["value"].to_dict()

    ax = axes[0]
    sub = nulls[nulls["design"] == "source_role_shuffle"]
    vals = sub["source_equal_mi_bits"].dropna().to_numpy()
    q975 = float(np.quantile(vals, 0.975))
    xmax = max(0.08, q975 * 1.25)
    ax.hist(vals, bins=32, color=COL["light"], edgecolor="white")
    ax.axvline(q975, color=COL["grey"], lw=1.0, ls="--")
    ax.set_xlim(0, xmax)
    ax.annotate(
        f"observed\n{obs_lookup['source_equal_mi_bits']:.3f}",
        xy=(xmax, ax.get_ylim()[1] * 0.78),
        xytext=(xmax * 0.68, ax.get_ylim()[1] * 0.78),
        arrowprops={"arrowstyle": "-|>", "lw": 0.9, "color": COL["blue"]},
        color=COL["blue"],
        fontsize=6.0,
        va="center",
    )
    ax.text(q975, ax.get_ylim()[1] * 0.92, "q97.5", ha="center", fontsize=5.6, color=COL["grey"])
    ax.set_xlabel("null source-equal MI (bits)")
    ax.set_ylabel("null draws")
    ax.text(-0.18, 1.08, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    ax = axes[1]
    sub = null_summary[
        (null_summary["design"] == "source_role_shuffle")
        & (null_summary["metric"].isin(["source_equal_mi_bits", "source_equal_normalized_mi", "source_equal_prespecified_route_support"]))
    ].copy()
    labels = ["MI", "NMI", "route"]
    y = np.arange(len(sub))
    for yi, (_, row) in zip(y, sub.iterrows()):
        ax.hlines(yi, row["null_q025"], row["null_q975"], color=COL["grey"], lw=1.0)
    ax.scatter(sub["null_mean"], y, s=18, color=COL["light"], edgecolor=COL["grey"], zorder=3, label="null mean")
    ax.scatter(sub["observed"], y, s=22, color=COL["blue"], zorder=4, label="observed")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("observed vs null interval")
    ax.text(-0.18, 1.08, "b", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    ax = axes[2]
    sub = null_summary[
        null_summary["metric"].isin(["source_equal_mi_bits", "source_layer_equal_mi_bits", "source_equal_prespecified_route_support"])
    ].copy()
    sub["effect"] = sub["observed"] - sub["null_q975"]
    sub["label"] = sub["design"].map({"source_role_shuffle": "source", "source_layer_role_shuffle": "source+layer"}) + "\n" + sub["metric"].map(
        {
            "source_equal_mi_bits": "MI",
            "source_layer_equal_mi_bits": "layer MI",
            "source_equal_prespecified_route_support": "route",
        }
    )
    colors = [COL["blue"] if v > 0 else COL["red"] for v in sub["effect"]]
    ax.barh(np.arange(len(sub)), sub["effect"], color=colors)
    ax.axvline(0, color=COL["grey"], lw=0.7)
    ax.set_yticks(np.arange(len(sub)), sub["label"])
    ax.invert_yaxis()
    ax.set_xlabel("observed - null q97.5")
    ax.text(-0.18, 1.08, "c", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    fig.savefig(OUT / "Fig_information_theory_evidence.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Fig_information_theory_evidence.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Fig_information_theory_evidence.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    obs_csv = OUT / "information_observed_metrics.csv"
    summary_csv = OUT / "information_permutation_summary.csv"
    nulls_csv = OUT / "information_permutation_draws.csv"
    if obs_csv.exists() and summary_csv.exists() and nulls_csv.exists():
        obs = pd.read_csv(obs_csv)
        summary = pd.read_csv(summary_csv)
        nulls = pd.read_csv(nulls_csv)
    else:
        df = load_full_with_audit_overlay()
        obs = observed_metrics(df)
        nulls = permutation_null(df)
        summary = summarize_null(obs, nulls)

        with pd.ExcelWriter(OUT / "routeB_information_theory_evidence.xlsx", engine="openpyxl") as writer:
            obs.to_excel(writer, sheet_name="observed_metrics", index=False)
            summary.to_excel(writer, sheet_name="permutation_summary", index=False)
            nulls.to_excel(writer, sheet_name="permutation_draws", index=False)

        obs.to_csv(obs_csv, index=False)
        summary.to_csv(summary_csv, index=False)
        nulls.to_csv(nulls_csv, index=False)
    make_figure(obs, summary, nulls)
    (OUT / "README.md").write_text(
        "# Route-B information-theory evidence\n\n"
        "This analysis asks whether donor identity carries information about organizational role after source-aware controls. "
        "The main statistic is within-source donor-role mutual information, supplemented by normalized mutual information, "
        "source-layer MI and prespecified route support. Nulls shuffle role labels within source or source-layer blocks. "
        "The analysis is a dependence diagnostic, not an independent-observation significance test.\n",
        encoding="utf-8",
    )
    print(OUT)
    print(obs.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
