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
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_evidence_accumulation_forecast_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
DONOR_CODE = {d: i for i, d in enumerate(DONORS)}
ROLE_CODE = {r: i for i, r in enumerate(ROLES)}
ROLE_LABEL = {"scaffold": "scaffold", "energetic_incorporation": "energy", "transition": "transition"}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}
RNG_SEED = 20260701
N_LAYER_ACCUM = 4000

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
    rows = []
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        arr = np.array([ROLE_CODE[mapping[d]] for d in DONORS], dtype=int)
        short = "_".join(ROLE_LABEL[mapping[d]] for d in DONORS)
        label = "; ".join(f"{d}->{ROLE_LABEL[mapping[d]]}" for d in DONORS)
        rows.append((short, label, arr, mapping == PRESPEC))
    return rows


ROUTE_MAPS = route_maps()
PRESPEC_IDX = [i for i, (_, _, _, isp) in enumerate(ROUTE_MAPS) if isp][0]


def load_full_with_audit_overlay() -> pd.DataFrame:
    df = pd.read_csv(FULL_TABLE)
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
    return work.reset_index(drop=True)


def count_tensor(df: pd.DataFrame):
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
    return counts, layers, sources


def route_supports_from_counts(summed_by_source: np.ndarray) -> np.ndarray:
    total_by_source = summed_by_source.sum(axis=(1, 2))
    active = total_by_source > 0
    supports = []
    if not active.any():
        return np.full(len(ROUTE_MAPS), np.nan)
    for _, _, arr, _ in ROUTE_MAPS:
        matched = np.zeros(summed_by_source.shape[0], dtype=float)
        for d in range(len(DONORS)):
            matched += summed_by_source[:, d, arr[d]]
        supports.append(float(np.mean(matched[active] / total_by_source[active])))
    return np.asarray(supports)


def route_margin_from_counts(summed_by_source: np.ndarray):
    supports = route_supports_from_counts(summed_by_source)
    if np.isnan(supports).all():
        return np.nan, np.nan, np.nan, 99
    best_alt = np.nanmax(np.delete(supports, PRESPEC_IDX))
    rank = int(np.sum(supports > supports[PRESPEC_IDX]) + 1)
    return float(supports[PRESPEC_IDX]), float(best_alt), float(supports[PRESPEC_IDX] - best_alt), rank


def layer_accumulation(df: pd.DataFrame):
    rng = np.random.default_rng(RNG_SEED)
    counts, layers, sources = count_tensor(df)
    rows = []
    for iteration in range(1, N_LAYER_ACCUM + 1):
        order = rng.permutation(np.arange(len(layers)))
        active = np.zeros(len(layers), dtype=bool)
        for k, idx in enumerate(order, start=1):
            active[idx] = True
            pres, alt, margin, rank = route_margin_from_counts(counts[active].sum(axis=0))
            rows.append(
                {
                    "iteration": iteration,
                    "n_layers": k,
                    "prespecified_support": pres,
                    "best_alternative_support": alt,
                    "margin_vs_best_alternative": margin,
                    "prespecified_rank": rank,
                }
            )
    draws = pd.DataFrame(rows)
    summary = (
        draws.groupby("n_layers")
        .agg(
            iterations=("iteration", "count"),
            rank1_probability=("prespecified_rank", lambda x: float(np.mean(np.asarray(x) == 1))),
            margin_mean=("margin_vs_best_alternative", "mean"),
            margin_q025=("margin_vs_best_alternative", lambda x: float(np.nanquantile(x, 0.025))),
            margin_q50=("margin_vs_best_alternative", lambda x: float(np.nanquantile(x, 0.5))),
            margin_q975=("margin_vs_best_alternative", lambda x: float(np.nanquantile(x, 0.975))),
        )
        .reset_index()
    )
    return summary, draws.sample(min(20000, len(draws)), random_state=RNG_SEED), len(layers)


def current_observed(df: pd.DataFrame):
    counts, layers, sources = count_tensor(df)
    pres, alt, margin, rank = route_margin_from_counts(counts.sum(axis=0))
    return {
        "n_sources": len(sources),
        "n_layers": len(layers),
        "prespecified_support": pres,
        "best_alternative_support": alt,
        "margin_vs_best_alternative": margin,
        "prespecified_rank": rank,
    }


def future_adversarial_boundary(obs: dict):
    n_sources = obs["n_sources"]
    margin = obs["margin_vs_best_alternative"]
    rows = []
    for new_sources in range(1, 11):
        required_new_source_margin = -n_sources * margin / new_sources
        feasible = required_new_source_margin >= -1
        final_if_fully_adversarial = (n_sources * margin - new_sources) / (n_sources + new_sources)
        rows.append(
            {
                "new_equal_weight_sources": new_sources,
                "required_mean_new_source_margin_to_flip": required_new_source_margin,
                "feasible_if_each_new_source_margin_ge_minus1": feasible,
                "final_margin_if_each_new_source_fully_adversarial": final_if_fully_adversarial,
                "interpretation": (
                    "A new source margin is prespecified support minus nearest-alternative support. "
                    "The most adversarial possible value is -1."
                ),
            }
        )
    return pd.DataFrame(rows)


def make_figure(accum: pd.DataFrame, future: pd.DataFrame, obs: dict):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), constrained_layout=True)

    ax = axs[0, 0]
    ax.fill_between(accum["n_layers"], accum["margin_q025"], accum["margin_q975"], color="#DCE8F2", alpha=0.9)
    ax.plot(accum["n_layers"], accum["margin_q50"], color=COL["blue"], lw=1.8)
    ax.axhline(0, color="#111827", lw=0.8)
    ax.set_xlabel("evidence layers accumulated")
    ax.set_ylabel("margin vs nearest alternative")
    ax.text(-0.16, 1.05, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[0, 1]
    ax.plot(accum["n_layers"], accum["rank1_probability"], color=COL["green"], lw=1.8)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("evidence layers accumulated")
    ax.set_ylabel("rank-1 probability")
    ax.text(-0.16, 1.05, "b", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 0]
    ax.plot(
        future["new_equal_weight_sources"],
        future["required_mean_new_source_margin_to_flip"],
        color=COL["red"],
        marker="o",
        lw=1.4,
        ms=3.5,
    )
    ax.axhline(-1, color="#111827", lw=0.8, ls="--")
    ax.set_xlabel("new equal-weight sources")
    ax.set_ylabel("required new-source margin")
    ax.text(0.03, 0.10, "below -1 is impossible", transform=ax.transAxes, fontsize=6.7, color=COL["grey"])
    ax.text(-0.16, 1.05, "c", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 1]
    ax.plot(
        future["new_equal_weight_sources"],
        future["final_margin_if_each_new_source_fully_adversarial"],
        color=COL["purple"],
        marker="o",
        lw=1.4,
        ms=3.5,
    )
    ax.axhline(0, color="#111827", lw=0.8)
    ax.set_xlabel("new fully adversarial sources")
    ax.set_ylabel("projected final margin")
    ax.text(
        0.03,
        0.10,
        f"current margin={obs['margin_vs_best_alternative']:.3f}\ncurrent sources={obs['n_sources']}",
        transform=ax.transAxes,
        fontsize=6.7,
        color=COL["grey"],
    )
    ax.text(-0.16, 1.05, "d", transform=ax.transAxes, fontsize=12, weight="bold")

    fig.savefig(OUT / "evidence_accumulation_forecast_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "evidence_accumulation_forecast_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "evidence_accumulation_forecast_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    obs = current_observed(df)
    accum_summary, accum_draws, n_layers = layer_accumulation(df)
    future = future_adversarial_boundary(obs)
    obs_df = pd.DataFrame([obs])
    with pd.ExcelWriter(OUT / "evidence_accumulation_forecast_v1.xlsx", engine="openpyxl") as writer:
        obs_df.to_excel(writer, sheet_name="observed", index=False)
        accum_summary.to_excel(writer, sheet_name="layer_accumulation_summary", index=False)
        accum_draws.to_excel(writer, sheet_name="layer_accumulation_draws", index=False)
        future.to_excel(writer, sheet_name="future_adversarial_boundary", index=False)
    obs_df.to_csv(OUT / "observed.csv", index=False, encoding="utf-8-sig")
    accum_summary.to_csv(OUT / "layer_accumulation_summary.csv", index=False, encoding="utf-8-sig")
    accum_draws.to_csv(OUT / "layer_accumulation_draws_sample.csv", index=False, encoding="utf-8-sig")
    future.to_csv(OUT / "future_adversarial_boundary.csv", index=False, encoding="utf-8-sig")
    (OUT / "README_evidence_accumulation_forecast_v1.md").write_text(
        "# Evidence accumulation and future-source stress forecast\n\n"
        "This analysis asks how route-margin stability changes as evidence layers accumulate, "
        "and how many future equal-weight sources with adversarial donor-role composition would "
        "be needed to overturn the current source-equal route margin. It is a planning and "
        "sensitivity analysis, not a claim about the probability distribution of future studies.\n",
        encoding="utf-8",
    )
    make_figure(accum_summary, future, obs)
    print("Wrote", OUT)
    print(obs_df.to_string(index=False))
    print(accum_summary.tail(5).to_string(index=False))
    print(future.to_string(index=False))


if __name__ == "__main__":
    main()
