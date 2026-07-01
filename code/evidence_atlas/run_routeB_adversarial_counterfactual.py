from __future__ import annotations

import math
from itertools import permutations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
SOURCE_XLSX = ROOT / "10_WORKING_COPY_FROM_NC12" / "04_Source_Data" / "Source_Data_complete.xlsx"
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_adversarial_counterfactual_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
ROLE_LABELS = ["scaffold", "energy", "transition"]
ROLE_SHORT = dict(zip(ROLES, ROLE_LABELS))
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}

COL = {
    "host": "#356B8A",
    "symbiont": "#E0942A",
    "other": "#7A58A6",
    "grey": "#6B7280",
    "light": "#E8EEF3",
    "line": "#9AA7B3",
    "red": "#B95A55",
    "green": "#3B8D5A",
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
    return str(x).strip().lower().replace("injection", "energetic_incorporation")


def load_matrices():
    mat = pd.read_excel(SOURCE_XLSX, sheet_name="continuous_role_matrix")
    mat["donor_class"] = mat["donor_class"].map(norm)
    mat["functional_class"] = mat["functional_class"].map(norm)
    mat = mat[mat["donor_class"].isin(DONORS) & mat["functional_class"].isin(ROLES)].copy()
    P = mat.pivot(index="donor_class", columns="functional_class", values="source_equal_probability")
    C = mat.pivot(index="donor_class", columns="functional_class", values="raw_role_space_count")
    P = P.reindex(index=DONORS, columns=ROLES).astype(float)
    C = C.reindex(index=DONORS, columns=ROLES).fillna(0).astype(float)
    return P, C


def route_weight(M, mapping):
    return float(np.mean([M.loc[d, mapping[d]] for d in DONORS]))


def route_weights(M):
    rows = []
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        rows.append({
            "route_map": "; ".join([f"{d}->{ROLE_SHORT[mapping[d]]}" for d in DONORS]),
            "host_route": mapping["host"],
            "symbiont_route": mapping["symbiont"],
            "other_route": mapping["other"],
            "is_prespecified": mapping == PRESPEC,
            "weight": route_weight(M, mapping),
        })
    df = pd.DataFrame(rows).sort_values("weight", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


def targeted_shift_curve(P, best_alt):
    pres = route_weight(P, PRESPEC)
    alt = route_weight(P, best_alt)
    margin = pres - alt
    changed = [d for d in DONORS if PRESPEC[d] != best_alt[d]]
    total_needed = margin * len(DONORS) / 2.0
    max_available = sum(float(P.loc[d, PRESPEC[d]]) for d in changed)
    grid = np.linspace(0, min(max_available, max(total_needed * 1.35, total_needed + 0.05)), 160)
    rows = []
    for mass in grid:
        P2 = P.copy()
        remaining = mass
        for d in changed:
            available = float(P2.loc[d, PRESPEC[d]])
            take = min(available, remaining / max(1, len([x for x in changed if P2.loc[x, PRESPEC[x]] > 0])))
            P2.loc[d, PRESPEC[d]] -= take
            P2.loc[d, best_alt[d]] += take
            remaining -= take
        if remaining > 1e-9:
            for d in changed:
                available = float(P2.loc[d, PRESPEC[d]])
                take = min(available, remaining)
                P2.loc[d, PRESPEC[d]] -= take
                P2.loc[d, best_alt[d]] += take
                remaining -= take
                if remaining <= 1e-9:
                    break
        pres2 = route_weight(P2, PRESPEC)
        alt2 = route_weight(P2, best_alt)
        rows.append({
            "targeted_shift_mass": float(mass),
            "prespecified_weight": pres2,
            "best_alternative_weight": alt2,
            "margin": pres2 - alt2,
        })
    return pd.DataFrame(rows), total_needed, max_available


def uniform_noise_curve(P, best_alt):
    U = pd.DataFrame(np.ones((3, 3)) / 3, index=DONORS, columns=ROLES)
    rows = []
    for lam in np.linspace(0, 1, 101):
        Pn = (1 - lam) * P + lam * U
        rw = route_weights(Pn)
        pres = rw.loc[rw["is_prespecified"], "weight"].iloc[0]
        alt = route_weight(Pn, best_alt)
        rows.append({
            "uniform_diffusion_lambda": float(lam),
            "prespecified_weight": float(pres),
            "best_alternative_weight": float(alt),
            "margin": float(pres - alt),
            "prespecified_rank": int(rw.loc[rw["is_prespecified"], "rank"].iloc[0]),
        })
    return pd.DataFrame(rows)


def raw_row_counterfactual(C, best_alt):
    pres_count = sum(float(C.loc[d, PRESPEC[d]]) for d in DONORS)
    alt_count = sum(float(C.loc[d, best_alt[d]]) for d in DONORS)
    total = float(C.values.sum())
    margin_count = pres_count - alt_count
    flips_needed = math.floor(margin_count / 2) + 1
    pres_route_count = pres_count
    changed = [d for d in DONORS if PRESPEC[d] != best_alt[d]]
    available_changed_route_rows = sum(float(C.loc[d, PRESPEC[d]]) for d in changed)
    return pd.DataFrame([{
        "raw_count_interpretation": "descriptive_not_independent",
        "route_eligible_raw_rows": total,
        "prespecified_route_raw_rows": pres_count,
        "best_alternative_raw_rows": alt_count,
        "raw_margin_rows": margin_count,
        "minimum_targeted_row_flips_to_make_best_alt_win": flips_needed,
        "flips_as_fraction_of_route_eligible_rows": flips_needed / total,
        "flips_as_fraction_of_prespecified_route_rows": flips_needed / pres_route_count,
        "changed_donor_rows_available_for_targeted_flip": available_changed_route_rows,
    }])


def main():
    P, C = load_matrices()
    routes = route_weights(P)
    best_alt_row = routes[~routes["is_prespecified"]].iloc[0]
    best_alt = {
        "host": best_alt_row["host_route"],
        "symbiont": best_alt_row["symbiont_route"],
        "other": best_alt_row["other_route"],
    }
    pres_weight = float(routes.loc[routes["is_prespecified"], "weight"].iloc[0])
    alt_weight = float(best_alt_row["weight"])
    shift_curve, shift_needed, max_available = targeted_shift_curve(P, best_alt)
    noise_curve = uniform_noise_curve(P, best_alt)
    raw_cf = raw_row_counterfactual(C, best_alt)

    route_mass_sum = sum(float(P.loc[d, PRESPEC[d]]) for d in DONORS)
    summary = pd.DataFrame([{
        "test": "source_equal_matrix_adversarial_erasure",
        "prespecified_weight": pres_weight,
        "best_alternative_weight": alt_weight,
        "initial_margin": pres_weight - alt_weight,
        "best_alternative_map": best_alt_row["route_map"],
        "donor_rows_that_differ_from_best_alt": "; ".join([d for d in DONORS if PRESPEC[d] != best_alt[d]]),
        "minimum_targeted_probability_mass_shift": shift_needed,
        "shift_as_fraction_of_total_three_row_probability_mass": shift_needed / len(DONORS),
        "shift_as_fraction_of_prespecified_route_mass": shift_needed / route_mass_sum,
        "maximum_available_prespecified_mass_in_changed_rows": max_available,
        "uniform_diffusion_rank1_until_lambda": float(noise_curve.loc[noise_curve["prespecified_rank"].eq(1), "uniform_diffusion_lambda"].max()),
        "interpretation": (
            "A targeted adversary must transfer this much source-equal donor-role probability mass "
            "from prespecified cells to the nearest alternative cells before the prespecified route loses. "
            "This is a matrix-level robustness calculation, not an independent row-level test."
        ),
    }])

    with pd.ExcelWriter(OUT / "adversarial_counterfactual_v1.xlsx", engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="summary", index=False)
        routes.to_excel(xw, sheet_name="route_weights", index=False)
        P.reset_index().rename(columns={"index": "donor_class"}).to_excel(xw, sheet_name="source_equal_matrix", index=False)
        C.reset_index().rename(columns={"index": "donor_class"}).to_excel(xw, sheet_name="raw_role_counts", index=False)
        shift_curve.to_excel(xw, sheet_name="targeted_shift_curve", index=False)
        noise_curve.to_excel(xw, sheet_name="uniform_noise_curve", index=False)
        raw_cf.to_excel(xw, sheet_name="raw_row_counterfactual", index=False)

    summary.to_csv(OUT / "adversarial_summary.csv", index=False)
    shift_curve.to_csv(OUT / "targeted_shift_curve.csv", index=False)
    noise_curve.to_csv(OUT / "uniform_noise_curve.csv", index=False)
    raw_cf.to_csv(OUT / "raw_row_counterfactual.csv", index=False)

    make_figure(P, routes, shift_curve, shift_needed, noise_curve, raw_cf)
    (OUT / "README_adversarial_counterfactual_v1.md").write_text(
        "# Route-B adversarial counterfactual robustness\n\n"
        "This analysis asks how much source-equal donor-role probability mass must be "
        "moved from the prespecified route cells to the nearest alternative route cells "
        "before the prespecified eukaryogenesis route loses rank 1.\n\n"
        "Key outputs:\n\n"
        f"- Prespecified route weight: {pres_weight:.3f}\n"
        f"- Nearest alternative weight: {alt_weight:.3f}\n"
        f"- Minimum targeted probability-mass shift: {shift_needed:.3f}\n"
        f"- Shift as fraction of total three-row probability mass: {shift_needed / len(DONORS):.3f}\n"
        f"- Shift as fraction of prespecified route mass: {shift_needed / route_mass_sum:.3f}\n\n"
        "The raw-row counterfactual is descriptive because row-level entries are traceable "
        "but not independent.\n",
        encoding="utf-8",
    )
    print(OUT)


def make_figure(P, routes, shift_curve, shift_needed, noise_curve, raw_cf):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.7), constrained_layout=True)
    ax = axs[0, 0]
    im = ax.imshow(P.values, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(3), ROLE_LABELS)
    ax.set_yticks(range(3), DONORS)
    ax.set_title("source-equal donor-role matrix", fontsize=9, weight="bold")
    for i, d in enumerate(DONORS):
        for j, r in enumerate(ROLES):
            val = P.loc[d, r]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=7, color="white" if val > 0.55 else "#1F2937",
                    weight="bold" if PRESPEC[d] == r else "normal")
    ax.tick_params(length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("probability", fontsize=7)
    ax.text(-0.16, 1.06, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[0, 1]
    plot_routes = routes.sort_values("weight", ascending=True)
    colors = [COL["host"] if x else "#D6DEE6" for x in plot_routes["is_prespecified"]]
    ax.barh(np.arange(len(plot_routes)), plot_routes["weight"], color=colors, edgecolor="white")
    ax.set_yticks(np.arange(len(plot_routes)), [f"R{int(r)}" for r in plot_routes["rank"]])
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("source-equal route weight")
    ax.set_title("nearest alternative gap", fontsize=9, weight="bold")
    for y, (_, row) in enumerate(plot_routes.iterrows()):
        ax.text(row["weight"] + 0.01, y, f"{row['weight']:.2f}", va="center", fontsize=6.5)
    ax.text(-0.16, 1.06, "b", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 0]
    ax.plot(shift_curve["targeted_shift_mass"], shift_curve["margin"],
            color=COL["red"], lw=1.8)
    ax.axhline(0, color="#111827", lw=0.8)
    ax.axvline(shift_needed, color=COL["red"], ls="--", lw=1.0)
    ax.fill_between(shift_curve["targeted_shift_mass"], shift_curve["margin"], 0,
                    where=shift_curve["margin"] >= 0, color="#F4D9D7", alpha=0.55)
    ax.set_xlabel("targeted shifted probability mass")
    ax.set_ylabel("margin over nearest alternative")
    ax.set_title("adversarial route erasure", fontsize=9, weight="bold")
    ax.text(shift_needed, max(shift_curve["margin"]) * 0.55,
            f"threshold\n{shift_needed:.3f}", ha="left", va="center",
            fontsize=6.6, color=COL["red"])
    ax.text(-0.16, 1.06, "c", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axs[1, 1]
    ax.plot(noise_curve["uniform_diffusion_lambda"], noise_curve["margin"],
            color=COL["green"], lw=1.8)
    ax.axhline(0, color="#111827", lw=0.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("uniform diffusion noise")
    ax.set_ylabel("margin over nearest alternative")
    ax.set_title("diffuse-noise degradation", fontsize=9, weight="bold")
    flips = int(raw_cf["minimum_targeted_row_flips_to_make_best_alt_win"].iloc[0])
    frac = float(raw_cf["flips_as_fraction_of_route_eligible_rows"].iloc[0])
    ax.text(0.06, 0.18,
            f"descriptive raw-row flip\nthreshold: {flips} rows\n({frac:.1%} of route-space rows)",
            transform=ax.transAxes, fontsize=6.5, color=COL["grey"],
            bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#D6DEE6"})
    ax.text(-0.16, 1.06, "d", transform=ax.transAxes, fontsize=12, weight="bold")

    fig.savefig(OUT / "adversarial_counterfactual_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "adversarial_counterfactual_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "adversarial_counterfactual_v1.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
