from __future__ import annotations

import math
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
SOURCE_XLSX = ROOT / "10_WORKING_COPY_FROM_NC12" / "04_Source_Data" / "Source_Data_complete.xlsx"
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_expanded_computational_evidence_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
ROLE_SHORT = {"scaffold": "scaffold", "energetic_incorporation": "energy", "transition": "transition"}
PRESPEC = {"host": "scaffold", "symbiont": "energetic_incorporation", "other": "transition"}


def norm(x):
    if pd.isna(x):
        return ""
    return str(x).strip().lower().replace("injection", "energetic_incorporation")


def load_probability_and_counts():
    mat = pd.read_excel(SOURCE_XLSX, sheet_name="continuous_role_matrix")
    mat["donor_class"] = mat["donor_class"].map(norm)
    mat["functional_class"] = mat["functional_class"].map(norm)
    mat = mat[mat["donor_class"].isin(DONORS) & mat["functional_class"].isin(ROLES)].copy()
    prob = mat.pivot(index="donor_class", columns="functional_class", values="source_equal_probability").reindex(index=DONORS, columns=ROLES).astype(float)
    counts = mat.pivot(index="donor_class", columns="functional_class", values="raw_role_space_count").reindex(index=DONORS, columns=ROLES).fillna(0).astype(float)
    return prob, counts


def route_weight(P, mapping):
    return float(np.mean([P.loc[d, mapping[d]] for d in DONORS]))


def nmi_from_matrix(P):
    P = np.asarray(P, dtype=float)
    P = P / P.sum(axis=1, keepdims=True)
    donor_prior = np.ones(P.shape[0]) / P.shape[0]
    joint = P * donor_prior[:, None]
    role_marg = joint.sum(axis=0)
    donor_marg = joint.sum(axis=1)
    mi = 0.0
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if joint[i, j] > 0:
                mi += joint[i, j] * math.log2(joint[i, j] / (donor_marg[i] * role_marg[j]))
    return mi, mi / math.log2(P.shape[1])


def fitted_model_matrices(P, C):
    # Candidate probability models Q(role | donor).
    models = {}
    # Diffuse mixture.
    models["diffuse_uniform"] = np.ones((3, 3)) / 3
    # Global role composition independent of donor.
    global_role = np.asarray(C.sum(axis=0), dtype=float)
    global_role = (global_role + 1) / (global_role.sum() + len(ROLES))
    models["donor_independent_global_role"] = np.tile(global_role, (3, 1))
    # Shared-theta route.
    theta = np.mean([P.loc[d, PRESPEC[d]] for d in DONORS])
    Q = np.zeros((3, 3))
    for i, d in enumerate(DONORS):
        for j, r in enumerate(ROLES):
            Q[i, j] = theta if PRESPEC[d] == r else (1 - theta) / 2
    models["route_shared_theta"] = Q
    # Donor-specific route.
    Q = np.zeros((3, 3))
    for i, d in enumerate(DONORS):
        theta_d = float(P.loc[d, PRESPEC[d]])
        for j, r in enumerate(ROLES):
            Q[i, j] = theta_d if PRESPEC[d] == r else (1 - theta_d) / 2
    models["route_donor_specific_theta"] = Q
    # Best hard alternatives, shared theta by mapping.
    for perm in permutations(ROLES):
        mapping = dict(zip(DONORS, perm))
        label = "map_" + "_".join([ROLE_SHORT[mapping[d]] for d in DONORS])
        theta_m = np.mean([P.loc[d, mapping[d]] for d in DONORS])
        Q = np.zeros((3, 3))
        for i, d in enumerate(DONORS):
            for j, r in enumerate(ROLES):
                Q[i, j] = theta_m if mapping[d] == r else (1 - theta_m) / 2
        models[label] = Q
    # Single-role concentration models.
    for role in ROLES:
        theta_r = np.mean([P.loc[d, role] for d in DONORS])
        Q = np.zeros((3, 3))
        for j, r in enumerate(ROLES):
            Q[:, j] = theta_r if r == role else (1 - theta_r) / 2
        models["single_role_" + ROLE_SHORT[role]] = Q
    return models


def model_comparison(P, C):
    models = fitted_model_matrices(P, C)
    count_mat = C.values.astype(float)
    n = count_mat.sum()
    rows = []
    for name, Q in models.items():
        Q = np.clip(Q, 1e-12, 1)
        Q = Q / Q.sum(axis=1, keepdims=True)
        loglik_raw = float((count_mat * np.log(Q)).sum())
        # Source-equal cross entropy uses each donor row of P equally.
        Pmat = P.values.astype(float)
        Pmat = Pmat / Pmat.sum(axis=1, keepdims=True)
        source_equal_ce = float(-(Pmat * np.log2(Q)).sum(axis=1).mean())
        # Simple k bookkeeping for information criteria.
        if name == "diffuse_uniform":
            k = 0
        elif name == "donor_independent_global_role":
            k = 2
        elif name == "route_shared_theta" or name.startswith("map_") or name.startswith("single_role_"):
            k = 1
        elif name == "route_donor_specific_theta":
            k = 3
        else:
            k = 1
        aic = 2 * k - 2 * loglik_raw
        bic = k * math.log(n) - 2 * loglik_raw
        rows.append({
            "model": name,
            "k_parameters": k,
            "raw_count_log_likelihood_descriptive": loglik_raw,
            "AIC_descriptive": aic,
            "BIC_descriptive": bic,
            "source_equal_cross_entropy_bits": source_equal_ce,
            "route_family": (
                "prespecified_route" if name in {"route_shared_theta", "route_donor_specific_theta", "map_scaffold_energy_transition"}
                else "alternative"
            ),
        })
    df = pd.DataFrame(rows).sort_values(["source_equal_cross_entropy_bits", "AIC_descriptive"]).reset_index(drop=True)
    df["rank_by_source_equal_CE"] = np.arange(1, len(df) + 1)
    df = df.sort_values("AIC_descriptive").reset_index(drop=True)
    min_aic = df["AIC_descriptive"].min()
    df["delta_AIC_descriptive"] = df["AIC_descriptive"] - min_aic
    w = np.exp(-0.5 * df["delta_AIC_descriptive"])
    df["AIC_weight_descriptive"] = w / w.sum()
    df = df.sort_values("source_equal_cross_entropy_bits").reset_index(drop=True)
    return df


def posterior_route_uncertainty(C, n_iter=20000, seed=20260630):
    rng = np.random.default_rng(seed)
    counts = C.values.astype(float)
    routes = []
    route_maps = []
    for perm in permutations(ROLES):
        route_maps.append(dict(zip(DONORS, perm)))
    pres_idx = [i for i, m in enumerate(route_maps) if m == PRESPEC][0]
    samples = []
    rank1_count = 0
    margin_samples = []
    nmi_samples = []
    conc_samples = []
    for _ in range(n_iter):
        P = np.vstack([rng.dirichlet(counts[i] + 1.0) for i in range(3)])
        weights = np.array([np.mean([P[j, ROLES.index(m[DONORS[j]])] for j in range(3)]) for m in route_maps])
        rank = int((-weights).argsort().tolist().index(pres_idx) + 1)
        if rank == 1:
            rank1_count += 1
        best_alt = np.max(np.delete(weights, pres_idx))
        margin_samples.append(weights[pres_idx] - best_alt)
        mi, nmi = nmi_from_matrix(P)
        nmi_samples.append(nmi)
        conc_samples.append(weights[pres_idx] / (1 / 3))
        samples.append(weights[pres_idx])
    arr = np.array(samples)
    margin = np.array(margin_samples)
    nmi_arr = np.array(nmi_samples)
    conc = np.array(conc_samples)
    summary = pd.DataFrame([{
        "posterior_type": "Dirichlet_row_count_descriptive",
        "iterations": n_iter,
        "prespecified_rank1_probability": rank1_count / n_iter,
        "prespecified_weight_mean": float(arr.mean()),
        "prespecified_weight_q025": float(np.quantile(arr, 0.025)),
        "prespecified_weight_q50": float(np.quantile(arr, 0.5)),
        "prespecified_weight_q975": float(np.quantile(arr, 0.975)),
        "margin_vs_best_alt_mean": float(margin.mean()),
        "margin_vs_best_alt_q025": float(np.quantile(margin, 0.025)),
        "margin_vs_best_alt_q50": float(np.quantile(margin, 0.5)),
        "margin_vs_best_alt_q975": float(np.quantile(margin, 0.975)),
        "route_concentration_over_diffuse_mean": float(conc.mean()),
        "normalized_mutual_information_mean": float(nmi_arr.mean()),
        "interpretation": "Posterior is descriptive because row-level units are not independent; it provides an uncertainty stress test around the observed donor-role count structure.",
    }])
    dist = pd.DataFrame({
        "prespecified_weight": arr,
        "margin_vs_best_alternative": margin,
        "route_concentration_over_diffuse": conc,
        "normalized_mutual_information": nmi_arr,
    })
    return summary, dist


def sensitivity_envelope():
    sheets = []
    for sheet, label in [
        ("leave_source_audit", "leave_one_source"),
        ("leave_layer_audit", "leave_one_layer"),
        ("subset_robustness", "subset"),
    ]:
        df = pd.read_excel(SOURCE_XLSX, sheet_name=sheet)
        df["stress_test_family"] = label
        sheets.append(df)
    all_df = pd.concat(sheets, ignore_index=True, sort=False)
    key_cols = ["stress_test_family", "subset", "n_units", "n_sources", "prespecified_rank", "prespecified_source_equal", "best_alternative_source_equal", "margin_vs_best_alternative"]
    present = [c for c in key_cols if c in all_df.columns]
    env = all_df[present].copy()
    env["passes_rank1"] = env["prespecified_rank"].eq(1)
    summary = env.groupby("stress_test_family").agg(
        n_tests=("subset", "count"),
        rank1_fraction=("passes_rank1", "mean"),
        min_margin=("margin_vs_best_alternative", "min"),
        median_margin=("margin_vs_best_alternative", "median"),
        max_margin=("margin_vs_best_alternative", "max"),
        min_prespecified_support=("prespecified_source_equal", "min"),
        median_prespecified_support=("prespecified_source_equal", "median"),
    ).reset_index()
    worst = env.sort_values("margin_vs_best_alternative").head(20)
    return summary, env, worst


def donor_role_network_diagnostics(P):
    rows = []
    for d in DONORS:
        probs = P.loc[d].values.astype(float)
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        max_entropy = math.log2(len(ROLES))
        pres = float(P.loc[d, PRESPEC[d]])
        off = 1 - pres
        top_role = P.loc[d].idxmax()
        rows.append({
            "donor_class": d,
            "prespecified_role": PRESPEC[d],
            "prespecified_probability": pres,
            "off_route_leakage": off,
            "top_role": top_role,
            "top_is_prespecified": top_role == PRESPEC[d],
            "entropy_bits": entropy,
            "entropy_fraction_of_diffuse": entropy / max_entropy,
            "dominance_index_1_minus_entropy_fraction": 1 - entropy / max_entropy,
        })
    df = pd.DataFrame(rows)
    route_weight = np.mean(df["prespecified_probability"])
    off_weight = 1 - route_weight
    mi, nmi = nmi_from_matrix(P.values)
    summary = pd.DataFrame([{
        "route_weight": route_weight,
        "off_route_weight": off_weight,
        "route_to_off_route_ratio": route_weight / off_weight,
        "normalized_mutual_information": nmi,
        "mean_donor_dominance_index": df["dominance_index_1_minus_entropy_fraction"].mean(),
        "main_leakage_edge": "host->transition",
        "main_leakage_probability": float(P.loc["host", "transition"]),
        "interpretation": "The donor-role graph is structured, but not perfectly diagonal; host-to-transition leakage marks the main boundary rather than invalidating the route architecture.",
    }])
    return summary, df


def plot_outputs(model_df, post_summary, post_dist, sens_summary, sens_env, network_summary, network_donor):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7.0,
        "axes.spines.right": False,
        "axes.spines.top": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    ax = axes[0, 0]
    plot = model_df.sort_values("source_equal_cross_entropy_bits").head(6).iloc[::-1]
    colors = ["#2F5D7C" if "route" in m and ("shared" in m or "donor_specific" in m or m == "map_scaffold_energy_transition") else "#D8DEE6" for m in plot["model"]]
    ax.barh(np.arange(len(plot)), plot["source_equal_cross_entropy_bits"], color=colors, edgecolor="white")
    ax.set_yticks(np.arange(len(plot)))
    label_map = {
        "route_donor_specific_theta": "route\ndonor-specific",
        "route_shared_theta": "route\nshared",
        "map_scaffold_energy_transition": "prespecified\nmap",
        "map_energy_transition_scaffold": "best alt.\nmap",
        "map_transition_energy_scaffold": "alt.\nmap",
        "map_transition_scaffold_energy": "alt.\nmap",
        "single_role_scaffold": "single-role\nscaffold",
    }
    ax.set_yticklabels([label_map.get(m, m.replace("_", "\n")) for m in plot["model"]], fontsize=6)
    ax.set_xlabel("source-equal cross entropy (bits)")
    ax.set_title("Model-family comparison", fontsize=8, weight="bold")

    ax = axes[0, 1]
    ax.hist(post_dist["margin_vs_best_alternative"], bins=45, color="#4F7896", alpha=0.9)
    ax.axvline(0, color="#B95A55", lw=1.0, ls="--")
    ax.axvline(post_summary["margin_vs_best_alt_q025"].iloc[0], color="#20262E", lw=0.8, ls=":")
    ax.axvline(post_summary["margin_vs_best_alt_q975"].iloc[0], color="#20262E", lw=0.8, ls=":")
    ax.set_xlabel("posterior margin vs best alternative")
    ax.set_ylabel("draws")
    ax.set_title("Posterior rank stability", fontsize=8, weight="bold")
    ax.text(0.03, 0.95, f"P(rank 1)={post_summary['prespecified_rank1_probability'].iloc[0]:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7)

    ax = axes[1, 0]
    env = sens_env.copy()
    fam_order = ["leave_one_source", "leave_one_layer", "subset"]
    xpos = {f: i for i, f in enumerate(fam_order)}
    for f in fam_order:
        vals = env.loc[env["stress_test_family"].eq(f), "margin_vs_best_alternative"].dropna()
        x = np.full(len(vals), xpos[f]) + np.linspace(-0.12, 0.12, len(vals)) if len(vals) else []
        ax.scatter(x, vals, s=18, color="#7A58A6" if f == "subset" else "#6B8F71", alpha=0.75, edgecolor="white", linewidth=0.3)
    ax.axhline(0, color="#B95A55", lw=1, ls="--")
    ax.set_xticks(range(len(fam_order)))
    ax.set_xticklabels(["leave\nsource", "leave\nlayer", "subsets"])
    ax.set_ylabel("margin vs best alternative")
    ax.set_title("Stress-test envelope", fontsize=8, weight="bold")

    ax = axes[1, 1]
    colors = ["#356B8A", "#E0942A", "#7A58A6"]
    ax.bar(np.arange(3), network_donor["prespecified_probability"], color=colors)
    ax.bar(np.arange(3), network_donor["off_route_leakage"], bottom=network_donor["prespecified_probability"], color="#D8DEE6")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["host", "symbiont", "other"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("probability")
    ax.set_title("Route concentration and leakage", fontsize=8, weight="bold")
    ax.text(0.98, 0.95, "colour = on-route role\ngrey = off-route leakage", transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#59636E")

    for label, ax in zip(["a", "b", "c", "d"], axes.ravel()):
        ax.text(-0.16, 1.10, label, transform=ax.transAxes, fontsize=11, weight="bold", va="top")
    fig.tight_layout(w_pad=2.2, h_pad=2.3)
    fig.savefig(OUT / "expanded_computational_evidence_v1.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "expanded_computational_evidence_v1.pdf", bbox_inches="tight")
    fig.savefig(OUT / "expanded_computational_evidence_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    P, C = load_probability_and_counts()
    model_df = model_comparison(P, C)
    post_summary, post_dist = posterior_route_uncertainty(C)
    sens_summary, sens_env, sens_worst = sensitivity_envelope()
    net_summary, net_donor = donor_role_network_diagnostics(P)

    xlsx = OUT / "expanded_computational_evidence_v1.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        model_df.to_excel(writer, sheet_name="model_family_comparison", index=False)
        post_summary.to_excel(writer, sheet_name="posterior_summary", index=False)
        post_dist.sample(min(5000, len(post_dist)), random_state=7).to_excel(writer, sheet_name="posterior_draws_sample", index=False)
        sens_summary.to_excel(writer, sheet_name="stress_test_summary", index=False)
        sens_env.to_excel(writer, sheet_name="stress_test_all", index=False)
        sens_worst.to_excel(writer, sheet_name="stress_test_worst", index=False)
        net_summary.to_excel(writer, sheet_name="network_summary", index=False)
        net_donor.to_excel(writer, sheet_name="network_by_donor", index=False)
    for name, df in [
        ("model_family_comparison.csv", model_df),
        ("posterior_summary.csv", post_summary),
        ("stress_test_summary.csv", sens_summary),
        ("stress_test_worst.csv", sens_worst),
        ("network_summary.csv", net_summary),
        ("network_by_donor.csv", net_donor),
    ]:
        df.to_csv(OUT / name, index=False, encoding="utf-8-sig")
    plot_outputs(model_df, post_summary, post_dist, sens_summary, sens_env, net_summary, net_donor)
    readme = OUT / "README_expanded_computational_evidence_v1.md"
    readme.write_text(
        "# Expanded computational evidence v1\n\n"
        "This folder adds four atlas-level diagnostics for the route-B manuscript:\n\n"
        "1. Model-family comparison against diffuse, donor-independent, single-role and alternative-route models.\n"
        "2. Dirichlet posterior rank-stability stress test based on donor-role row-count structure.\n"
        "3. Leave-source, leave-layer and subset stress-test envelope.\n"
        "4. Donor-role network concentration, leakage and information diagnostics.\n\n"
        "These analyses are designed to strengthen the paper as a computational evidence atlas. "
        "Row-count posterior and likelihood outputs are labelled descriptive because row-level units are traceable but not independent.\n",
        encoding="utf-8",
    )
    print(xlsx)
    print(model_df.head(8).to_string(index=False))
    print(post_summary.to_string(index=False))
    print(sens_summary.to_string(index=False))
    print(net_summary.to_string(index=False))


if __name__ == "__main__":
    main()

