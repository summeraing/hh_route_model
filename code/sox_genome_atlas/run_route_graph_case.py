from __future__ import annotations

import argparse
import math
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "case_id",
    "unit_id",
    "source_id",
    "evidence_layer",
    "donor_class",
    "functional_class",
    "route_eligible",
    "dependency_group",
}


def split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_map(value: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in split_csv(value):
        if "=" not in part:
            raise ValueError(f"Invalid map item: {part!r}; expected donor=role")
        donor, role = part.split("=", 1)
        out[donor.strip()] = role.strip()
    return out


def read_table(path: Path, sheet: str | None) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet or 0)
    return pd.read_csv(path)


def normalize_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})


def clean_label(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def validate(df: pd.DataFrame) -> None:
    missing = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing:
        raise ValueError(f"Input table missing required columns: {missing}")


def route_indicator(df: pd.DataFrame, mapping: dict[str, str]) -> pd.Series:
    return df.apply(lambda r: mapping.get(r["donor_class"]) == r["functional_class"], axis=1).astype(float)


def route_scores(
    df: pd.DataFrame,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
) -> pd.DataFrame:
    rows = []
    maps = []
    if len(donors) == len(roles):
        for perm in permutations(roles):
            maps.append(dict(zip(donors, perm)))
    else:
        maps.append(prespec)

    for mapping in maps:
        work = df[df["donor_class"].isin(donors) & df["functional_class"].isin(roles)].copy()
        if work.empty:
            raw = np.nan
            source_equal = np.nan
            n_units = 0
            n_sources = 0
        else:
            work["route_match"] = route_indicator(work, mapping)
            raw = float(work["route_match"].mean())
            per_source = work.groupby("source_id")["route_match"].mean()
            source_equal = float(per_source.mean())
            n_units = int(len(work))
            n_sources = int(per_source.shape[0])
        rows.append(
            {
                "route_map": "; ".join(f"{d}->{mapping[d]}" for d in donors if d in mapping),
                "is_prespecified_route": mapping == prespec,
                "n_units": n_units,
                "n_sources": n_sources,
                "raw_route_rate": raw,
                "source_equal_route_rate": source_equal,
            }
        )
    out = pd.DataFrame(rows).sort_values("source_equal_route_rate", ascending=False).reset_index(drop=True)
    out["rank_by_source_equal"] = np.arange(1, len(out) + 1)
    best_alt = out.loc[~out["is_prespecified_route"], "source_equal_route_rate"]
    pres = out.loc[out["is_prespecified_route"], "source_equal_route_rate"]
    if len(pres) and len(best_alt):
        out["prespecified_margin_vs_best_alternative"] = float(pres.iloc[0] - best_alt.max())
    else:
        out["prespecified_margin_vs_best_alternative"] = np.nan
    return out


def source_rates(df: pd.DataFrame, prespec: dict[str, str]) -> pd.DataFrame:
    work = df.copy()
    work["route_match"] = route_indicator(work, prespec)
    rows = []
    for source, g in work.groupby("source_id"):
        rows.append(
            {
                "source_id": source,
                "n_units": int(len(g)),
                "route_rate": float(g["route_match"].mean()),
                "n_layers": int(g["evidence_layer"].nunique()),
                "n_dependency_groups": int(g["dependency_group"].nunique()),
            }
        )
    return pd.DataFrame(rows).sort_values("n_units", ascending=False)


def donor_role_matrix(df: pd.DataFrame, donors: list[str], roles: list[str]) -> pd.DataFrame:
    rows = []
    for source, g in df.groupby("source_id"):
        tab = pd.crosstab(g["donor_class"], g["functional_class"]).reindex(index=donors, columns=roles, fill_value=0)
        total = tab.values.sum()
        if total == 0:
            continue
        prob = tab / total
        for donor in donors:
            for role in roles:
                rows.append(
                    {
                        "source_id": source,
                        "donor_class": donor,
                        "functional_class": role,
                        "raw_count": int(tab.loc[donor, role]),
                        "within_source_probability": float(prob.loc[donor, role]),
                    }
                )
    per_source = pd.DataFrame(rows)
    if per_source.empty:
        return per_source
    source_equal = (
        per_source.groupby(["donor_class", "functional_class"], as_index=False)
        .agg(
            source_equal_probability=("within_source_probability", "mean"),
            raw_count=("raw_count", "sum"),
            n_sources_observed=("source_id", "nunique"),
        )
    )
    return source_equal.sort_values(["donor_class", "functional_class"])


def leave_one_source(
    df: pd.DataFrame,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
) -> pd.DataFrame:
    rows = []
    for source in sorted(df["source_id"].unique()):
        sub = df[df["source_id"] != source]
        scores = route_scores(sub, donors, roles, prespec)
        pres = scores[scores["is_prespecified_route"]].iloc[0]
        best_alt = scores[~scores["is_prespecified_route"]].iloc[0] if (~scores["is_prespecified_route"]).any() else None
        rows.append(
            {
                "removed_source": source,
                "n_units": int(len(sub)),
                "n_sources": int(sub["source_id"].nunique()),
                "prespecified_rank": int(pres["rank_by_source_equal"]),
                "prespecified_source_equal": float(pres["source_equal_route_rate"]),
                "best_alternative_source_equal": float(best_alt["source_equal_route_rate"]) if best_alt is not None else np.nan,
                "margin_vs_best_alternative": float(pres["source_equal_route_rate"] - best_alt["source_equal_route_rate"]) if best_alt is not None else np.nan,
            }
        )
    return pd.DataFrame(rows)


def leave_one_layer(
    df: pd.DataFrame,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
) -> pd.DataFrame:
    rows = []
    for layer in sorted(df["evidence_layer"].unique()):
        sub = df[df["evidence_layer"] != layer]
        scores = route_scores(sub, donors, roles, prespec)
        pres = scores[scores["is_prespecified_route"]].iloc[0]
        best_alt = scores[~scores["is_prespecified_route"]].iloc[0] if (~scores["is_prespecified_route"]).any() else None
        rows.append(
            {
                "removed_layer": layer,
                "n_units": int(len(sub)),
                "n_sources": int(sub["source_id"].nunique()),
                "n_layers": int(sub["evidence_layer"].nunique()),
                "prespecified_rank": int(pres["rank_by_source_equal"]),
                "prespecified_source_equal": float(pres["source_equal_route_rate"]),
                "best_alternative_source_equal": float(best_alt["source_equal_route_rate"]) if best_alt is not None else np.nan,
                "margin_vs_best_alternative": float(pres["source_equal_route_rate"] - best_alt["source_equal_route_rate"]) if best_alt is not None else np.nan,
            }
        )
    return pd.DataFrame(rows)


def leave_one_group(
    df: pd.DataFrame,
    group_col: str,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
) -> pd.DataFrame:
    rows = []
    values = sorted(value for value in df[group_col].dropna().astype(str).unique() if value)
    for value in values:
        sub = df[df[group_col].astype(str) != value]
        if sub.empty:
            continue
        scores = route_scores(sub, donors, roles, prespec)
        pres = scores[scores["is_prespecified_route"]].iloc[0]
        best_alt = scores[~scores["is_prespecified_route"]].iloc[0]
        rows.append(
            {
                "group_column": group_col,
                "removed_group": value,
                "n_units": int(len(sub)),
                "n_dependency_groups": int(sub["dependency_group"].nunique()),
                "prespecified_rank": int(pres["rank_by_source_equal"]),
                "prespecified_source_equal": float(pres["source_equal_route_rate"]),
                "best_alternative_source_equal": float(best_alt["source_equal_route_rate"]),
                "margin_vs_best_alternative": float(
                    pres["source_equal_route_rate"] - best_alt["source_equal_route_rate"]
                ),
            }
        )
    return pd.DataFrame(rows)


def collapse_dependencies(df: pd.DataFrame, dependency_col: str = "dependency_group") -> pd.DataFrame:
    rows = []
    group_cols = ["source_id", dependency_col]
    for (source, dep), g in df.groupby(group_cols):
        donor_counts = g["donor_class"].value_counts()
        role_counts = g["functional_class"].value_counts()
        donor_winners = donor_counts.index[donor_counts.eq(donor_counts.iloc[0])].tolist()
        role_winners = role_counts.index[role_counts.eq(role_counts.iloc[0])].tolist()
        donor_tie = len(donor_winners) != 1
        role_tie = len(role_winners) != 1
        donor = "uncertain" if donor_tie else donor_winners[0]
        role = "uncertain" if role_tie else role_winners[0]
        rows.append(
            {
                "case_id": g["case_id"].iloc[0],
                "unit_id": f"{source}::{dep}",
                "source_id": source,
                "evidence_layer": ";".join(sorted(map(str, g["evidence_layer"].unique()))),
                "donor_class": donor,
                "functional_class": role,
                "route_eligible": not donor_tie and not role_tie,
                "dependency_group": dep,
                "dependency_scheme": dependency_col,
                "collapsed_n_rows": int(len(g)),
                "donor_tie": donor_tie,
                "role_tie": role_tie,
            }
        )
    return pd.DataFrame(rows)


def permutation_null(
    df: pd.DataFrame,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
    strata: list[str],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    observed = route_scores(df, donors, roles, prespec)
    pres = observed[observed["is_prespecified_route"]].iloc[0]
    best_alt = observed[~observed["is_prespecified_route"]].iloc[0]
    observed_margin = float(pres["source_equal_route_rate"] - best_alt["source_equal_route_rate"])

    donor_index = {label: index for index, label in enumerate(donors)}
    role_index = {label: index for index, label in enumerate(roles)}
    source_labels = sorted(df["source_id"].unique())
    source_index = {label: index for index, label in enumerate(source_labels)}
    donor_codes = df["donor_class"].map(donor_index).to_numpy(dtype=np.int16)
    role_codes = df["functional_class"].map(role_index).to_numpy(dtype=np.int16)
    source_codes = df["source_id"].map(source_index).to_numpy(dtype=np.int16)
    source_counts = np.bincount(source_codes, minlength=len(source_labels)).astype(float)

    map_codes = np.array(list(permutations(range(len(roles)))), dtype=np.int16)
    prespec_codes = np.array([role_index[prespec[donor]] for donor in donors], dtype=np.int16)
    prespec_matches = np.flatnonzero(np.all(map_codes == prespec_codes, axis=1))
    if len(prespec_matches) != 1:
        raise ValueError("Prespecified route is not represented exactly once in permutation maps.")
    prespec_index = int(prespec_matches[0])
    alternative_indices = np.array(
        [index for index in range(len(map_codes)) if index != prespec_index], dtype=int
    )
    donor_axis = np.arange(len(donors))

    grouper = strata[0] if len(strata) == 1 else strata
    group_indices = [
        np.asarray(index, dtype=int)
        for index in df.groupby(grouper, sort=True).indices.values()
    ]
    margins_arr = np.empty(iterations, dtype=float)
    contingency_size = len(source_labels) * len(donors) * len(roles)
    for iteration in range(iterations):
        permuted_roles = role_codes.copy()
        for index in group_indices:
            permuted_roles[index] = rng.permutation(permuted_roles[index])
        flat = (
            source_codes * (len(donors) * len(roles))
            + donor_codes * len(roles)
            + permuted_roles
        )
        counts = np.bincount(flat, minlength=contingency_size).reshape(
            len(source_labels), len(donors), len(roles)
        )
        map_scores = np.empty(len(map_codes), dtype=float)
        for map_index, mapping in enumerate(map_codes):
            matched_by_source = counts[:, donor_axis, mapping].sum(axis=1)
            map_scores[map_index] = np.mean(matched_by_source / source_counts)
        margins_arr[iteration] = map_scores[prespec_index] - np.max(map_scores[alternative_indices])
    null_sd = float(np.std(margins_arr, ddof=1)) if len(margins_arr) > 1 else 0.0
    unique_margins = int(np.unique(np.round(margins_arr, 12)).size)
    return pd.DataFrame(
        [
            {
                "strata": "+".join(strata),
                "iterations": iterations,
                "observed_margin": observed_margin,
                "null_mean_margin": float(np.mean(margins_arr)),
                "null_q025_margin": float(np.quantile(margins_arr, 0.025)),
                "null_q975_margin": float(np.quantile(margins_arr, 0.975)),
                "null_min_margin": float(np.min(margins_arr)),
                "null_max_margin": float(np.max(margins_arr)),
                "null_sd_margin": null_sd,
                "null_unique_margins": unique_margins,
                "null_non_degenerate": bool(null_sd > 0 and unique_margins > 1),
                "empirical_p_ge_observed": float((np.sum(margins_arr >= observed_margin) + 1) / (iterations + 1)),
            }
        ]
    )


def preflight_report(
    work: pd.DataFrame,
    min_units: int,
    min_sources: int,
    min_dependency_groups: int,
    max_source_fraction: float,
) -> pd.DataFrame:
    source_counts = work["source_id"].value_counts()
    donor_counts = work["donor_class"].value_counts()
    role_counts = work["functional_class"].value_counts()
    n_units = int(len(work))
    n_sources = int(work["source_id"].nunique())
    n_dependency_groups = int(work["dependency_group"].nunique())
    largest_source_fraction = float(source_counts.iloc[0] / n_units) if n_units else np.nan
    checks = [
        {
            "gate": "minimum_route_eligible_units",
            "observed": n_units,
            "threshold": min_units,
            "pass": n_units >= min_units,
        },
        {
            "gate": "minimum_sources",
            "observed": n_sources,
            "threshold": min_sources,
            "pass": n_sources >= min_sources,
        },
        {
            "gate": "minimum_dependency_groups",
            "observed": n_dependency_groups,
            "threshold": min_dependency_groups,
            "pass": n_dependency_groups >= min_dependency_groups,
        },
        {
            "gate": "maximum_single_source_fraction",
            "observed": largest_source_fraction,
            "threshold": max_source_fraction,
            "pass": largest_source_fraction <= max_source_fraction,
        },
        {
            "gate": "all_donor_classes_represented",
            "observed": "; ".join(f"{k}:{v}" for k, v in donor_counts.items()),
            "threshold": "all requested donor labels present",
            "pass": True,
        },
        {
            "gate": "all_role_classes_represented",
            "observed": "; ".join(f"{k}:{v}" for k, v in role_counts.items()),
            "threshold": "all requested role labels present",
            "pass": True,
        },
    ]
    return pd.DataFrame(checks)


def main() -> None:
    ap = argparse.ArgumentParser(description="Run a route-resolved evidence-graph analysis for one case.")
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--sheet", default=None)
    ap.add_argument("--case-id", required=True)
    ap.add_argument("--donors", required=True, help="Comma-separated donor labels.")
    ap.add_argument("--roles", required=True, help="Comma-separated role labels.")
    ap.add_argument("--prespec", required=True, help="Comma-separated donor=role mapping.")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--iterations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260701)
    ap.add_argument("--min-units", type=int, default=60)
    ap.add_argument("--min-sources", type=int, default=3)
    ap.add_argument("--min-dependency-groups", type=int, default=30)
    ap.add_argument("--max-source-fraction", type=float, default=0.70)
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--enforce-gates", action="store_true")
    args = ap.parse_args()

    donors = split_csv(args.donors)
    roles = split_csv(args.roles)
    prespec = parse_map(args.prespec)
    args.out.mkdir(parents=True, exist_ok=True)

    df = read_table(args.input, args.sheet)
    validate(df)
    df = df.copy()
    df["case_id"] = df["case_id"].map(clean_label)
    df = df[df["case_id"] == args.case_id].copy()
    if df.empty:
        raise ValueError(f"No rows found for case_id={args.case_id!r}")
    df["route_eligible"] = normalize_bool(df["route_eligible"])
    for col in ["source_id", "evidence_layer", "donor_class", "functional_class", "dependency_group", "unit_id"]:
        df[col] = df[col].map(clean_label)
    for col in ["genome_dependency_group", "operon_dependency_group", "homolog_dependency_group", "gtdb_class", "gtdb_order"]:
        if col in df.columns:
            df[col] = df[col].map(clean_label)
    work = df[
        df["route_eligible"]
        & df["donor_class"].isin(donors)
        & df["functional_class"].isin(roles)
    ].copy().reset_index(drop=True)
    if work.empty:
        raise ValueError("No route-eligible rows remain after donor/role filtering.")

    preflight = preflight_report(
        work,
        min_units=args.min_units,
        min_sources=args.min_sources,
        min_dependency_groups=args.min_dependency_groups,
        max_source_fraction=args.max_source_fraction,
    )
    requested_donors = set(donors)
    requested_roles = set(roles)
    present_donors = set(work["donor_class"].unique())
    present_roles = set(work["functional_class"].unique())
    preflight.loc[preflight["gate"] == "all_donor_classes_represented", "pass"] = requested_donors <= present_donors
    preflight.loc[preflight["gate"] == "all_role_classes_represented", "pass"] = requested_roles <= present_roles
    preflight.to_csv(args.out / "preflight_report.csv", index=False)
    failed = preflight[~preflight["pass"].astype(bool)]
    if args.enforce_gates and not failed.empty:
        raise SystemExit("Preflight gates failed:\n" + failed.to_string(index=False))
    if args.preflight_only:
        print(preflight.to_string(index=False))
        return

    scores = route_scores(work, donors, roles, prespec)
    matrix = donor_role_matrix(work, donors, roles)
    source = source_rates(work, prespec)
    leave_source = leave_one_source(work, donors, roles, prespec)
    leave_layer = leave_one_layer(work, donors, roles, prespec)
    taxon_holdouts: dict[str, pd.DataFrame] = {}
    for group_col in ["gtdb_class", "gtdb_order"]:
        if group_col in work.columns and work[group_col].replace("", np.nan).nunique() >= 2:
            taxon_holdouts[group_col] = leave_one_group(work, group_col, donors, roles, prespec)
    collapsed = collapse_dependencies(work)
    collapsed_eligible = collapsed[
        collapsed["route_eligible"]
        & collapsed["donor_class"].isin(donors)
        & collapsed["functional_class"].isin(roles)
    ].copy()
    if collapsed_eligible.empty:
        raise ValueError("No dependency-collapsed units remain after excluding donor/role ties.")
    collapsed_scores = route_scores(collapsed_eligible, donors, roles, prespec)
    dependency_summary_rows = [
        {
            "dependency_scheme": "dependency_group",
            "collapsed_units": int(len(collapsed)),
            "route_eligible_units": int(len(collapsed_eligible)),
            "tie_exclusions": int((~collapsed["route_eligible"]).sum()),
            "prespecified_rank": int(
                collapsed_scores.loc[collapsed_scores["is_prespecified_route"], "rank_by_source_equal"].iloc[0]
            ),
            "margin_vs_best_alternative": float(
                collapsed_scores["prespecified_margin_vs_best_alternative"].iloc[0]
            ),
        }
    ]
    additional_dependency_tables: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}
    for dependency_col in ["genome_dependency_group", "homolog_dependency_group"]:
        if dependency_col not in work.columns:
            continue
        collapsed_alt = collapse_dependencies(work, dependency_col=dependency_col)
        eligible_alt = collapsed_alt[
            collapsed_alt["route_eligible"]
            & collapsed_alt["donor_class"].isin(donors)
            & collapsed_alt["functional_class"].isin(roles)
        ].copy()
        if eligible_alt.empty:
            continue
        scores_alt = route_scores(eligible_alt, donors, roles, prespec)
        dependency_summary_rows.append(
            {
                "dependency_scheme": dependency_col,
                "collapsed_units": int(len(collapsed_alt)),
                "route_eligible_units": int(len(eligible_alt)),
                "tie_exclusions": int((~collapsed_alt["route_eligible"]).sum()),
                "prespecified_rank": int(
                    scores_alt.loc[scores_alt["is_prespecified_route"], "rank_by_source_equal"].iloc[0]
                ),
                "margin_vs_best_alternative": float(
                    scores_alt["prespecified_margin_vs_best_alternative"].iloc[0]
                ),
            }
        )
        additional_dependency_tables[dependency_col] = (collapsed_alt, eligible_alt, scores_alt)
    dependency_sensitivity = pd.DataFrame(dependency_summary_rows)
    perm_source = permutation_null(work, donors, roles, prespec, ["source_id"], args.iterations, args.seed)
    perm_layer = permutation_null(work, donors, roles, prespec, ["source_id", "evidence_layer"], args.iterations, args.seed + 1)

    summary = pd.DataFrame(
        [
            {
                "case_id": args.case_id,
                "n_route_eligible_units": int(len(work)),
                "n_sources": int(work["source_id"].nunique()),
                "n_layers": int(work["evidence_layer"].nunique()),
                "n_dependency_groups": int(work["dependency_group"].nunique()),
                "n_collapsed_units": int(len(collapsed)),
                "n_collapsed_route_eligible_units": int(len(collapsed_eligible)),
                "n_collapsed_tie_exclusions": int((~collapsed["route_eligible"]).sum()),
                "prespecified_rank": int(scores[scores["is_prespecified_route"]]["rank_by_source_equal"].iloc[0]),
                "prespecified_source_equal": float(scores[scores["is_prespecified_route"]]["source_equal_route_rate"].iloc[0]),
                "best_alternative_source_equal": float(scores[~scores["is_prespecified_route"]]["source_equal_route_rate"].iloc[0]),
                "margin_vs_best_alternative": float(scores["prespecified_margin_vs_best_alternative"].iloc[0]),
                "collapsed_prespecified_rank": int(collapsed_scores[collapsed_scores["is_prespecified_route"]]["rank_by_source_equal"].iloc[0]),
                "collapsed_margin_vs_best_alternative": float(collapsed_scores["prespecified_margin_vs_best_alternative"].iloc[0]),
                "source_null_non_degenerate": bool(perm_source["null_non_degenerate"].iloc[0]),
                "source_layer_null_non_degenerate": bool(perm_layer["null_non_degenerate"].iloc[0]),
            }
        ]
    )

    out_xlsx = args.out / f"{args.case_id}_route_graph_results.xlsx"
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        preflight.to_excel(writer, sheet_name="preflight_report", index=False)
        scores.to_excel(writer, sheet_name="route_scores", index=False)
        matrix.to_excel(writer, sheet_name="donor_role_matrix", index=False)
        source.to_excel(writer, sheet_name="source_rates", index=False)
        leave_source.to_excel(writer, sheet_name="leave_one_source", index=False)
        leave_layer.to_excel(writer, sheet_name="leave_one_layer", index=False)
        dependency_sensitivity.to_excel(writer, sheet_name="dependency_sensitivity", index=False)
        for group_col, table in taxon_holdouts.items():
            table.to_excel(writer, sheet_name=f"leave_one_{group_col}"[:31], index=False)
        collapsed.to_excel(writer, sheet_name="dependency_collapsed_units", index=False)
        collapsed_eligible.to_excel(writer, sheet_name="collapsed_eligible", index=False)
        collapsed_scores.to_excel(writer, sheet_name="dependency_collapsed_scores", index=False)
        pd.concat([perm_source, perm_layer], ignore_index=True).to_excel(writer, sheet_name="permutation_nulls", index=False)

    for name, table in {
        "summary": summary,
        "preflight_report": preflight,
        "route_scores": scores,
        "donor_role_matrix": matrix,
        "source_rates": source,
        "leave_one_source": leave_source,
        "leave_one_layer": leave_layer,
        "dependency_sensitivity": dependency_sensitivity,
        "dependency_collapsed_units": collapsed,
        "collapsed_eligible": collapsed_eligible,
        "dependency_collapsed_scores": collapsed_scores,
        "permutation_nulls": pd.concat([perm_source, perm_layer], ignore_index=True),
    }.items():
        table.to_csv(args.out / f"{name}.csv", index=False)

    for group_col, table in taxon_holdouts.items():
        table.to_csv(args.out / f"leave_one_{group_col}.csv", index=False)
    for dependency_col, (collapsed_alt, eligible_alt, scores_alt) in additional_dependency_tables.items():
        collapsed_alt.to_csv(args.out / f"collapsed_{dependency_col}.csv", index=False)
        eligible_alt.to_csv(args.out / f"collapsed_{dependency_col}_eligible.csv", index=False)
        scores_alt.to_csv(args.out / f"collapsed_{dependency_col}_scores.csv", index=False)

    readme = args.out / "README_results.md"
    readme.write_text(
        "\n".join(
            [
                f"# Route-graph results for {args.case_id}",
                "",
                f"Input: `{args.input}`",
                f"Route-eligible units: {len(work)}",
                f"Sources: {work['source_id'].nunique()}",
                f"Evidence layers: {work['evidence_layer'].nunique()}",
                f"Dependency groups: {work['dependency_group'].nunique()}",
                "",
                "Primary output workbook:",
                f"- `{out_xlsx.name}`",
                "",
                "Primary metrics are source-equal route rank, margin over the nearest alternative, dependency-collapse retention and source/source-layer permutation nulls.",
                "",
                "Preflight gates are included to distinguish a manuscript-ready transfer case from a small illustrative fixture.",
            ]
        ),
        encoding="utf-8",
    )

    print(summary.to_string(index=False))
    print(f"Wrote {out_xlsx}")


if __name__ == "__main__":
    main()
