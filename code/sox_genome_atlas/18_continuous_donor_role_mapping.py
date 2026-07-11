from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def normalize_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def conditional_matrix(
    counts: np.ndarray,
) -> np.ndarray:
    """Average P(role|donor, source) across sources where the donor is observed."""
    output = np.full((counts.shape[1], counts.shape[2]), np.nan, dtype=float)
    for donor_index in range(counts.shape[1]):
        donor_counts = counts[:, donor_index, :]
        denominator = donor_counts.sum(axis=1)
        observed_sources = denominator > 0
        if observed_sources.any():
            output[donor_index] = np.mean(
                donor_counts[observed_sources] / denominator[observed_sources, None], axis=0
            )
    return output


def source_equal_joint(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=(1, 2), keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        joint = counts / totals
    joint[~np.isfinite(joint)] = np.nan
    return np.nanmean(joint, axis=0)


def mutual_information(joint: np.ndarray) -> float:
    joint = np.asarray(joint, dtype=float)
    joint = joint / joint.sum()
    row = joint.sum(axis=1, keepdims=True)
    col = joint.sum(axis=0, keepdims=True)
    expected = row @ col
    mask = (joint > 0) & (expected > 0)
    return float(np.sum(joint[mask] * np.log(joint[mask] / expected[mask])))


def metrics_from_counts(counts: np.ndarray, expected_role_index: np.ndarray) -> dict[str, object]:
    conditional = conditional_matrix(counts)
    joint = source_equal_joint(counts)
    expected = conditional[np.arange(conditional.shape[0]), expected_role_index]
    valid_rows = np.isfinite(conditional).any(axis=1)
    best = np.full(conditional.shape[0], np.nan, dtype=float)
    best_index = np.full(conditional.shape[0], -1, dtype=np.int16)
    if valid_rows.any():
        safe = np.where(np.isfinite(conditional[valid_rows]), conditional[valid_rows], -np.inf)
        best[valid_rows] = np.max(safe, axis=1)
        best_index[valid_rows] = np.argmax(safe, axis=1)
    entropy = np.full(conditional.shape[0], np.nan, dtype=float)
    for index, row in enumerate(conditional):
        positive = row[row > 0]
        if positive.size:
            entropy[index] = -np.sum(positive * np.log(positive)) / np.log(len(row))
    return {
        "conditional": conditional,
        "joint": joint,
        "mean_prespecified_role_probability": float(np.nanmean(expected)),
        "mean_unconstrained_best_role_probability": float(np.nanmean(best)),
        "prespecified_to_unconstrained_gap": float(np.nanmean(best - expected)),
        "number_of_prespecified_rowwise_argmax_roles": int(
            np.sum(valid_rows & (best_index == expected_role_index))
        ),
        "source_equal_mutual_information": mutual_information(joint),
        "row_entropy": entropy,
        "best_role_index": best_index,
    }


def counts_from_codes(
    source_codes: np.ndarray,
    donor_codes: np.ndarray,
    role_codes: np.ndarray,
    n_sources: int,
    n_donors: int,
    n_roles: int,
) -> np.ndarray:
    flat = source_codes * (n_donors * n_roles) + donor_codes * n_roles + role_codes
    return np.bincount(flat, minlength=n_sources * n_donors * n_roles).reshape(
        n_sources, n_donors, n_roles
    )


def quantile_record(values: np.ndarray, prefix: str) -> dict[str, float]:
    return {
        f"{prefix}_q025": float(np.nanquantile(values, 0.025)),
        f"{prefix}_median": float(np.nanquantile(values, 0.5)),
        f"{prefix}_q975": float(np.nanquantile(values, 0.975)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--permutations", type=int)
    parser.add_argument("--bootstraps", type=int)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    spec_bytes = args.spec.read_bytes()
    spec = json.loads(spec_bytes)
    spec_hash = hashlib.sha256(spec_bytes).hexdigest()
    permutations = args.permutations or int(spec["permutation_iterations"])
    bootstraps = args.bootstraps or int(spec["bootstrap_iterations"])
    seed = int(spec["random_seed"])
    rng = np.random.default_rng(seed)

    donors = list(spec["donors"])
    roles = list(spec["roles"])
    donor_index = {value: index for index, value in enumerate(donors)}
    role_index = {value: index for index, value in enumerate(roles)}
    expected_role_index = np.array(
        [role_index[spec["prespecified_mapping"][donor]] for donor in donors], dtype=np.int16
    )

    frame = pd.read_csv(args.input, dtype=str).fillna("")
    frame = frame.loc[frame["case_id"].eq(spec["case_id"])].copy()
    frame = frame.loc[normalize_bool(frame["route_eligible"])].copy()
    frame = frame.loc[
        frame["donor_class"].isin(donors) & frame["functional_class"].isin(roles)
    ].copy()
    if frame.empty:
        raise ValueError("No route-eligible rows remain for continuous mapping.")

    source_labels = sorted(frame["source_id"].unique())
    source_index = {value: index for index, value in enumerate(source_labels)}
    source_codes = frame["source_id"].map(source_index).to_numpy(dtype=np.int16)
    donor_codes = frame["donor_class"].map(donor_index).to_numpy(dtype=np.int16)
    role_codes = frame["functional_class"].map(role_index).to_numpy(dtype=np.int16)
    observed_counts = counts_from_codes(
        source_codes, donor_codes, role_codes, len(source_labels), len(donors), len(roles)
    )
    observed = metrics_from_counts(observed_counts, expected_role_index)

    strata_columns = list(spec["permutation_strata"])
    grouper = strata_columns[0] if len(strata_columns) == 1 else strata_columns
    strata_indices = [
        np.asarray(index, dtype=int)
        for index in frame.groupby(grouper, sort=True).indices.values()
    ]
    permutation_support = np.empty(permutations, dtype=float)
    permutation_mi = np.empty(permutations, dtype=float)
    permutation_row_argmax = np.empty(permutations, dtype=np.int16)
    for iteration in range(permutations):
        permuted_roles = role_codes.copy()
        for index in strata_indices:
            permuted_roles[index] = rng.permutation(permuted_roles[index])
        counts = counts_from_codes(
            source_codes,
            donor_codes,
            permuted_roles,
            len(source_labels),
            len(donors),
            len(roles),
        )
        metric = metrics_from_counts(counts, expected_role_index)
        permutation_support[iteration] = metric["mean_prespecified_role_probability"]
        permutation_mi[iteration] = metric["source_equal_mutual_information"]
        permutation_row_argmax[iteration] = metric[
            "number_of_prespecified_rowwise_argmax_roles"
        ]

    dependency_column = str(spec["bootstrap_dependency_column"])
    if dependency_column not in frame.columns:
        raise ValueError(f"Missing bootstrap dependency column: {dependency_column}")
    group_tensors: dict[int, list[np.ndarray]] = {}
    for source, source_frame in frame.groupby("source_id", sort=True):
        tensors: list[np.ndarray] = []
        for _, group in source_frame.groupby(dependency_column, sort=True):
            tensor = np.zeros((len(donors), len(roles)), dtype=np.int64)
            d = group["donor_class"].map(donor_index).to_numpy(dtype=np.int16)
            r = group["functional_class"].map(role_index).to_numpy(dtype=np.int16)
            np.add.at(tensor, (d, r), 1)
            tensors.append(tensor)
        group_tensors[source_index[source]] = tensors

    bootstrap_conditional = np.empty(
        (bootstraps, len(donors), len(roles)), dtype=float
    )
    bootstrap_support = np.empty(bootstraps, dtype=float)
    bootstrap_gap = np.empty(bootstraps, dtype=float)
    bootstrap_mi = np.empty(bootstraps, dtype=float)
    for iteration in range(bootstraps):
        counts = np.zeros_like(observed_counts)
        for source_code, tensors in group_tensors.items():
            tensor_stack = np.stack(tensors)
            multiplicity = rng.multinomial(
                len(tensors), np.full(len(tensors), 1.0 / len(tensors))
            )
            counts[source_code] = np.tensordot(multiplicity, tensor_stack, axes=(0, 0))
        metric = metrics_from_counts(counts, expected_role_index)
        bootstrap_conditional[iteration] = metric["conditional"]
        bootstrap_support[iteration] = metric["mean_prespecified_role_probability"]
        bootstrap_gap[iteration] = metric["prespecified_to_unconstrained_gap"]
        bootstrap_mi[iteration] = metric["source_equal_mutual_information"]

    probability_rows: list[dict[str, object]] = []
    for donor_code, donor in enumerate(donors):
        for role_code, role in enumerate(roles):
            draws = bootstrap_conditional[:, donor_code, role_code]
            probability_rows.append(
                {
                    "donor_class": donor,
                    "functional_class": role,
                    "source_equal_probability": float(observed["conditional"][donor_code, role_code]),
                    "bootstrap_q025": float(np.nanquantile(draws, 0.025)),
                    "bootstrap_median": float(np.nanquantile(draws, 0.5)),
                    "bootstrap_q975": float(np.nanquantile(draws, 0.975)),
                    "is_prespecified_cell": role_code == expected_role_index[donor_code],
                }
            )

    profile_rows: list[dict[str, object]] = []
    for donor_code, donor in enumerate(donors):
        probabilities = observed["conditional"][donor_code]
        expected_code = int(expected_role_index[donor_code])
        order = np.argsort(-probabilities)
        profile_rows.append(
            {
                "donor_class": donor,
                "prespecified_role": roles[expected_code],
                "prespecified_probability": float(probabilities[expected_code]),
                "prespecified_rank_unconstrained": int(np.flatnonzero(order == expected_code)[0] + 1),
                "best_unconstrained_role": roles[int(observed["best_role_index"][donor_code])],
                "best_unconstrained_probability": float(np.nanmax(probabilities)),
                "normalized_role_entropy": float(observed["row_entropy"][donor_code]),
            }
        )

    support_q975 = float(np.quantile(permutation_support, 0.975))
    mi_q975 = float(np.quantile(permutation_mi, 0.975))
    summary = {
        "case_id": spec["case_id"],
        "n_route_eligible_units": int(len(frame)),
        "n_sources": int(len(source_labels)),
        "n_evidence_layers": int(frame["evidence_layer"].nunique()),
        "n_bootstrap_dependency_groups": int(frame[dependency_column].nunique()),
        "analysis_spec_sha256": spec_hash,
        "mean_prespecified_role_probability": observed[
            "mean_prespecified_role_probability"
        ],
        "mean_unconstrained_best_role_probability": observed[
            "mean_unconstrained_best_role_probability"
        ],
        "prespecified_to_unconstrained_gap": observed[
            "prespecified_to_unconstrained_gap"
        ],
        "number_of_prespecified_rowwise_argmax_roles": observed[
            "number_of_prespecified_rowwise_argmax_roles"
        ],
        "source_equal_mutual_information": observed["source_equal_mutual_information"],
        "permutation_support_q975": support_q975,
        "permutation_support_p": float(
            (np.sum(permutation_support >= observed["mean_prespecified_role_probability"]) + 1)
            / (permutations + 1)
        ),
        "permutation_mi_q975": mi_q975,
        "permutation_mi_p": float(
            (np.sum(permutation_mi >= observed["source_equal_mutual_information"]) + 1)
            / (permutations + 1)
        ),
        "all_prespecified_roles_are_rowwise_argmax": bool(
            observed["number_of_prespecified_rowwise_argmax_roles"] == len(donors)
        ),
        "prespecified_support_exceeds_null_q975": bool(
            observed["mean_prespecified_role_probability"] > support_q975
        ),
    }
    summary.update(quantile_record(bootstrap_support, "bootstrap_prespecified_support"))
    summary.update(quantile_record(bootstrap_gap, "bootstrap_unconstrained_gap"))
    summary.update(quantile_record(bootstrap_mi, "bootstrap_mutual_information"))

    pd.DataFrame(probability_rows).to_csv(
        args.outdir / "continuous_donor_role_probabilities.tsv", sep="\t", index=False
    )
    pd.DataFrame(profile_rows).to_csv(
        args.outdir / "continuous_donor_profiles.tsv", sep="\t", index=False
    )
    pd.DataFrame([summary]).to_csv(
        args.outdir / "continuous_mapping_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "iteration": np.arange(1, permutations + 1),
            "mean_prespecified_role_probability": permutation_support,
            "source_equal_mutual_information": permutation_mi,
            "prespecified_rowwise_argmax_roles": permutation_row_argmax,
        }
    ).to_csv(args.outdir / "continuous_permutation_draws.tsv.gz", sep="\t", index=False)
    args.outdir.joinpath("analysis_spec_used.json").write_bytes(spec_bytes)
    args.outdir.joinpath("analysis_spec_sha256.txt").write_text(spec_hash + "\n", encoding="utf-8")
    args.outdir.joinpath("CONTINUOUS_MAPPING_COMPLETE").write_text("complete\n", encoding="utf-8")
    print(pd.DataFrame([summary]).to_string(index=False))
    print(pd.DataFrame(profile_rows).to_string(index=False))


if __name__ == "__main__":
    main()
