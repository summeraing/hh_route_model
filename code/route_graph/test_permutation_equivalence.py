from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_route_graph_case as route_graph  # noqa: E402


def reference_null(
    frame: pd.DataFrame,
    donors: list[str],
    roles: list[str],
    prespec: dict[str, str],
    strata: list[str],
    iterations: int,
    seed: int,
) -> tuple[float, float, float, float, float, int]:
    rng = np.random.default_rng(seed)
    observed = route_graph.route_scores(frame, donors, roles, prespec)
    pres = observed.loc[observed["is_prespecified_route"]].iloc[0]
    alt = observed.loc[~observed["is_prespecified_route"]].iloc[0]
    observed_margin = float(pres["source_equal_route_rate"] - alt["source_equal_route_rate"])
    margins = []
    for _ in range(iterations):
        shuffled = frame.copy()
        role_values = shuffled["functional_class"].to_numpy().copy()
        group_key = strata[0] if len(strata) == 1 else strata
        for _, index in shuffled.groupby(group_key).groups.items():
            positions = np.asarray(list(index), dtype=int)
            role_values[positions] = rng.permutation(role_values[positions])
        shuffled["functional_class"] = role_values
        scores = route_graph.route_scores(shuffled, donors, roles, prespec)
        pres = scores.loc[scores["is_prespecified_route"]].iloc[0]
        alt = scores.loc[~scores["is_prespecified_route"]].iloc[0]
        margins.append(float(pres["source_equal_route_rate"] - alt["source_equal_route_rate"]))
    values = np.asarray(margins)
    return (
        observed_margin,
        float(values.mean()),
        float(values.std(ddof=1)),
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
        int(np.unique(np.round(values, 12)).size),
    )


def main() -> None:
    data = Path(__file__).resolve().parents[2] / "data" / "sox_transfer" / "sox_evidence_units_expanded_v2.csv"
    frame = pd.read_csv(data)
    frame["route_eligible"] = route_graph.normalize_bool(frame["route_eligible"])
    donors = ["lineage_core", "hgt_pathway", "mobile_context"]
    roles = ["conserved_metabolic_backbone", "variable_energy_module", "mobility_boundary"]
    prespec = {
        "lineage_core": "conserved_metabolic_backbone",
        "hgt_pathway": "variable_energy_module",
        "mobile_context": "mobility_boundary",
    }
    work = frame.loc[
        frame["route_eligible"]
        & frame["donor_class"].isin(donors)
        & frame["functional_class"].isin(roles)
    ].copy().reset_index(drop=True)

    for strata, seed in [(["source_id"], 123), (["source_id", "evidence_layer"], 124)]:
        vectorized = route_graph.permutation_null(
            work, donors, roles, prespec, strata, iterations=50, seed=seed
        ).iloc[0]
        observed = (
            vectorized["observed_margin"],
            vectorized["null_mean_margin"],
            vectorized["null_sd_margin"],
            vectorized["null_q025_margin"],
            vectorized["null_q975_margin"],
            vectorized["null_unique_margins"],
        )
        expected = reference_null(work, donors, roles, prespec, strata, 50, seed)
        if not np.allclose(observed[:5], expected[:5], atol=1e-12) or observed[5] != expected[5]:
            raise AssertionError(f"Permutation implementations differ for {strata}: {observed} != {expected}")
    print("VECTORIZED_PERMUTATION_EQUIVALENCE_PASS")


if __name__ == "__main__":
    main()
