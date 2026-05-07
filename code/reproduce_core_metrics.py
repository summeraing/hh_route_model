from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "core_metrics"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def as_float(value: str) -> float:
    return float(value)


def main() -> None:
    routes = read_csv(DATA / "route_map_specificity_adjudicated.csv")
    prespecified = [r for r in routes if r.get("is_prespecified_route", "").lower() == "true"][0]
    support = as_float(prespecified["source_equal_rate"])
    rank = int(float(prespecified["rank_by_source_equal"]))

    nulls = {r["metric"]: as_float(r["value"]) for r in read_csv(DATA / "frozen_null_boundary.csv")}
    boundary = nulls["frozen_q975_boundary"]
    margin = nulls["margin_over_q975"]

    dep = read_csv(DATA / "dependency_collapse.csv")
    dep_margins = [
        as_float(r["true_minus_best_alt_source_equal"])
        for r in dep
        if r.get("strategy") != "row_level_reference"
    ]
    min_dep_margin = min(dep_margins) if dep_margins else float("nan")

    agreement = {r["axis"]: r for r in read_csv(DATA / "agreement_summary.csv")}
    donor = agreement["donor_class"]
    functional = agreement["functional_class"]

    print("Core audit-completed metrics")
    print(f"prespecified_route_rank={rank}")
    print(f"adjudicated_source_equal_support={support:.6f}")
    print(f"frozen_q97_5_null_boundary={boundary:.6f}")
    print(f"margin_over_q97_5={margin:.6f}")
    print(f"minimum_dependency_collapse_margin={min_dep_margin:.6f}")
    print(f"donor_percent_agreement={as_float(donor['percent_agreement']):.6f}")
    print(f"donor_cohen_kappa={as_float(donor['cohen_kappa']):.6f}")
    print(f"functional_percent_agreement={as_float(functional['percent_agreement']):.6f}")
    print(f"functional_cohen_kappa={as_float(functional['cohen_kappa']):.6f}")

    assert rank == 1
    assert abs(support - 0.76945) < 1e-6
    assert support > boundary
    assert margin > 0
    assert min_dep_margin > 0


if __name__ == "__main__":
    main()
