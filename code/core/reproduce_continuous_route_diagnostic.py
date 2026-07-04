from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "eukaryogenesis_prior_release" / "core_metrics"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def as_float(value: str) -> float:
    return float(value)


def main() -> None:
    top_roles = read_csv(DATA / "continuous_route_assignment_diagnostic.csv")
    bootstrap = read_csv(DATA / "continuous_route_bootstrap_diagnostic.csv")

    print("Continuous donor-role diagnostic")
    for row in top_roles:
        donor = row["donor_class"]
        prespecified = row["prespecified_functional_class"]
        top = row["top_source_equal_functional_class"]
        probability = as_float(row["top_source_equal_probability"])
        margin = as_float(row["margin_top_minus_next"])
        ok = row["prespecified_is_top"].lower() == "yes"
        print(
            f"{donor}: prespecified={prespecified}; "
            f"top={top}; source_equal_probability={probability:.6f}; "
            f"margin_top_minus_next={margin:.6f}; prespecified_is_top={ok}"
        )
        assert ok
        assert margin > 0

    joint = [r for r in bootstrap if r["donor_class"] == "all_three_donors"][0]
    p_joint = as_float(joint["probability_prespecified_role_is_top"])
    print(f"joint_bootstrap_probability_all_three_prespecified_top={p_joint:.6f}")
    assert p_joint > 0.99


if __name__ == "__main__":
    main()
