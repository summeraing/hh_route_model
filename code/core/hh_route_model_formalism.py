"""
H-h route-partition formalism helpers for the manuscript analysis package.

This module is intentionally lightweight. It documents the equations used in
the manuscript and provides reproducible helper functions for route support,
source-equal averaging, H-h route costs and simple ABC-style distance scoring.

It is not a one-command rebuild of all public raw data.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import exp, sqrt
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


DONORS = ("host", "symbiont", "other")
ROUTES = ("scaffold", "injection", "transition")
PRESPECIFIED_ROUTE = {
    "host": "scaffold",
    "symbiont": "injection",
    "other": "transition",
}


@dataclass(frozen=True)
class EvidenceUnit:
    """Minimal row needed for route-support calculations."""

    source: str
    donor_class: str
    functional_class: str
    layer: str | None = None


@dataclass(frozen=True)
class NodeState:
    """State variables in the H-h energy-time-complexity formalism."""

    energy: float
    material: float
    info_efficiency: float
    local_complexity: float
    global_complexity: float
    h_mismatch: float


@dataclass(frozen=True)
class RouteParameters:
    """Route-specific changes entering the H-h cost equation."""

    delta_energy_burden: float
    delta_complexity: float
    mismatch_burden: float
    delta_info_efficiency: float
    long_term_yield: float = 0.0


@dataclass(frozen=True)
class CostWeights:
    """Weights for the route-cost function."""

    energy: float = 1.0
    complexity: float = 1.0
    mismatch: float = 1.0
    information: float = 1.0
    yield_: float = 1.0


def all_route_maps() -> List[Dict[str, str]]:
    """Return all six bijections from donor class to route class."""

    return [dict(zip(DONORS, routes)) for routes in permutations(ROUTES)]


def route_indicator(unit: EvidenceUnit, route_map: Mapping[str, str]) -> int:
    """Indicator I_u(pi)=1{f_u=pi(d_u)}."""

    return int(route_map.get(unit.donor_class) == unit.functional_class)


def raw_route_support(units: Sequence[EvidenceUnit], route_map: Mapping[str, str]) -> float:
    """Raw pooled route support."""

    if not units:
        raise ValueError("units must not be empty")
    return sum(route_indicator(u, route_map) for u in units) / len(units)


def source_equal_support(units: Sequence[EvidenceUnit], route_map: Mapping[str, str]) -> float:
    """Source-equal support: mean within-source route agreement."""

    by_source: Dict[str, List[EvidenceUnit]] = {}
    for unit in units:
        by_source.setdefault(unit.source, []).append(unit)
    if not by_source:
        raise ValueError("units must contain at least one source")
    per_source = [raw_route_support(rows, route_map) for rows in by_source.values()]
    return sum(per_source) / len(per_source)


def rank_route_maps(units: Sequence[EvidenceUnit]) -> List[Tuple[int, Dict[str, str], float]]:
    """Rank all donor-route bijections by source-equal support."""

    scored = [(route_map, source_equal_support(units, route_map)) for route_map in all_route_maps()]
    scored.sort(key=lambda x: x[1], reverse=True)
    return [(i + 1, route_map, score) for i, (route_map, score) in enumerate(scored)]


def net_information_throughput(
    information_flux: float,
    local_complexity: float,
    global_complexity: float,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> float:
    """Phi_net = Phi / (1 + alpha*C_local + beta*C_global)."""

    return information_flux / (1.0 + alpha * local_complexity + beta * global_complexity)


def node_energy_burden(
    baseline: float,
    maintenance_cost: float,
    local_complexity: float,
    edge_cost: float,
    total_edge_weight: float,
    yield_rate: float,
    energy: float,
    material: float,
) -> float:
    """Energetic/material burden for a node group."""

    return (
        baseline
        + maintenance_cost * local_complexity
        + edge_cost * total_edge_weight
        - yield_rate * energy * material
    )


def route_cost(params: RouteParameters, weights: CostWeights = CostWeights()) -> float:
    """Omega(r) for the coupled H-h model."""

    return (
        weights.energy * params.delta_energy_burden
        + weights.complexity * params.delta_complexity
        + weights.mismatch * params.mismatch_burden
        - weights.information * params.delta_info_efficiency
        - weights.yield_ * params.long_term_yield
    )


def route_probabilities(
    route_params: Mapping[str, RouteParameters],
    weights: CostWeights = CostWeights(),
    beta: float = 1.0,
) -> Dict[str, float]:
    """Softmax probability over route costs."""

    costs = {route: route_cost(params, weights) for route, params in route_params.items()}
    unnorm = {route: exp(-beta * cost) for route, cost in costs.items()}
    total = sum(unnorm.values())
    if total == 0:
        raise ValueError("softmax denominator is zero")
    return {route: value / total for route, value in unnorm.items()}


def abc_distance(
    observed: Mapping[str, float],
    simulated: Mapping[str, float],
    weights: Mapping[str, float] | None = None,
) -> float:
    """Weighted Euclidean ABC distance between observed and simulated summaries."""

    keys = sorted(observed.keys())
    if set(keys) != set(simulated.keys()):
        raise ValueError("observed and simulated summaries must have identical keys")
    weights = weights or {key: 1.0 for key in keys}
    return sqrt(sum(weights.get(key, 1.0) * (simulated[key] - observed[key]) ** 2 for key in keys))


def failure_boundary_margin(observed: float, null_q975: float) -> float:
    """Positive values indicate observed support exceeds the conservative null boundary."""

    return observed - null_q975


def example_observed_summary() -> Dict[str, float]:
    """Observed route summaries used in the analysis."""

    return {
        "strict_route": 0.704466,
        "host_scaffold": 0.801126,
        "symbiont_injection": 0.812500,
        "other_transition": 0.569444,
        "route_gap": 0.274957,
        "null_margin_min": 0.087490,
        "collapse_margin_min": 0.206667,
    }


def main() -> None:
    """Tiny smoke test and printed numeric lock."""

    observed = example_observed_summary()
    print("v76 observed summary")
    for key, value in observed.items():
        print(f"{key}: {value:.6f}")

    toy_params = {
        "scaffold": RouteParameters(0.10, 0.10, 0.05, 0.05, 0.00),
        "injection": RouteParameters(0.30, 0.20, 0.30, 0.60, 0.80),
        "transition": RouteParameters(0.25, 0.45, 0.40, 0.25, 0.20),
    }
    probs = route_probabilities(toy_params)
    print("\nexample coupled H-h route probabilities")
    for key, value in probs.items():
        print(f"{key}: {value:.3f}")


if __name__ == "__main__":
    main()

