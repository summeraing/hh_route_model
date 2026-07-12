#!/usr/bin/env python3
"""Known-truth calibration for the route-resolved evidence graph."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd


ROUTES = np.asarray(list(itertools.permutations(range(3))), dtype=np.int8)
TRUE_ROUTE = np.asarray([0, 1, 2], dtype=np.int8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--scenario-index", type=int, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()


def scenario_catalog(spec: dict) -> list[tuple[str, float]]:
    out: list[tuple[str, float]] = []
    for architecture in ("null", "full", "partial"):
        for signal in spec["simulation"]["architectures"][architecture]:
            out.append((architecture, float(signal)))
    return out


def source_weights(n_sources: int, dominant_fraction: float) -> np.ndarray:
    if n_sources == 1:
        return np.ones(1)
    rest = (1.0 - dominant_fraction) / (n_sources - 1)
    return np.asarray([dominant_fraction] + [rest] * (n_sources - 1))


def draw_role(
    rng: np.random.Generator,
    donor: int,
    architecture: str,
    signal: float,
    source: int,
    source_regime: str,
) -> int:
    if architecture == "null":
        return int(rng.integers(0, 3))

    active = architecture == "full" or donor in (0, 1)
    if not active:
        return int(rng.integers(0, 3))

    target = donor
    if source_regime == "dominant_adversarial" and source == 0:
        target = int((1, 0, 2)[donor])
    match_probability = 1.0 / 3.0 + (2.0 / 3.0) * signal
    if rng.random() < match_probability:
        return target
    alternatives = [x for x in range(3) if x != target]
    return int(alternatives[int(rng.integers(0, 2))])


def flip_labels(rng: np.random.Generator, labels: np.ndarray, rate: float) -> np.ndarray:
    if rate <= 0:
        return labels.copy()
    out = labels.copy()
    mask = rng.random(len(out)) < rate
    if not mask.any():
        return out
    offsets = rng.integers(1, 3, size=int(mask.sum()))
    out[mask] = (out[mask] + offsets) % 3
    return out


def source_equal_scores(
    donor: np.ndarray,
    role: np.ndarray,
    source: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    scores = []
    present_sources = np.unique(source)
    if weights is None:
        weights = np.ones(len(donor), dtype=float)
    for route in ROUTES:
        match = role == route[donor]
        source_rates = []
        for src in present_sources:
            idx = source == src
            denom = float(weights[idx].sum())
            if denom > 0:
                source_rates.append(float(np.dot(match[idx].astype(float), weights[idx]) / denom))
        scores.append(float(np.mean(source_rates)))
    return np.asarray(scores)


def raw_scores(donor: np.ndarray, role: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    if weights is None:
        weights = np.ones(len(donor), dtype=float)
    denom = float(weights.sum())
    return np.asarray([
        float(np.dot((role == route[donor]).astype(float), weights) / denom)
        for route in ROUTES
    ])


def donor_role_profile(
    donor: np.ndarray,
    role: np.ndarray,
    source: np.ndarray,
    source_equal: bool,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    if weights is None:
        weights = np.ones(len(donor), dtype=float)
    profile = np.zeros((3, 3), dtype=float)
    for d in range(3):
        if source_equal:
            per_source = []
            for src in np.unique(source):
                idx = (source == src) & (donor == d)
                denom = float(weights[idx].sum())
                if denom > 0:
                    per_source.append([
                        float(weights[idx & (role == r)].sum() / denom) for r in range(3)
                    ])
            if per_source:
                profile[d] = np.mean(np.asarray(per_source), axis=0)
        else:
            idx = donor == d
            denom = float(weights[idx].sum())
            if denom > 0:
                profile[d] = [float(weights[idx & (role == r)].sum() / denom) for r in range(3)]
    return profile


def summarize_method(
    method: str,
    donor: np.ndarray,
    role: np.ndarray,
    source: np.ndarray,
    group: np.ndarray,
) -> dict:
    if method == "raw_pooled":
        scores = raw_scores(donor, role)
        profile = donor_role_profile(donor, role, source, source_equal=False)
    elif method == "source_equal":
        scores = source_equal_scores(donor, role, source)
        profile = donor_role_profile(donor, role, source, source_equal=True)
    elif method == "source_equal_dependency_collapsed":
        frame = pd.DataFrame({"group": group, "source": source, "donor": donor, "role": role})
        collapsed = []
        for (_, src), chunk in frame.groupby(["group", "source"], sort=False):
            n = float(len(chunk))
            for d in range(3):
                for r in range(3):
                    count = float(((chunk["donor"] == d) & (chunk["role"] == r)).sum())
                    if count:
                        collapsed.append((int(src), d, r, count / n))
        c = pd.DataFrame(collapsed, columns=["source", "donor", "role", "weight"])
        scores = source_equal_scores(
            c["donor"].to_numpy(dtype=int),
            c["role"].to_numpy(dtype=int),
            c["source"].to_numpy(dtype=int),
            c["weight"].to_numpy(dtype=float),
        )
        profile = donor_role_profile(
            c["donor"].to_numpy(dtype=int),
            c["role"].to_numpy(dtype=int),
            c["source"].to_numpy(dtype=int),
            source_equal=True,
            weights=c["weight"].to_numpy(dtype=float),
        )
    else:
        raise ValueError(method)

    prespecified_index = int(np.where((ROUTES == TRUE_ROUTE).all(axis=1))[0][0])
    prespecified_score = float(scores[prespecified_index])
    alternatives = np.delete(scores, prespecified_index)
    best_alternative = float(alternatives.max())
    rank = 1 + int(np.sum(scores > prespecified_score + 1e-12))
    return {
        "method": method,
        "prespecified_score": prespecified_score,
        "best_alternative_score": best_alternative,
        "margin": prespecified_score - best_alternative,
        "prespecified_rank": rank,
        "full_components_correct": int(np.all(np.argmax(profile, axis=1) == TRUE_ROUTE)),
        "partial_components_correct": int(
            np.argmax(profile[0]) == 0 and np.argmax(profile[1]) == 1
        ),
        "profile_d0_r0": float(profile[0, 0]),
        "profile_d1_r1": float(profile[1, 1]),
        "profile_d2_r2": float(profile[2, 2]),
    }


def run_combo(payload: tuple) -> list[dict]:
    (
        architecture,
        signal,
        dominance,
        mean_group_size,
        noise,
        source_regime,
        n_sources,
        n_groups,
        replicates,
        methods,
        base_seed,
        combo_index,
    ) = payload
    rows: list[dict] = []
    weights = source_weights(n_sources, dominance)
    for rep in range(replicates):
        rng = np.random.default_rng(base_seed + combo_index * 100000 + rep)
        group_source = rng.choice(n_sources, size=n_groups, p=weights)
        group_donor = rng.integers(0, 3, size=n_groups)
        group_role = np.asarray([
            draw_role(rng, int(d), architecture, signal, int(s), source_regime)
            for d, s in zip(group_donor, group_source)
        ], dtype=np.int8)
        if mean_group_size <= 1:
            group_size = np.ones(n_groups, dtype=int)
        else:
            group_size = 1 + rng.poisson(mean_group_size - 1, size=n_groups)
            group_size = np.clip(group_size, 1, 50)
        group = np.repeat(np.arange(n_groups, dtype=int), group_size)
        source = np.repeat(group_source, group_size)
        donor = flip_labels(rng, np.repeat(group_donor, group_size), noise)
        role = flip_labels(rng, np.repeat(group_role, group_size), noise)

        for method in methods:
            summary = summarize_method(method, donor, role, source, group)
            summary.update({
                "architecture": architecture,
                "signal_strength": signal,
                "dominant_source_fraction": dominance,
                "mean_dependency_group_size": mean_group_size,
                "label_noise_fraction": noise,
                "source_regime": source_regime,
                "replicate": rep,
                "n_rows": int(len(group)),
                "n_dependency_groups": n_groups,
                "seed": base_seed + combo_index * 100000 + rep,
            })
            rows.append(summary)
    return rows


def main() -> None:
    args = parse_args()
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    sim = spec["simulation"]
    catalog = scenario_catalog(spec)
    if not 0 <= args.scenario_index < len(catalog):
        raise SystemExit(f"scenario-index must be 0..{len(catalog)-1}")
    architecture, signal = catalog[args.scenario_index]

    regimes = ["homogeneous"] if architecture == "null" else sim["source_regime"]
    combos = []
    combo_index = 0
    for dominance, group_size, noise, regime in itertools.product(
        sim["dominant_source_fraction"],
        sim["mean_dependency_group_size"],
        sim["label_noise_fraction"],
        regimes,
    ):
        combos.append((
            architecture,
            signal,
            float(dominance),
            int(group_size),
            float(noise),
            regime,
            int(sim["sources"]),
            int(sim["dependency_groups"]),
            int(sim["replicates_per_scenario"]),
            list(sim["methods"]),
            int(sim["seed"]) + args.scenario_index * 100000000,
            combo_index,
        ))
        combo_index += 1

    all_rows: list[dict] = []
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for result in pool.map(run_combo, combos):
                all_rows.extend(result)
    else:
        for combo in combos:
            all_rows.extend(run_combo(combo))

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(all_rows)
    frame.to_csv(output, index=False)
    metadata = {
        "spec_id": spec["spec_id"],
        "scenario_index": args.scenario_index,
        "architecture": architecture,
        "signal_strength": signal,
        "rows": len(frame),
        "output": str(output),
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
