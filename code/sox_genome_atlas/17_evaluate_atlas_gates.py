from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def assess(
    label: str, root: Path, continuous_root: Path | None = None
) -> tuple[dict[str, object], list[dict[str, object]]]:
    summary = pd.read_csv(root / "summary.csv").iloc[0]
    preflight = pd.read_csv(root / "preflight_report.csv")
    nulls = pd.read_csv(root / "permutation_nulls.csv")
    dependency = read_optional(root / "dependency_sensitivity.csv")
    leave_layer = read_optional(root / "leave_one_layer.csv")
    leave_class = read_optional(root / "leave_one_gtdb_class.csv")
    leave_order = read_optional(root / "leave_one_gtdb_order.csv")

    layer_null = nulls.loc[nulls["strata"].eq("source_id+evidence_layer")].iloc[0]
    preflight_pass = bool(preflight["pass"].map(as_bool).all())
    null_non_degenerate = as_bool(layer_null["null_non_degenerate"])
    null_exceeded = bool(
        summary["margin_vs_best_alternative"] > layer_null["null_q975_margin"]
        and layer_null["empirical_p_ge_observed"] <= 0.05
    )
    dependency_rank1_fraction = (
        float(dependency["prespecified_rank"].eq(1).mean()) if not dependency.empty else float("nan")
    )
    layer_rank1_fraction = (
        float(leave_layer["prespecified_rank"].eq(1).mean()) if not leave_layer.empty else float("nan")
    )
    class_rank1_fraction = (
        float(leave_class["prespecified_rank"].eq(1).mean()) if not leave_class.empty else float("nan")
    )
    order_rank1_fraction = (
        float(leave_order["prespecified_rank"].eq(1).mean()) if not leave_order.empty else float("nan")
    )
    dependency_all_rank1 = bool(
        not dependency.empty
        and dependency["prespecified_rank"].eq(1).all()
        and dependency["margin_vs_best_alternative"].gt(0).all()
    )
    continuous = pd.DataFrame()
    if continuous_root is not None:
        continuous = read_optional(continuous_root / "continuous_mapping_summary.tsv")
    continuous_all_argmax = False
    continuous_exceeds_null = False
    continuous_support = float("nan")
    continuous_gap = float("nan")
    if not continuous.empty:
        record = continuous.iloc[0]
        continuous_all_argmax = as_bool(record["all_prespecified_roles_are_rowwise_argmax"])
        continuous_exceeds_null = as_bool(record["prespecified_support_exceeds_null_q975"])
        continuous_support = float(record["mean_prespecified_role_probability"])
        continuous_gap = float(record["prespecified_to_unconstrained_gap"])

    metrics = {
        "model": label,
        "preflight_pass": preflight_pass,
        "prespecified_rank": int(summary["prespecified_rank"]),
        "margin": float(summary["margin_vs_best_alternative"]),
        "collapsed_rank": int(summary["collapsed_prespecified_rank"]),
        "collapsed_margin": float(summary["collapsed_margin_vs_best_alternative"]),
        "source_layer_null_non_degenerate": null_non_degenerate,
        "source_layer_null_q975": float(layer_null["null_q975_margin"]),
        "source_layer_null_p": float(layer_null["empirical_p_ge_observed"]),
        "observed_exceeds_source_layer_null": null_exceeded,
        "dependency_rank1_fraction": dependency_rank1_fraction,
        "dependency_all_rank1_positive_margin": dependency_all_rank1,
        "leave_one_layer_rank1_fraction": layer_rank1_fraction,
        "leave_one_class_rank1_fraction": class_rank1_fraction,
        "leave_one_order_rank1_fraction": order_rank1_fraction,
        "continuous_mean_prespecified_probability": continuous_support,
        "continuous_unconstrained_gap": continuous_gap,
        "continuous_all_prespecified_roles_rowwise_argmax": continuous_all_argmax,
        "continuous_support_exceeds_null_q975": continuous_exceeds_null,
    }
    gates = [
        {"model": label, "gate": "preflight", "pass": preflight_pass},
        {"model": label, "gate": "prespecified_rank_1", "pass": int(summary["prespecified_rank"]) == 1},
        {"model": label, "gate": "collapsed_rank_1", "pass": int(summary["collapsed_prespecified_rank"]) == 1},
        {"model": label, "gate": "source_layer_null_non_degenerate", "pass": null_non_degenerate},
        {"model": label, "gate": "observed_exceeds_source_layer_null", "pass": null_exceeded},
        {"model": label, "gate": "all_dependency_schemes_rank_1", "pass": dependency_all_rank1},
        {"model": label, "gate": "continuous_all_rows_argmax", "pass": continuous_all_argmax},
        {"model": label, "gate": "continuous_exceeds_null_q975", "pass": continuous_exceeds_null},
    ]
    return metrics, gates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary", type=Path)
    parser.add_argument("--fixed", type=Path)
    parser.add_argument("--primary-continuous", type=Path)
    parser.add_argument("--fixed-continuous", type=Path)
    parser.add_argument("--decision-spec", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    models: list[tuple[str, Path, Path | None]] = []
    if args.primary and (args.primary / "summary.csv").exists():
        models.append(("ModelFinder", args.primary, args.primary_continuous))
    if args.fixed and (args.fixed / "summary.csv").exists():
        models.append(("LG+F+R6", args.fixed, args.fixed_continuous))
    if not models:
        raise FileNotFoundError("No completed route-result directory was supplied.")

    metrics_rows: list[dict[str, object]] = []
    gate_rows: list[dict[str, object]] = []
    for label, root, continuous_root in models:
        metrics, gates = assess(label, root, continuous_root)
        metrics_rows.append(metrics)
        gate_rows.extend(gates)
    metrics_frame = pd.DataFrame(metrics_rows)
    gate_frame = pd.DataFrame(gate_rows)

    core_gate_names = {
        "preflight",
        "prespecified_rank_1",
        "collapsed_rank_1",
        "source_layer_null_non_degenerate",
        "observed_exceeds_source_layer_null",
    }
    core_by_model = (
        gate_frame.loc[gate_frame["gate"].isin(core_gate_names)]
        .groupby("model")["pass"]
        .all()
    )
    all_core_gates = bool(core_by_model.all())
    any_core_model = bool(core_by_model.any())
    model_concordance = bool(
        len(metrics_frame) == 2
        and metrics_frame["prespecified_rank"].eq(1).all()
        and metrics_frame["margin"].gt(0).all()
        and metrics_frame["collapsed_margin"].gt(0).all()
    )
    dependency_concordance = bool(
        metrics_frame["dependency_all_rank1_positive_margin"].all()
    )
    continuous_concordance = bool(
        metrics_frame["continuous_all_prespecified_roles_rowwise_argmax"].all()
        and metrics_frame["continuous_support_exceeds_null_q975"].all()
    )

    if len(metrics_frame) < 2:
        decision = "PRELIMINARY_SUPPORT" if any_core_model else "PRELIMINARY_BOUNDARY"
    elif all_core_gates and model_concordance and dependency_concordance and continuous_concordance:
        decision = "PRIMARY_SUPPORT"
    elif any_core_model:
        decision = "CONDITIONAL_SUPPORT"
    else:
        decision = "BOUNDARY_RESULT"

    decision_record = {
        "decision": decision,
        "models_available": metrics_frame["model"].tolist(),
        "all_core_gates_pass": all_core_gates,
        "any_tree_model_passes_core_gates": any_core_model,
        "model_concordance": model_concordance,
        "dependency_concordance": dependency_concordance,
        "continuous_mapping_concordance": continuous_concordance,
    }
    if args.decision_spec:
        import hashlib

        spec_bytes = args.decision_spec.read_bytes()
        decision_record["decision_spec_sha256"] = hashlib.sha256(spec_bytes).hexdigest()
        (args.outdir / "decision_spec_used.json").write_bytes(spec_bytes)
    metrics_frame.to_csv(args.outdir / "atlas_model_metrics.tsv", sep="\t", index=False)
    gate_frame.to_csv(args.outdir / "atlas_gate_results.tsv", sep="\t", index=False)
    (args.outdir / "atlas_decision.json").write_text(
        json.dumps(decision_record, indent=2), encoding="utf-8"
    )
    lines = [
        "# SOX atlas decision",
        "",
        f"Decision: **{decision}**",
        "",
        f"Model concordance: {model_concordance}",
        f"All core gates pass: {all_core_gates}",
        "",
        "```text",
        metrics_frame.to_string(index=False),
        "```",
        "",
        "```text",
        gate_frame.to_string(index=False),
        "```",
    ]
    (args.outdir / "README_atlas_decision.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(decision_record, indent=2))


if __name__ == "__main__":
    main()
