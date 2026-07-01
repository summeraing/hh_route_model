from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, log_loss
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


ROOT = Path.cwd() / "TOP_JOURNAL_REBUILD_ANALYSES_v1" / "iMETA_ROUTE_B_REBUILD_20260630"
FULL_TABLE = (
    Path.cwd()
    / "TOP_JOURNAL_REBUILD_ANALYSES_v1"
    / "v77_clean"
    / "01_MAIN_FIGURES"
    / "Fig1"
    / "source_package"
    / "Fig1B_selected_original_network"
    / "source_combined_v4_flat_export.csv"
)
SOURCE_XLSX = ROOT / "10_WORKING_COPY_FROM_NC12" / "04_Source_Data" / "Source_Data_complete.xlsx"
OUT = ROOT / "02_ANALYSIS_UPGRADES" / "routeB_predictive_feature_ablation_v1"
OUT.mkdir(parents=True, exist_ok=True)

DONORS = ["host", "symbiont", "other"]
ROLES = ["scaffold", "energetic_incorporation", "transition"]
ROLE_ORDER = ROLES
RNG_SEED = 20260701

COL = {
    "blue": "#356B8A",
    "orange": "#E0942A",
    "purple": "#7A58A6",
    "grey": "#6B7280",
    "light": "#D8DEE6",
    "red": "#B95A55",
    "green": "#3B8D5A",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 7,
    }
)


def norm(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = s.replace("injection", "energetic_incorporation")
    s = s.replace("energy", "energetic_incorporation")
    return s


def load_full_with_audit_overlay() -> pd.DataFrame:
    df = pd.read_csv(FULL_TABLE)
    df["row_index"] = np.arange(len(df))
    df["donor_original"] = df["donor_class"].map(norm)
    df["functional_original"] = df["functional_class"].map(norm)
    df["donor_final"] = df["donor_original"]
    df["functional_final"] = df["functional_original"]

    audit = pd.read_excel(SOURCE_XLSX, sheet_name="audit_overlay")
    for _, row in audit.iterrows():
        idx = int(row["combined_index"])
        if 0 <= idx < len(df):
            df.loc[idx, "donor_final"] = norm(row["adjudicated_donor_class"])
            df.loc[idx, "functional_final"] = norm(row["adjudicated_functional_class"])

    work = df[df["donor_final"].isin(DONORS) & df["functional_final"].isin(ROLES)].copy()
    text_cols = ["evidence_layer", "evidence_unit", "unit_id", "unit_label", "module_family", "compartment", "source_table"]
    for col in text_cols:
        if col not in work.columns:
            work[col] = ""
    work["text_features"] = work[text_cols].fillna("").astype(str).agg(" | ".join, axis=1)
    work["source_layer"] = work["study_id"].astype(str) + "::" + work["evidence_layer"].astype(str)
    return work.reset_index(drop=True)


def make_pipeline(feature_set: str) -> Pipeline:
    transformers = []
    if feature_set in {"source_layer", "source_layer_text", "source_layer_donor", "source_layer_text_donor"}:
        transformers.append(("source_layer", OneHotEncoder(handle_unknown="ignore"), ["study_id", "evidence_layer"]))
    if feature_set in {"donor", "source_layer_donor", "source_layer_text_donor", "text_donor"}:
        transformers.append(("donor", OneHotEncoder(handle_unknown="ignore"), ["donor_final"]))
    if feature_set in {"text", "source_layer_text", "source_layer_text_donor", "text_donor"}:
        transformers.append(
            (
                "text",
                TfidfVectorizer(
                    lowercase=True,
                    ngram_range=(1, 2),
                    min_df=2,
                    max_features=1200,
                    token_pattern=r"(?u)\b[\w:/.-]+\b",
                ),
                "text_features",
            )
        )
    if not transformers:
        raise ValueError(feature_set)
    pre = ColumnTransformer(transformers=transformers, sparse_threshold=0.3)
    clf = LogisticRegression(
        max_iter=3000,
        class_weight="balanced",
        solver="liblinear",
        random_state=RNG_SEED,
    )
    return Pipeline([("features", pre), ("clf", clf)])


def split_iterator(df: pd.DataFrame, cv_design: str):
    y = df["functional_final"].to_numpy()
    if cv_design == "stratified_row":
        splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG_SEED)
        yield from splitter.split(df, y)
    elif cv_design == "leave_source":
        groups = df["study_id"].to_numpy()
        splitter = GroupKFold(n_splits=df["study_id"].nunique())
        yield from splitter.split(df, y, groups)
    elif cv_design == "leave_layer":
        counts = df["evidence_layer"].value_counts()
        keep_layers = counts[counts >= 8].index
        sub_idx = np.where(df["evidence_layer"].isin(keep_layers))[0]
        sub = df.iloc[sub_idx].reset_index(drop=True)
        groups = sub["evidence_layer"].to_numpy()
        y_sub = sub["functional_final"].to_numpy()
        splitter = GroupKFold(n_splits=min(10, sub["evidence_layer"].nunique()))
        for train_rel, test_rel in splitter.split(sub, y_sub, groups):
            yield sub_idx[train_rel], sub_idx[test_rel]
    else:
        raise ValueError(cv_design)


def fit_and_score(df: pd.DataFrame, feature_set: str, cv_design: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    labels = ROLE_ORDER
    for fold, (train_idx, test_idx) in enumerate(split_iterator(df, cv_design), start=1):
        train = df.iloc[train_idx].copy()
        test = df.iloc[test_idx].copy()
        if train["functional_final"].nunique() < len(labels):
            continue
        pipe = make_pipeline(feature_set)
        pipe.fit(train, train["functional_final"])
        proba_raw = pipe.predict_proba(test)
        model_classes = list(pipe.named_steps["clf"].classes_)
        proba = np.zeros((len(test), len(labels)), dtype=float)
        for j, cls in enumerate(model_classes):
            proba[:, labels.index(cls)] = proba_raw[:, j]
        proba = np.clip(proba, 1e-9, 1.0)
        proba = proba / proba.sum(axis=1, keepdims=True)
        pred = np.array(labels)[np.argmax(proba, axis=1)]
        ll = log_loss(test["functional_final"], proba, labels=labels)
        acc = accuracy_score(test["functional_final"], pred)
        macro_f1 = f1_score(test["functional_final"], pred, labels=labels, average="macro", zero_division=0)
        rows.append(
            {
                "cv_design": cv_design,
                "feature_set": feature_set,
                "fold": fold,
                "n_train": len(train),
                "n_test": len(test),
                "test_sources": ";".join(sorted(test["study_id"].astype(str).unique())),
                "test_layers": ";".join(sorted(test["evidence_layer"].astype(str).unique()))[:500],
                "log_loss": ll,
                "bits_per_row": ll / np.log(2),
                "accuracy": acc,
                "macro_f1": macro_f1,
            }
        )
        for i, (_, r) in enumerate(test.iterrows()):
            pred_rows.append(
                {
                    "cv_design": cv_design,
                    "feature_set": feature_set,
                    "fold": fold,
                    "row_index": int(r["row_index"]),
                    "study_id": r["study_id"],
                    "evidence_layer": r["evidence_layer"],
                    "donor_final": r["donor_final"],
                    "true_role": r["functional_final"],
                    "pred_role": pred[i],
                    **{f"p_{lab}": proba[i, j] for j, lab in enumerate(labels)},
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cv, feat), sub in scores.groupby(["cv_design", "feature_set"]):
        weights = sub["n_test"].to_numpy(dtype=float)
        rows.append(
            {
                "cv_design": cv,
                "feature_set": feat,
                "folds": len(sub),
                "n_test_total": int(weights.sum()),
                "mean_bits_per_row": float(np.average(sub["bits_per_row"], weights=weights)),
                "mean_log_loss": float(np.average(sub["log_loss"], weights=weights)),
                "mean_accuracy": float(np.average(sub["accuracy"], weights=weights)),
                "mean_macro_f1": float(np.average(sub["macro_f1"], weights=weights)),
            }
        )
    summary = pd.DataFrame(rows)
    comparisons = []
    pairs = [
        ("source_layer", "source_layer_donor"),
        ("text", "text_donor"),
        ("source_layer_text", "source_layer_text_donor"),
        ("donor", "source_layer_donor"),
    ]
    for cv in summary["cv_design"].unique():
        ss = summary[summary["cv_design"] == cv].set_index("feature_set")
        for base, plus in pairs:
            if base in ss.index and plus in ss.index:
                comparisons.append(
                    {
                        "cv_design": cv,
                        "baseline_feature_set": base,
                        "augmented_feature_set": plus,
                        "bits_saved_by_augmented": float(ss.loc[base, "mean_bits_per_row"] - ss.loc[plus, "mean_bits_per_row"]),
                        "accuracy_gain": float(ss.loc[plus, "mean_accuracy"] - ss.loc[base, "mean_accuracy"]),
                        "macro_f1_gain": float(ss.loc[plus, "mean_macro_f1"] - ss.loc[base, "mean_macro_f1"]),
                    }
                )
    return summary, pd.DataFrame(comparisons)


def make_figure(summary: pd.DataFrame, comparisons: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.35), gridspec_kw={"wspace": 0.55})

    cv = "stratified_row"
    ss = summary[summary["cv_design"] == cv].set_index("feature_set")
    pairs = [
        ("source_layer", "source_layer_donor", "source+layer"),
        ("text", "text_donor", "text"),
        ("source_layer_text", "source_layer_text_donor", "source+layer+text"),
    ]
    ax = axes[0]
    x = np.arange(len(pairs))
    width = 0.34
    base_vals = [ss.loc[a, "mean_bits_per_row"] for a, _, _ in pairs]
    donor_vals = [ss.loc[b, "mean_bits_per_row"] for _, b, _ in pairs]
    ax.bar(x - width / 2, base_vals, width=width, color=COL["light"], label="without donor")
    ax.bar(x + width / 2, donor_vals, width=width, color=COL["blue"], label="+ donor")
    ax.set_xticks(x, [lab for _, _, lab in pairs], rotation=25, ha="right")
    ax.set_ylabel("cross-entropy (bits/row)")
    ax.legend(frameon=False, fontsize=5.8, loc="upper right")
    ax.text(-0.18, 1.08, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    ax = axes[1]
    comp = comparisons[comparisons["baseline_feature_set"].isin(["source_layer", "text", "source_layer_text"])].copy()
    cv_order = ["stratified_row", "leave_source", "leave_layer"]
    base_order = ["source_layer", "text", "source_layer_text"]
    marker = {"source_layer": "o", "text": "s", "source_layer_text": "^"}
    color = {"source_layer": COL["blue"], "text": COL["grey"], "source_layer_text": COL["purple"]}
    y_pos = {cv: i for i, cv in enumerate(cv_order)}
    for base in base_order:
        sub = comp[comp["baseline_feature_set"] == base]
        ax.scatter(
            sub["bits_saved_by_augmented"],
            [y_pos[x] for x in sub["cv_design"]],
            marker=marker[base],
            s=28,
            color=color[base],
            label={"source_layer": "source+layer", "text": "text", "source_layer_text": "source+layer+text"}[base],
        )
    ax.axvline(0, color=COL["grey"], lw=0.7)
    ax.set_yticks(np.arange(len(cv_order)), ["row CV", "leave-source", "leave-layer"])
    ax.invert_yaxis()
    ax.set_xlabel("bits saved by adding donor")
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    ax.text(-0.18, 1.08, "b", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    ax = axes[2]
    ss = summary.set_index(["cv_design", "feature_set"])
    y = np.arange(len(cv_order))
    base = [ss.loc[(cv, "source_layer_text"), "mean_macro_f1"] for cv in cv_order]
    donor = [ss.loc[(cv, "source_layer_text_donor"), "mean_macro_f1"] for cv in cv_order]
    for yi, b, d in zip(y, base, donor):
        ax.plot([b, d], [yi, yi], color=COL["grey"], lw=1.0)
    ax.scatter(base, y, color=COL["light"], edgecolor=COL["grey"], s=28, label="source+layer+text")
    ax.scatter(donor, y, color=COL["blue"], s=28, label="+ donor")
    ax.set_yticks(y, ["row CV", "leave-source", "leave-layer"])
    ax.invert_yaxis()
    ax.set_xlabel("macro-F1")
    ax.set_xlim(0, 1.02)
    ax.legend(frameon=False, fontsize=5.8, loc="lower right")
    ax.text(-0.18, 1.08, "c", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

    fig.savefig(OUT / "Fig_predictive_feature_ablation.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Fig_predictive_feature_ablation.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Fig_predictive_feature_ablation.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_full_with_audit_overlay()
    feature_sets = ["source_layer", "source_layer_donor", "donor", "text", "text_donor", "source_layer_text", "source_layer_text_donor"]
    cv_designs = ["stratified_row", "leave_source", "leave_layer"]
    score_parts = []
    pred_parts = []
    for cv in cv_designs:
        for feat in feature_sets:
            scores, preds = fit_and_score(df, feat, cv)
            score_parts.append(scores)
            pred_parts.append(preds)
    scores = pd.concat(score_parts, ignore_index=True)
    preds = pd.concat(pred_parts, ignore_index=True)
    summary, comparisons = summarize(scores)

    with pd.ExcelWriter(OUT / "routeB_predictive_feature_ablation.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="model_summary", index=False)
        comparisons.to_excel(writer, sheet_name="feature_comparisons", index=False)
        scores.to_excel(writer, sheet_name="fold_scores", index=False)
        preds.to_excel(writer, sheet_name="row_predictions", index=False)
    summary.to_csv(OUT / "predictive_feature_summary.csv", index=False)
    comparisons.to_csv(OUT / "predictive_feature_comparisons.csv", index=False)
    scores.to_csv(OUT / "predictive_feature_fold_scores.csv", index=False)
    preds.to_csv(OUT / "predictive_feature_row_predictions.csv", index=False)
    make_figure(summary, comparisons)
    (OUT / "README.md").write_text(
        "# Route-B predictive feature-ablation check\n\n"
        "This analysis asks whether donor class adds held-out predictive information for functional role beyond "
        "source, evidence layer and text/module features. It is a supervised diagnostic of label structure, not "
        "an independent biological validation. Models are multinomial logistic regressions with one-hot categorical "
        "features and TF-IDF text features. Designs include stratified row CV, leave-one-source and leave-one-layer CV.\n",
        encoding="utf-8",
    )
    print(OUT)
    print(summary.to_string(index=False))
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
