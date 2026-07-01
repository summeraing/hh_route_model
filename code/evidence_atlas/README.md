# Evidence-atlas analysis scripts

This folder contains the key calculation scripts for the route-resolved evolutionary evidence-atlas analysis.

The scripts are provided for transparent inspection and rerun planning. They operate on the processed evidence tables used in the manuscript package and are not a one-command reconstruction from every original public source file.

Included calculation layers:

- probabilistic donor-role mixture and rule-based audit-label cross-check
- model-family and posterior rank-stability diagnostics
- predictive cross-validation and source/layer permutation nulls
- hierarchical source/layer/row resampling and Shapley contribution analysis
- evidence-layer accumulation and future-source adversarial boundary
- source-aware information-theory dependence tests
- supervised predictive feature-ablation checks
- unsupervised correspondence, graph-modularity and clustering route-recovery diagnostics
- block-level MDL and held-out model-evidence diagnostics
- matrix-level and row-level adversarial counterfactual diagnostics

The companion processed data workbooks are in `data/evidence_atlas/`.
