# Route-resolved evolutionary evidence-atlas code and data

This repository contains the calculation code, processed numerical tables and structural-interactivity support package for a route-resolved evolutionary evidence-atlas analysis of eukaryogenesis.

The analysis separates contributor identity from organizational role and tests whether public evolutionary evidence supports localized energetic incorporation within a continuing host scaffold rather than diffuse ancestry mixture or wholesale replacement.

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Included

- `code/hh_route_model_formalism.py`: route-support and H-h formalism helper functions.
- `code/reproduce_core_metrics.py`: small validator that reads the derived calculation tables and prints the locked audit-completed metrics.
- `code/reproduce_continuous_route_diagnostic.py`: validator for the continuous donor-role probability diagnostic.
- `code/run_expanded_fes_validation.py`: UniProt/MAFFT/IQ-TREE workflow used for the expanded Fe-S marker-tree validation; this script requires network access and external phylogenetic tools.
- `code/evidence_atlas/`: route-resolved evidence-atlas calculation scripts for probabilistic mixture, model-family comparison, predictive cross-validation, permutation nulls, information theory, feature ablation, unsupervised route recovery, block-level MDL, resampling, Shapley and adversarial diagnostics.
- `data/core_metrics/`: compact derived tables for route specificity, frozen null boundary, dependency collapse, independent-coding agreement, continuous donor-role diagnostics and expanded Fe-S validation summaries.
- `data/evidence_atlas/`: processed workbooks for the full evidence-atlas computational evidence package, including a consolidated Source Data workbook and module-level result workbooks.
- `data/structural_interactivity/`: AF3 structural-interactivity support package containing representative model CIF files, AF3 confidence outputs, ranking-score files, integrated interface metrics and interpretation-boundary metadata.
- `metadata/`: manifest and citation metadata.

## Not included

This repository intentionally excludes manuscript files, cover letters, Supplementary Information files, final figures, figure-design working files, PPT files, image assets and editorial planning materials.

The journal upload files remain in the submission package. This repository is meant to be the public code, derived-table and structural-output companion, not a mirror of the submission folder.

## Quick check

```bash
python code/reproduce_core_metrics.py
python code/reproduce_continuous_route_diagnostic.py
```

The expanded evidence-atlas modules are provided as auditable calculation scripts and processed workbooks:

```text
code/evidence_atlas/
data/evidence_atlas/
```

The expanded Fe-S workflow is not a lightweight quick check because it retrieves homologues and calls external alignment/tree-building tools. The script is included for transparency and rerun planning:

```bash
python code/run_expanded_fes_validation.py --help
```

The AF3 support package is provided as a compact archive for inspection rather than as a rerunnable AF3 workflow:

```text
data/structural_interactivity/AF3_structural_interactivity_public_release_20260608.zip
```

This archive contains 19 representative model CIF files, 19 AF3 confidence JSON outputs, 19 ranking-score CSV files, integrated interface-metric tables and a README. It does not include proprietary AF3 model parameters and is not interpreted as direct ancestral interface reconstruction.

Expected core outputs:

- prespecified route rank: 1
- adjudicated source-equal support: 0.76945
- frozen q97.5 null boundary: 0.565984
- margin over q97.5 boundary: 0.203466
- donor-class agreement: 0.771104; kappa = 0.464655
- functional-class agreement: 0.793831; kappa = 0.466401

Expected continuous diagnostic outputs:

- host -> scaffold: source-equal probability = 0.540984; margin = 0.088525
- symbiont -> energetic incorporation: source-equal probability = 0.967320; margin = 0.934641
- other -> transition: source-equal probability = 0.761905; margin = 0.619048
- joint bootstrap probability that all three prespecified roles remain top-ranked: 0.999800

Additional evidence-atlas diagnostics include:

- source-layer MDL best-ranked model: shared-theta route model
- leave-source held-out loss for shared-theta route model: 1.220 bits
- leave-source-layer held-out loss for shared-theta route model: 1.163 bits
- source-layer block bootstrap rank-1 probability: 0.9819
- unsupervised correspondence residual rank: 1
- unsupervised graph-modularity rank: 1
