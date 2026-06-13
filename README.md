# Core calculation code for route-partitioned integration

This repository contains the core calculation code, compact derived numerical tables and structural-interactivity support package needed to inspect the main locked claims for:

**Route-partitioned integration reconciles scaffold continuity with energetic innovation in eukaryogenesis**

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Included

- `code/hh_route_model_formalism.py`: route-support and H-h formalism helper functions.
- `code/reproduce_core_metrics.py`: small validator that reads the derived calculation tables and prints the locked audit-completed metrics.
- `code/reproduce_continuous_route_diagnostic.py`: validator for the continuous donor-role probability diagnostic.
- `code/run_expanded_fes_validation.py`: UniProt/MAFFT/IQ-TREE workflow used for the expanded Fe-S marker-tree validation; this script requires network access and external phylogenetic tools.
- `data/core_metrics/`: compact derived tables for route specificity, frozen null boundary, dependency collapse, independent-coding agreement, continuous donor-role diagnostics and expanded Fe-S validation summaries.
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

The expanded Fe-S workflow is not a lightweight quick check because it retrieves homologues and calls external alignment/tree-building tools. The script is included for transparency and rerun planning:

```bash
python code/run_expanded_fes_validation.py --help
```

The AF3 support package is provided in two layers.

First, the compact 19-pair structural-interactivity archive is retained for inspection:

```text
data/structural_interactivity/AF3_structural_interactivity_public_release_20260608.zip
```

Second, the AF3-enhanced W2/W2C structural-interactome update is provided as processed tables and representative structures:

```bash
python code/reproduce_af3_structural_gate.py
```

This script recomputes the W2 high-confidence edge counts, route-label null result and W2C seed-stability calls from `data/af3/`. The AF3 outputs are interpreted as a calibrated structural-compatibility screen, not as direct ancestral in vivo complex reconstruction.

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
