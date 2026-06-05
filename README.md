# Core calculation code for route-partitioned integration

This minimal repository contains only the core calculation code and the small derived numerical tables needed to verify the main locked claims for:

**Route-partitioned integration reconciles scaffold continuity with energetic innovation in eukaryogenesis**

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Included

- `code/hh_route_model_formalism.py`: route-support and H-h formalism helper functions.
- `code/reproduce_core_metrics.py`: small validator that reads the derived calculation tables and prints the locked audit-completed metrics.
- `code/reproduce_continuous_route_diagnostic.py`: validator for the continuous donor-role probability diagnostic.
- `code/run_expanded_fes_validation.py`: UniProt/MAFFT/IQ-TREE workflow used for the expanded Fe-S marker-tree validation; this script requires network access and external phylogenetic tools.
- `data/core_metrics/`: compact derived tables for route specificity, frozen null boundary, dependency collapse, independent-coding agreement, continuous donor-role diagnostics and expanded Fe-S validation summaries.
- `metadata/`: manifest and citation metadata.

## Not included

This repository intentionally excludes manuscript files, cover letters, Supplementary Information files, final figures, figure-plotting scripts, PPT files, image assets, full source-data workbooks and editorial planning materials.

The full Source Data workbook and journal upload files remain in the submission package. This repository is meant to be the small public code/calculation companion, not a mirror of the submission folder.

## Quick check

```bash
python code/reproduce_core_metrics.py
python code/reproduce_continuous_route_diagnostic.py
```

The expanded Fe-S workflow is not a lightweight quick check because it retrieves homologues and calls external alignment/tree-building tools. The script is included for transparency and rerun planning:

```bash
python code/run_expanded_fes_validation.py --help
```

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
