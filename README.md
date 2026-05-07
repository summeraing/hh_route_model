# Core calculation code for route-partitioned integration

This minimal repository contains only the core calculation code and the small derived numerical tables needed to verify the main locked claims for:

**Route-partitioned integration reconciles scaffold continuity with energetic innovation in eukaryogenesis**

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Included

- `code/hh_route_model_formalism.py`: route-support and H-h formalism helper functions.
- `code/reproduce_core_metrics.py`: small validator that reads the derived calculation tables and prints the locked audit-completed metrics.
- `data/core_metrics/`: compact derived tables for route specificity, frozen null boundary, dependency collapse and independent-coding agreement.
- `metadata/`: manifest and citation metadata.

## Not included

This repository intentionally excludes manuscript files, cover letters, Supplementary Information files, final figures, figure-plotting scripts, PPT files, image assets, full source-data workbooks and editorial planning materials.

The full Source Data workbook and journal upload files remain in the submission package. This repository is meant to be the small public code/calculation companion, not a mirror of the submission folder.

## Quick check

```bash
python code/reproduce_core_metrics.py
```

Expected core outputs:

- prespecified route rank: 1
- adjudicated source-equal support: 0.76945
- frozen q97.5 null boundary: 0.565984
- margin over q97.5 boundary: 0.203466
- donor-class agreement: 0.771104; kappa = 0.464655
- functional-class agreement: 0.793831; kappa = 0.466401
