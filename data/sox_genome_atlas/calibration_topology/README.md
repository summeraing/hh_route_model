# Route-graph calibration and SOX tree-robustness upgrade

This directory contains the frozen analysis code and processed outputs for three additions to the manuscript:

1. known-truth calibration of raw, source-equal and source-equal-plus-dependency-collapse route scoring;
2. fixed-tree/ModelFinder recipient-support consensus and module-completeness tests;
3. SH, weighted-SH and approximately unbiased tests of GTDB-constrained SOX family trees.

The directory intentionally excludes manuscript files and figure-rendering code. The public analysis package is designed to expose computational inputs, specifications, statistical outputs and audit trails rather than journal-specific presentation files.

## Layout

- `config/frozen_analysis_spec.json`: frozen simulation and test specification.
- `FROZEN_ANALYSIS_PLAN.md`: hypotheses, estimands, boundaries and decision rules.
- `scripts/01_route_simulation.py`: generates known-truth simulation replicates.
- `scripts/02_aggregate_simulation.py`: calibrates matched null gates and recovery metrics.
- `scripts/03_cross_model_consensus.py`: builds fixed/ModelFinder consensus support and permutation tests.
- `scripts/04_prepare_topology_test.py`: prepares optimized GTDB-constrained family-tree candidates.
- `scripts/05_parse_topology_tests.py`: parses IQ-TREE SH, weighted-SH and AU outputs and applies Holm correction.
- `slurm/`: scheduler scripts used for all server computations; no analysis was run on the login node.
- `results/01_route_simulation/`: held-out null gates, full parameter-grid summaries and replicate-level results.
- `results/02_cross_model_consensus/`: family-genome and genome-level consensus tables.
- `results/03_topology_tests/`: candidate trees, IQ-TREE logs and family-level topology-test summaries.
- `results/UPGRADE_RESULTS_SUMMARY.md`: concise verified findings and interpretation boundaries.

## Interpretation boundary

The known-truth simulation measures operating characteristics under specified generators and is not evidence that either empirical dataset follows a simulation generator. The formal topology tests reject strict GTDB congruence but do not identify transfer direction or timing. AlphaFold 3 model weights and journal-specific figure scripts are not included here.
