# v2026.07.12-calibration-consensus-topology

This release adds three prespecified computational checks to the SOX genome-atlas benchmark.

## Added

- Known-truth simulations comparing raw, source-equal and source-equal-plus-dependency-collapse route scoring across source imbalance, label noise, dependence and partial-route conditions.
- Fixed-tree/ModelFinder cross-tree recipient-support consensus and module-completeness permutation tests.
- SH, weighted-SH and approximately unbiased tests of optimized GTDB-constrained topologies for six SOX families under two tree models, with Holm correction across families.
- Frozen analysis specifications, scheduler files, replicate-level outputs, IQ-TREE audit files and concise verified summaries.

## Main calibration results

- Diffuse-null false-positive rates remain 2.3-2.8% across the three route scores.
- Under a representative 85% dominant-source scenario, source-equal and collapsed scores recover the complete route in 79.75% and 81.50% of simulations, whereas raw pooling recovers 0%.
- Cross-tree support profiles agree above permutation expectation (mean rho = 0.202, P = 0.0104); the stable-in-both subset remains a stated boundary result.
- All twelve family-by-model GTDB constraints are rejected after Holm-adjusted AU testing.

## Interpretation boundary

Simulation operating characteristics apply to the specified generators and are not evidence that an empirical dataset follows a simulation generator. Topology tests reject strict GTDB congruence but do not identify horizontal-transfer direction or timing. Manuscript drafts, journal-formatted figures and editorial files remain excluded.
