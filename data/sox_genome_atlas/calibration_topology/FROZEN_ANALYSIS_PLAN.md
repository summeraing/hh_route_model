# Frozen three-day analysis plan

Spec ID: `imeta_3day_upgrade_20260712_v1`

The analysis specification was frozen before inspecting the new results. All
outcomes will be retained, including null results and results that weaken the
manuscript's current interpretation.

## Question 1: Is the route graph calibrated under known truth?

Synthetic evidence graphs contain five sources, three donors, three roles and
300 dependency groups. Full-route, partial-route and diffuse-null truth are
crossed with source dominance, repeated-row dependence and label noise. A
dominant-adversarial regime tests whether one row-rich source can overturn the
majority of independent sources. Raw pooling, source equalization and
source-equal dependency collapse are compared without changing their scoring
rules. Half of the null replicates define matched 97.5% margin thresholds; the
other half estimate false-positive rates.

## Question 2: Does the SOX completeness association survive both trees?

Terminal transfer-recipient support from the fixed LG+F+R6 and ModelFinder
GeneRax analyses is merged by genome and family. Conservative support is the
minimum across tree designs; mean support and a binary stable-in-both metric
are secondary. Association with module completeness is tested by shuffling
genome labels within GTDB order. Confidence intervals use within-order
bootstrap resampling. Analyses are repeated after omitting each SOX family and
Holm corrected.

## Question 3: Is species-tree congruence rejected by sequence likelihood?

For SoxA, SoxB, SoxC, SoxX, SoxY and SoxZ, IQ-TREE compares the unconstrained
maximum-likelihood topology with a maximum-likelihood tree inferred under the
corresponding pruned GTDB topology constraint. Tests are run under both the
family-specific ModelFinder model and fixed LG+F+R6. SH, weighted-SH and AU
tests use 10,000 RELL replicates. Family-wise P values are Holm corrected within
tree design and test.

## Reporting boundary

Simulation is a method-calibration analysis. Cross-model consensus is a
conservative sensitivity analysis. Topology tests establish sequence-level
discordance but do not alone prove HGT direction, timing or mechanism.
