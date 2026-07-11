# Validation log

## Public smoke-test contract

- Command: `python code/run_smoke_tests.py`
- Expected endpoint: `SMOKE_TESTS_PASS`
- Covered components: eukaryogenesis core metrics, continuous donor-role diagnostic, SOX-HGT route graph, and processed SOX-AF3 gate evaluation.
- Excluded from the smoke test: raw literature extraction, sequence retrieval, phylogenetic inference and AlphaFold 3 model generation.

## Third-system schema test

- Command: see `QUICKSTART.md`.
- Expected result: the toy prespecified route ranks first and all explicitly lowered toy preflight thresholds pass.
- Interpretation: software/schema test only; the toy rows are not biological evidence.

## Current validation

- Date: 2026-07-11 (Asia/Shanghai).
- Runtime: Python 3.11-compatible Anaconda environment on Windows.
- Public smoke tests: PASS (`SMOKE_TESTS_PASS`).
- Third-system template: PASS; 9 route-eligible rows, 3 sources, 9 dependency groups, prespecified rank 1 and margin 0.667.
