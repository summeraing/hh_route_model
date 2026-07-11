# Release notes: v2026.07.11-reproducible-route-graph

Date: 2026-07-11

This release makes the dual-case route-resolved evidence-graph analysis easier to audit, reproduce and transfer without changing the reported biological claims.

## Reproducibility additions

- Added a pinned Python 3.11 environment (`environment.yml` and `requirements.txt`).
- Added `QUICKSTART.md` with tested commands, expected outputs and scope limits.
- Added an executable nine-row third-system template under `examples/third_system_template/`.
- Added a repository-level MIT license and an explicit validation log.
- Added journal-neutral SOX-AF3 gate labels and a clean processed-summary filename.
- Expanded the public Source Data companion with coding agreement, adjudication and dependency-collapse outputs.

## Validation

`python code/run_smoke_tests.py` passed on 2026-07-11. The eukaryogenesis benchmark reproduced the reported adjudicated route support and the SOX-HGT case reproduced its route ranking, sensitivity outputs and structural-screen gate call. The third-system template also completed with the prespecified mapping ranked first.

## Interpretation limits

- Row-level evidence entries are traceable records, not independent biological replicates.
- The SOX case is a transfer validation, not a complete reconstruction of sulfur-oxidation evolution.
- AlphaFold 3 outputs are prediction-based structural screens and are not biochemical binding measurements.
