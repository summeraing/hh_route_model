# Objective-Level Reproducibility Audit

Date: 2026-07-03

Objective: establish a reproducible route-resolved evidence-graph framework for testing whether donor origin in complex evolutionary systems is randomly mixed or maps to organizational roles, with eukaryogenesis as the deep benchmark and a prokaryotic HGT metabolic pathway as transfer-validation case.

## Current Evidence

| Requirement | Evidence in this public packet | Status |
|---|---|---|
| Reusable route-resolved evidence-graph framework | `code/route_graph/run_route_graph_case.py`; `code/core/hh_route_model_formalism.py`; `README.md` | Present |
| Test donor-to-role mapping versus alternatives | SOX route scores in `data/sox_transfer/route_scores.csv`; smoke-test route rank check in `code/run_smoke_tests.py` | Present |
| Source-bias correction | Source-equal scoring in `run_route_graph_case.py`; source-rate and leave-one-source outputs in `data/sox_transfer/` | Present |
| Dependency evidence collapse | `dependency_collapsed_scores.csv`; dependency-collapse implementation in `run_route_graph_case.py` | Present |
| Eukaryogenesis deep benchmark | Prior-release eukaryogenesis code/data carried in `data/eukaryogenesis_prior_release/`; current Source Data workbook in `data/source_data/` | Present as processed public package |
| Prokaryotic HGT transfer-validation case | SOX-HGT evidence units, route outputs, annotation cross-check and stress rerun in `data/sox_transfer/` | Present |
| Structural cross-validation | SOX-AF3 summaries, gate calls and AF3 job inputs in `data/sox_af3/` and `af3_jobs/` | Present as prediction-based screen |
| Phylogenetic cross-validation | Fe-S validation scripts and prior-release eukaryogenesis resources in `code/core/` and `data/eukaryogenesis_prior_release/` | Present as retained benchmark layer |
| One-command public smoke test | `python code/run_smoke_tests.py` | Present |
| Checksum traceability | `metadata/FILE_MANIFEST_SHA256.csv` | Present |

## Verified Local Commands

From the packet root:

```bash
python code/run_smoke_tests.py
```

Expected terminal sentinel:

```text
SMOKE_TESTS_PASS
```

The smoke test verifies:

- Eukaryogenesis benchmark: adjudicated prespecified route rank 1, source-equal support 0.769450, q97.5 null boundary 0.565984, positive dependency-collapse margin, and joint continuous-route bootstrap probability 0.999800.
- SOX transfer case: 65 route-eligible rows, prespecified route rank 1, source-equal support 1.000, margin 0.308712.
- SOX-AF3 gate layer: processed summaries support the structural-gate separation decision string.

## Not Yet Proved by This Local Packet

- A fresh GitHub release has not been verified from this audit file alone.
- A fresh Zenodo version DOI has not been verified from this audit file alone.
- Raw AF3 model directories are not included; compact processed summaries and job inputs are included instead.
- External independent reproduction on a clean machine has not been performed here.

## Interpretation

The local public-release packet is ready for repository upload and archival. The full objective should be treated as not externally complete until the packet is pushed to GitHub, archived by Zenodo and the resulting DOI/version information is checked.
