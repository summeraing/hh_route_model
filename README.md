# Route-resolved evidence-graph framework

This repository contains public code and processed data for a route-resolved evidence-graph framework. The framework tests whether donor origin in a complex evolutionary system is randomly mixed across organizational roles or maps to stable roles after source balancing, dependency collapse and stress testing.

The release supports two analysis cases:

1. **Eukaryogenesis benchmark.** Public route-resolved evidence units, source-equal route support, Fe-S phylogenetic checks and AF3 structural-interactivity support screens.
2. **SOX-HGT transfer validation.** A prokaryotic sulfur-oxidation metabolic-pathway case testing the same donor-to-role mapping workflow on an independent HGT system.

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Repository layout

- `code/core/`: core eukaryogenesis reproduction scripts retained from prior releases.
- `code/route_graph/`: generic route-graph and SOX transfer-case scripts.
- `code/af3_postprocess/`: AF3 summary and gate-evaluation scripts.
- `data/source_data/`: consolidated Source Data workbook for the current dual-case manuscript.
- `data/sox_transfer/`: SOX evidence units, route scores, sensitivity outputs and annotation-stress outputs.
- `data/sox_af3/`: SOX-AF3 processed summaries, gate calls and target matrices.
- `af3_jobs/`: compact AF3 JSON and SLURM submission packets used to run the SOX structural screen on a compute node.
- `data/eukaryogenesis_prior_release/`: selected code/data carried forward from earlier public releases.
- `figures/final_png/`: final figure PNGs for reader orientation.
- `metadata/`: release notes, citation file, traceability reports and SHA256 manifest.

## Minimal quick checks

Install the minimal Python dependencies:

```bash
python -m pip install -r requirements.txt
```

Run both public smoke tests:

```bash
python code/run_smoke_tests.py
```

The smoke test now covers both cases: eukaryogenesis core metrics and
continuous donor-role diagnostics, followed by the SOX-HGT route graph and
SOX-AF3 structural gate check.

Run the SOX route graph:

```bash
python code/route_graph/run_route_graph_case.py \
  --input data/sox_transfer/sox_evidence_units_expanded_v2.csv \
  --case-id sox_hgt_transfer \
  --donors lineage_core,hgt_pathway,mobile_context \
  --roles conserved_metabolic_backbone,variable_energy_module,mobility_boundary \
  --prespec lineage_core=conserved_metabolic_backbone,hgt_pathway=variable_energy_module,mobile_context=mobility_boundary \
  --out out/sox_route_graph \
  --iterations 2000 --enforce-gates
```

Evaluate SOX-AF3 gates from processed summaries:

```bash
python code/af3_postprocess/evaluate_sox_af3_gates.py \
  --summary-csv data/sox_af3/sox_af3_combined_priority1_tier2.csv \
  --output-md out/sox_af3_gate_check.md \
  --output-csv out/sox_af3_gate_check.csv
```

Core eukaryogenesis checks from prior releases are retained where available:

```bash
python code/core/reproduce_core_metrics.py
python code/core/reproduce_continuous_route_diagnostic.py
```

## Expected case-level outputs

Eukaryogenesis benchmark:

- adjudicated prespecified route rank: 1
- adjudicated source-equal support: 0.769450
- frozen q97.5 null boundary: 0.565984
- margin over q97.5: 0.203466
- minimum dependency-collapse margin: 0.206667
- joint bootstrap probability that all three prespecified roles are top: 0.999800

SOX transfer case:

- route-eligible rows: 65
- prespecified route rank: 1
- source-equal support: 1.000
- best alternative: 0.691
- margin over best alternative: 0.309
- annotation-stress route-eligible rows: 64
- annotation-stress prespecified route rank: 1
- annotation-stress margin over best alternative: 0.315

SOX-AF3 layer:

- the public release contains processed AF3 summary tables and gate calls, not a claim of experimentally measured binding;
- AF3 predictions are used as structural-screen support for interface plausibility and decoy separation;
- raw AF3 model directories may be too large or environment-specific for repository storage, so compact summaries and job inputs are included.

## Scope

This repository intentionally excludes manuscript drafts, cover letters, Supplementary Information documents, peer-review correspondence, PPT files and planning notes. It is a public code and processed-data companion for the route-resolved framework.

## Release

Current release tag: `v2026.07.03-dual-case-af3-framework`  
Zenodo concept DOI: https://doi.org/10.5281/zenodo.20453582
