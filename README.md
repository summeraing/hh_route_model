# Route-resolved evidence-graph framework

This repository contains public code and processed data for a route-resolved evidence-graph framework. The framework tests whether donor origin in a complex evolutionary system is randomly mixed across organizational roles or maps to stable roles after source balancing, dependency collapse and stress testing.

The release supports two principal analysis cases and one retained pilot:

1. **Eukaryogenesis benchmark.** Public route-resolved evidence units, source-equal route support, Fe-S phylogenetic checks and AF3 structural-interactivity support screens.
2. **SOX genome atlas.** A comparative-genomics benchmark built from 7,453 RefSeq reference complete/chromosome assemblies, 427 strict SOX-positive genomes, GTDB R232 species trees, six core SOX family trees, matched genomic controls and dependency-aware route tests.

The earlier 65-row SOX-HGT literature-coded case is retained as a legacy smoke-test example, not as the primary SOX result.

Author: Zhiqiang Xia  
Correspondence: zqiangx@gmail.com

## Repository layout

- `code/core/`: core eukaryogenesis reproduction scripts retained from prior releases.
- `code/route_graph/`: generic route-graph and SOX transfer-case scripts.
- `code/sox_genome_atlas/`: RefSeq/GTDB comparative-genomics, reconciliation, matched-control and sensitivity scripts.
- `code/sox_genome_atlas/calibration_topology/`: known-truth calibration, cross-tree consensus and formal topology-test scripts.
- `code/af3_postprocess/`: AF3 summary and gate-evaluation scripts.
- `data/source_data/`: consolidated Source Data workbook for the current dual-case manuscript.
- `data/sox_transfer/`: SOX evidence units, route scores, sensitivity outputs and annotation-stress outputs.
- `data/sox_genome_atlas/`: processed atlas evidence, family trees, reconciliations, matched controls and model sensitivities.
- `data/sox_genome_atlas/calibration_topology/`: frozen specifications and complete processed outputs for route-score calibration, fixed/ModelFinder consensus and SH/AU topology tests.
- `data/sox_af3/`: SOX-AF3 processed summaries, gate calls and target matrices.
- `af3_jobs/`: compact AF3 JSON and SLURM submission packets used to run the SOX structural screen on a compute node.
- `data/eukaryogenesis_prior_release/`: selected code/data carried forward from earlier public releases.
- `figures/final_png/`: final figure PNGs for reader orientation.
- `metadata/`: release notes, citation file, traceability reports and SHA256 manifest.
- `slurm/calibration_topology/`: scheduler specifications used for the added server analyses; no production computation was run on the login node.

## Minimal quick checks

Install the minimal Python dependencies:

```bash
python -m pip install -r requirements.txt
```

For a pinned Conda environment, use `conda env create -f environment.yml`. A complete five-minute walkthrough, expected outputs and the third-system schema example are provided in `QUICKSTART.md`.

Run both public smoke tests:

```bash
python code/run_smoke_tests.py
```

The smoke test covers eukaryogenesis core metrics, the legacy SOX-HGT route graph,
the SOX-AF3 structural gate and invariant checks for the released SOX genome atlas.

Validate the genome-atlas release directly:

```bash
python code/sox_genome_atlas/validate_released_results.py \
  --data-root data/sox_genome_atlas
```

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
  --summary-csv data/sox_af3/sox_af3_combined_primary_boundary.csv \
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

SOX genome atlas:

- 7,453 RefSeq genomes scanned with zero GFF parse failures;
- 427 strict SOX-positive genomes, including 398 mapped to GTDB R232;
- six core families show 36.5-47.4 calibrated transfer events per 100 tips across ModelFinder-selected and fixed-tree analyses;
- the frozen symmetric three-route map ranks 5 of 6 and is rejected;
- the retained result is localized HGT-like evidence in the SOX module, with matched mobile and composition controls rejecting a general recent-island explanation.
- known-truth simulations keep diffuse-null false-positive rates at 2.3-2.8% and show that source-equal scoring recovers localized signals that raw pooling loses under severe source imbalance;
- fixed-tree and ModelFinder recipient-support profiles agree above permutation expectation (mean cross-tree rho = 0.202, P = 0.0104), while the stable-in-both subset is reported as a boundary result;
- strict GTDB-constrained topologies are rejected for all six tested SOX families under both tree models after Holm-adjusted AU tests; these tests reject strict congruence but do not infer transfer direction or timing.

Legacy 65-row SOX transfer case:

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

## License

The analysis code is released under the MIT License. Public-source data and derived tables remain subject to the provenance and reuse conditions of their original sources.

## Release

Current release tag: `v2026.07.12-calibration-consensus-topology`
Zenodo concept DOI: https://doi.org/10.5281/zenodo.20453582
