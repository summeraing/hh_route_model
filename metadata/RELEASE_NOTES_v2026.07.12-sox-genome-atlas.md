# v2026.07.12-sox-genome-atlas

This release upgrades the SOX application from a 65-row literature-coded pilot to an independent comparative-genomics benchmark.

## Added

- Core Python workflow for scanning 7,453 RefSeq reference complete/chromosome prokaryotic genomes and mapping strict candidates to GTDB R232.
- Processed QC tables for 427 strict SOX-positive genomes and 724 SOX clusters.
- Six ModelFinder-selected and fixed-model SOX family trees with GeneRax reconciliation summaries.
- Matched same-genome controls for mobile-element context and local composition.
- Frozen route-map, continuous donor-role, dependency-collapse and cross-layer sensitivity outputs.
- Leave-one-family sensitivity under fixed and ModelFinder-selected trees.
- A compact release validator integrated into the public smoke tests.

## Interpretation boundary

The frozen symmetric three-route SOX map ranks fifth of six and is not supported. The retained result is narrower: HGT-like evidence is concentrated in the SOX module, while matched controls do not support a general mobile-boundary or recent compositionally anomalous island signature. AF3 outputs remain predicted-interface support only.

## Full archive

The GitHub release asset contains the compact complete reproducibility archive, including processed tables and cluster job specifications. Raw RefSeq genomes, GTDB database files and AlphaFold 3 model weights are not redistributed.
