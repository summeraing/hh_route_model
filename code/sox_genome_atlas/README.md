# SOX genome-atlas workflow

This directory contains the core Python implementation used to build and test the genome-scale SOX benchmark. The numbered scripts follow the analysis order from RefSeq universe preparation and candidate detection through gene-tree/species-tree discordance, evidence construction, route testing, matched controls, GeneRax summaries and cross-layer sensitivity.

The full phylogenetic workflow additionally requires MAFFT, IQ-TREE 2, GeneRax and a local GTDB R232 species-tree release. Cluster computations must be submitted through the scheduler; the scripts do not require or recommend computation on a login node.

Run the compact release check from the repository root:

```bash
python code/sox_genome_atlas/validate_released_results.py \
  --data-root data/sox_genome_atlas
```

This check validates released results; it does not reconstruct 7,453 genomes, infer trees or rerun GeneRax.
