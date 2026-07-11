# Frozen SOX genome-atlas analysis contract

## Scope

The SOX genome atlas is an independent comparative-genomics benchmark based on 7,453 RefSeq reference complete/chromosome prokaryotic genomes and GTDB R232. It tests whether donor evidence maps to organizational roles without using role identity to define donor labels.

The primary biological dependency units are genomes, position-defined operon clusters and homolog families. Row-level source-by-layer entries remain traceable but are not treated as independent biological replicates.

## Independent evidence axes

Donor evidence is derived from gene-tree/species-tree discordance, calibrated transfer-recipient support, within-order distribution, composition and matched mobile context. Organizational role is derived separately from SOX feature identity, pathway position and operon architecture.

The prespecified symmetric map was:

- vertical-like evidence to sulfur backbone;
- HGT-like evidence to SOX module;
- mobile context to positional boundary.

All six donor-role bijections, the complete continuous donor-role matrix and dependency-collapsed analyses are released whether or not the prespecified map passes.

## Frozen quality gates

- Scan all 7,453 frozen RefSeq assemblies and report parse failures and GTDB mapping failures.
- Exclude pseudogene, partial CDS and known product-name collision records from strict candidate calls while retaining an audit ledger.
- Require at least 20 single-copy genomes for a SOX family tree.
- Preserve mixed on-route and off-route observations within database-derived evidence layers.
- Require at least two evidence modalities, at least 200 dependency groups and no single source contributing more than 80% of route-eligible rows.
- Require non-zero variance in the source and source-by-layer permutation nulls.
- Report raw, source-equal, dependency-collapsed, leave-one-layer and continuous-map outputs.
- Treat post-boundary analyses as exploratory sensitivity tests rather than replacements for the frozen rank test.
- Treat AlphaFold 3 outputs only as predicted-interface support for extant protein combinations.

## Frozen outcome

- 7,453/7,453 GFF files were parsed successfully.
- 427 strict SOX-positive genomes and 724 position-defined clusters were retained; 398 genomes mapped to GTDB R232.
- Six core families passed the tree gate and showed extensive gene-tree/species-tree discordance under ModelFinder-selected and fixed LG+F+R6 analyses.
- The source-dominance gate failed because RefSeq-derived rows contributed 85.4% of route-eligible entries.
- The symmetric route ranked fifth of six before and after dependency collapse and was rejected.
- Continuous mapping retained a localized HGT-like-to-SOX-module component.
- Same-genome matched controls did not support general mobile-boundary enrichment or recent compositionally anomalous SOX islands.
- The association between calibrated recipient support and module completeness was positive under both tree designs; corrected leave-one-family significance was obtained only under the fixed-tree sensitivity design.

This boundary outcome is the reported result. No thresholds, donor labels or role labels were changed to recover the prespecified symmetric map.
