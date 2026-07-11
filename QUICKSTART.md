# Quick start

## 1. Create the environment

Conda:

```bash
conda env create -f environment.yml
conda activate route-evidence-graph
```

or pip:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

Tested with Python 3.11. Randomized analyses use explicit seeds in the command-line scripts.

## 2. Run the public smoke tests

```bash
python code/run_smoke_tests.py
```

Expected terminal endpoint:

```text
SMOKE_TESTS_PASS
```

The smoke tests reproduce the eukaryogenesis core metrics, run the legacy SOX-HGT route graph with 200 permutation draws, evaluate the processed SOX-AF3 gate table and validate key invariants in the released SOX genome atlas. Typical desktop runtime is below five minutes; it does not rerun sequence searches, phylogenetic inference, GeneRax or AlphaFold 3.

## 3. Validate the SOX genome-atlas release

```bash
python code/sox_genome_atlas/validate_released_results.py \
  --data-root data/sox_genome_atlas
```

Expected endpoint:

```text
SOX_GENOME_ATLAS_RELEASE_VALIDATION_PASS
```

The complete reconstruction workflow is documented in `code/sox_genome_atlas/README.md` and `data/sox_genome_atlas/ANALYSIS_CONTRACT.md`. It requires RefSeq/GTDB source files plus MAFFT, IQ-TREE 2 and GeneRax, and should be submitted through a compute scheduler.

## 4. Run the legacy SOX-HGT pilot

```bash
python code/route_graph/run_route_graph_case.py \
  --input data/sox_transfer/sox_evidence_units_expanded_v2.csv \
  --case-id sox_hgt_transfer \
  --donors lineage_core,hgt_pathway,mobile_context \
  --roles conserved_metabolic_backbone,variable_energy_module,mobility_boundary \
  --prespec lineage_core=conserved_metabolic_backbone,hgt_pathway=variable_energy_module,mobile_context=mobility_boundary \
  --out out/sox_route_graph \
  --iterations 2000 --seed 20260701 --enforce-gates
```

Primary outputs are `summary.csv`, `route_scores.csv`, `donor_role_matrix.csv`, `leave_one_source.csv`, `dependency_collapsed_scores.csv`, `permutation_nulls.csv` and an Excel workbook containing the same tables.

## 5. Run the third-system template

The toy template checks file structure and command syntax; it is not manuscript evidence.

```bash
python code/route_graph/run_route_graph_case.py \
  --input examples/third_system_template/evidence_units_template.csv \
  --case-id toy_transfer \
  --donors inherited,foreign,mobile \
  --roles backbone,module,boundary \
  --prespec inherited=backbone,foreign=module,mobile=boundary \
  --out out/toy_transfer \
  --iterations 100 --seed 11 \
  --min-units 9 --min-sources 3 --min-dependency-groups 9 --max-source-fraction 0.5 --enforce-gates
```

## Reproducibility boundary

The repository reproduces route scoring, stress tests, genome-atlas result invariants and processed AF3 gate summaries from released processed inputs. Raw public genomes, GTDB databases and AF3 model weights are not redistributed. The complete compact atlas archive is provided as a release asset for inspection and staged reconstruction.
