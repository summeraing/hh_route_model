#!/usr/bin/env bash
#SBATCH --job-name=sox_topotest
#SBATCH --partition=base
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/04_topology_%A_%a.out
#SBATCH --error=logs/04_topology_%A_%a.err

set -euo pipefail
ROOT="${ROOT:?ROOT must point to the three-day upgrade directory}"
ATLAS="${ATLAS:?ATLAS must point to the frozen SOX atlas directory}"
FAMILIES=(soxa soxb soxc soxx soxy soxz)
DESIGNS=(fixed_LG_F_R6 ModelFinder)
FAMILY="${FAMILIES[$((SLURM_ARRAY_TASK_ID % 6))]}"
DESIGN="${DESIGNS[$((SLURM_ARRAY_TASK_ID / 6))]}"

if [[ "${DESIGN}" == "fixed_LG_F_R6" ]]; then
  INPUT_ROOT="${ATLAS}/results/22_generax_inputs_fixed"
  MODEL="LG+F+R6"
else
  INPUT_ROOT="${ATLAS}/results/31_generax_inputs_modelselected"
  MODEL="$(grep -m1 'Best-fit model' "${ATLAS}/results/07_gene_trees/${FAMILY}/${FAMILY}.iqtree" | awk -F': ' '{print $2}')"
fi

IN="${INPUT_ROOT}/${FAMILY}"
OUT="${ROOT}/results/03_topology_tests/${DESIGN}/${FAMILY}"
mkdir -p "${OUT}" "${ROOT}/logs"
ALIGNMENT="${IN}/${FAMILY}.gtdb_shared.trimmed.faa"
UNCONSTRAINED="${IN}/${FAMILY}.gtdb_shared.tree"
SPECIES_TREE="${IN}/${FAMILY}.gtdb_rooted_binary.tree"
MAPPING="${IN}/${FAMILY}.mapping.link"

python "${ROOT}/scripts/04_prepare_topology_test.py" \
  --mapping "${MAPPING}" \
  --species-tree "${SPECIES_TREE}" \
  --out "${OUT}/gtdb_constraint_gene_labels.tree"

iqtree -s "${ALIGNMENT}" -m "${MODEL}" \
  -g "${OUT}/gtdb_constraint_gene_labels.tree" \
  -nt AUTO -ntmax "${SLURM_CPUS_PER_TASK}" \
  -pre "${OUT}/constrained" -redo

cat "${UNCONSTRAINED}" "${OUT}/constrained.treefile" > "${OUT}/candidate_trees.tre"
iqtree -s "${ALIGNMENT}" -m "${MODEL}" \
  -z "${OUT}/candidate_trees.tre" -n 0 -zb 10000 -zw -au \
  -nt AUTO -ntmax "${SLURM_CPUS_PER_TASK}" \
  -pre "${OUT}/topology_test" -redo

printf '%s\n' "${MODEL}" > "${OUT}/model_used.txt"
date -Iseconds > "${OUT}/TOPOLOGY_TEST_COMPLETE"
