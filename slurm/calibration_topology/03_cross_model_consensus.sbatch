#!/usr/bin/env bash
#SBATCH --job-name=sox_consensus
#SBATCH --partition=base
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/03_cross_model_consensus_%j.out
#SBATCH --error=logs/03_cross_model_consensus_%j.err

set -euo pipefail
ROOT="${ROOT:?ROOT must point to the three-day upgrade directory}"
ATLAS="${ATLAS:?ATLAS must point to the frozen SOX atlas directory}"
mkdir -p "${ROOT}/results/02_cross_model_consensus" "${ROOT}/logs"
python "${ROOT}/scripts/03_cross_model_consensus.py" \
  --spec "${ROOT}/config/frozen_analysis_spec.json" \
  --fixed-support "${ATLAS}/results/39_generax_calibrated_summary/generax_per_branch_support.tsv" \
  --mfp-support "${ATLAS}/results/44_generax_modelselected_summary/generax_per_branch_support.tsv" \
  --genome-table "${ATLAS}/results/41_cross_layer_calibrated/sox_cross_layer_genome_table.tsv" \
  --outdir "${ROOT}/results/02_cross_model_consensus"
date -Iseconds > "${ROOT}/results/02_cross_model_consensus/CONSENSUS_COMPLETE"
