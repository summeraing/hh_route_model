#!/usr/bin/env bash
#SBATCH --job-name=route_sim_cal
#SBATCH --partition=base
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/01_route_sim_%A_%a.out
#SBATCH --error=logs/01_route_sim_%A_%a.err

set -euo pipefail
ROOT="${ROOT:?ROOT must point to the three-day upgrade directory}"
mkdir -p "${ROOT}/results/01_route_simulation/raw" "${ROOT}/logs"
python "${ROOT}/scripts/01_route_simulation.py" \
  --spec "${ROOT}/config/frozen_analysis_spec.json" \
  --scenario-index "${SLURM_ARRAY_TASK_ID}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --out "${ROOT}/results/01_route_simulation/raw/scenario_${SLURM_ARRAY_TASK_ID}.csv"
