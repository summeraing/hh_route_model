#!/usr/bin/env bash
#SBATCH --job-name=route_sim_sum
#SBATCH --partition=base
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/02_route_sim_aggregate_%j.out
#SBATCH --error=logs/02_route_sim_aggregate_%j.err

set -euo pipefail
ROOT="${ROOT:?ROOT must point to the three-day upgrade directory}"
python "${ROOT}/scripts/02_aggregate_simulation.py" \
  --spec "${ROOT}/config/frozen_analysis_spec.json" \
  --input-dir "${ROOT}/results/01_route_simulation/raw" \
  --outdir "${ROOT}/results/01_route_simulation/summary"
date -Iseconds > "${ROOT}/results/01_route_simulation/SIMULATION_COMPLETE"
