#!/usr/bin/env bash
#SBATCH --job-name=sox_toposum
#SBATCH --partition=base
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/05_topology_aggregate_%j.out
#SBATCH --error=logs/05_topology_aggregate_%j.err

set -euo pipefail
ROOT="${ROOT:?ROOT must point to the three-day upgrade directory}"
python "${ROOT}/scripts/05_parse_topology_tests.py" \
  --root "${ROOT}/results/03_topology_tests" \
  --outdir "${ROOT}/results/03_topology_tests/summary"
date -Iseconds > "${ROOT}/results/03_topology_tests/TOPOLOGY_COMPLETE"
