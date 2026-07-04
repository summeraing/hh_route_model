#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SBATCH_FILE="${ROOT_DIR}/server_jobs/sox_af3_slurm/collect.sbatch"

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch is not available. Run this on the cluster login node." >&2
  exit 1
fi

mkdir -p "${ROOT_DIR}/logs"
cd "${ROOT_DIR}"
sbatch --export=ALL,ROOT_DIR="${ROOT_DIR}" "${SBATCH_FILE}"
