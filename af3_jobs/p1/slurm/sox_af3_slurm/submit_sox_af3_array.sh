#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JOB_DIR="${ROOT_DIR}/json_screen"
OUT_DIR="${ROOT_DIR}/af3_output_screen"
JOB_LIST="${ROOT_DIR}/jobs.txt"
SBATCH_FILE="${ROOT_DIR}/server_jobs/sox_af3_slurm/sox_af3_array.sbatch"

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch is not available. Run this on the cluster login node." >&2
  exit 1
fi

if [[ ! -d "${JOB_DIR}" ]]; then
  echo "ERROR: Missing AF3 JSON directory: ${JOB_DIR}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"
find "${JOB_DIR}" -maxdepth 1 -type f -name "*.json" | sort > "${JOB_LIST}"

N_JOBS="$(wc -l < "${JOB_LIST}" | tr -d ' ')"
if [[ "${N_JOBS}" == "0" ]]; then
  echo "ERROR: No AF3 JSON jobs found in ${JOB_DIR}" >&2
  exit 1
fi

echo "Submitting ${N_JOBS} AF3 jobs as a Slurm array."
echo "Root: ${ROOT_DIR}"
echo "Job list: ${JOB_LIST}"

cd "${ROOT_DIR}"
sbatch --array="1-${N_JOBS}" \
  --export=ALL,ROOT_DIR="${ROOT_DIR}",JOB_LIST="${JOB_LIST}",OUT_DIR="${OUT_DIR}" \
  "${SBATCH_FILE}"
