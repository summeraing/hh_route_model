#!/usr/bin/env bash
#SBATCH --job-name=sox_af3_collect
#SBATCH --output=logs/sox_af3_collect_%j.out
#SBATCH --error=logs/sox_af3_collect_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

: "${ROOT_DIR:?ROOT_DIR must be exported by submit_collect_sox_af3.sh}"

AF3_OUTPUT_DIR="${ROOT_DIR}/af3_output_screen"
TARGET_MATRIX="${ROOT_DIR}/target_matrix/sox_af3_target_matrix_v1.csv"
SUMMARY_SCRIPT="${ROOT_DIR}/scripts/summarize_sox_af3_results.py"
OUTPUT_CSV="${AF3_OUTPUT_DIR}/sox_af3_summary_from_slurm.csv"

if [[ ! -d "${AF3_OUTPUT_DIR}" ]]; then
  echo "ERROR: Missing AF3 output directory: ${AF3_OUTPUT_DIR}" >&2
  exit 1
fi
if [[ ! -f "${TARGET_MATRIX}" ]]; then
  echo "ERROR: Missing target matrix: ${TARGET_MATRIX}" >&2
  exit 1
fi
if [[ ! -f "${SUMMARY_SCRIPT}" ]]; then
  echo "ERROR: Missing summary script: ${SUMMARY_SCRIPT}" >&2
  exit 1
fi

python "${SUMMARY_SCRIPT}" \
  --af3-output-dir "${AF3_OUTPUT_DIR}" \
  --target-matrix "${TARGET_MATRIX}" \
  --output-csv "${OUTPUT_CSV}"

echo "Wrote ${OUTPUT_CSV}"

