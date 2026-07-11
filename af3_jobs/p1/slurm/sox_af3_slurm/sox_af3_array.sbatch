#!/usr/bin/env bash
#SBATCH --job-name=sox_af3_array
#SBATCH --output=logs/sox_af3_%A_%a.out
#SBATCH --error=logs/sox_af3_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00

set -euo pipefail

: "${ROOT_DIR:?ROOT_DIR must be exported by submit_sox_af3_array.sh}"
: "${JOB_LIST:?JOB_LIST must be exported by submit_sox_af3_array.sh}"
: "${OUT_DIR:?OUT_DIR must be exported by submit_sox_af3_array.sh}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
: "${AF3_SIF:?Set AF3_SIF to the AlphaFold 3 container path}"
: "${AF3_MODEL_DIR:?Set AF3_MODEL_DIR to the AlphaFold 3 model directory}"
: "${AF3_DB_DIR:?Set AF3_DB_DIR to the AlphaFold 3 database directory}"

AF3_RUN_SCRIPT="${AF3_RUN_SCRIPT:-/app/alphafold/run_alphafold.py}"

JSON_PATH="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_LIST}")"
if [[ -z "${JSON_PATH}" || ! -f "${JSON_PATH}" ]]; then
  echo "ERROR: Invalid JSON for task ${SLURM_ARRAY_TASK_ID}: ${JSON_PATH}" >&2
  exit 1
fi

JSON_NAME="$(basename "${JSON_PATH}")"
JOB_NAME="${JSON_NAME%.json}"
JOB_OUT="${OUT_DIR}/${JOB_NAME}"
mkdir -p "${JOB_OUT}"

echo "[SOX-AF3] task=${SLURM_ARRAY_TASK_ID} json=${JSON_NAME}"
echo "[SOX-AF3] output=${JOB_OUT}"

singularity exec --nv \
  -B "$(dirname "${JSON_PATH}"):/af_input" \
  -B "${JOB_OUT}:/af_output" \
  -B "${AF3_MODEL_DIR}:/models" \
  -B "${AF3_DB_DIR}:/public_databases" \
  "${AF3_SIF}" \
  python "${AF3_RUN_SCRIPT}" \
    --json_path="/af_input/${JSON_NAME}" \
    --model_dir=/models \
    --db_dir=/public_databases \
    --output_dir=/af_output \
    --flash_attention_implementation=triton
