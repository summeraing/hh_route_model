# SOX-AF3 Slurm job packet

Purpose: run the SOX-AF3 predicted-interface screen as scheduled GPU jobs, not on the login node.

This directory is a submission template. It assumes the server packet contains:

- `json_screen/*.json`: AlphaFold 3 JSON jobs.
- `af3_output_screen/`: output directory created before or during the run.
- `scripts/summarize_sox_af3_results.py`: summary parser from the public release, only needed for the collection step.
- `target_matrix/sox_af3_target_matrix_v1.csv`: AF3 target matrix, only needed for the collection step.

The current server AF3 asset paths are hard-coded in `sox_af3_array.sbatch`:

```bash
AF3_SIF=/jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/06.PPI/af3_assets/alphafold3/alphafold3.sif
AF3_MODEL_DIR=/jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/06.PPI/af3_assets/models
AF3_DB_DIR=/jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/06.PPI/af3_assets/public_databases
AF3_RUN_SCRIPT=/app/alphafold/run_alphafold.py
```

## Server-side run

From the uploaded packet root on the server:

```bash
cd /jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/01.network_model/SOX_AF3_slurm_packet
bash server_jobs/sox_af3_slurm/submit_sox_af3_array.sh
```

The wrapper writes `jobs.txt` and submits one Slurm array task per JSON file.

Check status:

```bash
squeue -u "$USER"
```

After all AF3 array tasks finish, summarize the outputs with:

```bash
bash server_jobs/sox_af3_slurm/submit_collect_sox_af3.sh
```

This submits `collect.sbatch`, a lightweight Slurm job that parses AF3 confidence JSONs into:

```text
af3_output_screen/sox_af3_summary_from_slurm.csv
```

## Rules

- Do not run `singularity exec --nv ... run_alphafold.py` directly on the login node.
- The login node is only for file transfer, `sbatch`, `squeue` and lightweight checks.
- If a cluster requires a partition or account, add it to the `#SBATCH` header in `sox_af3_array.sbatch` before submission.
