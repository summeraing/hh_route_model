# Expanded SOX-AF3 job inputs

This directory contains the additional JSON inputs and manifests used for the expanded SOX structural screen. The JSON files are computational inputs, not experimental observations.

Run GPU predictions through the local scheduler rather than on a login node. The Slurm array template under `af3_jobs/p1/slurm/sox_af3_slurm/` accepts cluster-specific `AF3_SIF`, `AF3_MODEL_DIR` and `AF3_DB_DIR` environment variables and can be adapted to this JSON directory.

Processed confidence summaries used by the manuscript are provided under `data/sox_af3/`; rerunning AlphaFold 3 requires licensed model parameters and local database assets that are not redistributed here.
