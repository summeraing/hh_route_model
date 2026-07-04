# SOX-AF3 Tier 2 Run Status

Date: 2026-07-02

## Purpose

This Tier 2 screen is an exploratory strengthening run for the SOX/HGT transfer-validation case. It expands the structural-interactivity check beyond the priority-1 SOX-AF3 design. The final Slurm summary has now been collected. The screen is suitable as a diagnostic structural-gate sensitivity layer, but it should not be presented as a clean confirmatory validation that all boundary or mobile-element pairings are structurally incompatible.

## Local packet

- Rebuilt short-path packet source: `%TEMP%\SOX_AF3_TIER2_SERVER_PACKET_20260702`
- Local upload zip: `TOP_JOURNAL_REBUILD_ANALYSES_v1/SOX_AF3_TIER2_SLURM_20260702.zip`
- Zip SHA256: `dfdbeca2fa39127a2367d7eb35acf08d1552ea5dd4c99f4a84b878dfa2e124b9`
- Files in zip: 18
- AF3 JSON jobs in zip: 6
- Internal zip path check: no absolute paths or Windows backslashes

## Server packet

- Server directory: `/jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/01.network_model/SOX_AF3_tier2_slurm_packet_20260702`
- Server-side file count after unzip: 18
- Server-side JSON count after unzip: 6
- Slurm script syntax check: passed

## Submission

The first submission attempt failed because the sbatch template requested GPU resources without selecting the cluster `gpu` partition. The local template and first SOX-AF3 packet template were patched to include:

```bash
#SBATCH --partition=gpu
```

The patched sbatch file was copied to the server Tier 2 packet and submitted with:

```bash
cd /jinxianstor/home/xiazhiqiang1/02.work/02.2026_work/01.network_model/SOX_AF3_tier2_slurm_packet_20260702
bash server_jobs/sox_af3_slurm/submit_sox_af3_array.sh
```

Submitted Slurm array job:

- Job ID: `274850`
- Tasks: `274850_[1-6]`
- Partition: `gpu`
- Status at submission check: pending

Follow-up checks:

- The array entered running state on `gpunode1` for all six tasks.
- A dynamic attempt to extend the running job time limit to 24 hours with `scontrol update JobId=274850 TimeLimit=24:00:00` failed with `Access/permission denied`.
- Follow-up queue checks showed the tasks running normally within the 12-hour limit. Slurm `TIME` values such as `14:24` and `19:27` are minute:second values at this early stage, not hours.
- Logs show AF3 started correctly and entered standard MSA/template-processing stages; no collected summary is available yet.
- A dependent collection job was submitted with `afterany:274850` so that the summary step runs automatically after the AF3 array exits.
- Collection job ID: `274856`
- Collection job status at submission check: `PENDING (Dependency)` in the `base` partition

## Partial summary

A non-final partial summary was generated while two rpoB-decoy jobs were still running:

- Remote partial CSV: `af3_output_screen/sox_af3_tier2_partial_summary_current.csv`
- Local copy: `07_SOX_AF3_TIER2/sox_af3_tier2_partial_summary_current.csv`
- Parsed completed jobs: 4

Partial results:

| job_id | comparison_group | chain_set | median ipTM | observed gate |
|---|---|---|---:|---|
| SOX_AF3_FUNC_004 | functional_route | soxX+soxA+soxY+soxZ | 0.74 | high |
| SOX_AF3_NEG_013 | mobile_boundary_negative | hdrA+transposase | 0.57 | moderate |
| SOX_AF3_POS_008 | cognate_positive | hdrA+hdrB+hdrC | 0.49 | moderate |
| SOX_AF3_POS_009 | cognate_positive | dsrA+dsrB | 0.83 | high |

Interpretation of partial results:

- The functional SOX assembly and DsrAB positive are strong in the partial screen.
- HdrABC is moderate, not high.
- The hdrA+transposase negative control is also moderate. This means Tier 2 currently does not support a simple "all positives high, all mobile-boundary negatives low" interpretation.
- These partial results should not be promoted to manuscript claims. If included later, Tier 2 should be framed as a boundary/diagnostic screen, not as straightforward confirmatory evidence.

A later 5-job partial summary was generated after `SOX_AF3_DECOY_016` completed:

- Remote partial CSV: `af3_output_screen/sox_af3_tier2_partial_summary_5jobs.csv`
- Local copy: `07_SOX_AF3_TIER2/sox_af3_tier2_partial_summary_5jobs.csv`
- Parsed completed jobs: 5

New result:

| job_id | comparison_group | chain_set | median ipTM | observed gate |
|---|---|---|---:|---|
| SOX_AF3_DECOY_016 | housekeeping_decoy | soxY+rpoB_decoy | 0.28 | low |

Updated interpretation:

- The housekeeping decoy `soxY+rpoB` is low, which supports the decoy-gate logic for at least one unrelated pairing.
- The mobile-boundary negative `hdrA+transposase` remains moderate, so the Tier 2 result is mixed.
- `SOX_AF3_DECOY_015` was still running at this intermediate checkpoint and required final collection before interpretation.

## Final collected summary

The Slurm array and dependent collection job have completed. The final summary was retrieved locally:

- Remote CSV: `af3_output_screen/sox_af3_summary_from_slurm.csv`
- Local copy: `07_SOX_AF3_TIER2/sox_af3_summary_from_slurm.csv`
- Parsed completed jobs: 6

Final results:

| job_id | comparison_group | chain_set | expected gate | median ipTM | observed gate |
|---|---|---|---|---:|---|
| SOX_AF3_DECOY_015 | housekeeping_decoy | soxB+rpoB_decoy | low | 0.17 | low |
| SOX_AF3_DECOY_016 | housekeeping_decoy | soxY+rpoB_decoy | low | 0.28 | low |
| SOX_AF3_FUNC_004 | functional_route | soxX+soxA+soxY+soxZ | moderate_or_high | 0.74 | high |
| SOX_AF3_NEG_013 | mobile_boundary_negative | hdrA+transposase | low | 0.57 | moderate |
| SOX_AF3_POS_008 | cognate_positive | hdrA+hdrB+hdrC | moderate_or_high | 0.49 | moderate |
| SOX_AF3_POS_009 | cognate_positive | dsrA+dsrB | high | 0.83 | high |

Final interpretation:

- The SOX functional assembly is high-confidence by the Tier 2 gate.
- The DsrAB cognate positive is high-confidence.
- HdrABC is moderate, not high, consistent with a weaker or less completely captured structural signal.
- Both rpoB housekeeping decoys are low, supporting the unrelated-decoy gate.
- The hdrA+transposase mobile-boundary negative is moderate. This prevents a broad claim that mobile-boundary pairings are uniformly low-confidence.
- Tier 2 should therefore be used as a boundary-aware diagnostic layer: it supports separation between functional/cognate route assemblies and unrelated housekeeping decoys, but it does not justify a simple positive-versus-mobile-boundary dichotomy.

Recommended manuscript use:

- Use Tier 2 in Supplementary Results or Source Data as an exploratory sensitivity screen.
- In the main text, if mentioned, write that the expanded SOX-AF3 screen supported functional and cognate route assemblies over unrelated housekeeping decoys, while also identifying a mobile-boundary control with moderate predicted interface confidence.
- Do not use Tier 2 to claim direct structural proof of SOX-HGT history or universal mobile-element incompatibility.
