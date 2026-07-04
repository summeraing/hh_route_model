# SOX-AF3 Tier 2 target selection

Source target matrix: `TOP_JOURNAL_REBUILD_ANALYSES_v1\iMETA_ROUTE_B_REBUILD_20260630\17_SOX_AF3_VALIDATION_DESIGN_20260701\target_matrix\sox_af3_target_matrix_v1.csv`
Sequence manifest: `TOP_JOURNAL_REBUILD_ANALYSES_v1\iMETA_ROUTE_B_REBUILD_20260630\17_SOX_AF3_VALIDATION_DESIGN_20260701\sequence_candidates_v3\sox_sequence_manifest_candidate_review.csv`
Priority max: 3
Completed priority-1 jobs excluded: 7
Included Tier 2/3 jobs: 6
Excluded jobs: 10

## Included jobs

- `SOX_AF3_FUNC_004`: soxX+soxA+soxY+soxZ (functional_route, expected moderate_or_high)
- `SOX_AF3_POS_008`: hdrA+hdrB+hdrC (cognate_positive, expected moderate_or_high)
- `SOX_AF3_POS_009`: dsrA+dsrB (cognate_positive, expected high)
- `SOX_AF3_NEG_013`: hdrA+transposase (mobile_boundary_negative, expected low)
- `SOX_AF3_DECOY_015`: soxB+rpoB_decoy (housekeeping_decoy, expected low)
- `SOX_AF3_DECOY_016`: soxY+rpoB_decoy (housekeeping_decoy, expected low)

## Excluded jobs

- `SOX_AF3_POS_001`: already_completed_priority1
- `SOX_AF3_POS_002`: already_completed_priority1
- `SOX_AF3_FUNC_003`: already_completed_priority1
- `SOX_AF3_FUNC_005`: missing_sequences:soxC;soxD
- `SOX_AF3_POS_006`: already_completed_priority1
- `SOX_AF3_POS_007`: missing_sequences:soeA;soeB;soeC
- `SOX_AF3_NEG_010`: already_completed_priority1
- `SOX_AF3_NEG_011`: already_completed_priority1
- `SOX_AF3_NEG_012`: already_completed_priority1
- `SOX_AF3_NEG_014`: missing_sequences:soeB

## Interpretation boundary

Tier 2 is an exploratory strengthening screen. It extends the SOX predicted-interface layer beyond the priority-1 manuscript set but should not be interpreted as experimental binding evidence. Sequence rows marked `candidate_needs_review` require manual accession/provenance review before being promoted to primary manuscript evidence.
