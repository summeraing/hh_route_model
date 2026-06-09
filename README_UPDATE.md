# GitHub Update For Existing Repository

Repository: https://github.com/summeraing/hh_route_model

Do not replace the whole repository. Add only these files to the existing structure:

- `README.md` (overwrite existing file)
- `code/reproduce_continuous_route_diagnostic.py`
- `code/run_expanded_fes_validation.py`
- `data/core_metrics/README.md` (overwrite existing file)
- `data/core_metrics/continuous_donor_role_probability_matrix.csv`
- `data/core_metrics/continuous_route_assignment_diagnostic.csv`
- `data/core_metrics/continuous_route_bootstrap_diagnostic.csv`
- `data/core_metrics/expanded_fes_compartment_summary.csv`
- `data/core_metrics/expanded_fes_fixed_rerun_neighbor_summary.csv`
- `data/core_metrics/expanded_fes_full_neighbor_summary.csv`
- `data/core_metrics/expanded_fes_marker_interpretation.csv`
- `data/core_metrics/expanded_fes_marker_manifest.csv`
- `data/core_metrics/expanded_fes_retrieval_summary.csv`
- `metadata/CITATION.cff` (overwrite existing file)
- `metadata/FILE_MANIFEST_SHA256.csv` (overwrite existing file)

Optional quick check:

```bash
python code/reproduce_continuous_route_diagnostic.py
python code/run_expanded_fes_validation.py --help
```

Expected result: the continuous diagnostic script prints the three prespecified donor-role top assignments and the joint bootstrap probability, then exits without assertion errors. The Fe-S validation command should print command-line help.

After updating GitHub, create a new Zenodo version from the updated repository. If Zenodo gives a new version-specific DOI, replace the DOI in the manuscript Data Availability, Code Availability, cover letter and submission-system fields before journal upload.
