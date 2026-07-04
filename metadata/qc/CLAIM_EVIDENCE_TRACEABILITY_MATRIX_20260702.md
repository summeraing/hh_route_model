# Claim-to-evidence traceability matrix

Date: 2026-07-02

## Purpose

This matrix links each major manuscript claim to its figure, Source Data sheets, public-release scripts and current reproducibility status. It is intended to make the route-resolved framework auditable before submission and to prevent claims from drifting away from data, code or figures.

## Summary

- Claims tracked: 10.
- Claims with text and data present: 10.
- Claims with partial reproducibility risk: 0.
- Final main-figure assets exist for Figures 1-6 in the working figure folder.
- Fe-S and eukaryogenesis AF3 layers remain framed as support-layer screens: Source Data and render scripts are present, but the manuscript does not claim a full one-command raw-tree or raw-AF3 reconstruction for those layers.

## Matrix

| Claim | Main-text location | Evidence status | Primary figure | Source Data sheets | Public scripts | Current risk |
|---|---|---|---|---|---|---|
| C01. Route-resolved evidence graphs separate donor identity from organizational role and score prespecified donor-role maps under source-balanced support with dependency collapse. | Introduction/Methods lines 19-35 and 97-109 | Text and data present | Figure 1 | `Vita_index`; `unit_granularity`; `route_space_defense`; `source_profile`; `source_layer_profile` | `run_route_graph_case.py` | Final Figure 1 asset exists; author visual review pending. |
| C02. In eukaryogenesis, the prespecified host-scaffold / alphaproteobacterial-energetic / transition-mosaic map ranks first after source-equal scoring and audit, with support 0.769 and margin 0.355. | Results lines 39-45 | Text and data present | Figure 2 | `route_rank_audit`; `source_rates_audit`; `leave_source_audit`; `leave_layer_audit`; `subset_robustness`; `audit_overlay` | `run_routeB_expanded_computational_evidence.py`; `run_route_graph_case.py` | Final Figure 2 asset exists; method-dense but acceptable for current branch. |
| C03. Continuous donor-role modelling supports the same qualitative route organization without requiring a hard one-to-one bijection. | Results line 43 | Text and data present | Figure 2 | `continuous_role_matrix`; `continuous_top_roles`; `continuous_bootstrap`; `calc_prob_summary`; `calc_prob_routes`; `calc_prob_donor_diag` | `run_routeB_computational_upgrades.py` | Final Figure 2 asset exists; method-dense but acceptable for current branch. |
| C04. The eukaryogenesis route signal remains stable under hierarchical resampling, projected audit-label noise, Shapley/source contribution analysis and adversarial row redirection. | Results line 45 and Figure 4 legend | Text and data present | Figure 4 | `calc_hier_summary`; `calc_hier_draws`; `calc_source_shapley`; `calc_layer_shapley`; `calc_row_adv_summary`; `calc_row_adv_selected`; `calc_adversarial_summary`; `calc_targeted_shift`; `calc_uniform_noise` | `run_routeB_hierarchical_resampling_and_shapley.py`; `run_routeB_row_level_adversarial.py`; `run_routeB_adversarial_counterfactual.py` | Final Figure 4 asset exists; robustness figure remains method-dense by design. |
| C05. Expanded Fe-S phylogenies support localized energetic incorporation: mitochondrial-integration markers show alphaproteobacterial enrichment, whereas cytosolic CIA markers do not. | Results lines 49-51 | Text and data present | Figure 3 | `fes_marker_manifest`; `fes_retrieval`; `fes_full_trees`; `fes_fixed_rerun`; `fes_compartment`; `fes_marker_readout`; `FeS_readout` | Source Data tables and Figure 3 render script present; no one-command public raw-tree reconstruction claimed for this support layer | Support layer is traceable from Source Data and Figure 3 rendering; full raw tree reconstruction is not claimed as a one-command public rebuild. |
| C06. Focused eukaryogenesis AF3 screens support predicted structural plausibility for selected local energetic/Fe-S interfaces while unrelated scaffold-energy decoys remain low after confidence gating. | Results line 53 | Text and data present | Figure 3 | `AF3_FeS_pairs`; `AF3_ATP_pairs`; `AF3_controls`; `AF3_integrated_metrics`; `af3_confidence`; `af3_interface_metrics`; `af3_digest`; `AF3_integrated_metrics_tsv`; `AF3_negative_controls_tsv` | Source Data tables and Figure 3 render script present; eukaryogenesis AF3 raw generation/post-processing not claimed as one-command public rebuild | Support layer is traceable from Source Data and Figure 3 rendering; eukaryogenesis AF3 predictions remain prediction-limited. |
| C07. The SOX-HGT transfer case passes preflight gates and the prespecified lineage-core / HGT-pathway / mobile-context route ranks first with support 1.000 and margin 0.309. | Results lines 57-61 | Text and data present | Figure 5 | `SOX_evidence_v2`; `SOX_summary`; `SOX_routes`; `SOX_leave_source`; `SOX_perm_nulls` | `build_sox_evidence_units_v2_from_public_tables.py`; `run_route_graph_case.py` | Final Figure 5 asset exists; outdated provisional file removed from working figure folder. |
| C08. SOX annotation stress preserves the route result after excluding a low-confidence dsrABC rare-presence row, while source-by-layer null degeneracy defines a boundary. | Results lines 63-65 | Text and data present | Figure 5 | `SOX_ann_summary`; `SOX_ann_flagged`; `SOX_stress_summary`; `SOX_stress_leave`; `SOX_stress_perm` | `build_sox_annotation_crosscheck.py`; `run_route_graph_case.py` | Final Figure 5 asset exists; source-by-layer null degeneracy retained as boundary condition. |
| C09. SOX-AF3 screens support route-relevant assemblies over unrelated housekeeping decoys, while Tier 2 identifies a moderate mobile-boundary control and prevents a blanket incompatibility claim. | Results lines 69-75 | Text and data present | Figure 6 | `SOX_AF3_FINAL`; `SOX_AF3_FAMILY_GATE`; `SOX_AF3_TIER2`; `SOX_AF3_COMBINED`; `SOX-AF3 processed CSVs` | `build_sox_af3_json_inputs.py`; `make_sox_af3_server_packet.py`; `summarize_sox_af3_results.py`; `evaluate_sox_af3_gates.py`; `prepare_sox_af3_tier2_targets.py` | Final Figure 6 asset exists; keep claim prediction-limited and boundary-aware. |
| C10. The framework reports failure boundaries rather than forcing every subset to support the route model; eukaryogenesis narrow mechanistic subset and SOX source-by-layer null define limits. | Results lines 75-79 and Discussion lines 89-93 | Text and data present | Figure 4 and Figure 5 | `subset_robustness`; `calc_future_boundary`; `calc_accum_summary`; `SOX_stress_perm`; `SOX_perm_nulls` | `run_routeB_evidence_accumulation_forecast.py`; `run_route_graph_case.py` | Final Figure 4 and Figure 5 assets exist; failure-boundary framing retained. |

## Reproducibility implications

The route-graph core and SOX transfer case have the strongest public-script coverage. The eukaryogenesis Fe-S and eukaryogenesis AF3 layers are traceable through Source Data tables, methods text and Figure 3 rendering code, but they should remain described as support-layer screens rather than as fully automated raw-data reconstructions.

## Immediate actions before submission

1. Author visual review of final Figures 1-6.
2. Keep all AF3 claims prediction-limited.
3. Synchronize GitHub/Zenodo release if these final figure scripts or Source Data changes are meant to be public.
4. Run `python tools/run_all_readiness_checks.py` after any further text or figure change.
