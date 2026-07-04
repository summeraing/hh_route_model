# SOX transfer-validation case v2 results and limits

Date: 2026-07-01

## Input

Evidence table:

`15_iMETA_ROUTE_B_DUAL_CASE_TRANSFER_20260701/01_HGT_SOX_CASE/curated_seed/sox_evidence_units_expanded_v2.csv`

Main public sources:

- Gregersen et al. 2011, doi:10.3389/fmicb.2011.00116
- Berben et al. 2019, doi:10.3389/fmicb.2019.00160
- Ghosh et al. 2009, doi:10.1016/j.resmic.2009.07.003
- Meyer et al. 2007, doi:10.1111/j.1462-2920.2007.01407.x

The table combines:

- automated expansion of Gregersen 2011 Table 1 strain-by-gene presence data;
- Berben 2019 sulfur-gene module rows from table/text statements;
- seed rows from Ghosh 2009 and Meyer 2007;
- support rows retained but excluded from strict route scoring.

## SOX case-local route map

Donor classes:

- `lineage_core`
- `hgt_pathway`
- `mobile_context`

Role classes:

- `conserved_metabolic_backbone`
- `variable_energy_module`
- `mobility_boundary`

Prespecified route:

- `lineage_core -> conserved_metabolic_backbone`
- `hgt_pathway -> variable_energy_module`
- `mobile_context -> mobility_boundary`

## Preflight

The expanded v2 table passes the current minimum gate:

- route-eligible units: 65, threshold 60;
- sources: 4, threshold 3;
- dependency groups: 48, threshold 30;
- largest source fraction: 0.554, threshold 0.700;
- all donor classes represented;
- all role classes represented.

## Primary result

The prespecified SOX route ranks first:

- source-equal support: 1.000;
- best alternative support: 0.691;
- margin: 0.309.

Dependency collapse retains the result:

- collapsed route-eligible units: 48;
- collapsed prespecified rank: 1;
- collapsed margin: 0.307.

Leave-one-source tests retain the route as rank 1:

- removing Berben 2019: margin 0.106;
- removing Ghosh 2009: margin 0.412;
- removing Gregersen 2011: margin 0.306;
- removing Meyer 2007: margin 0.412.

## Null tests

Source-only role shuffling:

- observed margin: 0.309;
- q97.5 null margin: 0.261;
- empirical P = 0.0005.

Source-by-layer role shuffling:

- degenerate in v2;
- observed margin equals q97.5 null margin;
- empirical P = 1.000.

Interpretation: the SOX v2 table passes minimum preflight and demonstrates framework transfer, but the current source-by-layer null is not yet informative because several expanded source-layer blocks remain label-homogeneous. This must be disclosed if v2 is used in the manuscript.

## Recommended manuscript wording

Use:

> In a prokaryotic sulfur-oxidation transfer case, the same route-graph implementation passed preflight gates and ranked a case-local lineage-core, HGT-pathway and mobile-context map first under source-equal and dependency-collapsed scoring. This transfer case demonstrates workflow portability, while the homogeneous source-layer null identifies where finer event-level HGT curation is still needed.

Avoid:

> The SOX case proves a universal HGT route model.

Avoid:

> Source-layer null tests independently validate the SOX route.

## Next strengthening step

The highest-value next curation step is to add mixed source-layer blocks:

1. Extract event-level HGT rows from Meyer 2007 soxB phylogenies.
2. Extract component-level SoxXA/SoxYZ/SoxCD rows from Ghosh 2009 full text or figures.
3. Extract strain-level sulfur gene repertoire rows from Berben 2019 supplementary tables if available.

These additions should reduce label homogeneity inside source-layer strata and make source-by-layer permutation nulls meaningful.

