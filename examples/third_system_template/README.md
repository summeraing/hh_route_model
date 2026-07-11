# Third-system template

`evidence_units_template.csv` is a deliberately small toy dataset for checking the route-graph schema and adapting the workflow to another biological system. Replace the donor and role vocabularies only after defining biologically meaningful contributor classes, organizational roles, source identifiers and dependency groups.

Required columns are `case_id`, `unit_id`, `source_id`, `evidence_layer`, `donor_class`, `functional_class`, `route_eligible` and `dependency_group`. Recommended provenance columns include `source_type`, `evidence_unit`, `module_id`, `gene_symbol`, `taxon_or_genome`, `confidence`, `evidence_text`, `provenance` and `notes`.

The toy rows are not used in the manuscript and must not be interpreted biologically.
