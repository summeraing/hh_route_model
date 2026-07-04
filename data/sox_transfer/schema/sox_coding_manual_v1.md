# SOX HGT transfer-case coding manual v1

## Scientific question

The SOX transfer case asks whether donor-origin signals in a prokaryotic HGT metabolic island are diffusely mixed across organizational roles or preferentially routed into specific roles.

The prespecified map is:

- `lineage_core -> conserved_metabolic_backbone`
- `hgt_pathway -> variable_energy_module`
- `mobile_context -> mobility_boundary`

## Evidence-unit definition

An evidence unit is one traceable public-data entry linking a SOX-associated gene, homolog group, operon component, genomic context, phylogenetic placement or structural/interface candidate to both:

1. a donor-origin class; and
2. an organizational role class.

Rows are not assumed to be independent. Repeated genes from the same source, cluster or homolog family must share a `dependency_group` for collapse tests.

## Donor classes

### lineage_core

Use when evidence supports a vertically retained or lineage-conserved sulfur-metabolism signal. Examples:

- sulfur-metabolism genes present across most or all strains in a focal lineage;
- genes whose phylogeny follows the organismal background more closely than the SOX island;
- core sulfur-energy modules treated by the source as part of the stable repertoire.

### hgt_pathway

Use when evidence supports an acquired, lineage-variable or taxonomically discordant metabolic component. Examples:

- SOX catalytic genes whose phylogeny conflicts with host taxonomy;
- operon modules shared across distant taxa with patchy distribution;
- pathway components explicitly annotated as HGT-derived in the source.

### mobile_context

Use when evidence supports a mobile element, recombination boundary or island-carriage signal. Examples:

- transposase, integrase, resolvase or recombinase;
- plasmid-borne context;
- genomic island boundary marker;
- insertion sequence adjacent to SOX modules.

### uncertain

Use when donor origin is not traceable from the evidence unit. `uncertain` rows remain in the ledger but are excluded from strict route support.

## Organizational-role classes

### conserved_metabolic_backbone

Use for sulfur-energy genes or operons treated by a source as a lineage-conserved metabolic backbone.

### variable_energy_module

Use for SOX catalytic, electron-transfer or sulfur-energy modules treated by a source as acquired, variable, lineage-discordant or gained/lost.

### mobility_boundary

Use for island boundaries, plasmid context, recombination, transposition and mobility-associated roles.

### support

Use for descriptive rows that support the evidence table but do not define one of the three strict route roles.

### uncertain

Use when organizational role is ambiguous. Exclude from strict route support.

## Ambiguity rules

- If a row has a catalytic SOX gene and a mobile-context annotation, code the gene's functional role as `variable_energy_module` and use `mobile_context` only if donor origin is specifically mobility-associated.
- If a row describes a transposase adjacent to a SOX operon, code donor as `mobile_context` and role as `mobility_boundary`.
- If a source reports a whole operon without gene-level detail, use one row for the operon and assign `dependency_group` at the operon level.
- If multiple rows come from the same gene family and source, assign the same `dependency_group` so collapse tests do not treat them as independent.

## Minimum evidence required for analysis

For a first transfer-validation run:

- at least 3 independent public sources or database-derived evidence layers;
- at least 60 route-eligible rows;
- all three donor classes represented;
- all three role classes represented;
- no single source contributing more than 70% of route-eligible rows unless source-equal and leave-one-source tests are emphasized.
