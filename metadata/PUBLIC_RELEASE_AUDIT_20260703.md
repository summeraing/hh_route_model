# Public Release Packet Audit

Packet: `GH_PACKET_v2026_07_03_dual_case_af3_framework`

## Included

- analysis code required for route-graph scoring and SOX-AF3 post-processing;
- processed evidence tables and Source Data;
- compact AF3 job inputs, manifests and processed summaries;
- final figure PNGs for orientation;
- QA and traceability reports.

## Excluded

- manuscript drafts;
- cover letters;
- Supplementary Information documents;
- peer-review correspondence;
- PPT files;
- planning notes;
- raw bulky AF3 output directories.

## Upload use

This folder is intended as a repository-root overlay for `https://github.com/summeraing/hh_route_model`.

Recommended tag: `v2026.07.03-dual-case-af3-framework`

## Local verification

Run:

```bash
python code/run_smoke_tests.py
```

Expected sentinel:

```text
SMOKE_TESTS_PASS
```

The smoke test verifies the eukaryogenesis benchmark core metrics and
continuous-route diagnostics before running the SOX-HGT route graph and SOX-AF3
gate check.
