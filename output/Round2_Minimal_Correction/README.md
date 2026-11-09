# Corrected moment and distribution workflow

From the repository root:

```matlab
addpath('output/Round2_Minimal_Correction/Programs','-begin')
run_round2('check')          % formulas, GND inversion, processed FEM inputs
run_round2('tables')         % corrected predictions; FEM-referenced errors
run_round2('single')         % standalone moment panels
run_round2('distributions')  % PDFs/CDFs, JS/KL, exceedance, Q99
```

Inputs are the four measured-statistics workbooks in `output` and compact
field metrics in `Results/Cache/Field_Distributions`. The Summary sheet is
authoritative for the Vol20 and point-process inputs. No FEM solve is run.

To recompute metrics from raw fields, put them in `output/FieldValue`, or set
`POROUSMEDIA_OUTPUT_DIR` to the external output directory. Then call
`round2_collect_field_metrics(true)`. This reads existing fields, not ODBs.
Afterward, `run_round2_distributions('plots')` uses the resulting metrics.

Byte/time provenance checks apply to fresh raw-field processing. Bundled
aggregate metrics permit plotting without raw files. Original templates are
in `output/Style_Templates`. Generated figures/logs are ignored by Git.
Re-exporting CSVs can update paths without changing measured FEM numbers.
