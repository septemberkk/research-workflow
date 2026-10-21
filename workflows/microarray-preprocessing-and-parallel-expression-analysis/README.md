# Microarray Preprocessing & Parallel Expression Analysis

Independent preparation of mRNA and miRNA expression arrays: explicit sample
scope, reviewed probe annotations, arithmetic feature consolidation and
branch-specific normalization.

The two assays remain separate. Common group names do not establish biological
sample pairing, and the workflow does not infer miRNA targets, regulatory
relationships, differential expression or functional enrichment.

## Processing

Each layer retains the selected sample order, records every input probe's mapping
outcome and averages probes assigned to the same unambiguous feature. Missing
values are rejected rather than replaced by zero. Sample names and feature
identifiers are not shortened, renamed or updated through guessed aliases.

mRNA defaults to no normalization; miRNA defaults to quantile normalization.
Both choices are explicit and independently configurable. Averaging operates on
the declared supplied scale. No background correction or logarithmic transform
is performed, and unlogged intensities are never described as log2 fold changes.

## Inputs

Three UTF-8 tab-separated tables are supplied independently for each assay:

- Expression: `probe_id` and numeric sample columns.
- Annotation: `probe_id` and `feature_id`.
- Sample scope: `sample_id`, `group` and an explicit `include` decision of
  `yes` or `no`.

The declared scale is `intensity` or `log2`. Intensity values must be nonnegative.
All expression values must be finite; every sample column must have one metadata
record. At least two samples must be selected in each layer.
Quantile normalization requires at least two retained feature rows.

Blank feature assignments and compound annotations remain visible as unmapped or
ambiguous records. They are not resolved by taking the first symbol. Probe IDs,
annotation keys and sample keys must be unique. Identifiers are case-sensitive;
quoted fields, unsafe spreadsheet prefixes, invisible formatting and boundary
whitespace are rejected without silently rewriting their contents.

## Outputs

Each assay exports selected samples, annotation provenance, per-probe mapping
decisions, feature memberships, probe counts, collapsed and processed matrices,
and sample-level distribution summaries. Optional figures compare each assay's
collapsed and processed expression distributions independently.

Run metadata records declared scales, normalization choices, software versions
and completion state. Input and output manifests contain role-based checksums
without input names or locations. Exports retain supplied identifiers and should
be reviewed before sharing.

## Scope

Only illustrative inputs are included. The workflow prepares expression
profiles, not differential-expression results.

R is required; limma is needed only for quantile normalization. No package
installation, network access or source-input modification occurs during execution.
Licensed under the MIT License.
