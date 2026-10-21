# Methods and interpretation

## Sample scope

Expression columns are matched exactly to the supplied sample table. Selection
and output order follow that table, not assumed column positions or inferred
sample-name prefixes. Excluded samples remain documented but do not contribute
to the normalization reference. No cross-assay matching is performed.

## Probe annotation and consolidation

A probe key has at most one annotation row. Blank feature IDs are unmapped.
Compound feature descriptions containing separators or whitespace are treated
as ambiguous. No alias conversion, version update, identifier repair or
first-symbol selection is applied.

Unambiguous probe-to-feature memberships remain in a long-form export.
Within each selected sample, probes sharing a feature ID contribute equally to
its arithmetic mean. Features are then ordered by their exact identifier.
Averaging intensities and averaging log2 values are different operations; the
declared scale is recorded and is not transformed automatically.

Prepared feature matrices are descriptive outputs. No probe-level or
feature-level differential-expression model is fitted.

## Independent normalization

The no-normalization option preserves the collapsed matrix. Quantile
normalization uses limma's `normalizeQuantiles` with explicit tie handling and
the complete selected-sample feature matrix, separately for each assay.
It aligns empirical distributions; it does not infer whether that assumption
fits the experimental design, perform background correction or change intensity
values into log2 units.

The feature set, annotation source and normalization choice remain
experiment-specific decisions. No normalization reference is shared between mRNA
and miRNA, and no cross-layer correlation or targeting relationship is estimated.

## Inspection and provenance

Distribution summaries and optional boxplots describe the prepared values.
They do not provide a differential-expression test, a fold-change estimate,
evidence of sample pairing or a biological mechanism.

Manifests bind the provided annotation, sample scope and expression files to
their outputs using MD5 checksums. These are comparison records, not security
signatures or authenticated authorship. Input locations are deliberately omitted.

[limma documentation](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
describes the quantile-normalization function and its treatment of ties.
