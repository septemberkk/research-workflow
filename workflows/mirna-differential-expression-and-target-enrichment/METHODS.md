# Methods and interpretation

## Analysis unit and scale

One run represents a comparison between two independent groups. Sample labels
and sample order are explicit. A reference/comparison design is fitted without
an intercept; the contrast is comparison minus reference. Positive logFC means
higher miRNA expression in the comparison group.

Inputs are already processed log2 expression values. Negative log2 values are
permitted. No background correction or automatic log transform is performed.
Raw sequencing counts need an appropriate count-aware workflow instead.

Repeated miRNA IDs are rejected unless arithmetic averaging is explicitly
requested. Averaging occurs on the supplied log2 values before any requested
between-array normalization. Distinct IDs are never merged through alias or
case conversion.

## Differential expression

The analysis calls limma's `normalizeBetweenArrays(method = "quantile")` unless
normalization is explicitly disabled, then `lmFit`, `contrasts.fit` and
`eBayes(trend = FALSE, robust = FALSE)`. It exports all features with
`topTable(number = Inf, adjust.method = "BH", sort.by = "P")`.

The selection rule uses strict inequalities: BH-adjusted value below the chosen
threshold and absolute logFC above the chosen effect-size cutoff. Defaults are
0.05 and 1. This is an empirical-Bayes test against a zero-effect null followed
by a descriptive candidate screen, not a TREAT test of a minimum effect.
The final effect-filtered subset is not guaranteed to retain a nominal FDR bound.

The workflow does not check biological independence, fit nuisance covariates,
infer valid normalization from a value range or replace experiment-specific QC.
Versions and settings are recorded so results are not detached from their method.

## Target evidence and functional follow-up

A reviewed target table is supplied by the researcher. It records the database,
evidence type and source reference for each miRNA–gene relationship. The workflow
does not identify whether a database label represents experimental support,
prediction or relevance to the current tissue.

All supplied evidence rows remain in the review export. Inclusion requires both
an explicit include decision and selection of the source miRNA. Evidence rows
are not collapsed; gene-list exports alone are deduplicated. Direction always
refers to miRNA expression. A target shared by up- and down-selected miRNAs can
appear in both downstream gene sets.

An enrichment analysis still requires an appropriate background universe,
organism, identifier namespace, database version, correction method and a
record of the exact submitted foreground. No background is inferred from selected
targets, and no pathway significance, gene-expression direction or causal
mechanism is generated here.

## Figures and provenance

The volcano plot uses BH values, clipping zeros only for display. The heatmap
uses selected rows in statistical-table order, retains metadata sample order,
and displays row z-scores without clustering. Constant rows are omitted from
that figure only. Empty selections receive an explicit skipped status rather
than a substitute result.

MD5 manifests support file comparison, not authentication or a security claim.
Input manifests identify roles, sizes and checksums without names or locations.
Software manifests contain portable source-file labels. Exported identifiers
remain research content; no output is uploaded automatically.

## Reference

[limma package documentation](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
describes the normalization, linear-model, empirical-Bayes and ranking functions.
