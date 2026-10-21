# Single-Cell RNA-seq Analysis Workflow

Method notes on a Seurat-based exploratory workflow, from supplied gene-count
matrices through quality review, embedding, annotation and descriptive
cell-type composition.

This is a methods-only resource, not an executable analysis package or a
packaged reproduction of a study. No study data, cell-level exports, embeddings,
annotation results or biological findings are bundled.

## Method outline

Count matrices and sample provenance → quality review → normalization and PCA →
neighbors, clustering and UMAP → reviewed annotations → sample-level summaries.

[METHODS.md](METHODS.md) distinguishes the exploratory stages from quality-control
completeness requirements and sample-level inference. Thresholds, reference
choices and cluster labels are dataset-specific rather than transferable defaults.

Licensed under the [MIT License](LICENSE).
