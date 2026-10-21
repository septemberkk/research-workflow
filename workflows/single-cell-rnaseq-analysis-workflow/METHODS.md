# Method notes

## Inputs and observational units

The starting material is a sparse gene-by-cell count matrix for each capture,
with stable feature keys and cell barcodes. Species, tissue, assay, chemistry
and annotation provenance belong to the input record.

Cell identity combines capture identity with the original barcode. Biological
sample, capture channel, processing batch and comparison group are separate
metadata concepts; one must not be inferred from another. Pooling or technical
replication does not create additional independent biological samples.

Raw counts remain distinct from normalized expression and scaled feature values.

## Quality review

Detected genes, total counts and mitochondrial fraction are inspected by
capture as well as overall. Thresholds depend on observed distributions and
biology; no universal numerical cutoffs are supplied here. Mitochondrial
fractions require an appropriate annotation, not a silent zero when genes
cannot be identified. Retention summaries retain denominators and exclusion
reasons.

Simple metric filters do not establish singlet status or remove ambient RNA.
Full QC requires capture-aware doublet assessment, with the score, class and
partition retained; [scDblFinder](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)
documents these quantities. Ambient-RNA assessment depends on the assay and
available inputs. Missing assessments remain unassessed, not clean.

Doublet and ambient-RNA assessment are completeness requirements, not completed
analyses reported by this methods-only resource.

## Normalization, embedding and clustering

The described exploratory sequence uses Seurat LogNormalize, variable-feature
selection, scaling, PCA, neighbor construction, clustering and UMAP.
[LogNormalize](https://satijalab.org/seurat/reference/normalizedata)
uses library-size scaling followed by a natural-log transformation; it does
not turn the raw count assay into log2 fold changes.

Feature selection, dimensionality, clustering resolution and random seeds are
analysis records rather than fixed presets. Merging objects is not batch
correction. Batch correction or latent integration requires a separately
justified model; none is claimed here.

Split views use the same parent embedding. Independently recomputed embeddings
do not share coordinates. UMAP geometry is descriptive, not a significance
test, trajectory or measure of treatment effect.

## Annotation

Computational cluster IDs are separate from biological labels. Marker patterns
and a compatible reference support reviewed, coarse-to-fine annotation;
ambiguous groups retain unknown or mixed labels. Evidence, confidence and
manual exclusions remain attached to the annotation record.

A cluster-marker comparison supports exploration, not automatic validation
of a cell type, subtype or experimental effect. No study-specific marker
panel or cluster-to-label assignment is supplied.

## Composition and inference

Composition summaries report cell counts and the retained-cell denominator
within each biological sample, including valid zero counts for absent
categories. A pooled condition proportion is cell-weighted, not an estimate
based on independent biological replicates.

Cells from the same biological sample are not independent experimental
replicates. [Seurat's sample-level comparison guidance](https://satijalab.org/seurat/articles/de_vignette.html)
distinguishes cell-level comparisons from replicate-aware inference.
Condition-level differential expression, differential abundance, enrichment
and causal or regulatory conclusions are outside this resource.
