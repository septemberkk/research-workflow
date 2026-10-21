# Methods

## Expression and sample definition

The analysis matrix contains normalized, continuous expression with samples as
rows and stable gene identifiers as columns. Sample identifiers link expression
to trait records without relying on row order. Expression scale, feature
mapping, missing-value handling, low-expression or low-variance filtering, and
sample exclusions are part of the analysis definition.

Raw sequencing counts require suitable normalization and transformation before
correlation-network analysis. Trait-based differential-expression filtering is
not a neutral feature-selection step for an unsupervised network. Outliers,
batch structure and biological heterogeneity can influence correlations; any
adjustment has a stated rationale and can also remove relevant variation.
These considerations follow the
[authors' methodological FAQ](https://edo98811.github.io/WGCNA_official_documentation/faq.html).

## Network and modules

The network definition specifies the correlation measure, its missingness and
outlier treatment, signed or unsigned adjacency, soft-thresholding power and
topological-overlap formulation. Unsigned adjacency treats positive and
negative correlations of equal magnitude alike; signed adjacency distinguishes
their direction. Candidate powers are assessed alongside scale-free-fit and
connectivity diagnostics. A fit criterion is a selection aid, not proof of
biological validity.

Topological overlap summarizes network-neighborhood similarity. Hierarchical
clustering of the resulting dissimilarity supports module identification.
Tree-cut settings, minimum module size, block partitioning when used, and the
eigengene-based merging criterion define the reported partition. Unassigned
genes remain distinguishable from identified modules. Color labels identify
modules within an analysis, not biological identities shared across analyses.
The [WGCNA reference manual](https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf)
defines the relevant network and module functions.

## Eigengenes and trait associations

A module eigengene summarizes module expression through its first principal
component. Its orientation and explained variance qualify interpretation.
Module-trait associations use explicitly encoded traits and matched samples.
Effect estimates, contributing sample counts, missingness and the
multiple-comparison policy accompany significance measures. Covariates and
repeated observations require an appropriate statistical design; a simple
correlation does not resolve confounding or dependence.

Gene-level measures answer different questions:

| Measure | Meaning |
| --- | --- |
| Module membership (kME) | Correlation between a gene's expression and a module eigengene. |
| Intramodular connectivity | Sum of adjacency weights to other genes in the same module. |
| Gene significance (GS) | A declared measure of gene-trait association or external relevance. |

These distinctions follow the
[original WGCNA paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/).
Neither high connectivity nor high module membership alone establishes trait
association. A hub designation depends on the network and ranking definition;
it does not demonstrate causal or regulatory importance.

## Reporting and interpretation

A reporting contract links gene identifiers, module assignments, eigengenes and
association definitions to the expression preparation, parameter choices and
software versions. A Cytoscape node table preserves gene-to-module assignments;
an edge table identifies whether weights represent adjacency or topological
overlap, together with the edge-selection rule. Co-expression edges are not
physical-interaction or experimentally verified regulatory edges.

Functional enrichment and independent module-preservation assessment are
separate analyses. Module-trait associations remain exploratory evidence rather
than validated mechanisms. No fixed parameter set or study-specific findings
are supplied here.
