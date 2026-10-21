# Method notes

## Inputs and alignment

Each study contributes a separately prepared, continuous gene-level expression
matrix. Feature identifiers, annotation provenance, species and measurement
scales must be compatible across studies; matching text alone does not establish
that compatibility. Probe consolidation is part of within-study preparation.

The common matrix uses the exact feature-key intersection, joined by identifier
rather than row position. Study-specific rows are outside this intersection,
not filled with zeros. Feature-retention counts describe what the join excludes.

An explicit sample table associates every selected column with a sample ID,
study and biological group. Group scope and exclusions are documented separately
from column order. Sample columns and their metadata refer to the same ordered
samples.

## Batch model

The documented model uses study origin as the batch variable, includes the
biological group in the covariate matrix and uses parametric priors in
`sva::ComBat`. Its input is cleaned, normalized continuous array expression,
not raw RNA-seq counts. These roles are described in the
[sva reference manual](https://bioconductor.org/packages/release/bioc/manuals/sva/man/sva.pdf).

Batch adjustment cannot separate completely confounded study and group effects.
The [ComBat implementation](https://raw.githubusercontent.com/bioc/sva/devel/R/ComBat.R)
checks design rank for confounding. Including a group term describes the model;
it does not establish that every biological signal is preserved or that
incompatible studies become interchangeable.

## Descriptive diagnostics

Before/after values for selected features describe the supplied and adjusted
matrices without changing their declared scale.

Sample PCA uses samples as observations and shared genes as variables, with
gene centering and no unit-variance scaling. Each stage has its own fitted basis
and explained-variance labels. Separately fitted coordinates are not a
before/after trajectory. The
[R PCA reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html)
defines the centering, scaling and variance quantities.

Full-feature sample correlations can provide a complementary descriptive view,
with the same sample order and a common correlation scale. Neither stronger
correlations nor apparent PCA separation establishes successful biological
correction.

## Interpretation boundary

Aligned matrices, the sample model, adjusted values and stage-specific
diagnostics are distinct analytical materials. This documentation does not
supply computed versions of them.

Differential-expression coefficients, significance thresholds, enrichment
results and regulatory or causal conclusions are outside this resource.
Batch adjustment alone neither establishes a log2 measurement scale nor
justifies a log2 fold-change interpretation.
