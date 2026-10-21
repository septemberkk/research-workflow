# Cross-Study Integration Workflow

Method notes on combining independently prepared gene-level microarray matrices,
aligning shared features and inspecting study-aware batch adjustment.

This is a methods-only resource, not an executable analysis pipeline or a
packaged reproduction of a study. It contains no analysis scripts, expression
matrices, sample records or biological results.

## Method outline

Prepared study matrices → shared-feature alignment → sample, study and group
model → ComBat adjustment → separate before/after diagnostics.

The scope starts with supplied gene-level array expression. Probe annotation
and within-study preparation come before integration; differential expression
and enrichment are separate downstream analyses.

[METHODS.md](METHODS.md) describes the input assumptions, model and interpretation
boundaries.

Licensed under the [MIT License](LICENSE).
