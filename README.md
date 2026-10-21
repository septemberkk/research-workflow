# Research Workflows

Independent research-analysis workflows with explicit inputs, methods and outputs.

| Workflow | Purpose | Scope |
| --- | --- | --- |
| [RNA-seq Differential Expression & Functional Enrichment](workflows/rna-seq-differential-expression-and-functional-enrichment/) | Two-group edgeR differential expression, an optional heatmap, and gene lists for downstream enrichment. | Independent biological samples in two groups; identifier mapping and GO/KEGG analysis are separate steps. |
| [Functional Module Evidence Mapping](workflows/functional-module-evidence-mapping/) | Reviewed term-to-module mappings, traceable gene annotations, heatmap matrices and Cytoscape-ready node tables. | Uses reviewed enrichment exports; it does not discover modules, run enrichment analysis or build Cytoscape networks. |

Each workflow contains its own source, method notes, input contract and synthetic
examples. Follow its README for runtime requirements and usage.

Only synthetic example data are included. Analysis exports retain supplied
identifiers and should be reviewed before sharing.

Browser utilities are shared separately in [Research Tools](https://github.com/septemberkk/research-tools-public).

Licensed under the [MIT License](LICENSE).
