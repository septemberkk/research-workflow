# Research Workflows

Independent research-analysis workflows and method notes with explicit inputs,
methods and scope.

| Workflow | Purpose | Scope |
| --- | --- | --- |
| [RNA-seq Differential Expression & Functional Enrichment](workflows/rna-seq-differential-expression-and-functional-enrichment/) | Two-group edgeR differential expression, an optional heatmap, and gene lists for downstream enrichment. | Independent biological samples in two groups; identifier mapping and GO/KEGG analysis are separate steps. |
| [miRNA Differential Expression & Target Enrichment](workflows/mirna-differential-expression-and-target-enrichment/) | Log2-expression limma analysis, optional figures, reviewed target evidence and gene-list handoff. | Independent two-group comparisons; target queries, identifier conversion and GO/KEGG calculations are separate steps. |
| [Microarray Preprocessing & Parallel Expression Analysis](workflows/microarray-preprocessing-and-parallel-expression-analysis/) | Explicit sample selection, probe mapping, arithmetic feature consolidation and branch-specific normalization. | Independent mRNA and miRNA preparation; no differential-expression model, cross-assay pairing or target inference. |
| [Cross-Study Integration Workflow](workflows/cross-study-integration-workflow/) | Shared-feature alignment, study/group models, ComBat adjustment and descriptive diagnostics. | Methods-only documentation; no executable pipeline, study data or biological results. |
| [Single-Cell RNA-seq Analysis Workflow](workflows/single-cell-rnaseq-analysis-workflow/) | Count provenance, quality review, Seurat-based exploration, reviewed annotation and sample-level composition. | Methods-only documentation; no study data, runnable analysis package or biological findings. |
| [Functional Module Evidence Mapping](workflows/functional-module-evidence-mapping/) | Reviewed term-to-module mappings, traceable gene annotations, heatmap matrices and Cytoscape-ready node tables. | Uses reviewed enrichment exports; it does not discover modules, run enrichment analysis or build Cytoscape networks. |

Runnable workflows include their own source, input contracts and synthetic
examples. Methods-only entries describe an analysis approach without executable
code or a study-reproduction claim. Each entry states its scope in its README.

Only synthetic example data are included. Analysis exports retain supplied
identifiers and should be reviewed before sharing.

Browser utilities are shared separately in [Research Tools](https://github.com/septemberkk/research-tools-public).

Licensed under the [MIT License](LICENSE).
