# miRNA Differential Expression & Target Enrichment

An independent R workflow for two-group limma analysis of **log2-scale miRNA
expression**, with optional preparation of reviewed target-evidence exports.
Only invented example identifiers and values are bundled.

```text
log2 expression + explicit sample groups
    -> optional duplicate-ID averaging
    -> quantile normalization, or an explicit no-normalization mode
    -> limma comparison - reference
    -> complete statistics and screened miRNA lists
    -> optional reviewed target evidence and gene-list handoff
```

Target databases are not queried. Identifier conversion, target prediction and
GO/KEGG calculations are separate steps. The workflow does not implement raw
count analysis, array background correction, log transformation, paired designs,
batch adjustment, repeated measures or a time-course model.

## Runtime and execution

R 4.1 or newer and an existing installation of limma are required. Optional
figures use base R. The workflow never installs packages or accesses the network.

```sh
Rscript run_all.R --check-environment
Rscript run_all.R --example --output OUTPUT --plots
Rscript run_all.R --expression EXPRESSION --samples SAMPLES --reference CONTROL --comparison TREATMENT --input-scale log2 --output OUTPUT
```

Uppercase arguments are caller-supplied values. Output must be a new directory
with an existing parent. Input locations are runtime arguments, not embedded
defaults or entries in exported metadata.

| Option | Meaning |
| --- | --- |
| `--normalization quantile` | Quantile-normalize the supplied log2 matrix; the default. |
| `--normalization none` | Preserve already normalized log2 expression. |
| `--duplicates error` | Reject repeated miRNA IDs; the default. |
| `--duplicates mean` | Average repeated IDs arithmetically on the supplied log2 scale. |
| `--fdr`, `--logfc` | Candidate screen: BH value strictly below 0.05 and absolute log2FC strictly above 1 by default. |
| `--targets TARGETS` | An optional reviewed target-evidence TSV. |
| `--plots`, `--heatmap-rows` | Volcano plot and row-z-score heatmap, limited to 40 selected rows by default. |

## Input contract

All tables are UTF-8 TSV with headers, unquoted fields, consistent column counts
and no blank physical lines. Identifiers are preserved as case-sensitive strings.
Blank IDs, boundary whitespace, formula-like prefixes and invisible formatting
characters are rejected; aliases are not repaired or inferred.

- Expression: `miRNA_id` followed by numeric sample columns. Values must already
  be on the log2 expression scale; declaring the scale is the caller's responsibility.
- Samples: exactly `sample_id` and `group`. Every expression column must match
  one unique metadata row. Exactly the two named groups and at least two
  independent biological samples per group are required. Metadata order determines
  output sample order, not the original expression-column order.
- Reviewed targets: exactly `miRNA_id`, `target_gene_id`, `database`,
  `evidence_type`, `reference` and `review`. Fields must be nonblank. Review
  is explicitly `include` or `exclude`; every miRNA must occur in the assayed matrix.
  Adapt a database export to this contract without discarding its evidence provenance.

For separate experimental contrasts, prepare an explicit two-group input for
each run. Sharing a reference group between runs does not turn the outputs into
a longitudinal test or an interaction analysis.

## Outputs

| Output | Contents |
| --- | --- |
| `normalized_expression.tsv`, `samples_used.tsv`, `design.tsv` | Actual analyzed matrix, sample order and two-column group design. |
| `differential_all.tsv` | All unique tested miRNAs, limma statistics, selection flag and direction. |
| `differential_selected.tsv`, `miRNAs_up.txt`, `miRNAs_down.txt` | Candidate-screen results; empty selections remain valid. |
| `target_review.tsv` | Every supplied evidence row and an explicit inclusion/exclusion reason. |
| `target_evidence_selected.tsv`, `target_genes_all.txt` | Included evidence for selected miRNAs and unique target IDs. |
| `targets_of_up_miRNAs.txt`, `targets_of_down_miRNAs.txt` | Target sets grouped by the **source miRNA's** direction, not inferred target-gene expression. A gene can appear in both sets. |
| `volcano.pdf`, `heatmap.pdf`, `plot_status.tsv` | Optional figures and explicit status for skipped or unrequested plots. |
| `input_manifest.tsv`, `source_manifest.tsv`, `output_manifest.tsv`, `run_metadata.tsv` | Checksums, settings, group sizes, software versions and completion state, without input locations. |

Target files are created only when reviewed evidence is supplied. Gene lists are
not enrichment results and do not define an enrichment background. Database
support is not proof of a mechanism in the current experiment. Outputs retain
supplied identifiers and must be reviewed before sharing.

See [methods and interpretation](METHODS.md) and [synthetic inputs](examples/README.md).
Licensed under the [MIT License](LICENSE).
