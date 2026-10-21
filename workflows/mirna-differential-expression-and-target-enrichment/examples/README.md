# Synthetic example inputs

All values and identifiers are invented. `MIR_EXAMPLE_*` and `GENE_EXAMPLE_*`
are illustrative labels, not biological entities or database-supported findings.

`expression.tsv` contains deterministic log2-scale values for thirty features
and six independent-sample placeholders. `samples.tsv` defines Control and
Treatment. No random sampling is used.

`reviewed_targets.tsv` illustrates included and excluded evidence records and
shared targets. Its synthetic database, evidence and reference fields make no
claim about real interactions. The example run includes this table; callers
with real data supply their own reviewed exports.
