# ReNU splicing signature

A command-line tool accompanying [Dawes et al., A cross-tissue splicing
signature as a quantitative biomarker for ReNU
syndrome](https://doi.org/10.64898/2026.09.26.26364077).

Quantifies the six-event RNU4-2 n.64_65insT splicing signature from
[rMATS](https://doi.org/10.1038/s41596-023-00944-2) junction counts. Reports an
unweighted signature score and LASSO `P(case)` for every case and control,
using the supplied matched controls as a baseline.

## Install

Install with [uv](https://docs.astral.sh/uv/guides/tools/):

```bash
uv tool install --python 3.11 'git+https://github.com/Computational-Rare-Disease-Genomics-WHG/renu_splicing_signature.git'
renu-score --help
```

## Input: GRCh38 rMATS A5SS junction counts

**This tool requires GRCh38 coordinates. GRCh37/hg19 input is
incompatible. No coordinate conversion or liftover is performed.**

Provide the tab-separated `A5SS.MATS.JC.txt` from rMATS, with both cases and
matched controls. Use the complete, unfiltered output.

Required columns (exact spelling):

| Columns                              | Meaning                                                          |
| ------------------------------------ | ---------------------------------------------------------------- |
| `chr`, `strand`                      | Chromosome and strand                                            |
| `longExonStart_0base`, `longExonEnd` | Long exon start/end                                              |
| `shortES`, `shortEE`                 | Short exon start/end                                             |
| `flankingES`, `flankingEE`           | Flanking exon start/end                                          |
| `IJC_SAMPLE_1`, `SJC_SAMPLE_1`       | Group 1 inclusion/skipping junction counts (cases by default)    |
| `IJC_SAMPLE_2`, `SJC_SAMPLE_2`       | Group 2 inclusion/skipping junction counts (controls by default) |

The other columns are ignored.

> [!NOTE]
>
> - rMATS uses zero-based, half-open coordinates, as in BED format.
>   Contigs with or without a `chr` prefix are accepted. Events are matched on
>   chromosome, strand, and all six coordinate columns above.
> - Counts are comma-separated integers in a consistent sample order. Zero
>   counts are allowed. IJC/SJC vectors must have matching lengths, and each
>   group's size must be constant across rows.
> - `NA` and empty entries represent missing data.

A stub example input is included in
[examples/A5SS.MATS.JC.txt](examples/A5SS.MATS.JC.txt):

```bash
renu-score examples/A5SS.MATS.JC.txt --output scores.tsv --details events.tsv
```

## Usage

```bash
renu-score A5SS.MATS.JC.txt --output scores.tsv

# Group 1 contains controls instead of cases:
renu-score A5SS.MATS.JC.txt --control-group 1 --output scores.tsv

# Also save the underlying measurements and event contributions:
renu-score A5SS.MATS.JC.txt --output scores.tsv --details events.tsv

# Supply sample names, one per line in the corresponding rMATS vector order:
renu-score A5SS.MATS.JC.txt --group1-samples case_ids.txt --group2-samples control_ids.txt
```

The scores TSV contains:

| Column                               | Meaning                                      |
| ------------------------------------ | -------------------------------------------- |
| `sample`, `group`, `role`            | Sample name, rMATS group, case/control role  |
| `uniform_sum`                        | Sum of oriented, control-centered log ratios |
| `weighted_sum`                       | LASSO-weighted sum, excluding the intercept  |
| `lasso_linpred`                      | Weighted sum plus intercept                  |
| `lasso_prob`                         | LASSO `P(case)`                              |
| `n_events_used`, `n_events_expected` | Usable and expected events per sample        |
| `n_zero_coverage_events`             | Events with IJC=SJC=0 in this sample         |
| `status`                             | `complete`, `partial`, or `unavailable`      |

The optional `--details` TSV includes one row for every sample/event pair:
counts, log ratio, control mean, direction, weight, both contributions, and
event status. Missing values appear as `NA`.

## Configuration

The six splice events and fitted model parameters are bundled with the tool.
[examples/insT.toml](examples/insT.toml) is an editable copy of the default
configuration. This tool does not train a new splicing signature model.

To change the configuration, copy the example TOML, edit it, and select it:

```bash
renu-score A5SS.MATS.JC.txt --config my-config.toml --output scores.tsv
```

## References

- [Companion analysis pipeline](https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker).
- [Dawes et al., A cross-tissue splicing signature as a quantitative biomarker
  for ReNU syndrome](https://doi.org/10.64898/2026.09.26.26364077).
