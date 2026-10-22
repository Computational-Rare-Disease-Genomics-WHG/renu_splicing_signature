# ReNU splicing signature

Quantify the six-event RNU4-2 n.64_65insT splicing signature from rMATS junction
counts. Reports an unweighted signature score and LASSO `P(case)` for every
case and control, using the supplied matched controls as the baseline.


## Install

Install with [uv](https://docs.astral.sh/uv/guides/tools/):

```bash
# Replace OWNER/REPOSITORY with the GitHub repository containing this package.
uv tool install --python 3.11 'git+https://github.com/OWNER/REPOSITORY.git'
renu-score --help
```

## Input: GRCh38 rMATS A5SS junction counts

**The bundled signature requires GRCh38 coordinates. GRCh37/hg19 input is
incompatible. No coordinate conversion or liftover is performed.**

Provide the tab-separated `A5SS.MATS.JC.txt` from rMATS, with both cases and
matched controls. Use the complete, unfiltered output.

Required columns (exact spelling):

| Columns | Meaning |
| --- | --- |
| `chr`, `strand` | Chromosome and strand |
| `longExonStart_0base`, `longExonEnd` | Long exon start/end |
| `shortES`, `shortEE` | Short exon start/end |
| `flankingES`, `flankingEE` | Flanking exon start/end |
| `IJC_SAMPLE_1`, `SJC_SAMPLE_1` | Group 1 (default "cases") inclusion/skipping junction counts |
| `IJC_SAMPLE_2`, `SJC_SAMPLE_2` | Group 2 (default "controls") inclusion/skipping junction counts |

The other columns are ignored.

>[!NOTE]
> rMATS outputs are 0-based half-open, as in BED format.
> Outputs with or without a 'chr' prefix are accepted.
> Splicing events are matched on `chr`, `strand`, `longExonStart_0base`, `longExonEnd`
> `shortES`, `shortEE`, `flankingES`, and `flankingEE`.
> Extra columns are allowed, but ignored.
> Counts are comma-separated nonnegative integers in a consistent sample order.
> IJC/SJC vectors must have matching lengths, and each group's size must stay
> constant between rows.
> `NA` and empty entries represent missing data. Preserve comma positions for missing
> samples. For example, `10,NA,20` represents three samples.

A complete synthetic input is included in
[examples/A5SS.MATS.JC.txt](examples/A5SS.MATS.JC.txt):

```bash
renu-score examples/A5SS.MATS.JC.txt --output scores.tsv --details events.tsv
```

## Usage

```bash
renu-score A5SS.MATS.JC.txt --output scores.tsv

# Group 1 contains controls instead of cases:
renu-score A5SS.MATS.JC.txt --controls-group 1 --output scores.tsv

# Also save the underlying measurements and event contributions:
renu-score A5SS.MATS.JC.txt --output scores.tsv --details events.tsv

# Supply sample names, one per line in the corresponding rMATS vector order:
renu-score A5SS.MATS.JC.txt --group1-names case_ids.txt --group2-names control_ids.txt
```

The scores TSV contains:

| Column | Meaning |
| --- | --- |
| `sample`, `group`, `role` | Sample name, rMATS group, case/control role |
| `uniform_sum` | Sum of oriented, control-centered log ratios |
| `weighted_sum` | LASSO-weighted sum, excluding the intercept |
| `lasso_linpred` | Weighted sum plus intercept |
| `lasso_prob` | LASSO `P(case)` |
| `n_events_used`, `n_events_expected` | Usable and expected events per sample |
| `n_zero_coverage_events` | Events with IJC=SJC=0 in this sample |
| `status` | `complete`, `partial`, or `unavailable` |

The optional `--details` TSV includes one row for every sample/event pair:
counts, log ratio, control mean, direction, weight, both contributions, and
event status. Missing values appear as `NA`.

## Configuration

Copy [examples/insT.toml](examples/insT.toml), edit it, and select it:

```bash
renu-score A5SS.MATS.JC.txt --signature my-signature.toml --output scores.tsv
```

The TOML is a complete model: schema version, name, reference genome,
intercept, and `[[events]]` entries with `id`, `strand`, `direction` (+1/-1),
and nonnegative `weight`. Event IDs encode the gene, chromosome, and six
coordinates. A different number of A5SS events is supported. Directions and
weights must be supplied for new events; this tool does not train a model.

The bundled default is `intercept = -2.927847423`. Change it in TOML;
there is no CLI intercept override. The example TOML and installed model
contain the same six Table S2 parameters and work independently of your
current directory.

Set the minimum total junction read depth in TOML, before any `[[events]]`:

```toml
min_read_depth = 0
```

An event is usable only when `IJC + SJC >= min_read_depth` for that sample.
The setting must be a nonnegative integer and defaults to `0` if omitted.
Zero-coverage events remain missing even when the threshold is `0`.
Below-threshold events are excluded only from the affected sample's score
and, for controls, from that event's baseline. They generate warnings and
appear as `below_min_read_depth` in `--details`. Partial scores retain the
original weights and intercept.

## Provenance

Based on the scoring mathematics in
[rnu42_splicing_biomarker](https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker),
commit `1d1a5dbd9a80d04f84d8187cf84231acbb1f8eee`, particularly
`score_ipsc_against_insT_signature.R` and `export_insT_signature.R`.
Coordinates, directions, and weights come from Table S2 of
[Dawes et al., A cross-tissue splicing signature as a quantitative biomarker
for ReNU syndrome](https://doi.org/10.64898/2026.09.26.26364077).
The fitted `intercept_at_1se` was supplied separately by Alexander Blakes.
The [Nava et al. supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02184-4/MediaObjects/41588_2025_2184_MOESM4_ESM.xlsx)
was used to confirm the rMATS input layout; its Table 12 contains only three
of the selected events. The included example is synthetic.

## Development

```bash
uv run --python 3.11 --extra dev pytest
uv build
```
