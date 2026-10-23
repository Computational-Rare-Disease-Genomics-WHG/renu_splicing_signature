# ReNU splicing signature

A command-line tool accompanying [Dawes et al., A cross-tissue splicing
signature as a quantitative biomarker for ReNU
syndrome](https://doi.org/10.64898/2026.09.26.26364077).

Quantifies the six-event RNU4-2 n.64_65insT splicing signature or the provisional
six-event Stem III signature from
[rMATS](https://doi.org/10.1038/s41596-023-00944-2) junction counts.

## Install

Install with [uv](https://docs.astral.sh/uv/guides/tools/):

```bash
uv tool install --python '>=3.11' 'git+https://github.com/Computational-Rare-Disease-Genomics-WHG/renu_splicing_signature.git'

renu-score --help
```

## Input

Requires the unfiltered `A5SS.MATS.JC.txt` output from rMATS.

A stub example input is included in
[examples/A5SS.MATS.JC.txt](examples/A5SS.MATS.JC.txt):

> [!NOTE]
> - This tool requires GRCh38 coordinates. GRCh37/hg19 input is
> incompatible. No coordinate conversion or liftover is performed.
> - This tool accepts contigs with or without a `chr` prefix.

## Usage

Run on the bundled example without providing an input file:

```bash
# Test on the bundled example data
renu-score --example

# Provide the rMATS data for your samples; write output to a file
renu-score A5SS.MATS.JC.txt --output scores.tsv

# Use the provisional Stem III signature:
renu-score A5SS.MATS.JC.txt --stem3 --output stem3-scores.tsv

# rMATS group 1 contains controls instead of cases:
renu-score A5SS.MATS.JC.txt --control-group 1 --output scores.tsv

# Also report the underlying event contributions:
renu-score A5SS.MATS.JC.txt --output scores.tsv --details events.tsv

# Supply sample names, one per line. The order should match the order of samples in the rMATS input:
renu-score A5SS.MATS.JC.txt --group1-samples case_ids.txt --group2-samples control_ids.txt
```

## Output
The default output contains:

| Column                               | Meaning                                      |
| ------------------------------------ | -------------------------------------------- |
| `sample`, `group`, `role`            | Sample name, rMATS group, case/control role  |
| `uniform_sum`                        | Sum of oriented, control-centered log ratios |
| `weighted_sum`                       | LASSO-weighted sum, excluding the intercept  |
| `lasso_linpred`                      | Weighted sum plus intercept                  |
| `lasso_prob`                         | LASSO `P(case)`                              |
| `n_events_used`, `n_events_expected` | Usable and expected events per sample        |
| `n_zero_coverage_events`             | Events with IJC=SJC=0 in this sample         |
| `status`                             | Status of splicing signature: `complete`, `partial`, or `unavailable`      |

The optional `--details` output offers more granular information for every sample/event pair.

## Configuration

This tool does not train a new splicing signature model.

The six splicing events and fitted model parameters are bundled with the tool.
[examples/insT.toml](examples/insT.toml) is an editable copy of the default
configuration.

Use `--stem3` to select the separate provisional Stem III signature, with the
six events and fitted weights from supplementary Table S3 and its own intercept.
[examples/stem3.toml](examples/stem3.toml) is an editable copy. `--stem3` and
`--config` are mutually exclusive; to customize Stem III, edit this copy and pass
it with `--config`. The synthetic example input includes all events for both
signatures, so `--example --stem3` produces complete Stem III scores.

To change the configuration, copy and edit the example TOML, and use it with the `--config` CLI option:

```bash
renu-score A5SS.MATS.JC.txt --config my-config.toml --output scores.tsv
```

## References

- [Companion analysis pipeline](https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker).
- [Dawes et al., A cross-tissue splicing signature as a quantitative biomarker
  for ReNU syndrome](https://doi.org/10.64898/2026.09.26.26364077).
