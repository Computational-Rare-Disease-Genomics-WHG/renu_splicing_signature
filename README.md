# ReNU splicing signature

A command-line tool accompanying [Dawes et al., A cross-tissue splicing
signature as a quantitative biomarker for ReNU
syndrome](https://doi.org/10.64898/2026.09.26.26364077).

Quantifies the six-event RNU4-2 n.64_65insT splicing signature from
[rMATS](https://doi.org/10.1038/s41596-023-00944-2) junction counts.

## Install

Install with [uv](https://docs.astral.sh/uv/guides/tools/):

```bash
uv tool install --python '>=3.11' 'git+https://github.com/Computational-Rare-Disease-Genomics-WHG/renu_splicing_signature.git'
```

Test the installation:

```bash
renu-score --help
```

Test run on bundled example data:

```bash
renu-score --example
```

## Input

Requires the unfiltered `A5SS.MATS.JC.txt` output from rMATS.

A stub example input is included in
[examples/A5SS.MATS.JC.txt](examples/A5SS.MATS.JC.txt).

> [!IMPORTANT]
> [Dawes et al.](https://doi.org/10.64898/2026.09.26.26364077) defined the splicing
> signature with rMATS-turbo (v4.3.0) with the following arguments:
> --libType fr-firststrand, --anchorLength 1, and --novelSS. The signature was developed
> on the GENCODE v49 and GENCODE v46 annotations. We suggest that users match these
> settings to enable direct comparisons between samples.

> [!NOTE]
> - This tool requires GRCh38 coordinates. GRCh37/hg19 input is
> incompatible.

> [!NOTE]
> The tool takes no account of rMATS statistical outputs, but does
> require the output from two rMATS “sample groups”
> (`--s1` and `--s2` or `--b1` and `--b2`). **The minimum
> number of samples required is therefore one “case” and one “control”.**
> The "case" and "control" labels are arbitrary - a ReNU syndrome
> sample could be included among the "controls", for example.

## Usage

```bash
# Provide the rMATS input data; write output to a file
renu-score A5SS.MATS.JC.txt --output scores.tsv

# Use the provisional Stem III signature
renu-score A5SS.MATS.JC.txt --stem3 --output stem3-scores.tsv

# rMATS group 1 contains controls instead of cases
renu-score A5SS.MATS.JC.txt --control-group 1 --output scores.tsv

# Also report the underlying event contributions
renu-score A5SS.MATS.JC.txt --output scores.tsv --details events.tsv

# Supply sample names, one per line. The order should match the order of samples in the rMATS input
renu-score A5SS.MATS.JC.txt --group1-samples case_ids.txt --group2-samples control_ids.txt
```

## Output

The default output contains:

| Column                               | Meaning                                                               |
| ------------------------------------ | --------------------------------------------------------------------- |
| `sample`, `group`, `role`            | Sample name, rMATS group, case/control role                           |
| `uniform_sum`                        | Sum of oriented, control-centered log ratios                          |
| `weighted_sum`                       | LASSO-weighted sum, excluding the intercept                           |
| `lasso_linpred`                      | Weighted sum plus intercept                                           |
| `lasso_prob`                         | LASSO `P(case)`                                                       |
| `n_events_used`, `n_events_expected` | Usable and expected events per sample                                 |
| `n_zero_coverage_events`             | Events with IJC=SJC=0 in this sample                                  |
| `status`                             | Status of splicing signature: `complete`, `partial`, or `unavailable` |

The optional `--details` output offers more granular information for every sample/event pair.

| Column | Meaning |
| ------ | ------- |
| `logit_psi` | Natural log of `(IJC + 0.5) / (SJC + 0.5)`; same as `log_ratio` |
| `delta_logit_psi` | `logit_psi` minus the mean logit across usable controls (`control_mean`) |
| `delta_logit_psi_weighted` | `delta_logit_psi` multiplied by the event's model weight |
| `delta_logit_psi_direction` | `direction` multiplied by `delta_logit_psi`; same as `contribution` |
| `delta_logit_psi_direction_weighted` | `delta_logit_psi_direction` multiplied by the event's model weight |

These statistics use pseudocount-adjusted junction counts rather than rMATS's
reported PSI. The `delta_logit_psi` and `delta_logit_psi_weighted` columns retain
the inclusion/skipping direction. The two direction-adjusted columns apply the
signature's `direction`, matching `contribution` and `weighted_contribution`.
Unusable measurements are reported as `NA`.
If an event has no usable controls, its logit is still reported for usable
samples, but all delta statistics are `NA`.

## Configuration

This tool does not train a new splicing signature model.

The six splicing events and fitted model parameters for the n.64_65insT signature are bundled with the tool in
[config/insT.toml](config/insT.toml).

Use `--stem3` to test the provisional Stem III signature, bundled in
[config/stem3.toml](config/stem3.toml).

To customise the signature, copy and edit the example TOML and use the `--config` CLI option:

```bash
renu-score A5SS.MATS.JC.txt --config my-config.toml --output scores.tsv
```

## References

- [Companion analysis pipeline](https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker).
- [Dawes et al., A cross-tissue splicing signature as a quantitative biomarker
  for ReNU syndrome](https://doi.org/10.64898/2026.09.26.26364077).
