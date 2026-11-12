# ReNU splicing signature

A command-line tool accompanying [Dawes et al., A cross-tissue splicing
signature as a quantitative biomarker for ReNU
syndrome](https://doi.org/10.64898/2026.09.26.26364077).

Quantifies the eleven-event *RNU4-2* n.64_65insT splicing signature or the
provisional seven-event Stem III signature from
[`rMATS`](https://doi.org/10.1038/s41596-023-00944-2) junction counts.

## Install

Install with [`uv`](https://docs.astral.sh/uv/guides/tools/):

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

Requires the unfiltered `A5SS.MATS.JC.txt` output from `rMATS`.

A stub example input is included in
[`examples/A5SS.MATS.JC.txt`](examples/A5SS.MATS.JC.txt).

> [!IMPORTANT]
> [Dawes et al.](https://doi.org/10.64898/2026.09.26.26364077) defined the splicing
> signature with `rMATS-turbo` (`v4.3.0`) with the following arguments:
> `--libType fr-firststrand`, `--anchorLength 1`, and `--novelSS`. The signature was developed
> on the GENCODE v49 and GENCODE v46 annotations. We suggest that users match these
> settings to enable direct comparisons between samples.

> [!NOTE]
> - `renu-score` requires the `rMATS` output to contain two “sample groups”
> (run `rMATS` with `--s1` and `--s2` or `--b1` and `--b2`).
> - We suggest running `renu-score` with a bare minimum of two controls in the `rMATS` output. The more cases and controls available, the better.
> - This tool requires GRCh38 coordinates. GRCh37/hg19 input is incompatible.

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
| `sample`, `group`, `role`            | Sample name, `rMATS` group, case/control role                         |
| `uniform_sum`                        | Sum of oriented, control-centered log ratios                          |
| `weighted_sum`                       | LASSO-weighted sum, excluding the intercept                           |
| `lasso_linpred`                      | Weighted sum plus intercept                                           |
| `lasso_prob`                         | LASSO `P(case)`                                                       |
| `n_events_used`, `n_events_expected` | Usable and expected events per sample                                 |
| `n_zero_coverage_events`             | Events with `IJC=SJC=0` in this sample                                |
| `status`                             | Status of splicing signature: `complete`, `partial`, or `unavailable` |

The optional `--details` output offers more granular information for every
sample/event pair.

## Configuration

This tool does not train a new splicing signature model.

The eleven splicing events and fitted model parameters for the n.64_65insT signature are bundled with the tool in
[`config/insT.toml`](config/insT.toml).

Use `--stem3` to test the provisional seven-event Stem III signature, bundled in
[`config/stem3.toml`](config/stem3.toml).

Commonly, not every splicing event is available for each sample. In this case, the tool
uses custom weights calculated for every possible combination of events.
These custom weights are bundled here:

- [insT: 2,047 combinations of 11 unique events](src/renu_splicing_signature/data/subset_refit_weights_insT.csv)
- [Stem III: 127 combinations of 7 unique events](src/renu_splicing_signature/data/subset_refit_weights_stemIII.csv)

## References

- [Dawes et al., A cross-tissue splicing signature as a quantitative biomarker
  for ReNU syndrome](https://doi.org/10.64898/2026.09.26.26364077).
- [Companion analysis pipeline](https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker).
