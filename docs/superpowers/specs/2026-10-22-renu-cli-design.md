# RNU4-2 signature CLI design

Status: requirements settled after user review, including sample-specific
zero-coverage handling.

## Intent and agreed requirements

Create a small Python CLI for quantifying the six-event RNU4-2 n.64_65insT
splicing signature in new samples. It adapts the scoring mathematics in
Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker, rather than
recreating the analysis and model training pipeline.

The user confirmed:

- Input is an rMATS output containing both cases and matched controls.
- Group 1 defaults to cases and group 2 to controls; an option can reverse this.
- Report the unweighted signature score and original LASSO P(case).
- Bundle the six selected events, with a plaintext configuration override.
- Default intercept is -2.927847423, supplied by the user; configure it in TOML
  rather than through a CLI option.
- Allow a partial signature with a warning when selected events are absent.
  Missingness is handled per sample, without dropping an event from other
  samples that have data for it.
- Use a text-only CLI, installable from GitHub using uv tool install. Standard
  Python data science dependencies are permitted.
- Generated sample names use cases_1, cases_2, etc. and controls_1, controls_2,
  etc., according to the selected control group.
- Do not add a --print-signature flag.
- Explicitly list required rMATS columns and the GRCh38 requirement in README.
- Configure a minimum total junction read depth with TOML min_read_depth,
  default 0. Exclude events below the threshold separately per sample and
  from control baselines; zero coverage remains missing at any threshold.

## Approach

CLI with python standard-library and standard python data science tools. Accepts rMATS tab-separated output, emits TSV, and
uses TOML configuration. Small installation and direct pipeline use.

## Scope and constraints

Implement scoring, configuration, input validation, useful diagnostics,
packaging, documentation, and meaningful automated checks. No fitting,
retraining, alignment, rMATS execution, cohort statistics, plots, or clinical
classification threshold is required. Input should be the complete A5SS JC
output, without filtering events by significance first.

Python >=3.11 provides tomllib. Allow standard Python data science tools. A build backend is permitted as an isolated build
dependency. Package the default model as a resource so installation does not
depend on the working directory. The local environment currently has Python
3.9 and no uv executable; acquiring a supported interpreter/tooling is an
implementation environment task.

Work is created in /home/ec2-user/renu_splicing_signature. Git was initialized
at the user's request; Git metadata operations require elevated filesystem
access in this environment. Prepare the installable package and GitHub
installation instructions; a live GitHub installation requires publication
in a repository whose destination has not yet been supplied.

## Interface and output

Builder proposals, not additional user requirements:

- Distribution: renu-splicing-signature. Executable: renu-score.
- Main usage: `renu-score A5SS.MATS.JC.txt --output scores.tsv`.
- `--control-group {1,2}` defaults to 2.
- `--config PATH` reads a complete replacement model from TOML.
- `--output PATH` writes TSV; absent this option, write TSV to stdout.
- `--details PATH` optionally writes per-sample, per-event intermediate values.
- `--group1-samples PATH` and `--group2-samples PATH` optionally supply one sample
  name per line in the corresponding rMATS vector order. Otherwise use
  cases_1, cases_2, etc. and controls_1, controls_2, etc. according to role.
  Names must be nonempty, unique, and match group sizes.
- `--version` and standard `--help` are available.

Emit one row per sample from both groups, in group/vector order, with columns
sample, group, role, uniform_sum, weighted_sum, lasso_linpred, lasso_prob,
n_events_used, n_events_expected, n_zero_coverage_events, and status.
The weighted sum excludes the intercept; lasso_linpred includes it.
Status is complete, partial, or unavailable, based on event availability. Coverage has its
own count and warning. Emit diagnostics to stderr, preserving clean TSV stdout.

The details TSV identifies the event and sample, IJC, SJC, log ratio,
control mean, direction, oriented contribution, weight, weighted contribution,
and event status. Include absent events with NA numeric fields so omissions
are visible. Neither output rounds intermediate calculations.

## Input contract

Read plain UTF-8 tab-separated A5SS.MATS.JC.txt. Accept a UTF-8 BOM. Required
columns are chr, strand, longExonStart_0base, longExonEnd, shortES, shortEE,
flankingES, flankingEE, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2,
SJC_SAMPLE_2. Other columns are ignored, including duplicate ID columns
present in ordinary rMATS output. P values and IncLevel columns do not affect
scoring.

Identify events by chromosome, strand, and all six coordinates, including the
flanking exon. Do not use rMATS numeric IDs or require gene symbols to match.
Normalize a missing chr prefix consistently between input and configuration;
otherwise use exact coordinates and preserve rMATS start/end conventions.
Bundled coordinates are GRCh38; no liftover is performed.

Counts for each matching row are comma-separated nonnegative integers.
NA and empty entries represent missing counts. A missing member of an IJC/SJC
pair makes that sample's event unavailable; both-zero counts also make it
unavailable, as confirmed by the user.
Require consistent group sizes across matching rows and paired IJC/SJC
vectors. Both groups must contain at least one sample. Reject malformed count
tokens rather than interpreting them as zero; recognize agreed missing-value
tokens separately. Reject duplicate
rows for a selected event because choosing one would be ambiguous. Invalid
configuration, missing required columns, and unusable inputs yield a concise
error and nonzero exit status, without a normal Python traceback.

## Model configuration and provenance

One TOML file contains a schema version, model name, reference genome,
intercept, and an array of events. Each event has an id in the exact
A5SS:gene:chr:six-coordinates format, strand, direction (+1 or -1), and weight
(the beta_at_lambda1se value). A custom file fully replaces the bundled model
and includes its own intercept. Intercept changes are made in TOML; there is
no CLI intercept override. Ship a visible example/default TOML in the source
repository and link to it from README so users can copy and edit it.
Permit a different number of configured A5SS events; expected event count is
the number in the selected model. Reject empty models, duplicate coordinate
identities, invalid geometry, invalid directions, and nonfinite parameters.
Weights must be nonnegative, matching the original constrained LASSO model.
Adding different events requires supplying their direction and fitted weight;
the CLI does not infer them from the new cases.

Default event parameters from the preprint's Table S2:

| Gene | Strand | Direction | Weight |
| --- | --- | --- | --- |
| CTDSP2 | - | -1 | 0.825615979283261 |
| PTBP3 | - | -1 | 0.783806558002269 |
| MAP4K4 | + | +1 | 0.737422747142889 |
| TOP1 | + | -1 | 0.560782708354878 |
| STX16 | + | +1 | 0.401105638015341 |
| RPS15A | - | +1 | 0.290736572720099 |

Coordinates, provided by the user and independently matched to Table S2:

```text
A5SS:CTDSP2:chr12:57827533:57827590:57827551:57827590:57826995:57827097
A5SS:PTBP3:chr9:112297831:112297916:112297840:112297916:112275843:112276013
A5SS:MAP4K4:chr2:101887087:101887261:101887087:101887237:101887777:101887937
A5SS:TOP1:chr20:41113969:41114155:41113969:41114137:41115370:41115439
A5SS:STX16:chr20:58667489:58667608:58667489:58667597:58667986:58668127
A5SS:RPS15A:chr16:18790243:18790334:18790309:18790334:18788980:18789118
```

Upstream commit inspected: 1d1a5dbd9a80d04f84d8187cf84231acbb1f8eee.
Relevant source scripts are scripts/score_ipsc_against_insT_signature.R and
scripts/gel/export_insT_signature.R.

Sources:

- https://github.com/Computational-Rare-Disease-Genomics-WHG/rnu42_splicing_biomarker
- https://www.medrxiv.org/content/10.64898/2026.09.26.26364077v2
- https://www.medrxiv.org/content/medrxiv/early/2026/09/29/2026.09.26.26364077/DC1/embed/media-1.xlsx?download=true
- Nava example: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02184-4/MediaObjects/41588_2025_2184_MOESM4_ESM.xlsx

## Scoring and partial signatures

For each available event j and sample i:

```text
log_ratio[i,j] = log((IJC[i,j] + 0.5) / (SJC[i,j] + 0.5))
control_mean[j] = mean(log_ratio[c,j] for controls c with data for event j)
contribution[i,j] = direction[j] * (log_ratio[i,j] - control_mean[j])
uniform_sum[i] = sum(contribution[i,j] for events j available in sample i)
weighted_sum[i] = sum(weight[j] * contribution[i,j] for events j available in sample i)
lasso_linpred[i] = intercept + weighted_sum[i]
lasso_prob[i] = sigmoid(lasso_linpred[i])
```

Use a numerically stable log ratio and sigmoid. The formula uses raw JC counts,
without effective length correction or the rounded rMATS IncLevel fields.

If a selected event has no data for an individual sample, omit it only from
that sample's sum. Retain the original weights and intercept and do not
rescale the remaining events. Warn with sample and event identities, mark
affected scores partial, and report each sample's event count. Keep the event
for other samples with usable data. A row absent from the entire rMATS file
necessarily makes the event unavailable for every sample.

Compute each event's control baseline from controls with usable data for that
event. An event with no usable controls cannot be centered and is unavailable
to score in this run; warn explicitly. A sample with no
usable events is NA scores and status unavailable. If no sample can be scored,
exit with an error. Explain in documentation that a partial-model probability
differs from applying the complete published model.

IJC=SJC=0 represents an unavailable event in that sample. Exclude it from that
sample's sums and, for a control sample, from that event's control mean.
Report zero coverage and warn. A single zero count paired with a positive
count remains usable, with the original 0.5 pseudocount. Do not pool counts
before the log transform or weight controls by depth. Require IJC+SJC to be
at least the configured nonnegative integer min_read_depth, default 0. An
omitted setting also defaults to 0 for existing custom models. Below-threshold
events have detail status below_min_read_depth and generate warnings. A sample
exactly at the threshold is usable when it has positive coverage.

## Packaging and documentation

Use a src-layout Python package with a console entry point and bundled TOML.
The README explains installation from a GitHub URL using uv tool install,
usage, group conventions, sample order/naming, copying/editing the model TOML,
the exact formulas, partial and zero-coverage behavior, and source attribution.
List every required rMATS column explicitly. State prominently that bundled
event coordinates require GRCh38 input, that GRCh37 input is incompatible,
and that no liftover is performed. Include a small synthetic complete rMATS example so
users can see the file and expected output without downloading patient data.
Do not invent a destination repository URL or claim it has been uploaded.

Document the user-provided provenance of the intercept separately from the
published Table S2 weights. Implement the mathematics independently with
attribution; do not assume an upstream license that is not present.

## Validation and completion criteria

Meaningful automated checks cover independently calculated numeric scores,
per-control log centering, both control-group assignments, partial and empty
signatures, zero coverage, duplicate/ambiguous matches, full coordinate and
strand matching, malformed count vectors, custom models and TOML intercepts,
sample names, and clean CLI TSV output. Use only synthetic examples plus a
small, attributed subset of the public Nava table if useful; the Nava A5SS
table has only three selected events and is not a complete-model reference.

Build the package and install it into an isolated uv tool directory, then
invoke the installed executable outside the source checkout to demonstrate
that bundled configuration is included and readable. Verify source install
from a local Git URL if a real Git checkout can be created in a writable
temporary location. Live GitHub installation remains contingent on publishing
the completed source. Report checks actually run and any remaining limitation.

## Decisions left to the builder

Choose the internal module split, isolated build backend, test runner,
formatting of floating-point TSV values without intermediate rounding, and
wording of diagnostics. Keep the user-facing interface small and follow the
decisions above.

## Minimum coverage follow-up resolved

The user requested a configurable minimum read depth, default 0. It applies
to total IJC+SJC per event and sample, including controls used for the baseline.
