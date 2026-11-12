# Subset false-positive checks implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Bundle subset false-positive statistics, warn on nonzero false positives, and suppress high-FPR scores without a failing exit status.

**Architecture:** Extend each `SubsetFit` with `n_fp: int` and `fpr_pct: float` from the source exports. Check the selected fit independently per sample after availability and centering, before score arithmetic. Emit one sample-specific warning identifying all usable junctions when `n_fp > 0` or `fpr_pct >= 5`.

**Tech Stack:** Existing Python standard-library package, pytest, hatchling, ruff.

**Spec:** User follow-up: bundle `n_fp`/`fpr_pct`; warn with sample and usable junctions for `n_fp > 0`; do not provide scores at `fpr_pct >= 5`. Latest clarification: high-FPR samples only warn and must not cause a nonzero exit status.

## Global constraints

- Preserve coefficient/intercept precision and exact subset membership, including zero weights.
- `fpr_pct` is already percentage units: compare directly to 5, inclusively.
- High-FPR rows retain identity and usable event counts, set all four scores to `NA`, and use status `high_fpr`.
- Retain raw count/log-ratio/baseline diagnostics; suppress signature contribution fields for rejected samples and mark their usable detail rows `high_fpr`.
- Other samples and control centering remain independent of the FPR guard. Empty samples and fixed custom models retain existing behavior.
- Continue with exit code 0 for high FPR, including if all samples have high FPR.

## Review focus

- Exact 5% and values just below it; percentage units must not be confused with proportions.
- Cases and controls with different availability; reject only the offending subset.
- Usable zero-weight events appear in warning junction lists.
- All-rejected runs still produce identifiable NA rows and exit 0.
- Validate metadata and preserve it in installed packages without changing model weights.

### Task 1: Bundle and load metadata

**Files:** Both bundled subset CSVs, `src/renu_splicing_signature/subsets.py`, `tests/test_fpr.py`.
**Interfaces:** Add `SubsetFit.n_fp: int` and `SubsetFit.fpr_pct: float`; existing mappings and loading signatures remain unchanged.

- [x] Write failing tests for real subset statistics and malformed metadata.
- [x] Run `.venv/bin/python -m pytest tests/test_fpr.py -q`; observe missing metadata failures.
- [x] Add the two original columns verbatim; validate nonnegative integer counts and finite percentages in [0, 100].
- [x] Rerun targeted tests and independently compare every retained CSV cell to source exports.

### Task 2: Warn and suppress scores per sample

**Files:** `scoring.py`, `tests/test_fpr.py`, `tests/test_subsets.py`, `tests/test_cli.py`, `tests/test_distribution.py`, `README.md`.
**Interfaces:** Keep `ScoreResult` and CLI flags unchanged; use its warning list and introduce `high_fpr` row/detail status.

- [x] Add failing boundary, mixed-case/control, zero-coefficient, all-rejected CLI, detail-suppression, and custom-model tests. Extend all-subset expectations with the new rejection rule.
- [x] Run targeted tests; observe absent warnings and unexpected scores.
- [x] Select the fit before weighting; warn once with name, gene names for usable junctions and statistics; suppress score/contribution fields at >=5%, preserving other diagnostics and samples.
- [x] Verify installed-wheel warning/rejection behavior and update concise README guidance.
- [x] Run full pytest suite, ruff lint/format, and git diff checks; request one independent review and resolve material findings.

## Verification result

182 tests passed, including installed-wheel warning and high-FPR runs for both
signatures. Every retained cell in all 2,174 compact rows matches the original
exports. Ruff lint/format and git diff checks pass. Independent review found
no critical, important, or minor issues. Changes remain uncommitted.

## Warning wording follow-up

Warnings lead with whether scores are provided, explain false positives as
controls incorrectly flagged as cases during model evaluation, and list genes
for usable junctions. Percentages use two decimals, with explicit upper bounds
for tiny positive rates and rates just below 5% to avoid misleading rounding.
The selection, suppression, and exit-code rules are unchanged.

Warning wording verification: 79 targeted tests and the complete 182-test
suite pass, including installed-wheel checks. Ruff lint/format and diff
checks pass. Changes remain uncommitted.
