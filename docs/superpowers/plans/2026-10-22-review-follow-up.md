# Code Review Follow-up Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan inline, with one independent review at the end.

**Goal:** Implement the six recommendations accepted by the user after the code review.

**Architecture:** Retain the four existing modules. Preserve raw counts in parsing, classify measurement availability once in scoring, and keep scoring diagnostics with scoring results or failures.

**Tech Stack:** Python >=3.11, standard library, hatchling, pytest.

**Spec:** The six numbered recommendations in the preceding review, accepted with “Implement these changes”; existing scoring requirements remain in `docs/superpowers/specs/2026-10-22-renu-cli-design.md`.

## Global Constraints

- Keep numerical scoring, CLI options, TSV columns, and sample ordering compatible.
- Continue treating either missing count as an unavailable measurement, preserving the observed partner in details.
- Unknown model and event keys now cause a useful configuration error.
- Work in the supplied checkout on `feat/renu-cli`; preserve the existing staged `.gitignore` edit separately from this follow-up.

## Review Focus

- Either IJC or SJC missing: preserve the other count and omit the contribution.
- No usable controls: report event/sample exclusion reasons, fail cleanly, preserve existing output files.
- Absent events: emit one warning per absent event with affected samples.
- Misspelled root or event keys: fail with the unknown key named.
- Installed wheel: console entry point reads bundled configuration outside the checkout without source imports or network access.

## Tasks

### 1. Preserve partially observed counts

- [x] Add failing parser and CLI detail regressions for missing IJC and SJC.
- [x] Make `Counts.ijc` and `.sjc` optional, preserve pairs, and exclude incomplete pairs in scoring.
- [x] Run parser, scoring, and CLI tests (60 passed).

### 2. Simplify measurement classification and diagnostics

- [x] Add failing regressions for failure diagnostics and duplicate absent-event warnings.
- [x] Classify log ratio and exclusion status once per measurement and reuse after control centering.
- [x] Remove `Dataset.warnings`; scoring owns availability warnings.
- [x] Carry collected warnings through an unusable-score error and print them on stderr.
- [x] Run parser, scoring, and CLI tests (64 passed).

### 3. Validate configuration keys

- [x] Add failing root/event unknown-key regressions.
- [x] Reject keys outside the version 1 schema and document this behavior.
- [x] Run model and CLI tests (56 passed).

### 4. Verify installed execution and review

- [x] Build and install a wheel into a fresh virtual environment with no network or source-path injection.
- [x] Invoke its console entry point outside the checkout, checking complete scores and detail output.
- [x] Run the full suite and whitespace checks, then request one independent review of the working diff (95 passed; whitespace check clean).
- [x] Address material findings and report results (independent review found no substantive findings or deferred polish).

## Results

- Original suite: 85 passed.
- New behavior regressions failed before implementation, covering preserved counts, failure diagnostics, duplicate warnings, and unknown configuration fields.
- Final suite: 95 passed, including an offline wheel install and installed execution outside the checkout.
- `git diff --check` passed.
- Independent final review confirmed all six recommendations were implemented.
- The pre-existing staged `.gitignore` edit is outside this follow-up and was preserved.
