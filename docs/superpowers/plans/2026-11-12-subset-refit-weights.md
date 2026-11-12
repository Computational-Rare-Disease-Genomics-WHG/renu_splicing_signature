# Subset refit weights implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Score every sample with the supplied refit for its available junctions.

**Architecture:** Keep availability and control centering separate from applying weights. Bundle intercept/weight-only CSVs, indexed by sets of event IDs. Bundled TOMLs opt in by signature name; custom models without that setting retain fixed parameters.

**Tech Stack:** Python >=3.11 standard library, pytest, hatchling, ruff.

**Spec:** `docs/superpowers/specs/2026-11-12-subset-refit-weights-design.md`

## Global constraints

- Preserve all supplied coefficient and intercept precision and use the full-set rows too.
- No fitting, classification thresholds, new runtime dependencies, or network access.
- Keep uniform scores, missingness statuses, centering, and sample order unchanged.
- Work on the existing `subset_weights` feature branch and leave the changes available for review.

## Review focus

- A fitted zero weight still identifies an available junction and counts toward completeness.
- Different subsets of the same size must select different refits, independently for controls too.
- Exclusions for missing counts, absent rows, low depth, and absent control baselines affect selection.
- Edited custom models must not silently inherit unrelated bundled subset parameters.
- Installed tools must find and protect bundled data outside the checkout.

### Task 1: Bundle and load subset refits

**Files:** Create `src/renu_splicing_signature/subsets.py` and two compact CSVs in `src/renu_splicing_signature/data/`; modify `model.py`, both bundled TOMLs and `config/` copies; test `tests/test_model.py` and `tests/test_subsets.py`.

**Interfaces:** `SubsetFit(intercept: float, weights: dict[str, float])`; `Model.subset_fits: dict[frozenset[str], SubsetFit] | None`; `load_subset_fits(signature: str, events: tuple[Event, ...])` maps event IDs to weights and validates the supplied signature identity.

- [x] Write failing tests asserting 2,047/127 loaded fits, full-set coefficients/intercepts, singleton/zero coefficients, matching editable copies, rejection of invalid subset identifiers and incompatible event definitions.
- [x] Run `.venv/bin/python -m pytest tests/test_model.py tests/test_subsets.py -q`; confirm failures for missing feature/new parameters.
- [x] Generate compact CSVs preserving original intercept/weight strings. Add strict loading, explicit optional TOML `subset_weights = "insT"` / `"stem3"`, and updated full-set parameters.
- [x] Run those tests again; expect all pass. Compare every compact cell to both Downloads exports with a one-off independent script.

### Task 2: Apply sample-specific refits

**Files:** Modify `src/renu_splicing_signature/scoring.py`; extend `tests/test_scoring.py` and `tests/test_subsets.py`.

**Interfaces:** `score_samples(model: Model, dataset: Dataset) -> ScoreResult` retains its public signature. `DetailRow.weight` allows `None` for events absent from the selected refit. Existing fixed custom models retain their behavior.

- [x] Add failing numeric tests using real insT/Stem III subsets with different case/control availability, all exclusion reasons, full and empty samples, and zero coefficients. Assert weighted detail contributions sum to weighted scores and selected intercepts reconcile to linear predictors.
- [x] Run `.venv/bin/python -m pytest tests/test_scoring.py tests/test_subsets.py -q`; confirm feature failures.
- [x] Produce centered, unweighted details first; select the refit after each sample's usable event set is known; apply its weights consistently to details and sums, preserving event-major detail order and errors for nonfinite scores.
- [x] Run the targeted tests; expect all pass, including existing fixed-model calculations.

### Task 3: CLI, packaging, and documentation

**Files:** Modify `cli.py`, `README.md`, `tests/test_cli.py`, and `tests/test_distribution.py`; add data provenance in `src/renu_splicing_signature/data/`.

**Interfaces:** CLI flags and score columns remain unchanged. Details show actual per-sample coefficients. Both wheel and sdist include compact tables.

- [x] Add failing checks for CLI subset scores/details, output paths aliasing bundled tables, sdist resources, and installed wheel partial scoring for both signatures.
- [x] Run the targeted new CLI/distribution tests and confirm the expected failures.
- [x] Protect bundled tables from output overwrites and document subset selection, full insT update, precision, provenance, and custom-model opt-out.
- [x] Run `.venv/bin/python -m pytest -q`, `ruff check src tests`, `ruff format --check src tests`, and `git diff --check`; expect all pass.
- [x] Request one independent whole-change review, resolve material findings, and report the verified result.

## Verification result

157 tests passed, including all subset scoring checks and installed wheel
execution for full and partial signatures. All compact parameter strings were
verified against the supplied exports. Ruff lint/format and git diff checks
passed. An independent review identified one implicit model-input overwrite
gap; its regression failed before the fix and passes afterward. No remaining
review findings. Changes are left uncommitted on `subset_weights`.
