"""Apply signature refits after centering on matched controls."""

import math
from dataclasses import dataclass, replace

from .model import Model, SignatureError
from .rmats import Counts, Dataset
from .subsets import SubsetFit


@dataclass(frozen=True)
class ScoreRow:
    sample: str
    group: int
    role: str
    uniform_sum: float | None
    weighted_sum: float | None
    lasso_linpred: float | None
    lasso_prob: float | None
    n_events_used: int
    n_events_expected: int
    n_zero_coverage_events: int
    status: str


@dataclass(frozen=True)
class DetailRow:
    sample: str
    group: int
    role: str
    event: str
    ijc: int | None
    sjc: int | None
    log_ratio: float | None
    control_mean: float | None
    direction: int
    weight: float | None
    status: str
    logit_psi: float | None
    delta_logit_psi: float | None
    delta_logit_psi_weighted: float | None
    delta_logit_psi_direction: float | None
    delta_logit_psi_direction_weighted: float | None


@dataclass
class ScoreResult:
    rows: list[ScoreRow]
    details: list[DetailRow]
    warnings: list[str]


def _sigmoid(value: float) -> float:
    if value >= 0:
        return 1 / (1 + math.exp(-value))
    exp_value = math.exp(value)
    return exp_value / (1 + exp_value)


def _measurement(pair: Counts | None, min_read_depth: int) -> tuple[float | None, str]:
    """Calculate the log ratio or identify why the measurement is unavailable."""
    if pair is None:
        return None, "absent_from_file"
    if pair.ijc is None or pair.sjc is None:
        return None, "missing_counts"
    if pair.ijc == pair.sjc == 0:
        return None, "zero_coverage"
    if pair.ijc + pair.sjc < min_read_depth:
        return None, "below_min_read_depth"
    # log(2*count+1) differs from log(count+0.5) by log(2), which cancels;
    # integer arithmetic avoids overflow at conversion.
    return math.log(2 * pair.ijc + 1) - math.log(2 * pair.sjc + 1), "used"


def _select_fit(model: Model, details: list[DetailRow]) -> SubsetFit | None:
    available = frozenset(d.event for d in details if d.status == "used")
    if model.subset_fits is None or not available:
        return None
    fit = model.subset_fits.get(available)
    if fit is None:
        raise SignatureError(f"No exact subset refit for {details[0].sample}")
    return fit


def _apply_weights(
    model: Model, details: list[DetailRow], fit: SubsetFit | None
) -> tuple[float, list[DetailRow]]:
    intercept = model.intercept
    if fit is not None:
        intercept, weights = fit.intercept, fit.weights
    elif model.subset_fits is None:
        weights = {e.id: e.weight for e in model.events}
    else:
        weights = {}
    weighted_details = []
    for detail in details:
        weight = weights.get(detail.event)
        weighted = oriented_weighted = None
        if detail.status == "used":
            weighted = weight * detail.delta_logit_psi
            oriented_weighted = weight * detail.delta_logit_psi_direction
            if not math.isfinite(oriented_weighted):
                raise SignatureError(
                    f"Nonfinite weighted contribution for {detail.sample}, {detail.event}"
                )
        weighted_details.append(
            replace(
                detail,
                weight=weight,
                delta_logit_psi_weighted=weighted,
                delta_logit_psi_direction_weighted=oriented_weighted,
            )
        )
    return intercept, weighted_details


def _refit_warning(sample: str, fit: SubsetFit, usable: list[DetailRow]) -> str:
    high_fpr = fit.fpr_pct >= 5
    outcome = "scores not provided" if high_fpr else "scores provided with caution"
    rate = f"{fit.fpr_pct:.2f}%"
    if 0 < fit.fpr_pct < 0.01:
        rate = "<0.01%"
    elif fit.fpr_pct < 5 and rate == "5.00%":
        rate = "<5%"
    controls = "control sample" if fit.n_fp == 1 else "control samples"
    cases = "a case" if fit.n_fp == 1 else "cases"
    lines = [
        f"Sample {sample}: {outcome}.",
        f"The model for these junctions incorrectly flagged {fit.n_fp:,} "
        f"{controls} as {cases} during evaluation ({rate}).",
    ]
    if high_fpr:
        lines.append("This false-positive rate is at or above the 5% limit.")
    genes = ", ".join(d.event.split(":")[1] for d in usable)
    lines.append(f"Usable junctions (by gene): {genes}.")
    return "\n  ".join(lines)


def score_samples(model: Model, dataset: Dataset) -> ScoreResult:
    """Score each sample with its own available subset; never impute missingness."""
    n_samples = len(dataset.samples)
    zero_counts = [0] * n_samples
    details = []
    warnings = []
    for event in model.events:
        counts = dataset.event_counts.get(event.key)
        measurements = [
            _measurement(
                counts[index] if counts is not None else None, model.min_read_depth
            )
            for index in range(n_samples)
        ]
        control_values = [
            value
            for sample, (value, _) in zip(dataset.samples, measurements)
            if sample.role == "control" and value is not None
        ]
        baseline = (
            math.fsum(control_values) / len(control_values) if control_values else None
        )
        if baseline is None and counts is not None:
            warnings.append(
                f"No usable matched controls for event {event.id}; event cannot be centered"
            )
        missing_by_status: dict[str, list[str]] = {}
        for index, (sample, (log_ratio, status)) in enumerate(
            zip(dataset.samples, measurements)
        ):
            pair = counts[index] if counts is not None else None
            contribution = delta_logit_psi = None
            if status == "zero_coverage":
                zero_counts[index] += 1
            if status == "used" and baseline is None:
                status = "no_control_baseline"
            if status == "used":
                delta_logit_psi = log_ratio - baseline
                contribution = event.direction * delta_logit_psi
            if status != "used":
                missing_by_status.setdefault(status, []).append(sample.name)
            details.append(
                DetailRow(
                    sample.name,
                    sample.group,
                    sample.role,
                    event.id,
                    pair.ijc if pair is not None else None,
                    pair.sjc if pair is not None else None,
                    log_ratio,
                    baseline,
                    event.direction,
                    None,
                    status,
                    log_ratio,
                    delta_logit_psi,
                    None,
                    contribution,
                    None,
                )
            )
        for status, names in missing_by_status.items():
            warnings.append(
                f"Event {event.id} unavailable for {', '.join(names)} ({status})"
            )
    rows = []
    for index, sample in enumerate(dataset.samples):
        sample_details = details[index::n_samples]
        usable = [d for d in sample_details if d.status == "used"]
        used = len(usable)
        fit = _select_fit(model, sample_details)
        high_fpr = fit is not None and fit.fpr_pct >= 5
        if fit is not None and (fit.n_fp > 0 or high_fpr):
            warnings.append(_refit_warning(sample.name, fit, usable))
        if high_fpr:
            sample_details = [
                replace(
                    d,
                    weight=fit.weights.get(d.event),
                    status="high_fpr" if d.status == "used" else d.status,
                    delta_logit_psi_direction=None,
                )
                for d in sample_details
            ]
        else:
            intercept, sample_details = _apply_weights(model, sample_details, fit)
        details[index::n_samples] = sample_details
        if used and not high_fpr:
            usable = [d for d in sample_details if d.status == "used"]
            try:
                uniform_sum = math.fsum(d.delta_logit_psi_direction for d in usable)
                weighted_sum = math.fsum(
                    d.delta_logit_psi_direction_weighted for d in usable
                )
                linpred = intercept + weighted_sum
            except OverflowError as exc:
                raise SignatureError(
                    f"Score overflow for {sample.name}; check model parameters"
                ) from exc
            if not math.isfinite(linpred):
                raise SignatureError(
                    f"Nonfinite linear predictor for {sample.name}; check model parameters"
                )
            probability = _sigmoid(linpred)
            status = "complete" if used == len(model.events) else "partial"
        else:
            uniform_sum = weighted_sum = linpred = probability = None
            status = "high_fpr" if high_fpr else "unavailable"
        rows.append(
            ScoreRow(
                sample.name,
                sample.group,
                sample.role,
                uniform_sum,
                weighted_sum,
                linpred,
                probability,
                used,
                len(model.events),
                zero_counts[index],
                status,
            )
        )
    if not any(row.n_events_used for row in rows):
        raise SignatureError(
            "No samples have usable signature events with a matched-control baseline",
            warnings=tuple(warnings),
        )
    return ScoreResult(rows, details, warnings)
