"""Apply a fixed signature after centering on matched controls."""
from dataclasses import dataclass
import math

from .model import Model, SignatureError
from .rmats import Dataset


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
    weight: float
    contribution: float | None
    weighted_contribution: float | None
    status: str


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


def score_samples(model: Model, dataset: Dataset) -> ScoreResult:
    """Score each sample with its own available subset; never impute missingness."""
    n_samples = len(dataset.samples)
    contributions: list[list[float]] = [[] for _ in dataset.samples]
    weighted: list[list[float]] = [[] for _ in dataset.samples]
    zero_counts = [0] * n_samples
    details = []
    warnings = list(dataset.warnings)
    for event in model.events:
        counts = dataset.event_counts.get(event.key)
        log_ratios: list[float | None] = []
        for index in range(n_samples):
            pair = counts[index] if counts is not None else None
            if pair is None:
                log_ratios.append(None)
            elif pair.ijc == pair.sjc == 0:
                zero_counts[index] += 1
                log_ratios.append(None)
            else:
                # log(2*count+1) differs from log(count+0.5) by log(2),
                # which cancels; integer arithmetic avoids overflow at conversion.
                log_ratios.append(math.log(2 * pair.ijc + 1) - math.log(2 * pair.sjc + 1))
        control_values = [
            value for sample, value in zip(dataset.samples, log_ratios)
            if sample.role == "control" and value is not None
        ]
        baseline = math.fsum(control_values) / len(control_values) if control_values else None
        if baseline is None and counts is not None:
            warnings.append(f"No usable matched controls for event {event.id}; event cannot be centered")
        missing_by_status: dict[str, list[str]] = {}
        for index, (sample, log_ratio) in enumerate(zip(dataset.samples, log_ratios)):
            pair = counts[index] if counts is not None else None
            contribution = weighted_contribution = None
            if counts is None:
                status = "absent_from_file"
            elif pair is None:
                status = "missing_counts"
            elif pair.ijc == pair.sjc == 0:
                status = "zero_coverage"
            elif baseline is None:
                status = "no_control_baseline"
            else:
                status = "used"
                contribution = event.direction * (log_ratio - baseline)
                weighted_contribution = event.weight * contribution
                if not math.isfinite(weighted_contribution):
                    raise SignatureError(f"Nonfinite weighted contribution for {sample.name}, {event.id}")
                contributions[index].append(contribution)
                weighted[index].append(weighted_contribution)
            if status != "used":
                missing_by_status.setdefault(status, []).append(sample.name)
            details.append(DetailRow(
                sample.name, sample.group, sample.role, event.id,
                pair.ijc if pair is not None else None,
                pair.sjc if pair is not None else None,
                log_ratio, baseline, event.direction, event.weight,
                contribution, weighted_contribution, status,
            ))
        for status, names in missing_by_status.items():
            warnings.append(f"Event {event.id} unavailable for {', '.join(names)} ({status}); scores use available events only")
    rows = []
    for index, sample in enumerate(dataset.samples):
        used = len(contributions[index])
        if used:
            try:
                uniform_sum = math.fsum(contributions[index])
                weighted_sum = math.fsum(weighted[index])
                linpred = model.intercept + weighted_sum
            except OverflowError as exc:
                raise SignatureError(f"Score overflow for {sample.name}; check model parameters") from exc
            if not math.isfinite(linpred):
                raise SignatureError(f"Nonfinite linear predictor for {sample.name}; check model parameters")
            probability = _sigmoid(linpred)
            status = "complete" if used == len(model.events) else "partial"
        else:
            uniform_sum = weighted_sum = linpred = probability = None
            status = "unavailable"
        rows.append(ScoreRow(
            sample.name, sample.group, sample.role, uniform_sum, weighted_sum,
            linpred, probability, used, len(model.events), zero_counts[index], status,
        ))
    if not any(row.n_events_used for row in rows):
        raise SignatureError("No samples have usable signature events with a matched-control baseline")
    return ScoreResult(rows, details, warnings)
