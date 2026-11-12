"""Load the bundled refits, keyed by the exact set of available event IDs."""

import csv
import math
import tomllib
from dataclasses import dataclass
from importlib.resources import files

from .model import Event, SignatureError

SUBSET_TABLES = {
    "insT": "subset_refit_weights_insT.csv",
    "stem3": "subset_refit_weights_stemIII.csv",
}


@dataclass(frozen=True)
class SubsetFit:
    intercept: float
    weights: dict[str, float]


def _fit(row: dict[str, str], columns: dict[str, str], line: int) -> SubsetFit:
    try:
        intercept = float(row["intercept"])
        weights = {
            event_id: float(row[column])
            for column, event_id in columns.items()
            if row[column] != "NA"
        }
        valid = math.isfinite(intercept) and all(
            math.isfinite(weight) and weight >= 0 for weight in weights.values()
        )
    except (TypeError, ValueError) as exc:
        raise SignatureError(f"Invalid subset_weights numbers at row {line}") from exc
    if not valid or not weights:
        raise SignatureError(f"Invalid subset_weights fit at row {line}")
    return SubsetFit(intercept, weights)


def load_subset_fits(
    signature: str, events: tuple[Event, ...]
) -> dict[frozenset[str], SubsetFit]:
    if not isinstance(signature, str) or signature not in SUBSET_TABLES:
        raise SignatureError("Model subset_weights must be 'insT' or 'stem3'")
    data = files("renu_splicing_signature").joinpath("data")
    with data.joinpath(f"{signature}.toml").open("rb") as handle:
        canonical = tomllib.load(handle)
    expected = {(e["id"], e["strand"], e["direction"]) for e in canonical["events"]}
    if {(e.id, e.strand, e.direction) for e in events} != expected:
        raise SignatureError(
            "Model events differ from subset_weights; remove subset_weights "
            "to use a custom fixed model"
        )
    columns = {"w_" + e.id.split(":")[1]: e.id for e in events}
    fits = {}
    with data.joinpath(SUBSET_TABLES[signature]).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        if sorted(reader.fieldnames or []) != sorted(["intercept", *columns]):
            raise SignatureError("Invalid subset_weights CSV header")
        for row in reader:
            if set(row) != {"intercept", *columns}:
                raise SignatureError(f"Invalid subset_weights row {reader.line_num}")
            fit = _fit(row, columns, reader.line_num)
            subset = frozenset(fit.weights)
            if subset in fits:
                raise SignatureError(
                    f"Duplicate subset_weights at row {reader.line_num}"
                )
            fits[subset] = fit
    if len(fits) != 2 ** len(events) - 1:
        raise SignatureError("subset_weights must cover every nonempty event subset")
    return fits
