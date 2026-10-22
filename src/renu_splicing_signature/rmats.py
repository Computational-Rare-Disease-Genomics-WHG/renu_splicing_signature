"""Read sample-ordered junction counts without using rMATS statistical results."""
import csv
from dataclasses import dataclass
from pathlib import Path
import re

from .model import COORDINATE_COLUMNS, EventKey, Model, SignatureError, event_key

COUNT_COLUMNS = ("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
REQUIRED_COLUMNS = ("chr", "strand", *COORDINATE_COLUMNS, *COUNT_COLUMNS)


@dataclass(frozen=True)
class Counts:
    ijc: int | None
    sjc: int | None


@dataclass(frozen=True)
class Sample:
    name: str
    group: int
    role: str


@dataclass
class Dataset:
    samples: tuple[Sample, ...]
    event_counts: dict[EventKey, tuple[Counts, ...]]


def _parse_counts(value: str, context: str) -> list[int | None]:
    counts = []
    for token in value.split(","):
        token = token.strip()
        if token.upper() == "NA" or not token:
            counts.append(None)
        elif re.fullmatch(r"[0-9]+", token):
            counts.append(int(token))
        else:
            raise SignatureError(f"Invalid count {token!r} in {context}; expected a nonnegative integer or NA")
    return counts


def read_rmats(
    path: Path,
    model: Model,
    controls_group: int = 2,
    group_names: dict[int, list[str]] | None = None,
) -> Dataset:
    if controls_group not in (1, 2):
        raise SignatureError("controls_group must be 1 or 2")
    selected = {event.key: event for event in model.events}
    event_counts = {}
    group_sizes = None
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        missing = [name for name in REQUIRED_COLUMNS if name not in header]
        if missing:
            raise SignatureError("Missing required columns: " + ", ".join(missing))
        if any(header.count(name) != 1 for name in REQUIRED_COLUMNS):
            raise SignatureError("Duplicate required columns in rMATS header")
        indices = {name: header.index(name) for name in REQUIRED_COLUMNS}
        for row in reader:
            if not row:
                continue
            line = reader.line_num
            if len(row) != len(header):
                raise SignatureError(f"Row {line} has {len(row)} fields; expected {len(header)}")
            try:
                coords = tuple(int(row[indices[name]]) for name in COORDINATE_COLUMNS)
            except ValueError as exc:
                raise SignatureError(f"Invalid exon coordinates at row {line}") from exc
            key = event_key(row[indices["chr"]], row[indices["strand"]], coords)
            if key not in selected:
                continue
            if key in event_counts:
                raise SignatureError(f"Duplicate selected event at row {line}: {selected[key].id}")
            vectors = [_parse_counts(row[indices[name]], f"{name}, row {line}") for name in COUNT_COLUMNS]
            if len(vectors[0]) != len(vectors[1]) or len(vectors[2]) != len(vectors[3]):
                raise SignatureError(f"IJC/SJC vector length mismatch at row {line}")
            sizes = (len(vectors[0]), len(vectors[2]))
            if group_sizes is not None and group_sizes != sizes:
                raise SignatureError(f"Inconsistent group sizes at row {line}: {sizes}; expected {group_sizes}")
            group_sizes = sizes
            pairs = []
            for ijcs, sjcs in ((vectors[0], vectors[1]), (vectors[2], vectors[3])):
                pairs.extend(Counts(i, s) for i, s in zip(ijcs, sjcs))
            event_counts[key] = tuple(pairs)
    if group_sizes is None:
        raise SignatureError("No selected A5SS events found; check the input and reference genome (bundled model: GRCh38)")
    samples = []
    for group, size in enumerate(group_sizes, start=1):
        role = "control" if group == controls_group else "case"
        names = (group_names or {}).get(group, [f"{role}s_{i + 1}" for i in range(size)])
        if len(names) != size or any(not name.strip() or any(c in name for c in "\t\r\n") for name in names):
            raise SignatureError(f"Sample names for group {group} must contain {size} nonempty names")
        samples.extend(Sample(name, group, role) for name in names)
    if len({sample.name for sample in samples}) != len(samples):
        raise SignatureError("Sample names must be unique across both groups")
    return Dataset(tuple(samples), event_counts)
