"""Load and validate portable A5SS model parameters."""
from dataclasses import dataclass
from importlib.resources import files
import math
from pathlib import Path
import re
import tomllib

EventKey = tuple[str, str, tuple[int, ...]]
COORDINATE_COLUMNS = (
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE",
)


class SignatureError(ValueError):
    """An input cannot be used to quantify the signature."""


def event_key(chromosome: str, strand: str, coordinates: tuple[int, ...]) -> EventKey:
    return (chromosome.removeprefix("chr"), strand, coordinates)


def _finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SignatureError(f"{label} must be a finite number")
    try:
        number = float(value)
    except OverflowError as exc:
        raise SignatureError(f"{label} must be a finite number") from exc
    if not math.isfinite(number):
        raise SignatureError(f"{label} must be a finite number")
    return number


@dataclass(frozen=True)
class Event:
    id: str
    strand: str
    direction: int
    weight: float

    def __post_init__(self):
        parts = self.id.split(":") if isinstance(self.id, str) else []
        if (len(parts) != 9 or parts[0] != "A5SS" or not parts[1] or not parts[2]
                or any(re.fullmatch(r"[0-9]+", c) is None for c in parts[3:])):
            raise SignatureError(f"Invalid A5SS event id: {self.id!r}")
        if self.strand not in ("+", "-"):
            raise SignatureError(f"Invalid strand for {self.id}")
        if type(self.direction) is not int or self.direction not in (-1, 1):
            raise SignatureError(f"Direction must be +1 or -1 for {self.id}")
        weight = _finite_number(self.weight, f"Weight for {self.id}")
        if weight < 0:
            raise SignatureError(f"Weight must be nonnegative for {self.id}")
        ls, le, ss, se, fs, fe = self.coordinates
        valid_exons = ls < le and ss < se and fs < fe
        if self.strand == "+":
            valid_a5ss = ls == ss and se < le and le <= fs
        else:
            valid_a5ss = ls < ss and le == se and fe <= ls
        if not valid_exons or not valid_a5ss:
            raise SignatureError(f"Invalid A5SS exon geometry for {self.id}")

    @property
    def chromosome(self) -> str:
        return self.id.split(":")[2]

    @property
    def coordinates(self) -> tuple[int, ...]:
        return tuple(int(c) for c in self.id.split(":")[3:])

    @property
    def key(self) -> EventKey:
        return event_key(self.chromosome, self.strand, self.coordinates)


@dataclass(frozen=True)
class Model:
    name: str
    reference_genome: str
    intercept: float
    events: tuple[Event, ...]
    min_read_depth: int = 0


def load_model(path: Path | None = None) -> Model:
    """Load a complete custom model or the bundled Table S2 model."""
    resource = path if path is not None else files("renu_splicing_signature").joinpath("data/insT.toml")
    try:
        with resource.open("rb") as handle:
            config = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise SignatureError(f"Invalid TOML in {resource}: {exc}") from exc
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise SignatureError("Model schema_version must be 1")
    for field in ("name", "reference_genome"):
        if not isinstance(config.get(field), str) or not config[field].strip():
            raise SignatureError(f"Model {field} must be a nonempty string")
    intercept = _finite_number(config.get("intercept"), "Model intercept")
    min_read_depth = config.get("min_read_depth", 0)
    if type(min_read_depth) is not int or min_read_depth < 0:
        raise SignatureError("Model min_read_depth must be a nonnegative integer")
    definitions = config.get("events")
    if not isinstance(definitions, list) or not definitions:
        raise SignatureError("Model must contain at least one event")
    events = []
    for definition in definitions:
        if not isinstance(definition, dict):
            raise SignatureError("Each model event must be a TOML table")
        try:
            event = Event(**{key: definition[key] for key in ("id", "strand", "direction", "weight")})
        except KeyError as exc:
            raise SignatureError(f"Model event is missing {exc.args[0]}") from exc
        if any(event.key == previous.key for previous in events):
            raise SignatureError(f"Duplicate model event coordinates: {event.id}")
        events.append(event)
    return Model(config["name"], config["reference_genome"], intercept, tuple(events), min_read_depth)
