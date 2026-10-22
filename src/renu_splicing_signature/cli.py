"""Command-line interface and tabular output."""
import argparse
import csv
from dataclasses import asdict, fields
import io
import os
from pathlib import Path
import sys
import tempfile

from . import __version__
from .model import SignatureError, load_model
from .rmats import read_rmats
from .scoring import DetailRow, ScoreRow, score_samples


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="renu-score", allow_abbrev=False,
        description="Quantify the RNU4-2 splicing signature from GRCh38 A5SS rMATS JC counts.",
    )
    parser.add_argument("rmats", type=Path, help="A5SS.MATS.JC.txt containing both cases and matched controls")
    parser.add_argument("--controls-group", type=int, choices=(1, 2), default=2,
                        help="rMATS group containing matched controls (default: 2)")
    parser.add_argument("--signature", type=Path, help="replacement model TOML, including its intercept")
    parser.add_argument("--output", type=Path, help="scores TSV (default: stdout)")
    parser.add_argument("--details", type=Path, help="optional per-sample, per-event details TSV")
    parser.add_argument("--group1-names", type=Path, help="one group 1 sample name per line, in vector order")
    parser.add_argument("--group2-names", type=Path, help="one group 2 sample name per line, in vector order")
    parser.add_argument("--version", action="version", version=f"renu-score {__version__}")
    return parser


def _same_path(first: Path, second: Path) -> bool:
    return first.resolve() == second.resolve() or (
        first.exists() and second.exists() and first.samefile(second)
    )


def _check_output_paths(outputs: list[Path], inputs: list[Path]) -> None:
    for index, output in enumerate(outputs):
        if output.exists() and not output.is_file():
            raise SignatureError(f"Output path must be a file, not a directory or special device: {output}")
        for other in [*inputs, *outputs[:index]]:
            if _same_path(output, other):
                raise SignatureError(f"Output path {output} would overwrite an input or another output: {other}")


def _tsv(rows: list, row_type: type) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, [field.name for field in fields(row_type)], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: "NA" if value is None else value for key, value in asdict(row).items()})
    return buffer.getvalue()


def _write_files(outputs: list[tuple[Path, str]]) -> None:
    # Stage both files before replacing either destination. A bad details path
    # must not truncate a previous scores file.
    staged: list[tuple[Path, Path]] = []
    try:
        for destination, content in outputs:
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="", dir=destination.parent,
                prefix=f".{destination.name}.", delete=False,
            ) as handle:
                temporary = Path(handle.name)
                staged.append((temporary, destination))
                handle.write(content)
        for temporary, destination in staged:
            os.replace(temporary, destination)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        names_paths = {1: args.group1_names, 2: args.group2_names}
        outputs = [path for path in (args.output, args.details) if path is not None]
        inputs = [args.rmats, Path(__file__).parent / "data" / "insT.toml"]
        inputs.extend(path for path in (args.signature, *names_paths.values()) if path is not None)
        _check_output_paths(outputs, inputs)
        model = load_model(args.signature)
        names = {group: path.read_text(encoding="utf-8-sig").splitlines()
                 for group, path in names_paths.items() if path is not None}
        dataset = read_rmats(args.rmats, model, args.controls_group, names)
        result = score_samples(model, dataset)
        scores_text = _tsv(result.rows, ScoreRow)
        destinations = []
        if args.output is not None:
            destinations.append((args.output, scores_text))
        if args.details is not None:
            destinations.append((args.details, _tsv(result.details, DetailRow)))
        _write_files(destinations)
        for message in result.warnings:
            print(f"warning: {message}", file=sys.stderr)
        if args.output is None:
            sys.stdout.write(scores_text)
    except (SignatureError, OSError, UnicodeError, csv.Error) as exc:
        print(f"renu-score: error: {exc}", file=sys.stderr)
        return 2
    return 0
