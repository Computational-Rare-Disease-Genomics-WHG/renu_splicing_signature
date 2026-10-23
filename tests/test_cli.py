import csv
import io
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    return subprocess.run(
        [sys.executable, "-m", "renu_splicing_signature", *map(str, args)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
        cwd=ROOT,
    )


def parse_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def test_cli_emits_numeric_tsv(model_file, write_rmats):
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode == 0, result.stderr
    assert not result.stderr
    rows = parse_tsv(result.stdout)
    assert [r["sample"] for r in rows] == ["cases_1", "controls_1", "controls_2"]
    assert float(rows[0]["uniform_sum"]) == pytest.approx(3 * math.log(9))
    assert float(rows[0]["lasso_prob"]) == pytest.approx(
        1 / (1 + math.exp(2.927847423 - 3.75 * math.log(9)))
    )


def test_cli_uses_packaged_model_for_complete_example():
    result = run_cli(ROOT / "examples" / "A5SS.MATS.JC.txt")
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert all(r["n_events_used"] == "6" and r["status"] == "complete" for r in rows)
    assert float(rows[0]["lasso_prob"]) > float(rows[-1]["lasso_prob"])


def test_stem3_selects_separate_signature_on_example(tmp_path):
    details = tmp_path / "details.tsv"
    result = run_cli("--example", "--stem3", "--details", details)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert all(row["n_events_expected"] == "6" for row in rows)
    assert not result.stderr
    assert all(row["n_events_used"] == "6" for row in rows)
    assert all(row["status"] == "complete" for row in rows)
    events = {row["event"].split(":")[1] for row in parse_tsv(details.read_text())}
    assert events == {"GNB1", "KMT2C", "MAP4K4", "CBFB", "TOP1", "CINP"}
    reference = run_cli("--example", "--config", ROOT / "config" / "stem3.toml")
    assert reference.returncode == 0, reference.stderr
    assert result.stdout == reference.stdout


def test_stem3_and_config_are_mutually_exclusive(model_file, tmp_path):
    output = tmp_path / "scores.tsv"
    output.write_text("previous scores")
    result = run_cli(
        "--example", "--stem3", "--config", model_file, "--output", output
    )
    assert result.returncode == 2
    assert "not allowed" in result.stderr
    assert result.stdout == ""
    assert output.read_text() == "previous scores"


@pytest.mark.parametrize("write_output", [False, True])
def test_example_flag_scores_bundled_input(tmp_path, write_output):
    output, details = tmp_path / "scores.tsv", tmp_path / "details.tsv"
    extra = ["--output", output] if write_output else []
    result = run_cli("--example", "--details", details, *extra)
    assert result.returncode == 0, result.stderr
    assert not result.stderr
    scores = output.read_text() if write_output else result.stdout
    reference = run_cli(ROOT / "examples" / "A5SS.MATS.JC.txt")
    assert reference.returncode == 0, reference.stderr
    assert scores == reference.stdout
    assert len(parse_tsv(scores)) == 4
    assert len(parse_tsv(details.read_text())) == 24
    if write_output:
        assert result.stdout == ""


def test_example_flag_rejects_explicit_input_without_writing_output(tmp_path):
    output = tmp_path / "scores.tsv"
    output.write_text("previous scores")
    result = run_cli(
        ROOT / "examples" / "A5SS.MATS.JC.txt", "--example", "--output", output
    )
    assert result.returncode == 2
    assert "not allowed" in result.stderr
    assert result.stdout == ""
    assert output.read_text() == "previous scores"


def test_cli_requires_input_or_example():
    result = run_cli()
    assert result.returncode == 2
    assert "required" in result.stderr
    assert "--example" in result.stderr
    assert result.stdout == ""


def test_partial_warnings_do_not_pollute_stdout(model_file, write_rmats, rmats_rows):
    result = run_cli(write_rmats(rmats_rows[:1]), "--config", model_file)
    assert result.returncode == 0
    assert "warning:" in result.stderr and "MAP4K4" in result.stderr
    assert all(r["status"] == "partial" for r in parse_tsv(result.stdout))
    assert "warning:" not in result.stdout


@pytest.mark.parametrize(
    "ijc,sjc,min_depth,reason",
    [
        ("0,0", "0,0", 0, "zero_coverage"),
        ("NA,NA", "4,4", 0, "missing_counts"),
        ("4,4", "4,4", 9, "below_min_read_depth"),
    ],
)
def test_unusable_controls_report_reasons_without_replacing_outputs(
    model_file,
    write_rmats,
    rmats_rows,
    tmp_path,
    ijc,
    sjc,
    min_depth,
    reason,
):
    model_file.write_text(f"min_read_depth = {min_depth}\n" + model_file.read_text())
    for row in rmats_rows:
        row[14:16] = [ijc, sjc]
    output, details = tmp_path / "scores.tsv", tmp_path / "details.tsv"
    output.write_text("previous scores")
    details.write_text("previous details")
    result = run_cli(
        write_rmats(), "--config", model_file, "--output", output, "--details", details
    )
    assert result.returncode == 2
    assert result.stdout == ""
    assert "Traceback" not in result.stderr
    assert "error:" in result.stderr
    assert "No usable matched controls" in result.stderr
    assert reason in result.stderr
    assert "CTDSP2" in result.stderr and "MAP4K4" in result.stderr
    assert "controls_1" in result.stderr and "controls_2" in result.stderr
    assert output.read_text() == "previous scores"
    assert details.read_text() == "previous details"


def test_files_names_and_details(model_file, write_rmats, tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("patient\n")
    control_names = tmp_path / "control_names.txt"
    control_names.write_text("control_a\ncontrol_b\n")
    output, details = tmp_path / "scores.tsv", tmp_path / "details.tsv"
    result = run_cli(
        write_rmats(),
        "--config",
        model_file,
        "--output",
        output,
        "--details",
        details,
        "--group1-samples",
        names,
        "--group2-samples",
        control_names,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert [row["sample"] for row in parse_tsv(output.read_text())] == [
        "patient",
        "control_a",
        "control_b",
    ]
    rows = parse_tsv(details.read_text())
    assert len(rows) == 6
    assert rows[0]["status"] == "used"
    assert float(rows[0]["control_mean"]) == pytest.approx(-math.log(9) / 2)


@pytest.mark.parametrize("ijc,sjc", [("10", "NA"), ("NA", "10")])
def test_details_preserve_observed_partner_of_missing_count(
    model_file,
    write_rmats,
    rmats_rows,
    tmp_path,
    ijc,
    sjc,
):
    rmats_rows[0][12:14] = [ijc, sjc]
    details = tmp_path / "details.tsv"
    result = run_cli(write_rmats(), "--config", model_file, "--details", details)
    assert result.returncode == 0, result.stderr
    detail = parse_tsv(details.read_text())[0]
    assert (detail["ijc"], detail["sjc"]) == (ijc, sjc)
    assert detail["status"] == "missing_counts"
    assert detail["contribution"] == "NA"
    score = parse_tsv(result.stdout)[0]
    assert score["n_events_used"] == "1"
    assert float(score["weighted_sum"]) == pytest.approx(0.75 * math.log(9))


def test_intercept_changes_only_through_toml(model_file, write_rmats):
    model_file.write_text(model_file.read_text().replace("-2.927847423", "0.0"))
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode == 0
    row = parse_tsv(result.stdout)[0]
    assert row["lasso_linpred"] == row["weighted_sum"]


def test_reverse_groups_cli(model_file, write_rmats):
    result = run_cli(write_rmats(), "--config", model_file, "--control-group", "1")
    assert result.returncode == 0
    row = parse_tsv(result.stdout)[0]
    assert row["sample"] == "controls_1" and row["uniform_sum"] == "0.0"


def test_empty_sample_is_explicitly_na(model_file, write_rmats, rmats_rows):
    for row in rmats_rows:
        row[12:14] = ["0", "0"]
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode == 0
    assert parse_tsv(result.stdout)[0]["lasso_prob"] == "NA"
    assert parse_tsv(result.stdout)[0]["status"] == "unavailable"


@pytest.mark.parametrize(
    "args", [["--intercept", "0"], ["--print-signature"], ["--control-group", "3"]]
)
def test_removed_or_invalid_options_are_rejected(write_rmats, args):
    result = run_cli(write_rmats(), *args)
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("extra", [[], ["--group1-samples", "/nonexistent/names"]])
def test_missing_files_have_clear_errors(model_file, extra):
    result = run_cli("/nonexistent/rmats", "--config", model_file, *extra)
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert "Traceback" not in result.stderr


def test_malformed_input_has_no_traceback(model_file, write_rmats, rmats_rows):
    rmats_rows[0][12] = "broken"
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode != 0
    assert "count" in result.stderr and "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "column,index",
    [
        ("IJC_SAMPLE_1", 12),
        ("SJC_SAMPLE_1", 13),
        ("IJC_SAMPLE_2", 14),
        ("SJC_SAMPLE_2", 15),
    ],
)
def test_oversized_count_has_context_without_traceback(
    model_file,
    write_rmats,
    rmats_rows,
    monkeypatch,
    column,
    index,
):
    monkeypatch.setenv("PYTHONINTMAXSTRDIGITS", "640")
    tokens = rmats_rows[0][index].split(",")
    tokens[0] = "9" * 641
    rmats_rows[0][index] = ",".join(tokens)
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode == 2
    assert result.stdout == ""
    assert "error:" in result.stderr
    assert column in result.stderr and "row 2" in result.stderr
    assert "Traceback" not in result.stderr
    assert len(result.stderr) < 500


def test_unknown_configuration_key_has_clear_error(model_file, write_rmats):
    model_file.write_text("min_read_dept = 100\n" + model_file.read_text())
    result = run_cli(write_rmats(), "--config", model_file)
    assert result.returncode == 2
    assert result.stdout == ""
    assert "min_read_dept" in result.stderr
    assert "error:" in result.stderr and "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "target", ["input", "model", "names", "same_output", "hardlink"]
)
def test_cannot_overwrite_inputs_or_alias_outputs(
    model_file, write_rmats, tmp_path, target
):
    input_path = write_rmats()
    names = tmp_path / "names.txt"
    names.write_text("patient\n")
    if target == "hardlink":
        output = tmp_path / "alias.tsv"
        os.link(input_path, output)
    else:
        output = {
            "input": input_path,
            "model": model_file,
            "names": names,
            "same_output": tmp_path / "same.tsv",
        }[target]
    extra = ["--details", output] if target == "same_output" else []
    before = input_path.read_bytes(), model_file.read_bytes(), names.read_bytes()
    result = run_cli(
        input_path,
        "--config",
        model_file,
        "--group1-samples",
        names,
        "--output",
        output,
        *extra,
    )
    assert result.returncode != 0
    assert "error:" in result.stderr and "output" in result.stderr.lower()
    assert (
        input_path.read_bytes(),
        model_file.read_bytes(),
        names.read_bytes(),
    ) == before


def test_invalid_details_destination_does_not_truncate_existing_scores(
    model_file, write_rmats, tmp_path
):
    output = tmp_path / "scores.tsv"
    output.write_text("previous results")
    result = run_cli(
        write_rmats(),
        "--config",
        model_file,
        "--output",
        output,
        "--details",
        tmp_path / "missing" / "details.tsv",
    )
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert output.read_text() == "previous results"


def test_directory_details_destination_preserves_previous_scores(
    model_file, write_rmats, tmp_path
):
    output = tmp_path / "scores.tsv"
    output.write_text("previous results")
    details = tmp_path / "details.tsv"
    details.mkdir()
    result = run_cli(
        write_rmats(), "--config", model_file, "--output", output, "--details", details
    )
    assert result.returncode != 0
    assert output.read_text() == "previous results"
    assert details.is_dir()


def test_help_and_version():
    help_result = run_cli("--help")
    assert help_result.returncode == 0
    for option in (
        "--example",
        "--control-group",
        "--group1-samples",
        "--group2-samples",
        "--config",
        "--stem3",
    ):
        assert option in help_result.stdout
    assert "--intercept" not in help_result.stdout
    assert run_cli("--version").stdout.strip() == "renu-score 0.1.0"


def test_toml_minimum_depth_controls_scores_and_details(
    model_file, write_rmats, rmats_rows, tmp_path
):
    model_file.write_text("min_read_depth = 5\n" + model_file.read_text())
    rmats_rows[0][12:16] = ["3", "1", "1,3", "0,2"]
    rmats_rows[1][12:16] = ["0", "10", "10,5", "0,5"]
    details = tmp_path / "details.tsv"
    result = run_cli(write_rmats(), "--config", model_file, "--details", details)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert rows[0]["n_events_used"] == "1" and rows[0]["status"] == "partial"
    assert "below_min_read_depth" in result.stderr
    assert parse_tsv(details.read_text())[0]["status"] == "below_min_read_depth"
