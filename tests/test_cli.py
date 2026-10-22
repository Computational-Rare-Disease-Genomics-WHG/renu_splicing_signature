import csv
import io
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    return subprocess.run(
        [sys.executable, "-m", "renu_splicing_signature", *map(str, args)],
        capture_output=True, text=True, env=env, cwd=ROOT,
    )


def parse_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def test_cli_emits_numeric_tsv(model_file, write_rmats):
    result = run_cli(write_rmats(), "--signature", model_file)
    assert result.returncode == 0, result.stderr
    assert not result.stderr
    rows = parse_tsv(result.stdout)
    assert [r["sample"] for r in rows] == ["cases_1", "controls_1", "controls_2"]
    assert float(rows[0]["uniform_sum"]) == pytest.approx(3 * math.log(9))
    assert float(rows[0]["lasso_prob"]) == pytest.approx(1 / (1 + math.exp(2.927847423 - 3.75 * math.log(9))))


def test_cli_uses_packaged_model_for_complete_example():
    result = run_cli(ROOT / "examples" / "A5SS.MATS.JC.txt")
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert all(r["n_events_used"] == "6" and r["status"] == "complete" for r in rows)
    assert float(rows[0]["lasso_prob"]) > float(rows[-1]["lasso_prob"])


def test_partial_warnings_do_not_pollute_stdout(model_file, write_rmats, rmats_rows):
    result = run_cli(write_rmats(rmats_rows[:1]), "--signature", model_file)
    assert result.returncode == 0
    assert "warning:" in result.stderr and "MAP4K4" in result.stderr
    assert all(r["status"] == "partial" for r in parse_tsv(result.stdout))
    assert "warning:" not in result.stdout


def test_files_names_and_details(model_file, write_rmats, tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("patient\n")
    output, details = tmp_path / "scores.tsv", tmp_path / "details.tsv"
    result = run_cli(write_rmats(), "--signature", model_file, "--output", output,
                     "--details", details, "--group1-names", names)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert parse_tsv(output.read_text())[0]["sample"] == "patient"
    rows = parse_tsv(details.read_text())
    assert len(rows) == 6
    assert rows[0]["status"] == "used"
    assert float(rows[0]["control_mean"]) == pytest.approx(-math.log(9) / 2)


def test_intercept_changes_only_through_toml(model_file, write_rmats):
    model_file.write_text(model_file.read_text().replace("-2.927847423", "0.0"))
    result = run_cli(write_rmats(), "--signature", model_file)
    assert result.returncode == 0
    row = parse_tsv(result.stdout)[0]
    assert row["lasso_linpred"] == row["weighted_sum"]


def test_reverse_groups_cli(model_file, write_rmats):
    result = run_cli(write_rmats(), "--signature", model_file, "--controls-group", "1")
    assert result.returncode == 0
    row = parse_tsv(result.stdout)[0]
    assert row["sample"] == "controls_1" and row["uniform_sum"] == "0.0"


def test_empty_sample_is_explicitly_na(model_file, write_rmats, rmats_rows):
    for row in rmats_rows:
        row[12:14] = ["0", "0"]
    result = run_cli(write_rmats(), "--signature", model_file)
    assert result.returncode == 0
    assert parse_tsv(result.stdout)[0]["lasso_prob"] == "NA"
    assert parse_tsv(result.stdout)[0]["status"] == "unavailable"


@pytest.mark.parametrize("args", [["--intercept", "0"], ["--print-signature"], ["--controls-group", "3"]])
def test_removed_or_invalid_options_are_rejected(write_rmats, args):
    result = run_cli(write_rmats(), *args)
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("extra", [[], ["--group1-names", "/nonexistent/names"]])
def test_missing_files_have_clear_errors(model_file, extra):
    result = run_cli("/nonexistent/rmats", "--signature", model_file, *extra)
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert "Traceback" not in result.stderr


def test_malformed_input_has_no_traceback(model_file, write_rmats, rmats_rows):
    rmats_rows[0][12] = "broken"
    result = run_cli(write_rmats(), "--signature", model_file)
    assert result.returncode != 0
    assert "count" in result.stderr and "Traceback" not in result.stderr


@pytest.mark.parametrize("target", ["input", "model", "names", "same_output", "hardlink"])
def test_cannot_overwrite_inputs_or_alias_outputs(model_file, write_rmats, tmp_path, target):
    input_path = write_rmats()
    names = tmp_path / "names.txt"
    names.write_text("patient\n")
    if target == "hardlink":
        output = tmp_path / "alias.tsv"
        os.link(input_path, output)
    else:
        output = {"input": input_path, "model": model_file, "names": names,
                  "same_output": tmp_path / "same.tsv"}[target]
    extra = ["--details", output] if target == "same_output" else []
    before = input_path.read_bytes(), model_file.read_bytes(), names.read_bytes()
    result = run_cli(input_path, "--signature", model_file, "--group1-names", names, "--output", output, *extra)
    assert result.returncode != 0
    assert "error:" in result.stderr and "output" in result.stderr.lower()
    assert (input_path.read_bytes(), model_file.read_bytes(), names.read_bytes()) == before


def test_invalid_details_destination_does_not_truncate_existing_scores(model_file, write_rmats, tmp_path):
    output = tmp_path / "scores.tsv"
    output.write_text("previous results")
    result = run_cli(write_rmats(), "--signature", model_file, "--output", output,
                     "--details", tmp_path / "missing" / "details.tsv")
    assert result.returncode != 0
    assert "error:" in result.stderr
    assert output.read_text() == "previous results"


def test_help_and_version():
    help_result = run_cli("--help")
    assert help_result.returncode == 0 and "--signature" in help_result.stdout
    assert "--intercept" not in help_result.stdout
    assert run_cli("--version").stdout.strip() == "renu-score 0.1.0"
