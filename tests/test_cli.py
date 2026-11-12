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
    assert all(r["n_events_used"] == "11" and r["status"] == "complete" for r in rows)
    assert float(rows[0]["lasso_prob"]) > float(rows[-1]["lasso_prob"])


def test_stem3_selects_separate_signature_on_example(tmp_path):
    details = tmp_path / "details.tsv"
    result = run_cli("--example", "--stem3", "--details", details)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert all(row["n_events_expected"] == "7" for row in rows)
    assert not result.stderr
    assert all(row["n_events_used"] == "7" for row in rows)
    assert all(row["status"] == "complete" for row in rows)
    events = {row["event"].split(":")[1] for row in parse_tsv(details.read_text())}
    assert events == {"SF1", "PRRC2C", "NUP98", "GNB1", "MAP4K4", "HNRNPC", "MAP3K4"}
    reference = run_cli("--example", "--config", ROOT / "config" / "stem3.toml")
    assert reference.returncode == 0, reference.stderr
    assert result.stdout == reference.stdout


def test_stem3_and_config_are_mutually_exclusive(model_file, tmp_path):
    output = tmp_path / "scores.tsv"
    output.write_text("previous scores")
    result = run_cli("--example", "--stem3", "--config", model_file, "--output", output)
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
    assert len(parse_tsv(details.read_text())) == 44
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
    assert detail["delta_logit_psi_direction"] == "NA"
    score = parse_tsv(result.stdout)[0]
    assert score["n_events_used"] == "1"
    assert float(score["weighted_sum"]) == pytest.approx(0.75 * math.log(9))


def test_details_report_signed_and_direction_adjusted_logit_psi(
    model_file, write_rmats, tmp_path
):
    details = tmp_path / "details.tsv"
    result = run_cli(write_rmats(), "--config", model_file, "--details", details)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(details.read_text())
    # Hand-calculated multiples of ln(9), including both event directions
    # and the individual controls (centering uses their mean logit).
    expected = [
        (1, 1.5, 3, 1.5, 3),
        (-1, -0.5, -1, -0.5, -1),
        (0, 0.5, 1, 0.5, 1),
        (-1, -1.5, -0.75, 1.5, 0.75),
        (1, 0.5, 0.25, -0.5, -0.25),
        (0, -0.5, -0.25, 0.5, 0.25),
    ]
    columns = (
        "logit_psi",
        "delta_logit_psi",
        "delta_logit_psi_weighted",
        "delta_logit_psi_direction",
        "delta_logit_psi_direction_weighted",
    )
    assert len(rows) == len(expected)
    for row, multiples in zip(rows, expected):
        assert "contribution" not in row
        assert "weighted_contribution" not in row
        for column, multiple in zip(columns, multiples):
            assert float(row[column]) == pytest.approx(multiple * math.log(9))


@pytest.mark.parametrize(
    "reason",
    [
        "missing_counts",
        "zero_coverage",
        "below_min_read_depth",
        "absent_from_file",
        "no_control_baseline",
    ],
)
def test_details_logit_statistics_respect_missingness(
    model_file, write_rmats, rmats_rows, tmp_path, reason
):
    if reason == "missing_counts":
        rmats_rows[0][12] = "NA"
    elif reason == "zero_coverage":
        rmats_rows[0][12:14] = ["0", "0"]
    elif reason == "below_min_read_depth":
        model_file.write_text("min_read_depth = 5\n" + model_file.read_text())
    elif reason == "absent_from_file":
        rmats_rows.pop(0)
    else:
        rmats_rows[0][14:16] = ["0,0", "0,0"]
    details = tmp_path / "details.tsv"
    result = run_cli(write_rmats(), "--config", model_file, "--details", details)
    assert result.returncode == 0, result.stderr
    row = parse_tsv(details.read_text())[0]
    assert row["status"] == reason
    if reason == "no_control_baseline":
        assert float(row["logit_psi"]) == pytest.approx(math.log(9))
    else:
        assert row["logit_psi"] == "NA"
    for column in (
        "delta_logit_psi",
        "delta_logit_psi_weighted",
        "delta_logit_psi_direction",
        "delta_logit_psi_direction_weighted",
    ):
        assert row[column] == "NA"


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
    assert run_cli("--version").stdout.strip() == "renu-score 0.2.1"


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


@pytest.mark.parametrize("signature", ["insT", "stemIII"])
@pytest.mark.parametrize("destination", ["--output", "--details"])
def test_cannot_overwrite_bundled_subset_table(tmp_path, signature, destination):
    import shutil

    package = tmp_path / "renu_splicing_signature"
    shutil.copytree(ROOT / "src/renu_splicing_signature", package)
    table = package / "data" / f"subset_refit_weights_{signature}.csv"
    before = table.read_bytes()
    env = dict(os.environ, PYTHONPATH=str(tmp_path))
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "renu_splicing_signature",
            "--example",
            destination,
            str(table),
        ],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "overwrite" in result.stderr
    assert result.stdout == ""
    assert table.read_bytes() == before


@pytest.mark.parametrize(
    "stem3,weight,intercept",
    [
        (False, 5.36662803091663, -3.2633453653868),
        (True, 4.59763537561303, -3.52994047485274),
    ],
)
def test_cli_uses_singleton_refit_for_partial_samples(
    write_rmats, rmats_rows, tmp_path, stem3, weight, intercept
):
    details = tmp_path / "details.tsv"
    args = ["--stem3"] if stem3 else []
    result = run_cli(write_rmats(rmats_rows[1:]), *args, "--details", details)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    # MAP4K4 is direction +1 in both signatures. Controls have log ratios
    # ln(9) and zero, so the case is centered at -1.5*ln(9).
    expected_weighted = -1.5 * math.log(9) * weight
    assert float(rows[0]["weighted_sum"]) == pytest.approx(expected_weighted)
    assert float(rows[0]["lasso_linpred"]) == pytest.approx(
        intercept + expected_weighted
    )
    assert all(
        row["status"] == "partial" and row["n_events_used"] == "1" for row in rows
    )
    event_rows = parse_tsv(details.read_text())
    used = [row for row in event_rows if row["status"] == "used"]
    assert len(used) == 3
    assert all(float(row["weight"]) == weight for row in used)
    assert all(row["weight"] == "NA" for row in event_rows if row["status"] != "used")


@pytest.mark.parametrize("signature", ["insT", "stem3"])
@pytest.mark.parametrize("destination", ["--output", "--details"])
def test_copied_config_cannot_overwrite_implicit_bundled_model(
    tmp_path, signature, destination
):
    import shutil

    package = tmp_path / "renu_splicing_signature"
    shutil.copytree(ROOT / "src/renu_splicing_signature", package)
    bundled = package / "data" / f"{signature}.toml"
    before = bundled.read_bytes()
    copied_config = tmp_path / "copy.toml"
    copied_config.write_bytes(before)
    env = dict(os.environ, PYTHONPATH=str(tmp_path))
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "renu_splicing_signature",
            "--example",
            "--config",
            str(copied_config),
            destination,
            str(bundled),
        ],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "overwrite" in result.stderr
    assert result.stdout == ""
    assert bundled.read_bytes() == before


@pytest.mark.parametrize("stem3,risky_gene", [(False, "TOP1"), (True, "MAP3K4")])
@pytest.mark.parametrize("write_output", [False, True])
def test_high_fpr_warns_and_withholds_only_affected_sample(
    tmp_path, write_rmats, stem3, risky_gene, write_output
):
    with (ROOT / "examples/A5SS.MATS.JC.txt").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        event_rows = list(reader)
    gene_index = header.index("geneSymbol")
    for row in event_rows:
        if row[gene_index] != risky_gene:
            for column in ("IJC_SAMPLE_1", "SJC_SAMPLE_1"):
                index = header.index(column)
                row[index] = "0," + row[index].split(",")[1]
    details = tmp_path / "details.tsv"
    output = tmp_path / "scores.tsv"
    args = ["--stem3"] if stem3 else []
    if write_output:
        args += ["--output", output]
    result = run_cli(write_rmats(event_rows, header), *args, "--details", details)
    assert result.returncode == 0, result.stderr
    assert "warning:" in result.stderr and "scores not provided" in result.stderr
    assert "error:" not in result.stderr and "Traceback" not in result.stderr
    rows = parse_tsv(output.read_text() if write_output else result.stdout)
    assert [r["sample"] for r in rows] == [
        "cases_1",
        "cases_2",
        "controls_1",
        "controls_2",
    ]
    assert rows[0]["status"] == "high_fpr"
    assert rows[0]["n_events_used"] == "1"
    score_columns = ("uniform_sum", "weighted_sum", "lasso_linpred", "lasso_prob")
    assert all(rows[0][column] == "NA" for column in score_columns)
    assert all(
        r["status"] == "complete" and float(r["lasso_prob"]) >= 0 for r in rows[1:]
    )
    subset_warnings = [
        line for line in result.stderr.splitlines() if "scores not provided" in line
    ]
    assert len(subset_warnings) == 1
    assert f"Usable junctions (by gene): {risky_gene}." in result.stderr
    assert "cases_1" in subset_warnings[0]
    detail_rows = [
        r for r in parse_tsv(details.read_text()) if r["sample"] == "cases_1"
    ]
    assert all(
        r["delta_logit_psi_direction"]
        == r["delta_logit_psi_weighted"]
        == r["delta_logit_psi_direction_weighted"]
        == "NA"
        for r in detail_rows
    )
    selected = next(r for r in detail_rows if f"A5SS:{risky_gene}:" in r["event"])
    assert selected["status"] == "high_fpr"
    assert selected["log_ratio"] != "NA"
    if write_output:
        assert not result.stdout
    else:
        assert "warning:" not in result.stdout


@pytest.mark.parametrize("stem3,risky_gene", [(False, "TOP1"), (True, "MAP3K4")])
def test_all_high_fpr_samples_still_write_na_rows_and_exit_zero(
    write_rmats, stem3, risky_gene
):
    with (ROOT / "examples/A5SS.MATS.JC.txt").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        event_rows = [
            row for row in reader if row[header.index("geneSymbol")] == risky_gene
        ]
    args = ["--stem3"] if stem3 else []
    result = run_cli(write_rmats(event_rows, header), *args)
    assert result.returncode == 0, result.stderr
    rows = parse_tsv(result.stdout)
    assert len(rows) == 4
    assert all(
        r["status"] == "high_fpr"
        and r["uniform_sum"]
        == r["weighted_sum"]
        == r["lasso_linpred"]
        == r["lasso_prob"]
        == "NA"
        for r in rows
    )
    assert (
        len(
            [
                line
                for line in result.stderr.splitlines()
                if "scores not provided" in line
            ]
        )
        == 4
    )
