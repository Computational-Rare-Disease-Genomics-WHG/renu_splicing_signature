"""Verify release contents and execution of the installed wheel."""

import csv
import io
import os
import shutil
import subprocess
import sys
import tarfile
import venv
from pathlib import Path

import pytest


def test_sdist_ships_model_and_examples_without_workspace_artifacts(tmp_path):
    root = Path(__file__).resolve().parents[1]
    project = tmp_path / "source"
    project.mkdir()
    for name in ("pyproject.toml", "README.md", ".gitignore"):
        shutil.copy2(root / name, project / name)
    for name in ("src", "config", "examples"):
        shutil.copytree(root / name, project / name)
    for name in (".superpowers/sdd", ".venv", ".aws"):
        directory = project / name
        directory.mkdir(parents=True)
        (directory / "local-only.txt").write_text("must not be distributed")
    destination = tmp_path / "dist"
    destination.mkdir()
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from hatchling.build import build_sdist; build_sdist(sys.argv[1])",
            str(destination),
        ],
        cwd=project,
        check=True,
        capture_output=True,
        text=True,
    )
    with tarfile.open(next(destination.glob("*.tar.gz"))) as archive:
        paths = [Path(member.name) for member in archive.getmembers()]
    assert any(
        str(p).endswith("src/renu_splicing_signature/data/insT.toml") for p in paths
    )
    assert any(str(p).endswith("config/insT.toml") for p in paths)
    assert any(
        str(p).endswith("src/renu_splicing_signature/data/stem3.toml") for p in paths
    )
    assert any(str(p).endswith("config/stem3.toml") for p in paths)
    for signature in ("insT", "stemIII"):
        assert any(
            str(p).endswith(
                f"src/renu_splicing_signature/data/subset_refit_weights_{signature}.csv"
            )
            for p in paths
        )
    assert any(str(p).endswith("examples/A5SS.MATS.JC.txt") for p in paths)
    assert not any(
        {".superpowers", ".venv", ".aws"}.intersection(p.parts) for p in paths
    )


@pytest.mark.parametrize("example", [False, True])
@pytest.mark.parametrize("stem3", [False, True])
def test_installed_wheel_scores_outside_checkout_with_bundled_model(
    tmp_path, example, stem3
):
    root = Path(__file__).resolve().parents[1]
    wheel_dir = tmp_path / "dist"
    wheel_dir.mkdir()
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from hatchling.build import build_wheel; build_wheel(sys.argv[1])",
            str(wheel_dir),
        ],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    environment = tmp_path / "environment"
    venv.EnvBuilder(with_pip=True).create(environment)
    scripts = environment / ("Scripts" if os.name == "nt" else "bin")
    python = scripts / ("python.exe" if os.name == "nt" else "python")
    executable = scripts / ("renu-score.exe" if os.name == "nt" else "renu-score")
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    subprocess.run(
        [
            str(python),
            "-I",
            "-m",
            "pip",
            "install",
            "--no-index",
            "--no-deps",
            "--disable-pip-version-check",
            str(next(wheel_dir.glob("*.whl"))),
        ],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    if example:
        input_args = ["--example"]
    else:
        input_path = tmp_path / "A5SS.MATS.JC.txt"
        shutil.copy2(root / "examples" / "A5SS.MATS.JC.txt", input_path)
        input_args = [str(input_path)]
    details = tmp_path / "details.tsv"
    signature_args = ["--stem3"] if stem3 else []
    result = subprocess.run(
        [str(executable), *input_args, *signature_args, "--details", str(details)],
        check=False,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert not result.stderr
    rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
    assert [row["sample"] for row in rows] == [
        "cases_1",
        "cases_2",
        "controls_1",
        "controls_2",
    ]
    expected_events = 7 if stem3 else 11
    assert all(
        row["status"] == "complete"
        and int(row["n_events_used"]) == expected_events
        and int(row["n_events_expected"]) == expected_events
        for row in rows
    )
    assert float(rows[0]["lasso_prob"]) > float(rows[-1]["lasso_prob"])
    with details.open(newline="") as handle:
        detail_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(detail_rows) == 4 * expected_events
    assert all(row["status"] == "used" for row in detail_rows)
    if stem3:
        assert {row["event"].split(":")[1] for row in detail_rows} == {
            "SF1",
            "PRRC2C",
            "NUP98",
            "GNB1",
            "MAP4K4",
            "HNRNPC",
            "MAP3K4",
        }

    # Run a partial signature with the installed executable, outside the checkout.
    partial = tmp_path / "partial.txt"
    with (root / "examples/A5SS.MATS.JC.txt").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        event_rows = [
            row for row in reader if row[header.index("geneSymbol")] == "MAP4K4"
        ]
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(event_rows)
    result = subprocess.run(
        [str(executable), str(partial), *signature_args, "--details", str(details)],
        check=False,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
    expected_weight = 4.59763537561303 if stem3 else 5.36662803091663
    expected_intercept = -3.52994047485274 if stem3 else -3.2633453653868
    assert len(rows) == 4
    assert all(
        row["status"] == "partial" and row["n_events_used"] == "1" for row in rows
    )
    for row in rows:
        assert float(row["lasso_linpred"]) - float(
            row["weighted_sum"]
        ) == pytest.approx(expected_intercept)
    with details.open(newline="") as handle:
        detail_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(
        float(row["weight"]) == expected_weight
        for row in detail_rows
        if row["status"] == "used"
    )
    assert "scores provided with caution" in result.stderr
    assert "incorrectly flagged" in result.stderr
    assert "Usable junctions (by gene): MAP4K4." in result.stderr

    # The installed tool must use the bundled FPR metadata to suppress a
    # high-FPR singleton, with warning-only exit behavior for every sample.
    risky_gene = "MAP3K4" if stem3 else "TOP1"
    with (root / "examples/A5SS.MATS.JC.txt").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        event_rows = [
            row for row in reader if row[header.index("geneSymbol")] == risky_gene
        ]
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(event_rows)
    result = subprocess.run(
        [str(executable), str(partial), *signature_args, "--details", str(details)],
        check=False,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
    assert len(rows) == 4
    assert all(
        row["status"] == "high_fpr"
        and row["uniform_sum"]
        == row["weighted_sum"]
        == row["lasso_linpred"]
        == row["lasso_prob"]
        == "NA"
        for row in rows
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
    with details.open(newline="") as handle:
        detail_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["delta_logit_psi_direction_weighted"] == "NA" for row in detail_rows)
