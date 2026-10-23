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
    for name in ("src", "examples"):
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
    assert any(str(p).endswith("examples/insT.toml") for p in paths)
    assert any(str(p).endswith("examples/A5SS.MATS.JC.txt") for p in paths)
    assert not any(
        {".superpowers", ".venv", ".aws"}.intersection(p.parts) for p in paths
    )


@pytest.mark.parametrize("example", [False, True])
def test_installed_wheel_scores_outside_checkout_with_bundled_model(tmp_path, example):
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
    result = subprocess.run(
        [str(executable), *input_args, "--details", str(details)],
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
    assert all(
        row["status"] == "complete" and row["n_events_used"] == "6" for row in rows
    )
    assert float(rows[0]["lasso_prob"]) > float(rows[-1]["lasso_prob"])
    with details.open(newline="") as handle:
        detail_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(detail_rows) == 24
    assert all(row["status"] == "used" for row in detail_rows)
