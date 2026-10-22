"""Build a real source archive, including a deliberately present scratch file."""
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile


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
        [sys.executable, "-c", "import sys; from hatchling.build import build_sdist; build_sdist(sys.argv[1])", str(destination)],
        cwd=project, check=True, capture_output=True, text=True,
    )
    with tarfile.open(next(destination.glob("*.tar.gz"))) as archive:
        paths = [Path(member.name) for member in archive.getmembers()]
    assert any(str(p).endswith("src/renu_splicing_signature/data/insT.toml") for p in paths)
    assert any(str(p).endswith("examples/insT.toml") for p in paths)
    assert not any({".superpowers", ".venv", ".aws"}.intersection(p.parts) for p in paths)
