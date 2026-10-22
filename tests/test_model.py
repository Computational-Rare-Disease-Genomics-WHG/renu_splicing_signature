import tomllib
from pathlib import Path

import pytest

from renu_splicing_signature.model import SignatureError, load_model


def test_published_default_parameters():
    model = load_model()
    assert model.intercept == -2.927847423
    assert model.reference_genome == "GRCh38"
    assert [(e.id.split(":")[1], e.strand, e.direction, e.weight) for e in model.events] == [
        ("CTDSP2", "-", -1, 0.825615979283261),
        ("PTBP3", "-", -1, 0.783806558002269),
        ("MAP4K4", "+", 1, 0.737422747142889),
        ("TOP1", "+", -1, 0.560782708354878),
        ("STX16", "+", 1, 0.401105638015341),
        ("RPS15A", "-", 1, 0.290736572720099),
    ]
    assert model.events[-1].coordinates == (18790243, 18790334, 18790309, 18790334, 18788980, 18789118)


def test_visible_model_matches_packaged_default():
    example = Path(__file__).resolve().parents[1] / "examples" / "insT.toml"
    assert load_model(example) == load_model()
    assert tomllib.loads(example.read_text())["schema_version"] == 1


def test_toml_controls_intercept(model_file):
    model_file.write_text(model_file.read_text().replace("-2.927847423", "-7.5"))
    assert load_model(model_file).intercept == -7.5


@pytest.mark.parametrize("old,new", [
    ("intercept = -2.927847423", "intercept = nan"),
    ("intercept = -2.927847423", "intercept = true"),
    ("intercept = -2.927847423", ""),
    ("schema_version = 1", "schema_version = 2"),
    ("direction = 1", "direction = 0"),
    ("direction = 1", "direction = true"),
    ("weight = 2.0", "weight = -1.0"),
    ("weight = 2.0", "weight = inf"),
    ('strand = "-"', 'strand = "?"'),
    ("57827551:57827590", "57827500:57827590"),
    ("57826995:57827097", "57828000:57828100"),
    ("A5SS:CTDSP2", "SE:CTDSP2"),
])
def test_rejects_invalid_models(model_file, old, new):
    model_file.write_text(model_file.read_text().replace(old, new))
    with pytest.raises(SignatureError):
        load_model(model_file)


def test_rejects_duplicate_coordinate_identity(model_file):
    content = model_file.read_text()
    model_file.write_text(content + "[[events]]" + content.split("[[events]]")[1])
    with pytest.raises(SignatureError, match="[Dd]uplicate"):
        load_model(model_file)


def test_rejects_empty_model(model_file):
    model_file.write_text('schema_version=1\nname="empty"\nreference_genome="GRCh38"\nintercept=0\nevents=[]')
    with pytest.raises(SignatureError):
        load_model(model_file)


def test_invalid_toml_has_a_useful_error(model_file):
    model_file.write_text("[invalid")
    with pytest.raises(SignatureError, match="TOML"):
        load_model(model_file)
