import tomllib
from pathlib import Path

import pytest

from renu_splicing_signature.model import SignatureError, load_model


def test_exported_default_parameters():
    model = load_model()
    assert model.intercept == -8.33698995873104
    assert model.min_read_depth == 0
    assert model.reference_genome == "GRCh38"
    assert [(e.id, e.strand, e.direction, e.weight) for e in model.events] == [
        (
            "A5SS:TOP1:chr20:41113969:41114155:41113969:41114137:41115370:41115439",
            "+",
            -1,
            1.59583263727201,
        ),
        (
            "A5SS:PTBP3:chr9:112297831:112297916:112297840:112297916:112275843:112276013",
            "-",
            -1,
            1.66554096320677,
        ),
        (
            "A5SS:CTDSP2:chr12:57827533:57827590:57827551:57827590:57826995:57827097",
            "-",
            -1,
            1.35459404495048,
        ),
        (
            "A5SS:RAPGEF6:chr5:131446482:131446703:131446488:131446703:131442348:131442537",
            "-",
            -1,
            1.1151504974413,
        ),
        (
            "A5SS:MAP4K4:chr2:101887087:101887261:101887087:101887237:101887777:101887937",
            "+",
            1,
            1.10203778267566,
        ),
        (
            "A5SS:NOSIP:chr19:49557080:49557249:49557089:49557249:49556874:49556993",
            "-",
            1,
            1.32609242384374,
        ),
        (
            "A5SS:STX16:chr20:58667489:58667608:58667489:58667597:58667986:58668127",
            "+",
            1,
            0.468438646246804,
        ),
        (
            "A5SS:RPS15A:chr16:18790243:18790334:18790309:18790334:18788980:18789118",
            "-",
            1,
            0.879075351189685,
        ),
        (
            "A5SS:GMPR2:chr14:24233478:24233598:24233478:24233586:24235736:24235820",
            "+",
            -1,
            0.438133779049095,
        ),
        (
            "A5SS:KLHDC4:chr16:87756398:87756477:87756421:87756477:87755193:87755292",
            "-",
            1,
            0.710113836564136,
        ),
        (
            "A5SS:PHF3:chr6:63711166:63711362:63711166:63711286:63711585:63726011",
            "+",
            -1,
            0.275871869501411,
        ),
    ]


def test_visible_model_matches_packaged_default():
    example = Path(__file__).resolve().parents[1] / "config" / "insT.toml"
    assert load_model(example) == load_model()
    assert tomllib.loads(example.read_text())["schema_version"] == 1


def test_exported_stem3_parameters_and_editable_copy():
    root = Path(__file__).resolve().parents[1]
    model = load_model(root / "src/renu_splicing_signature/data/stem3.toml")
    assert model == load_model(root / "config/stem3.toml")
    assert model.intercept == -9.60267638133653
    assert model.reference_genome == "GRCh38"
    assert [(e.id, e.strand, e.direction, e.weight) for e in model.events] == [
        (
            "A5SS:SF1:chr11:64767566:64767844:64767570:64767844:64767191:64767251",
            "-",
            1,
            4.96907898342119,
        ),
        (
            "A5SS:PRRC2C:chr1:171584418:171584541:171584418:171584526:171587002:171587221",
            "+",
            1,
            3.84154874596638,
        ),
        (
            "A5SS:NUP98:chr11:3760487:3760626:3760538:3760626:3753315:3753408",
            "-",
            1,
            3.79622221015688,
        ),
        (
            "A5SS:GNB1:chr1:1787321:1787437:1787325:1787437:1785285:1787053",
            "-",
            -1,
            3.71763940311118,
        ),
        (
            "A5SS:MAP4K4:chr2:101887087:101887261:101887087:101887237:101887777:101887937",
            "+",
            1,
            2.33849486301692,
        ),
        (
            "A5SS:HNRNPC:chr14:21230957:21231072:21230996:21231072:21230318:21230366",
            "-",
            -1,
            0.633537728191756,
        ),
        (
            "A5SS:MAP3K4:chr6:161093772:161093851:161093772:161093839:161097079:161097176",
            "+",
            -1,
            0.165414905155582,
        ),
    ]


def test_toml_controls_intercept(model_file):
    model_file.write_text(model_file.read_text().replace("-2.927847423", "-7.5"))
    assert load_model(model_file).intercept == -7.5


def test_missing_minimum_read_depth_defaults_to_zero(model_file):
    assert load_model(model_file).min_read_depth == 0


def test_minimum_read_depth_is_configurable(model_file):
    model_file.write_text("min_read_depth = 10\n" + model_file.read_text())
    assert load_model(model_file).min_read_depth == 10


def test_rejects_unknown_model_keys(model_file):
    model_file.write_text("min_read_dept = 100\n" + model_file.read_text())
    with pytest.raises(SignatureError, match="min_read_dept"):
        load_model(model_file)


def test_rejects_unknown_event_keys(model_file):
    content = model_file.read_text().replace(
        "weight = 2.0", "weight = 2.0\nweigth = 3.0"
    )
    model_file.write_text(content)
    with pytest.raises(SignatureError, match="weigth"):
        load_model(model_file)


@pytest.mark.parametrize("value", ["-1", "1.5", "true", '"10"', "inf", "nan"])
def test_minimum_read_depth_requires_nonnegative_integer(model_file, value):
    model_file.write_text(f"min_read_depth = {value}\n" + model_file.read_text())
    with pytest.raises(SignatureError, match="min_read_depth"):
        load_model(model_file)


@pytest.mark.parametrize(
    "old,new",
    [
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
    ],
)
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
    model_file.write_text(
        'schema_version=1\nname="empty"\nreference_genome="GRCh38"\nintercept=0\nevents=[]'
    )
    with pytest.raises(SignatureError):
        load_model(model_file)


def test_invalid_toml_has_a_useful_error(model_file):
    model_file.write_text("[invalid")
    with pytest.raises(SignatureError, match="TOML"):
        load_model(model_file)
