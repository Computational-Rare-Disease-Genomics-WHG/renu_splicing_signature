import pytest

from renu_splicing_signature.model import SignatureError, load_model
from renu_splicing_signature.rmats import read_rmats


def test_reads_standard_rmats_with_duplicate_id_headers(model_file, write_rmats):
    data = read_rmats(write_rmats(bom=True), load_model(model_file))
    assert [(s.name, s.group, s.role) for s in data.samples] == [
        ("cases_1", 1, "case"),
        ("controls_1", 2, "control"),
        ("controls_2", 2, "control"),
    ]
    assert next(iter(data.event_counts.values()))[0].ijc == 4


def test_reversed_groups_name_samples_by_role(model_file, write_rmats):
    data = read_rmats(write_rmats(), load_model(model_file), controls_group=1)
    assert [s.name for s in data.samples] == ["controls_1", "cases_1", "cases_2"]


@pytest.mark.parametrize("column,value", [(4, "+"), (10, "57827098")])
def test_wrong_strand_or_flank_does_not_match(
    model_file, write_rmats, rmats_rows, column, value
):
    rmats_rows[0][column] = value
    model = load_model(model_file)
    data = read_rmats(write_rmats(), model)
    assert model.events[0].key not in data.event_counts
    assert model.events[1].key in data.event_counts


def test_chromosome_prefix_is_optional(model_file, write_rmats, rmats_rows):
    rmats_rows[0][3] = "12"
    assert len(read_rmats(write_rmats(), load_model(model_file)).event_counts) == 2


@pytest.mark.parametrize("token", ["NA", ""])
def test_explicit_missing_counts_are_sample_specific(
    model_file, write_rmats, rmats_rows, token
):
    rmats_rows[0][14] = token + ",4"
    data = read_rmats(write_rmats(), load_model(model_file))
    counts = next(iter(data.event_counts.values()))
    assert counts[1].ijc is None
    assert counts[1].sjc == 4
    assert counts[2].ijc == 4


@pytest.mark.parametrize("value", ["-1", "1.5", "x", "inf"])
def test_rejects_invalid_counts(model_file, write_rmats, rmats_rows, value):
    rmats_rows[0][12] = value
    with pytest.raises(SignatureError, match="count"):
        read_rmats(write_rmats(), load_model(model_file))


def test_rejects_mismatched_vectors(model_file, write_rmats, rmats_rows):
    rmats_rows[0][14] = "1"
    with pytest.raises(SignatureError, match="[Ll]ength|[Ss]ize"):
        read_rmats(write_rmats(), load_model(model_file))


def test_rejects_inconsistent_group_sizes(model_file, write_rmats, rmats_rows):
    rmats_rows[0][12:14] = ["1,2", "3,4"]
    with pytest.raises(SignatureError, match="[Ll]ength|[Ss]ize"):
        read_rmats(write_rmats(), load_model(model_file))


def test_rejects_duplicate_selected_rows(model_file, write_rmats, rmats_rows):
    with pytest.raises(SignatureError, match="[Dd]uplicate"):
        read_rmats(write_rmats(rmats_rows + [rmats_rows[0]]), load_model(model_file))


def test_reports_missing_columns(model_file, write_rmats):
    with pytest.raises(SignatureError, match="required columns"):
        read_rmats(write_rmats(rows=[], headers=["ID"]), load_model(model_file))


def test_absent_event_keeps_remaining_events(model_file, write_rmats, rmats_rows):
    model = load_model(model_file)
    data = read_rmats(write_rmats(rmats_rows[:1]), model)
    assert set(data.event_counts) == {model.events[0].key}


def test_optional_names_preserve_vector_order(model_file, write_rmats):
    data = read_rmats(
        write_rmats(),
        load_model(model_file),
        group_names={1: ["patient"], 2: ["c1", "c2"]},
    )
    assert [s.name for s in data.samples] == ["patient", "c1", "c2"]


@pytest.mark.parametrize(
    "names", [{1: ["a", "b"]}, {1: [""]}, {1: ["controls_1"]}, {2: ["same", "same"]}]
)
def test_rejects_invalid_names(model_file, write_rmats, names):
    with pytest.raises(SignatureError, match="name"):
        read_rmats(write_rmats(), load_model(model_file), group_names=names)
