from dataclasses import replace
import math

import pytest

from renu_splicing_signature.model import SignatureError, load_model
from renu_splicing_signature.rmats import read_rmats
from renu_splicing_signature.scoring import score_samples


def test_scores_against_independent_two_event_calculation(model_file, write_rmats):
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    row = result.rows[0]
    assert row.uniform_sum == pytest.approx(3 * math.log(9))
    assert row.weighted_sum == pytest.approx(3.75 * math.log(9))
    assert row.lasso_linpred == pytest.approx(-2.927847423 + 3.75 * math.log(9))
    assert row.lasso_prob == pytest.approx(1 / (1 + math.exp(-row.lasso_linpred)))
    assert row.n_events_used == row.n_events_expected == 2
    assert row.status == "complete"
    assert result.rows[1].uniform_sum == pytest.approx(-math.log(9))
    assert result.rows[2].uniform_sum == pytest.approx(math.log(9))
    assert len(result.details) == 6
    assert not result.warnings


def test_zero_coverage_removes_event_only_in_affected_sample(model_file, write_rmats, rmats_rows):
    rmats_rows[0][12:14] = ["4,0", "0,0"]
    rmats_rows[1][12:14] = ["0,0", "4,4"]
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    assert [r.n_events_used for r in result.rows] == [2, 1, 2, 2]
    partial = result.rows[1]
    assert partial.n_zero_coverage_events == 1
    assert partial.uniform_sum == pytest.approx(1.5 * math.log(9))
    assert partial.weighted_sum == pytest.approx(0.75 * math.log(9))
    assert partial.lasso_linpred == pytest.approx(model.intercept + partial.weighted_sum)
    assert partial.status == "partial"
    assert any("cases_2" in w and "CTDSP2" in w for w in result.warnings)
    detail = next(d for d in result.details if d.sample == "cases_2" and "CTDSP2" in d.event)
    assert detail.status == "zero_coverage"
    assert detail.ijc == detail.sjc == 0
    assert detail.contribution is None


@pytest.mark.parametrize("ijc,sjc", [("0,4", "0,4"), ("NA,4", "4,4"), (",4", "4,4")])
def test_baseline_excludes_missing_controls(model_file, write_rmats, rmats_rows, ijc, sjc):
    rmats_rows[0][14:16] = [ijc, sjc]
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    detail = next(d for d in result.details if d.sample == "cases_1" and "CTDSP2" in d.event)
    assert detail.control_mean == 0
    assert detail.contribution == pytest.approx(math.log(9))
    assert result.rows[0].uniform_sum == pytest.approx(2.5 * math.log(9))


def test_missing_case_count_only_affects_that_sample(model_file, write_rmats, rmats_rows):
    rmats_rows[0][12] = "NA"
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    assert [r.n_events_used for r in result.rows] == [1, 2, 2]
    assert result.rows[0].n_zero_coverage_events == 0


def test_event_without_controls_cannot_be_centered(model_file, write_rmats, rmats_rows):
    rmats_rows[0][14:16] = ["0,0", "0,0"]
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    assert all(r.n_events_used == 1 for r in result.rows)
    assert any("controls" in w and "CTDSP2" in w for w in result.warnings)
    assert next(d for d in result.details if d.sample == "cases_1" and "CTDSP2" in d.event).status == "no_control_baseline"


def test_empty_sample_has_na_scores_not_intercept_probability(model_file, write_rmats, rmats_rows):
    for row in rmats_rows:
        row[12:14] = ["0", "0"]
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    row = result.rows[0]
    assert row.status == "unavailable"
    assert row.n_events_used == 0
    assert row.lasso_prob is row.uniform_sum is row.weighted_sum is row.lasso_linpred is None


def test_all_events_without_usable_controls_fail(model_file, write_rmats, rmats_rows):
    for row in rmats_rows:
        row[14:16] = ["0,0", "0,0"]
    model = load_model(model_file)
    with pytest.raises(SignatureError, match="usable"):
        score_samples(model, read_rmats(write_rmats(), model))


def test_whole_file_missing_event_keeps_na_details(model_file, write_rmats, rmats_rows):
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(rmats_rows[:1]), model))
    assert all(r.status == "partial" for r in result.rows)
    absent = [d for d in result.details if "MAP4K4" in d.event]
    assert len(absent) == 3
    assert all(d.status == "absent_from_file" and d.ijc is None for d in absent)


def test_absent_event_has_one_warning_with_affected_samples(model_file, write_rmats, rmats_rows):
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(rmats_rows[:1]), model))
    warnings = [warning for warning in result.warnings if "MAP4K4" in warning]
    assert len(warnings) == 1
    assert "absent_from_file" in warnings[0]
    assert all(name in warnings[0] for name in ("cases_1", "controls_1", "controls_2"))


def test_swapped_control_group_changes_baseline(model_file, write_rmats):
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model, controls_group=1))
    assert result.rows[0].uniform_sum == 0
    assert result.rows[0].lasso_linpred == model.intercept
    assert result.rows[1].uniform_sum == pytest.approx(-4 * math.log(9))


@pytest.mark.parametrize("intercept,expected", [(1000, 1.0), (-1000, 0.0)])
def test_extreme_linear_predictor_is_stable(model_file, write_rmats, intercept, expected):
    model = replace(load_model(model_file), intercept=intercept)
    result = score_samples(model, read_rmats(write_rmats(), model))
    assert all(r.lasso_prob == expected for r in result.rows)


def test_minimum_depth_filters_each_sample_and_control_baseline(model_file, write_rmats, rmats_rows):
    model_file.write_text("min_read_depth = 5\n" + model_file.read_text())
    rmats_rows[0][12:16] = ["3", "1", "1,3", "0,2"]
    rmats_rows[1][12:16] = ["0", "10", "10,5", "0,5"]
    model = load_model(model_file)
    result = score_samples(model, read_rmats(write_rmats(), model))
    assert [row.n_events_used for row in result.rows] == [1, 1, 2]
    assert result.rows[0].status == "partial"
    assert result.rows[0].n_zero_coverage_events == 0
    assert result.rows[0].uniform_sum == pytest.approx(1.5 * math.log(21))
    assert result.rows[0].weighted_sum == pytest.approx(0.75 * math.log(21))
    event_details = [d for d in result.details if "CTDSP2" in d.event]
    assert [d.status for d in event_details] == ["below_min_read_depth", "below_min_read_depth", "used"]
    assert all(d.control_mean == pytest.approx(math.log(1.4)) for d in event_details)
    assert event_details[0].ijc == 3 and event_details[0].sjc == 1
    assert event_details[0].contribution is None
    assert any("below_min_read_depth" in warning and "cases_1" in warning for warning in result.warnings)


def test_minimum_depth_can_make_an_event_uncenterable(model_file, write_rmats):
    model_file.write_text("min_read_depth = 9\n" + model_file.read_text())
    model = load_model(model_file)
    with pytest.raises(SignatureError, match="usable"):
        score_samples(model, read_rmats(write_rmats(), model))
