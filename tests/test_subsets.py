"""Subset refits use availability, including available zero-weight events."""

from itertools import combinations
from pathlib import Path

import pytest

from renu_splicing_signature.model import SignatureError, load_model

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("signature,n_fits", [("insT", 2047), ("stem3", 127)])
def test_bundled_refits_cover_every_nonempty_event_subset(signature, n_fits):
    model = load_model(ROOT / "config" / f"{signature}.toml")
    ids = [e.id for e in model.events]
    expected = {
        frozenset(subset)
        for k in range(1, len(ids) + 1)
        for subset in combinations(ids, k)
    }
    assert model.subset_fits is not None
    assert len(model.subset_fits) == n_fits
    assert set(model.subset_fits) == expected
    for subset, fit in model.subset_fits.items():
        assert set(fit.weights) == subset
    full = model.subset_fits[frozenset(ids)]
    assert full.intercept == model.intercept
    assert full.weights == {e.id: e.weight for e in model.events}


def test_old_custom_model_keeps_fixed_parameters(model_file):
    model = load_model(model_file)
    assert model.subset_fits is None
    assert model.intercept == -2.927847423


@pytest.mark.parametrize("setting", ['"unknown"', "true", "42", '"../insT"'])
def test_invalid_subset_table_identifier_has_a_clear_error(model_file, setting):
    model_file.write_text(f"subset_weights = {setting}\n" + model_file.read_text())
    with pytest.raises(SignatureError, match="subset_weights"):
        load_model(model_file)


@pytest.mark.parametrize("signature", ["insT", "stem3"])
@pytest.mark.parametrize("change", ["coordinates", "direction", "weight", "intercept"])
def test_subset_table_rejects_incompatible_custom_model(tmp_path, signature, change):
    content = (ROOT / "config" / f"{signature}.toml").read_text()
    if change == "coordinates":
        content = content.replace("chr20:41113969:41114155", "chr21:41113969:41114155")
        content = content.replace("chr11:64767566:64767844", "chr12:64767566:64767844")
    elif change == "direction":
        content = content.replace("direction = -1", "direction = 1", 1)
    elif change == "weight":
        content = content.replace("weight = ", "weight = 9", 1)
    else:
        content = content.replace("intercept = -", "intercept = -9", 1)
    path = tmp_path / "changed.toml"
    path.write_text(content)
    with pytest.raises(SignatureError, match="subset_weights"):
        load_model(path)


def test_copied_config_can_opt_out_to_customize_parameters(tmp_path):
    content = (ROOT / "config/insT.toml").read_text()
    content = content.replace('subset_weights = "insT"\n', "")
    content = content.replace("intercept = -8.33698995873104", "intercept = -3.5")
    path = tmp_path / "fixed.toml"
    path.write_text(content)
    model = load_model(path)
    assert model.subset_fits is None
    assert model.intercept == -3.5


@pytest.mark.parametrize("signature", ["insT", "stem3"])
def test_every_subset_scores_with_its_own_refit(signature):
    import csv
    import math

    from renu_splicing_signature.rmats import Counts, Dataset, Sample
    from renu_splicing_signature.scoring import score_samples

    model = load_model(ROOT / "config" / f"{signature}.toml")
    export = "insT" if signature == "insT" else "stemIII"
    path = ROOT / f"src/renu_splicing_signature/data/subset_refit_weights_{export}.csv"
    with path.open(newline="") as handle:
        fits = list(csv.DictReader(handle))
    samples = tuple(Sample(f"subset_{i}", 1, "case") for i in range(len(fits)))
    samples += (Sample("baseline", 2, "control"),)
    counts = {}
    for event in model.events:
        column = "w_" + event.id.split(":")[1]
        signal = Counts(4, 0) if event.direction == 1 else Counts(0, 4)
        counts[event.key] = tuple(
            signal if row[column] != "NA" else Counts(0, 0) for row in fits
        ) + (Counts(4, 4),)
    result = score_samples(model, Dataset(samples, counts))
    for score, fit in zip(result.rows, fits):
        coefficients = [
            float(v) for k, v in fit.items() if k.startswith("w_") and v != "NA"
        ]
        expected_weighted = math.log(9) * math.fsum(coefficients)
        expected_lp = float(fit["intercept"]) + expected_weighted
        assert score.n_events_used == len(coefficients)
        if float(fit["fpr_pct"]) >= 5:
            assert score.status == "high_fpr"
            assert (
                score.uniform_sum
                is score.weighted_sum
                is score.lasso_linpred
                is score.lasso_prob
                is None
            )
            continue
        assert score.uniform_sum == pytest.approx(len(coefficients) * math.log(9))
        assert score.weighted_sum == pytest.approx(expected_weighted)
        assert score.lasso_linpred == pytest.approx(expected_lp)
        assert score.lasso_prob == pytest.approx(1 / (1 + math.exp(-expected_lp)))
    assert result.rows[-1].weighted_sum == 0
    assert result.rows[-1].lasso_linpred == model.intercept
    for detail in result.details:
        if detail.sample == "baseline":
            continue
        fit = fits[int(detail.sample.removeprefix("subset_"))]
        expected = fit["w_" + detail.event.split(":")[1]]
        if expected == "NA":
            assert detail.weight is None
            assert detail.delta_logit_psi_direction_weighted is None
        elif float(fit["fpr_pct"]) >= 5:
            assert detail.status == "high_fpr"
            assert (
                detail.delta_logit_psi_direction
                is detail.delta_logit_psi_weighted
                is detail.delta_logit_psi_direction_weighted
                is None
            )
        else:
            assert detail.weight == float(expected)
            assert detail.delta_logit_psi_direction_weighted == pytest.approx(
                float(expected) * math.log(9)
            )
            assert detail.delta_logit_psi_weighted == pytest.approx(
                detail.direction * float(expected) * math.log(9)
            )


@pytest.mark.parametrize("signature", ["insT", "stem3"])
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
def test_subset_selection_follows_all_availability_checks(signature, reason):
    from dataclasses import replace

    from renu_splicing_signature.rmats import Counts, Dataset, Sample
    from renu_splicing_signature.scoring import score_samples

    model = load_model(ROOT / "config" / f"{signature}.toml")
    samples = (
        Sample("partial", 1, "case"),
        Sample("full", 1, "case"),
        Sample("control", 2, "control"),
    )
    counts = {
        event.key: (Counts(5, 1), Counts(5, 1), Counts(1, 5)) for event in model.events
    }
    event = model.events[0]
    if reason == "absent_from_file":
        del counts[event.key]
    elif reason == "no_control_baseline":
        counts[event.key] = (Counts(5, 1), Counts(5, 1), Counts(0, 0))
    else:
        pair = {
            "missing_counts": Counts(None, 4),
            "zero_coverage": Counts(0, 0),
            "below_min_read_depth": Counts(1, 1),
        }[reason]
        counts[event.key] = (pair, Counts(5, 1), Counts(1, 5))
        if reason == "below_min_read_depth":
            model = replace(model, min_read_depth=3)
    result = score_samples(model, Dataset(samples, counts))
    partial = result.rows[0]
    subset = frozenset(e.id for e in model.events[1:])
    fit = model.subset_fits[subset]
    assert partial.n_events_used == len(model.events) - 1
    assert partial.lasso_linpred - partial.weighted_sum == pytest.approx(fit.intercept)
    detail = next(
        d for d in result.details if d.sample == "partial" and d.event == event.id
    )
    assert detail.status == reason
    assert detail.weight is None
    if reason not in ("absent_from_file", "no_control_baseline"):
        assert result.rows[1].status == result.rows[2].status == "complete"
        assert result.rows[1].lasso_linpred - result.rows[
            1
        ].weighted_sum == pytest.approx(model.intercept)


@pytest.mark.parametrize("signature", ["insT", "stem3"])
def test_subset_model_empty_sample_remains_unavailable(signature):
    from renu_splicing_signature.rmats import Counts, Dataset, Sample
    from renu_splicing_signature.scoring import score_samples

    model = load_model(ROOT / "config" / f"{signature}.toml")
    samples = (Sample("empty", 1, "case"), Sample("control", 2, "control"))
    dataset = Dataset(
        samples, {e.key: (Counts(0, 0), Counts(4, 4)) for e in model.events}
    )
    result = score_samples(model, dataset)
    assert result.rows[0].status == "unavailable"
    assert result.rows[0].n_events_used == 0
    assert (
        result.rows[0].lasso_prob
        is result.rows[0].lasso_linpred
        is result.rows[0].weighted_sum
        is result.rows[0].uniform_sum
        is None
    )
    assert all(d.weight is None for d in result.details if d.sample == "empty")


def test_partial_controls_select_their_own_fits(write_rmats, rmats_rows):
    import math

    from renu_splicing_signature.rmats import read_rmats
    from renu_splicing_signature.scoring import score_samples

    model = load_model()
    # Case has both events, each control has a different singleton. Each event
    # has a log-ratio-zero baseline from its one usable control.
    rmats_rows[0][14:16] = ["4,0", "4,0"]
    rmats_rows[1][14:16] = ["0,4", "0,4"]
    result = score_samples(model, read_rmats(write_rmats(), model))
    case, ctdsp2_control, map4k4_control = result.rows
    expected = -1 * math.log(9) * (6.09810705946417 + 12.9708828765886)
    assert case.uniform_sum == pytest.approx(-2 * math.log(9))
    assert case.weighted_sum == pytest.approx(expected)
    assert case.lasso_linpred == pytest.approx(-8.82581854676958 + expected)
    assert ctdsp2_control.status == "high_fpr"
    assert ctdsp2_control.lasso_linpred is None
    assert map4k4_control.lasso_linpred == -3.2633453653868
    used = [d for d in result.details if d.status == "used"]
    assert [d.weight for d in used] == [
        6.09810705946417,
        12.9708828765886,
        5.36662803091663,
    ]


def test_missing_refit_is_an_error_instead_of_falling_back(write_rmats):
    from dataclasses import replace

    from renu_splicing_signature.rmats import read_rmats
    from renu_splicing_signature.scoring import score_samples

    model = replace(load_model(), subset_fits={})
    with pytest.raises(SignatureError, match="refit"):
        score_samples(model, read_rmats(write_rmats(), model))
