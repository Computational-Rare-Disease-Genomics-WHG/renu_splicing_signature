"""False-positive checks apply to each sample's selected subset refit."""

import csv
from pathlib import Path

import pytest

from renu_splicing_signature.model import SignatureError, load_model

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "signature,gene,n_fp,fpr_pct",
    [
        ("insT", "STX16", 21, 0.35140562248996),
        ("insT", "TOP1", 5339, 94.3452906874006),
        ("stem3", "MAP3K4", 341, 5.70615796519411),
    ],
)
def test_selected_fit_retains_exported_false_positive_statistics(
    signature, gene, n_fp, fpr_pct
):
    model = load_model(ROOT / "config" / f"{signature}.toml")
    event = next(e for e in model.events if e.id.split(":")[1] == gene)
    fit = model.subset_fits[frozenset([event.id])]
    assert fit.n_fp == n_fp
    assert fit.fpr_pct == fpr_pct
    full = model.subset_fits[frozenset(e.id for e in model.events)]
    assert full.n_fp == full.fpr_pct == 0


@pytest.mark.parametrize(
    "column,value",
    [
        ("n_fp", "-1"),
        ("n_fp", "1.5"),
        ("n_fp", "NA"),
        ("fpr_pct", "-0.1"),
        ("fpr_pct", "101"),
        ("fpr_pct", "nan"),
        ("fpr_pct", "inf"),
        ("fpr_pct", "NA"),
    ],
)
def test_invalid_false_positive_statistics_are_rejected(
    tmp_path, monkeypatch, column, value
):
    from renu_splicing_signature import subsets

    data = tmp_path / "data"
    data.mkdir()
    (data / "stem3.toml").write_bytes((ROOT / "config/stem3.toml").read_bytes())
    table = ROOT / "src/renu_splicing_signature/data/subset_refit_weights_stemIII.csv"
    with table.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row.setdefault("n_fp", "0")
        row.setdefault("fpr_pct", "0")
    rows[0][column] = value
    with (data / table.name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    monkeypatch.setattr(subsets, "files", lambda package: tmp_path)
    with pytest.raises(SignatureError, match="n_fp|fpr_pct"):
        load_model(ROOT / "config/stem3.toml")


@pytest.fixture
def mixed_subset_data():
    from renu_splicing_signature.rmats import Counts, Dataset, Sample

    model = load_model()
    samples = (
        Sample("risky_case", 1, "case"),
        Sample("full_case", 1, "case"),
        Sample("risky_control", 2, "control"),
        Sample("full_control", 2, "control"),
    )
    counts = {}
    for event in model.events:
        signal = Counts(4, 0) if event.direction == 1 else Counts(0, 4)
        top1 = event.id.split(":")[1] == "TOP1"
        counts[event.key] = (
            signal if top1 else Counts(0, 0),
            signal,
            signal if top1 else Counts(0, 0),
            Counts(4, 4),
        )
    return model, Dataset(samples, counts)


@pytest.mark.parametrize(
    "fpr_pct,rejected,display_rate",
    [
        (0, False, "0.00%"),
        (0.001, False, "<0.01%"),
        (4.999999, False, "<5%"),
        (5, True, "5.00%"),
        (5.000001, True, "5.00%"),
    ],
)
def test_fpr_boundary_is_inclusive_and_applies_per_sample(
    mixed_subset_data, fpr_pct, rejected, display_rate
):
    from dataclasses import replace
    import math

    from renu_splicing_signature.scoring import score_samples

    model, dataset = mixed_subset_data
    top1 = next(e for e in model.events if e.id.split(":")[1] == "TOP1")
    key = frozenset([top1.id])
    fits = dict(model.subset_fits)
    fits[key] = replace(fits[key], n_fp=1, fpr_pct=fpr_pct)
    result = score_samples(replace(model, subset_fits=fits), dataset)
    for row in result.rows:
        if row.sample.startswith("risky_"):
            assert row.status == ("high_fpr" if rejected else "partial")
            assert row.n_events_used == 1
            if rejected:
                assert (
                    row.uniform_sum
                    is row.weighted_sum
                    is row.lasso_linpred
                    is row.lasso_prob
                    is None
                )
            else:
                assert row.lasso_prob is not None
        else:
            assert row.status == "complete"
            assert row.lasso_prob is not None
    # The rejected control still provides a usable TOP1 baseline: model quality
    # does not change its observed counts or other samples' centering.
    top1_details = [d for d in result.details if d.event == top1.id]
    assert all(d.control_mean == pytest.approx(-math.log(9) / 2) for d in top1_details)
    warnings = [w for w in result.warnings if "control sample" in w]
    assert len(warnings) == 2
    for name in ("risky_case", "risky_control"):
        message = next(w for w in warnings if name in w)
        assert "incorrectly flagged 1 control sample as a case" in message
        assert display_rate in message
        assert "Usable junctions (by gene): TOP1." in message
        assert "A5SS:" not in message
        if rejected:
            assert "scores not provided" in message
            assert "at or above the 5% limit" in message
        else:
            assert "scores provided with caution" in message
    for detail in top1_details:
        if detail.sample.startswith("risky_") and rejected:
            assert detail.status == "high_fpr"
            assert (
                detail.delta_logit_psi_weighted
                is detail.delta_logit_psi_direction
                is detail.delta_logit_psi_direction_weighted
                is None
            )
            assert detail.log_ratio is not None and detail.ijc is not None


def test_zero_false_positives_do_not_add_subset_warnings(mixed_subset_data):
    from dataclasses import replace

    from renu_splicing_signature.scoring import score_samples

    model, dataset = mixed_subset_data
    fits = {
        key: replace(fit, n_fp=0, fpr_pct=0) for key, fit in model.subset_fits.items()
    }
    result = score_samples(replace(model, subset_fits=fits), dataset)
    assert not any("control sample" in w for w in result.warnings)
    assert all(row.lasso_prob is not None for row in result.rows)


def test_high_fpr_warns_even_if_false_positive_count_is_zero(mixed_subset_data):
    from dataclasses import replace

    from renu_splicing_signature.scoring import score_samples

    model, dataset = mixed_subset_data
    fits = {
        key: replace(fit, n_fp=0, fpr_pct=5) for key, fit in model.subset_fits.items()
    }
    result = score_samples(replace(model, subset_fits=fits), dataset)
    assert all(
        row.status == "high_fpr" and row.lasso_prob is None for row in result.rows
    )
    assert len([w for w in result.warnings if "scores not provided" in w]) == 4


def test_subset_warning_lists_available_zero_weight_junctions():
    from renu_splicing_signature.rmats import Counts, Dataset, Sample
    from renu_splicing_signature.scoring import score_samples

    model = load_model()
    selected = [e for e in model.events if e.id.split(":")[1] in ("TOP1", "KLHDC4")]
    samples = (
        Sample("zero_coefficient_case", 1, "case"),
        Sample("control", 2, "control"),
    )
    dataset = Dataset(samples, {e.key: (Counts(4, 1), Counts(1, 4)) for e in selected})
    result = score_samples(model, dataset)
    message = next(
        w
        for w in result.warnings
        if "zero_coefficient_case" in w and "control sample" in w
    )
    assert "Usable junctions (by gene): TOP1, KLHDC4." in message
    assert "5,339 control samples as cases" in message
    assert "94.35%" in message
    assert result.rows[0].n_events_used == 2
    assert result.rows[0].status == "high_fpr"
