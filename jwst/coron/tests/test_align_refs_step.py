import numpy as np
import pytest

from jwst.coron.align_refs_step import AlignRefsStep


@pytest.mark.parametrize("dataset", ["target_model", "target_model_miri"])
def test_align_refs_no_shift(request, dataset):
    target = request.getfixturevalue(dataset)

    # Align the target to itself (no shift)
    result = AlignRefsStep.call(target, target)

    # Successful completion
    assert result.meta.cal_step.align_psfs == "COMPLETE"

    # Input is not modified
    assert result is not target
    assert target.meta.cal_step.align_psfs is None

    # No significant change to result
    np.testing.assert_allclose(result.data, target.data, atol=1e-5)


@pytest.mark.parametrize(
    "dataset", [("target_model", "psf_model"), ("target_model_miri", "psf_model_miri")]
)
def test_align_refs_with_shift(request, dataset):
    target = request.getfixturevalue(dataset[0])
    psf = request.getfixturevalue(dataset[1])
    result = AlignRefsStep.call(target, psf)

    # Successful completion
    assert result.meta.cal_step.align_psfs == "COMPLETE"

    # Input is not modified
    assert result is not target
    assert result is not psf
    assert target.meta.cal_step.align_psfs is None
    assert psf.meta.cal_step.align_psfs is None

    # Result is shifted
    assert not np.allclose(result.data, psf.data, atol=1e-5)


def test_no_psf_mask(monkeypatch, target_model, psf_model):
    step = AlignRefsStep()
    monkeypatch.setattr(step, "get_reference_file", lambda *args: "N/A")

    result = step.run(target_model, psf_model)

    # Step is skipped
    assert result.meta.cal_step.align_psfs == "SKIPPED"

    # Input is not modified
    assert result is not target_model
    assert result is not psf_model
    assert target_model.meta.cal_step.align_psfs is None
    assert psf_model.meta.cal_step.align_psfs is None

    # No change to result
    np.testing.assert_array_equal(result.data, psf_model.data)


def test_no_bad_bit(caplog, target_model, psf_model):
    AlignRefsStep.call(target_model, psf_model, bad_bits=None)
    # Note: this is a debug message.
    # The "self.log" logger for steps is always set to DEBUG level,
    # so it's available in the caplog, regardless of default logcfg.
    assert "No bad bits provided; treating all pixels as good" in caplog.text
