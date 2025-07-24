import numpy as np

from jwst.coron.align_refs_step import AlignRefsStep


def test_align_refs_no_shift(target_model):
    # Align the target to itself (no shift)
    result = AlignRefsStep.call(target_model, target_model)

    # Successful completion
    assert result.meta.cal_step.align_psfs == "COMPLETE"
    assert result is not target_model

    # No significant change to result
    np.testing.assert_allclose(result.data, target_model.data, atol=1e-5)


def test_align_refs_with_shift(target_model, psf_model):
    result = AlignRefsStep.call(target_model, psf_model)

    # Successful completion
    assert result.meta.cal_step.align_psfs == "COMPLETE"
    assert result is not target_model
    assert result is not psf_model

    # Result is shifted
    assert not np.allclose(result.data, target_model.data, atol=1e-5)


def test_no_psf_mask(monkeypatch, target_model, psf_model):
    step = AlignRefsStep()
    monkeypatch.setattr(step, "get_reference_file", lambda *args: "N/A")

    result = step.run(target_model, psf_model)

    # Step is skipped
    assert result.meta.cal_step.align_psfs == "SKIPPED"
    assert result is not target_model
    assert result is not psf_model

    # No change to result
    np.testing.assert_array_equal(result.data, psf_model.data)


def test_no_bad_bit(caplog, target_model, psf_model):
    AlignRefsStep.call(target_model, psf_model, bad_bits=None)
    # Note: this is a debug message.
    # The "self.log" logger for steps is always set to DEBUG level,
    # so it's available in the caplog, regardless of default logcfg.
    assert "No bad bits provided; treating all pixels as good" in caplog.text
