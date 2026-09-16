import pytest

from jwst.coron.hlsp_step import HlspStep


def test_hlsp_raises_deprecation_warning(tmp_path, target_image):
    """Test that calling HlspStep triggers a DeprecationWarning."""
    # pytest.warns captures the warning and checks the message
    with pytest.warns(DeprecationWarning, match="deprecated"):
        step = HlspStep()
        step.output_dir = str(tmp_path)
        step.run(target_image)


@pytest.mark.filterwarnings("ignore:'HlspStep' has been deprecated:DeprecationWarning")
def test_hlsp_step(tmp_path, target_image):
    # This step expects image models
    input_model = target_image
    input_model.meta.filename = "test.fits"

    result = HlspStep.call(input_model, output_dir=str(tmp_path))

    # This step does not return a product
    assert result is None

    # It creates two files, regardless of "save_results" setting
    assert (tmp_path / "test_snr.fits").exists()
    assert (tmp_path / "test_contrast.fits").exists()
