from jwst.coron.hlsp_step import HlspStep


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
