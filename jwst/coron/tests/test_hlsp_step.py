from stdatamodels.jwst import datamodels

from jwst.coron.hlsp_step import HlspStep


def test_hlsp_step(tmp_path, target_model):
    # This step expects image models
    input_model = datamodels.ImageModel()
    input_model.data = target_model.data[0]
    input_model.err = target_model.err[0]
    input_model.dq = target_model.dq[0]
    input_model.update(target_model)
    input_model.meta.filename = "test.fits"

    result = HlspStep.call(input_model, output_dir=str(tmp_path))

    # This step does not return a product
    assert result is None

    # It creates two files
    assert (tmp_path / "test_snr.fits").exists()
    assert (tmp_path / "test_contrast.fits").exists()
