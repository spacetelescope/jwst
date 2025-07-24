from stdatamodels.jwst import datamodels

from jwst.coron.klip_step import KlipStep


def test_klip_step(target_model, psf_model):
    result = KlipStep.call(target_model, psf_model)
    assert result is not target_model
    assert result is not psf_model
    assert result.meta.cal_step.klip == "COMPLETE"
    assert isinstance(result, datamodels.CubeModel)
