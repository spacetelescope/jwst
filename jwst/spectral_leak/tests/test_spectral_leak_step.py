import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.spectral_leak import SpectralLeakStep


def make_mrs_multispec(base_model="MRSMultiSpecModel"):
    # make a model
    if "MRS" in base_model:
        multispec = datamodels.MRSMultiSpecModel()
        one_spec = datamodels.MRSSpecModel((10,))
    else:
        multispec = datamodels.MultiSpecModel()
        one_spec = datamodels.SpecModel((10,))

    # Assign the spec_table to convert it to a FITS record
    one_spec.spec_table = one_spec.spec_table

    # Add some placeholder values
    one_spec.spec_table["WAVELENGTH"] = np.linspace(5, 28, 10)
    one_spec.spec_table["FLUX"] = 0.1
    one_spec.spec_table["FLUX_ERROR"] = 0.01

    multispec.spec.append(one_spec)

    # Add metadata
    multispec.meta.instrument.name = "MIRI"
    multispec.meta.instrument.detector = "MIRIFULONG"
    multispec.meta.observation.date = "2024-01-01"
    multispec.meta.observation.time = "00:00:00"
    multispec.meta.exposure.type = "MIR_MRS"
    multispec.meta.instrument.channel = "12"
    multispec.meta.instrument.band = "MEDIUM"

    return multispec


@pytest.fixture()
def mrs_multispec():
    model = make_mrs_multispec()
    yield model
    model.close()


@pytest.fixture()
def mrs_multispec_container():
    spec1 = make_mrs_multispec()
    spec2 = make_mrs_multispec()
    spec2.meta.instrument.channel = "34"
    spec2.meta.instrument.band = "SHORT"
    container = ModelContainer([spec1, spec2])
    yield container
    container.close()


@pytest.fixture()
def multispec_container():
    spec1 = make_mrs_multispec(base_model="MultiSpecModel")
    spec2 = make_mrs_multispec(base_model="MultiSpecModel")
    spec2.meta.instrument.channel = "34"
    spec2.meta.instrument.band = "SHORT"
    container = ModelContainer([spec1, spec2])
    yield container
    container.close()


@pytest.mark.parametrize("dataset", ["mrs_multispec_container", "multispec_container"])
@pytest.mark.parametrize("multiple", [True, False])
def test_step_complete(request, dataset, multiple):
    container = request.getfixturevalue(dataset)
    if multiple:
        for model in container:
            model.meta.instrument.band = "MULTIPLE"

    result = SpectralLeakStep.call(container)

    for model, input_model in zip(result, container, strict=True):
        # Step is complete, but marked only in the channel 3 data
        if "3" in model.meta.instrument.channel:
            assert model.meta.cal_step.spectral_leak == "COMPLETE"
            assert not np.allclose(
                model.spec[0].spec_table["FLUX"], input_model.spec[0].spec_table["FLUX"]
            )
        else:
            assert model.meta.cal_step.spectral_leak is None
            np.testing.assert_allclose(
                model.spec[0].spec_table["FLUX"], input_model.spec[0].spec_table["FLUX"]
            )

        # Input is not modified
        assert model is not input_model
        assert input_model.meta.cal_step.spectral_leak is None


def test_step_skip_wrong_channel(caplog, mrs_multispec):
    result = SpectralLeakStep.call([mrs_multispec])

    # Step is skipped
    assert result[0].meta.cal_step.spectral_leak == "SKIPPED"
    assert "CH1B and CH3A were not found" in caplog.text

    # Input is not modified
    assert result[0] is not mrs_multispec
    assert mrs_multispec.meta.cal_step.spectral_leak is None


def test_step_skip_not_container(caplog, mrs_multispec):
    result = SpectralLeakStep.call(mrs_multispec)

    # Step is skipped
    assert result.meta.cal_step.spectral_leak == "SKIPPED"
    assert "is not a ModelContainer" in caplog.text

    # Input is not modified
    assert result is not mrs_multispec
    assert mrs_multispec.meta.cal_step.spectral_leak is None


def test_step_skip_wrong_type(caplog, mrs_multispec):
    model = datamodels.ImageModel()
    model.update(mrs_multispec)
    result = SpectralLeakStep.call([model])

    # Step is skipped
    assert result[0].meta.cal_step.spectral_leak == "SKIPPED"
    assert "is not an extracted spectrum" in caplog.text

    # Input is not modified
    assert result[0] is not model
    assert model.meta.cal_step.spectral_leak is None


def test_step_skip_extended(caplog, mrs_multispec):
    mrs_multispec.spec[0].source_type = "EXTENDED"
    result = SpectralLeakStep.call([mrs_multispec])

    # Step is skipped
    assert result[0].meta.cal_step.spectral_leak == "SKIPPED"
    assert "No spectral leak correction for extended source data" in caplog.text

    # Input is not modified
    assert result[0] is not mrs_multispec
    assert mrs_multispec.meta.cal_step.spectral_leak is None
