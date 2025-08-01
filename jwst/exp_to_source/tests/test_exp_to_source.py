import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.datamodels import ModelContainer
from jwst.exp_to_source import exp_to_source, multislit_to_container


@pytest.fixture
def run_exp_to_source(mock_input):
    outputs = exp_to_source(mock_input)
    return mock_input, outputs


def test_model_structure(run_exp_to_source):
    inputs, outputs = run_exp_to_source
    assert len(outputs) == 5
    assert len(outputs["3"].exposures) == 3
    for in_model in inputs:
        for slit in in_model.slits:
            expected_n_exposures = 3 - abs(3 - slit.source_id)
            exposures = outputs[str(slit.source_id)].exposures
            assert len(exposures) == expected_n_exposures
            fnames = {exp.meta.filename for exp in exposures}
            assert in_model.meta.filename in fnames
            for exposure in exposures:
                assert (exposure.data == slit.data).all()
                assert len(exposure.meta._instance) >= len(in_model.meta._instance)

            assert outputs[str(slit.source_id)].meta.filename != in_model.meta.filename


def test_model_roundtrip(tmp_path, run_exp_to_source):
    _, outputs = run_exp_to_source
    files = []
    for output in outputs:
        file_path = tmp_path / (output + ".fits")
        outputs[output].save(file_path)
        files.append(file_path)
    for file_path in files:
        expected_n_exposures = 3 - abs(3 - int(file_path.stem.split(".")[0]))
        with dm.MultiExposureModel(file_path) as multiexposure_model:
            assert len(multiexposure_model.exposures) == expected_n_exposures
            exp_files = set()
            for exposure in multiexposure_model.exposures:
                exp_files.add(exposure.meta.filename)
            assert len(exp_files) == len(multiexposure_model.exposures)
            assert multiexposure_model.meta.filename not in exp_files


def test_container_structure(mock_input):
    """Test for container usage."""

    container = ModelContainer(mock_input)

    # Make the source-based containers
    outputs = multislit_to_container(container)

    # See what we got.
    assert len(container) == 3
    assert len(outputs) == 5
    for model in container:
        for slit in model.slits:
            expected_n_exposures = 3 - abs(3 - slit.source_id)
            exposures = outputs[str(slit.source_id)]
            assert len(exposures) == expected_n_exposures
            fnames = {exp.meta.filename for exp in exposures}
            for exposure in exposures:
                assert (exposure.data == slit.data).all()
                assert np.array_equal(exposure.data, slit.data)
                assert model.meta.filename in fnames

    # Closeout
    container.close()
    for model in mock_input:
        model.close()
    for model in outputs.values():
        model.close()


def test_slit_exptype(mock_input):
    """Test for slit exposure type handling."""

    container = ModelContainer(mock_input)

    # Add a slit exposure type to each input
    for model in container:
        for slit in model.slits:
            if slit.source_id == 1:
                slit.meta.exposure = {"type": "NRS_MSASPEC"}
            else:
                slit.meta.exposure = {"type": "NRS_FIXEDSLIT"}

    # Make the source-based containers
    outputs = multislit_to_container(container)

    # Check that exposure type was passed from input to output
    assert len(container) == 3
    assert len(outputs) == 5
    for model in container:
        for slit in model.slits:
            exposures = outputs[str(slit.source_id)]
            for exposure in exposures:
                assert exposure.meta.exposure.type == slit.meta.exposure.type
                if slit.source_id == 1:
                    assert exposure.meta.exposure.type == "NRS_MSASPEC"
                else:
                    assert exposure.meta.exposure.type == "NRS_FIXEDSLIT"

    # Closeout
    container.close()
    for model in mock_input:
        model.close()
    for model in outputs.values():
        model.close()
