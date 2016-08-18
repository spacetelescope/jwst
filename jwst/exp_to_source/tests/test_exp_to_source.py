from glob import iglob
import os
import pytest
import tempfile

from ...datamodels import (MultiExposureModel, MultiSlitModel)
from ..exp_to_source import exp_to_source

INPUT_FILES = 'data/jwst_nod*_cal.fits'


@pytest.fixture
def run_exp_to_source(scope='module'):
    inputs = [
        MultiSlitModel(f)
        for f in iglob(t_path(INPUT_FILES))
    ]
    outputs = exp_to_source(inputs)
    return inputs, outputs


def test_model_structure(run_exp_to_source):
    inputs, outputs = run_exp_to_source
    assert len(outputs) == 5
    assert len(outputs[next(iter(outputs))].exposures) == 3
    for in_idx, in_model in enumerate(inputs):
        for slit in in_model.slits:
            exposure = outputs[slit.name].exposures[in_idx]
            assert (exposure.data == slit.data).all()
            assert len(exposure.meta._instance) == len(in_model.meta._instance)
            assert exposure.meta.filename == in_model.meta.filename
            assert outputs[slit.name].meta.filename != in_model.meta.filename


def test_model_roundtrip(run_exp_to_source):
    inputs, outputs = run_exp_to_source
    files = []
    with tempfile.TemporaryDirectory() as path:
        for output in outputs:
            file_path = os.path.join(path, output) + '.fits'
            outputs[output].save(file_path)
            files.append(file_path)
        for file_path in files:
            multiexposure_model = MultiExposureModel(file_path)
            assert len(multiexposure_model.exposures) == 3
            exp_files = set()
            for exposure in multiexposure_model.exposures:
                exp_files.add(exposure.meta.filename)
            assert len(exp_files) == len(multiexposure_model.exposures)
            assert multiexposure_model.meta.filename not in exp_files


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)
