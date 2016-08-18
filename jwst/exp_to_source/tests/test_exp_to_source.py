from glob import iglob
import os

from ...datamodels import MultiSlitModel
from ..exp_to_source import exp_to_source

INPUT_FILES = 'data/jwst_nod*_cal.fits'


def test_exp_to_source():
    inputs = [
        MultiSlitModel(f)
        for f in iglob(t_path(INPUT_FILES))
    ]
    outputs = exp_to_source(inputs)
    assert len(outputs) == 5
    assert len(outputs[next(iter(outputs))].exposures) == 3
    for in_idx, in_model in enumerate(inputs):
        for slit in in_model.slits:
            exposure = outputs[slit.name].exposures[in_idx]
            assert (exposure.data == slit.data).all()
            assert len(exposure.meta._instance) == len(in_model.meta._instance)
            assert exposure.meta.filename == in_model.meta.filename
            assert outputs[slit.name].meta.filename != in_model.meta.filename


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)
