#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pytest

from stpipe.config_parser import ValidationError

from jwst.datamodels import RampModel
from jwst.datamodels.mask import MaskModel

from jwst.dq_init import DQInitStep


# Tests derived from example code from Jira JP-345
# Tests are for override functionality only,  not DQ init.


def create_models():
    """Returns  RampModel(), MaskModel()."""
    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rawramp(nints, ngroups, ysize, xsize)

    # create a MaskModel to create baseline arrays for
    # dq and dq_def to be modified afterwards to create
    # a test case mask model.
    dq, dq_def = make_maskmodel(ysize, xsize)

    # Construct a MaskModel on the fly for use as override.
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = 'MIRI'
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = ysize

    return dm_ramp, ref_data


def make_maskmodel(ysize, xsize):
    # create a mask model for the dq_init step
    csize = (ysize, xsize)
    dq = np.zeros(csize, dtype=int)
    # how do we define a dq_def extension?
    mask = MaskModel()

    dqdef = [(0, 1, 'DO_NOT_USE', 'Bad Pixel do not use'),
             (1, 2, 'DEAD', 'Dead Pixel'),
             (2, 4, 'HOT', 'Hot pixel'),
             (3, 8, 'UNRELIABLE_SLOPE', 'Large slope variance'),
             (4, 16, 'RC', 'RC pixel'),
             (5, 32, 'REFERENCE_PIXEL', 'Reference Pixel')]

    dq_def = np.array((dqdef), dtype=mask.dq_def.dtype)

    return dq, dq_def


def make_rawramp(nints, ngroups, ysize, xsize):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data)
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


def test_invalid_override():
    """Test that a bogus override type is caught."""
    dm_ramp, ref_data = create_models()

    with pytest.raises(ValidationError):
        DQInitStep(override_mask=DQInitStep)


def test_valid_model_override():
    dm_ramp, ref_data = create_models()

    step = DQInitStep(override_mask=ref_data)

    # Verify get_reference_file() returns an override model.
    fetched_reference = step.get_reference_file(dm_ramp, 'mask')
    assert isinstance(fetched_reference, MaskModel), \
        "get_reference_file() should return a model for this override."

    # Verify no exceptions occur during DQ processing.
    step.process(dm_ramp)


def test_string_override():
    dm_ramp, ref_data = create_models()

    step = DQInitStep(override_mask="some_file.fits")

    # Verify stpipe treats string as filename and attempts to open
    with pytest.raises(FileNotFoundError):
        step.get_reference_file(dm_ramp, 'mask')
