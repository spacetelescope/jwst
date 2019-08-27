"""
Unit tests for group_scale correction
"""

from jwst.datamodels import RampModel
from jwst.group_scale.group_scale import do_correction
from jwst.group_scale import GroupScaleStep
import numpy as np
import pytest


def test_nframes_or_frame_divisor_is_none(make_rampmodel):
    """If nframes or frame_divisor is None, skip correction.
    Here I am just setting nframes=None
    """
    datmod = make_rampmodel(2, None, 4, 2048, 2048)
    output = do_correction(datmod)

    assert(output.meta.cal_step.group_scale == 'SKIPPED')


def test_nframes_equal_frame_divisor(make_rampmodel):
    """If nframes and frame_divisor are equal, skip correction
    """
    datmod = make_rampmodel(2, 4, 4, 2048, 2048)
    output = GroupScaleStep.call(datmod)
    print(output.meta.exposure.frame_divisor,
          output.meta.exposure.nframes)
    assert(output.meta.cal_step.group_scale == 'SKIPPED')


def test_nframes_not_equal_frame_divisor(make_rampmodel):
    """If nframes and frame_divisor are not equal, do correction
    """
    datmod = make_rampmodel(2, 2, 4, 2048, 2048)
    output = GroupScaleStep.call(datmod)

    # Assert that the step completed
    assert(output.meta.cal_step.group_scale == 'COMPLETE')

    # This assertion doesn't verify for correct output,
    # it just checks that the correction ran and that the data array
    # outputs are different than the inputs as requested in the document.
    assert not np.array_equal(output.data, datmod.data)


def test_nframes_is_none(make_rampmodel):
    """Make sure step is skipped if nframes is None
    """
    datmod = make_rampmodel(2, None, 4, 2048, 2048)
    output = GroupScaleStep.call(datmod)

    assert(output.meta.cal_step.group_scale == 'SKIPPED')


def test_nframes_is_power_of_two(make_rampmodel):
    """When frame_divisor is None, the correction will skip if
    nframes is a power of 2.
    """
    datmod = make_rampmodel(2, 4, None, 2048, 2048)
    output = GroupScaleStep.call(datmod)

    assert(output.meta.cal_step.group_scale == 'SKIPPED')


def test_nframes_is_not_power_of_two(make_rampmodel):
    """When frame_divisor is None, then do_correction will be applied if
    nframes is not a power of two but will be skip because if nframes or
    frame_divisor is none.
    """
    datmod = make_rampmodel(2, 3, None, 2048, 2048)
    output = GroupScaleStep.call(datmod)

    assert(output.meta.cal_step.group_scale == 'SKIPPED')

def test_scale_value(make_rampmodel):
    """Compare the ratio of the FRMDIVSR/NFRAMES from the data model input and
    compare to the output of the pipeline.
    """

    datmod = make_rampmodel(2, 2, 4, 2048, 2048)

    # Calculate the scale based off of the input.
    scale = datmod.meta.exposure.frame_divisor/datmod.meta.exposure.nframes

    output = GroupScaleStep.call(datmod)

    scale_from_data = np.unique(output.data/datmod.data)

    # Since the scale value is applied uniformly to the array, if we divide the output
    # by the input then we should get a single unique value (ie the scale) calculated
    # by the pipeline.
    assert(len(scale_from_data) == 1)

    # Make sure the scale calculated manually from the data model aboved matched what the
    # pipeline calculated.
    assert(scale == scale_from_data[0])


@pytest.fixture(scope='function')
def make_rampmodel():
    '''Make NIRSPEC IRS2 model for testing'''

    # NRS1 and NRS2 are size  2048x2048 pixels
    def _dm(ngroups, nframes, frame_divisor, ysize, xsize):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, ysize, xsize)
        data = np.random.randint(low=1, high=50, size=csize)

        # create a JWST datamodel for NIRSPEC data
        dm = RampModel(data=data)

        dm.meta.instrument.name = 'NIRSPEC'
        dm.meta.date = '2018-01-01'
        dm.meta.description = 'Fake data'
        dm.meta.reftype = 'RampModel'
        dm.meta.author = 'Mees'
        dm.meta.pedigree = 'Dummy'
        dm.meta.useafter = '2015-10-01T00:00:00'
        dm.meta.readpatt = 'NRS_IRS2'

        # Begin NIRSPEC specific for tests.
        dm.meta.exposure.frame_divisor = frame_divisor
        dm.meta.exposure.nframes = nframes

        return dm

    return _dm
