import numpy as np
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_array_equal
from stdatamodels.jwst.datamodels import ImageModel, dqflags
from stdatamodels.jwst.transforms.models import Slit

from jwst.assign_wcs import AssignWcsStep
from jwst.msaflagopen.msaflag_open import (
    boundingbox_to_indices,
    create_slitlets,
    get_failed_open_shutters,
    wcs_to_dq,
)
from jwst.msaflagopen import MSAFlagOpenStep
from jwst.stpipe import Step


MSA_FAILED_OPEN = dqflags.pixel["MSA_FAILED_OPEN"]


def make_nirspec_mos_model():
    im = ImageModel((2048, 2048))
    im.meta.wcsinfo = {
        "dec_ref": -0.00601415671349804,
        "ra_ref": -0.02073605215697509,
        "roll_ref": -0.0,
        "v2_ref": -453.5134,
        "v3_ref": -373.4826,
        "v3yangle": 0.0,
        "vparity": -1,
    }

    im.meta.instrument = {
        "detector": "NRS1",
        "filter": "F100LP",
        "grating": "G140M",
        "name": "NIRSPEC",
        "gwa_tilt": 37.0610,
        "gwa_xtilt": 0.0001,
        "gwa_ytilt": 0.0001,
        "msa_metadata_id": 12,
    }

    im.meta.observation = {"program_number": "1234", "date": "2016-09-05", "time": "8:59:37"}

    im.meta.exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 100,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "readpatt": "NRSRAPID",
        "sample_time": 10.0,
        "start_time": 58119.8333,
        "type": "NRS_MSASPEC",
        "zero_frame": False,
    }
    im.meta.instrument.msa_metadata_file = get_pkg_data_filename(
        "data/msa_configuration.fits", package="jwst.assign_wcs.tests"
    )
    im.meta.dither.position_number = 1
    return im


def test_get_failed_open_shutters():
    """test that failed open shutters are returned from reference file"""

    # Set up data model to retrieve reference file
    dm = ImageModel()
    dm.meta.instrument.name = "NIRSPEC"
    dm.meta.observation.date = "2016-09-05"
    dm.meta.observation.time = "8:59:37"

    # Get reference file and return all failed open shutters
    msa_oper = Step().get_reference_file(dm, "msaoper")
    result = get_failed_open_shutters(msa_oper)

    # get_failed_open_shutters returns 3 flaggable states
    # state, Internal state, and TA state.
    for shutter in result:
        assert (
            shutter["state"] == "open"
            or shutter["Internal state"] == "open"
            or shutter["TA state"] == "open"
        )


def test_create_slitlets():
    """Test that slitlets are Slit type and have all the necessary fields"""

    dm = ImageModel()
    dm.meta.instrument.name = "NIRSPEC"
    dm.meta.observation.date = "2016-09-05"
    dm.meta.observation.time = "8:59:37"
    msa_oper = Step().get_reference_file(dm, "msaoper")
    result = create_slitlets(msa_oper)

    for slit in result:
        # Test the returned data type
        assert isinstance(slit, Slit)


def test_wcs_to_dq():
    """Test that non nan values are assigned the values of flags"""

    # Make data array
    grid = np.zeros((10, 10))
    wcs_array = np.array((grid, grid))

    # Put in some nans randomly in the data.
    wcs_array.ravel()[np.random.choice(wcs_array.size, 10, replace=False)] = np.nan
    nans = np.isnan(wcs_array[0])

    result = wcs_to_dq(wcs_array, MSA_FAILED_OPEN)

    # wcs_to_dq create an array of zeros and if nans are present they
    # will have value zero where are non nan elements will have the value
    # of FLAG.
    assert_array_equal(result[nans], 0)
    assert_array_equal(result[~nans], MSA_FAILED_OPEN)


def test_boundingbox_from_indices():
    dm = ImageModel((10, 10))
    bbox = ((1, 2), (3, 4))

    result = boundingbox_to_indices(dm, bbox)

    assert result == (1, 3, 3, 5)


def test_msaflagopen_step():
    im = make_nirspec_mos_model()
    im = AssignWcsStep.call(im)
    result = MSAFlagOpenStep.call(im)

    nonzero = np.nonzero(result.dq)
    assert_array_equal(result.dq[nonzero], MSA_FAILED_OPEN)


def test_no_ref_file():
    im = ImageModel((10, 10))
    im.meta.instrument.name = "NIRCAM"

    # the step is skipped if there is no msaoper reference file
    result = MSAFlagOpenStep.call(im)
    assert result.meta.cal_step.msa_flagging == "SKIPPED"


def test_custom_ref_file():
    im = make_nirspec_mos_model()
    wavelength_range = get_pkg_data_filename("data/waverange.asdf", package="jwst.assign_wcs.tests")
    im = AssignWcsStep.call(im, override_wavelengthrange=wavelength_range)
    result = MSAFlagOpenStep.call(im)

    nonzero = np.nonzero(result.dq)
    assert_array_equal(result.dq[nonzero], MSA_FAILED_OPEN)
