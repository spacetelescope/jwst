import numpy as np
import os

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel
from jwst.msaflagopen.msaflag_open import (boundingbox_to_indices,
                                           create_slitlets,
                                           flag,
                                           get_failed_open_shutters,
                                           id_from_xy,
                                           or_subarray_with_array,
                                           wcs_to_dq)
from jwst.msaflagopen.msaflagopen_step import create_reference_filename_dictionary
from jwst.assign_wcs.tests import data
from jwst.transforms.models import Slit
from jwst.stpipe.step import Step

data_path = os.path.split(os.path.abspath(data.__file__))[0]

def get_file_path(filename):
    """Construct an absolute path."""
    return os.path.join(data_path, filename)

def test_or_subarray_with_array():
    """test bitwise or with array and subarray."""

    dq_array = np.ones((20, 20), dtype=int) * 1
    dq_subarray = np.ones((10, 10), dtype=int) * 2

    xmin, xmax = 5, 15
    ymin, ymax = 5, 15

    # Uses numpy.bitwise_or to compute the array
    result = or_subarray_with_array(dq_array, dq_subarray, xmin, xmax, ymin, ymax)

    # Use python built in bitwise or | to verify. Just overwrite the
    # array since we wont be using it anymore in this test.
    dq_array[ymin:ymax, xmin:xmax] = dq_array[ymin:ymax, xmin:xmax] | dq_subarray

    assert np.array_equal(result, dq_array)

def test_id_from_xy():
    """Test id from x y location of shutter"""

    # First row of msaoper.json
    data = {"Q": 1,"x": 1,"y": 1,"state": "closed","TA state": "closed",
            "Internal state": "normal",
            "Vignetted": "yes"}

    shutters_per_row = 365

    x = data['x']
    y = data['y']

    result = id_from_xy(x, y)

    assert (x + (y-1)*shutters_per_row == result)

def test_get_failed_open_shutters():
    """test that failed open shutters are returned from reference file"""

    # Set up data model to retrieve reference file
    dm = ImageModel()
    dm.meta.instrument.name = 'NIRSPEC'
    dm.meta.observation.date = '2016-09-05'
    dm.meta.observation.time = '8:59:37'
    
    # Get reference file and return all failed open shutters
    msa_oper = Step().get_reference_file(dm, 'msaoper')
    result = get_failed_open_shutters(msa_oper)

    # get_failed_open_shutters returns 3 flaggable states
    # state, Internal state, and TA state.
    for shutter in result:
        assert (shutter['state'] == 'open' or 
                shutter['Internal state'] == 'open' or
                shutter['TA state'] == 'open')

def test_create_slitlets():
    """Test that slitlets are Slit type and have all the necessary fields"""

    dm = ImageModel()
    dm.meta.instrument.name = 'NIRSPEC'
    dm.meta.observation.date = '2016-09-05'
    dm.meta.observation.time = '8:59:37'
    msa_oper = Step().get_reference_file(dm, 'msaoper')
    result = create_slitlets(dm, msa_oper)

    slit_fields = ('name','shutter_id','dither_position','xcen',
                   'ycen','ymin','ymax','quadrant','source_id',
                   'shutter_state','source_name','source_alias',
                   'stellarity','source_xpos','source_ypos')

    for slit in result:
        # Test the returned data type and fields.
        assert type(slit) == Slit
        assert slit._fields == slit_fields

def test_wcs_to_dq():
    """Test that non nan values are assigned the values of flags"""

    # Make data array
    wcs_array = np.random.randn(100,100)

    # Put in some nans randomly in the data.
    wcs_array.ravel()[np.random.choice(wcs_array.size, 100, replace=False)] = np.nan

    FLAG = 1024
    result = wcs_to_dq(wcs_array, FLAG)

    nans = np.where(np.isnan(wcs_array[0]))
    non_nan = np.where(~np.isnan(wcs_array[0]))
    # wcs_to_dq create an array of zeros and if nans are present they
    # will have value zero where are non nan elements will have the value
    # of FLAG.
    assert (all(x == 0 for x in result[nans]))
    assert (all(x == FLAG for x in result[non_nan]))

def test_boundingbox_from_indices():
    dm = ImageModel()
    dm.data = np.ones((10,10))

    bbox = ((1,2), (3,4))

    result = boundingbox_to_indices(dm, bbox)

    assert (result == (1, 3, 3, 5))

def test_flag():
    wcsinfo = {
        'dec_ref': -0.00601415671349804,
        'ra_ref': -0.02073605215697509,
        'roll_ref': -0.0,
        'v2_ref': -453.5134,
        'v3_ref': -373.4826,
        'v3yangle': 0.0,
        'vparity': -1}

    instrument = {
        'detector': 'NRS1',
        'filter': 'F100LP',
        'grating': 'G140M',
        'name': 'NIRSPEC',
        'gwa_tilt': 37.0610,
        'gwa_xtilt': 0.0001,
        'gwa_ytilt': 0.0001,
        'msa_metadata_id':12}

    observation = {
        'date': '2016-09-05',
        'time': '8:59:37'}

    exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'NRSRAPID',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'NRS_MSASPEC',
        'zero_frame': False}

    im = ImageModel()
    im.data = np.random.rand(2048, 2048)
    im.error = np.random.rand(2048, 2048)
    im.dq = np.zeros((2048, 2048))

    im.meta.wcsinfo._instance.update(wcsinfo)
    im.meta.instrument._instance.update(instrument)
    im.meta.observation._instance.update(observation)
    im.meta.exposure._instance.update(exposure)

    metafl = get_file_path('msa_configuration.fits')
    im.meta.instrument.msa_metadata_file = metafl
    im.meta.dither.position_number = 1

    step = AssignWcsStep()
    msa_oper = step.get_reference_file(im, 'msaoper')

    im = step.call(im)

    wcs_ref_files = create_reference_filename_dictionary(im)

    failed_slitlets = create_slitlets(im, msa_oper)

    result = flag(im, failed_slitlets, wcs_ref_files)

    msa_open_dq = 536870912

    # Get all dqflags that are nonzero
    nonzero = np.nonzero(result.dq)

    # Make sure that all nonzero dqs are equal to the msa_open_dq
    assert (all(x == msa_open_dq for x in result.dq[nonzero]))
