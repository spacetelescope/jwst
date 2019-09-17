import numpy as np

from jwst.datamodels import MSAModel
from jwst.msaflagopen.msaflag_open import (or_subarray_with_array,
                                           id_from_xy,
                                           get_failed_open_shutters,
                                           create_slitlets)

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

    shutters_per_row = 365

    # First row of msaoper.json
    data = {"Q": 1,"x": 1,"y": 1,"state": "closed","TA state": "closed",
            "Internal state": "normal",
            "Vignetted": "yes"}

    x = data['x']
    y = data['y']

    result = id_from_xy(x, y)

    assert (x + (y-1)*shutters_per_row == result)

def test_get_failed_open_shutters():
    """test that failed open shutters are returned from reference file"""

    result = get_failed_open_shutters('msa_oper.json')

    for shutter in result:
        assert (shutter['state'] == 'open')

def test_create_slitlets():
    
    dm = MSAModel()

    result = create_slitlets(dm, 'msa_oper.json')

    return result