"""
Test for flat_field.find_dispaxis
"""
import pytest

import numpy as np

from jwst.flatfield import flat_field

def test_find_dispaxis():

    # Create a dummy array of wavelengths.  The array is mostly zero,
    # but the non-zero values are clearly increasing more rapidly
    # in the horizontal direction than in the vertical direction.
    wl = np.zeros((30, 200), dtype=np.float32)
    wl[16, :] = np.arange(200, dtype=np.float32)
    wl[17, :] = np.arange(200, dtype=np.float32) + 0.01
    wl[18, :] = np.arange(200, dtype=np.float32) + 0.02
    wl[19, :] = np.arange(200, dtype=np.float32) + 0.03

    dispaxis = flat_field.find_dispaxis(wl)
    assert dispaxis == 1

    # This is the code that was in combine_fast_slow and that find_dispaxis
    # replaced.
    (ny, nx) = wl.shape
    mid_x = nx // 2
    mid_y = ny // 2
    dwlx = abs(wl[mid_y, mid_x+1] - wl[mid_y, mid_x-1])
    dwly = abs(wl[mid_y+1, mid_x] - wl[mid_y-1, mid_x])
    if dwlx >= dwly:
        dispaxis = 1
    else:
        dispaxis = 2

    assert dispaxis == 2
