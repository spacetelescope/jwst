"""
Test for flat_field.find_dispaxis
"""
import numpy as np

from jwst.flatfield import flat_field

def test_find_dispaxis():

    # Create two dummy arrays of wavelengths.  The arrays are mostly zero,
    # but the non-zero values are clearly increasing more rapidly in one
    # direction rather than the other.

    wl_h = np.zeros((30, 200), dtype=np.float32)
    wl_h[16, :] = np.arange(200, dtype=np.float32)
    wl_h[17, :] = np.arange(200, dtype=np.float32) + 0.01
    wl_h[18, :] = np.arange(200, dtype=np.float32) + 0.02
    wl_h[19, :] = np.arange(200, dtype=np.float32) + 0.03

    wl_v = np.zeros((200, 30), dtype=np.float32)
    wl_v[:, 17] = np.arange(200, dtype=np.float32)
    wl_v[:, 18] = np.arange(200, dtype=np.float32) + 0.01
    wl_v[:, 19] = np.arange(200, dtype=np.float32) + 0.02
    wl_v[:, 20] = np.arange(200, dtype=np.float32) + 0.03

    dispaxis = flat_field.find_dispaxis(wl_h)
    assert dispaxis == 1

    dispaxis = flat_field.find_dispaxis(wl_v)
    assert dispaxis == 2
