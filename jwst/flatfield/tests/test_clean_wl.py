"""
Test for flat_field.clean_wl
"""
import numpy as np

from jwst.flatfield import flat_field


def test_clean_wl_1():

    # The dispersion direction is horizontal.
    wl_h = np.zeros((5, 9), dtype=np.float32)
    wl_h[1, 3:7] = np.arange(3, 7, dtype=np.float32)
    wl_h[2, 2:8] = np.arange(2, 8, dtype=np.float32) + 0.01
    wl_h[3, 4] = 4.02

    # These are the "cleaned" arrays that we expect.
    wl_h_clean = np.array(
        [[2.01, 2.01, 2.01, 3.005, 4.01, 5.005, 6.005, 7.01, 7.01],
         [2.01, 2.01, 2.01, 3., 4., 5., 6., 7.01, 7.01],
         [2.01, 2.01, 2.01, 3.01, 4.01, 5.01, 6.01, 7.01, 7.01],
         [2.01, 2.01, 2.01, 3.005, 4.02, 5.005, 6.005, 7.01, 7.01],
         [2.01, 2.01, 2.01, 3.005, 4.01, 5.005, 6.005, 7.01, 7.01]],
        dtype=np.float64)

    # 1 means horizontal dispersion
    wl_hc = flat_field.clean_wl(wl_h, 1)
    assert np.allclose(wl_hc, wl_h_clean, atol=0., rtol=1.e-6)


def test_clean_wl_2():

    wl_h = np.zeros((5, 9), dtype=np.float32)
    wl_h[1, 3:7] = np.arange(3, 7, dtype=np.float32)
    wl_h[2, 2:8] = np.arange(2, 8, dtype=np.float32) + 0.01
    wl_h[3, 4] = 4.02

    # The dispersion direction is vertical.
    wl_v = np.rot90(wl_h)

    wl_v_clean = np.array([[7.01, 7.01, 7.01, 7.01, 7.01],
                           [7.01, 7.01, 7.01, 7.01, 7.01],
                           [6.005, 6., 6.01, 6.005, 6.005],
                           [5.005, 5., 5.01, 5.005, 5.005],
                           [4.01, 4., 4.01, 4.02, 4.01],
                           [3.005, 3., 3.01, 3.005, 3.005],
                           [2.01, 2.01, 2.01, 2.01, 2.01],
                           [2.01, 2.01, 2.01, 2.01, 2.01],
                           [2.01, 2.01, 2.01, 2.01, 2.01]],
                          dtype=np.float64)

    # 2 means vertical dispersion
    wl_vc = flat_field.clean_wl(wl_v, 2)
    assert np.allclose(wl_vc, wl_v_clean, atol=0., rtol=1.e-6)
