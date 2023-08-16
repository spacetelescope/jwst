import numpy as np

from jwst.refpix.irs2_subtract_reference import rm_intermittent_badpix


def test_rm_intermittent_badpix():
    data = np.ones((2, 3, 5, 3200), dtype=np.float32)
    # set the bad reference pixels flagged in the ref file
    data[..., 688] = 0.
    data[..., 2110] = 0.
    # set the intermittent bad pixels
    data[..., 648] = 10.
    data[..., 988] = 11.
    data[..., 1369] = 7.
    data[..., 2150] = 13.
    data[..., 3128] = 15.

    scipix_n, refpix_r = 16, 4
    ovr_corr_mitigation_ftr = 3.0
    rm_intermittent_badpix(data, scipix_n, refpix_r, ovr_corr_mitigation_ftr)

    compare = np.ones((2, 3, 5, 3200), dtype=np.float32)
    compare[..., 688] = 0.
    compare[..., 2110] = 0.

    assert np.allclose(data, compare)
