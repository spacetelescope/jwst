import numpy as np

from jwst import datamodels
from jwst.refpix.irs2_subtract_reference import flag_bad_refpix


def test_flag_bad_refpix():
    data = np.ones((2, 3, 3200, 5), dtype=np.float32)
    pixeldq = np.full((3200, 5), 0)

    # set some previously marked bad reference pixels
    data[:, :, 670, :] = 0.
    data[:, :, 688, :] = 0.
    data[:, :, 2110, :] = 0.
    pixeldq[670, :] = datamodels.dqflags.pixel['BAD_REF_PIXEL']
    pixeldq[688, :] = datamodels.dqflags.pixel['BAD_REF_PIXEL']
    pixeldq[2110, :] = datamodels.dqflags.pixel['BAD_REF_PIXEL']

    # set the intermittent bad pixels
    data[:, :, 648, :] = 10.
    data[:, :, 668, :] = 20.
    data[:, :, 988, :] = 11.
    data[:, :, 1369, :] = 7.
    data[:, :, 2150, :] = 13.
    data[:, :, 3128, :] = 15.

    scipix_n, refpix_r = 16, 4
    ovr_corr_mitigation_ftr = 3.0

    input_model = datamodels.RampModel(data)
    input_model.pixeldq = pixeldq
    input_model.meta.exposure.nrs_normal = scipix_n
    input_model.meta.exposure.nrs_reference = refpix_r

    flag_bad_refpix(input_model, n_sigma=ovr_corr_mitigation_ftr)

    compare = np.ones((2, 3, 3200, 5), dtype=np.float32)

    # previously marked pixel, lower pix bad, upper pix good
    compare[:, :, 688, :] = 1.
    # previously marked pixel, lower and upper pix good
    compare[:, :, 670, :] = 1.
    compare[:, :, 2110, :] = 1.

    # no lower pix, upper pix bad, neighbor okay
    compare[:, :, 648:650, :] = 1.
    # lower and upper pix bad, neighbor bad
    compare[:, :, 668, :] = 0.
    # lower pix bad, upper pix good
    compare[:, :, 669, :] = 1.
    # lower and upper pix good
    compare[:, :, 988:990, :] = 1.
    compare[:, :, 1368:1370, :] = 1.
    compare[:, :, 2150:2152, :] = 1.
    compare[:, :, 3128:3130, :] = 1.

    assert np.allclose(data, compare)
