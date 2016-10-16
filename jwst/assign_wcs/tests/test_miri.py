"""
Test MIRI MRS WCS transformation against IDT team data.

Notes:

1. New inverse polynomials were delivered in CDP5. However the test data table was not
   updated. The inverse transform tests against the IDT test data are commented out
   until we get new test data.

2. The center of the first pixel is (0, 0) so the pixel is from -.5 to .5.
   The test inputs (x, y data in the table provided by the IDT) include pixels just outside
   a slice. For example, x=60.995683, y=1014 is within pixel (61, 1014) and this happens
   to be the first pixel outside slice 21 in "1C". Because of the way they round the numbers,
   the IDT team considers this to be pixel (60, 1014), i.e. the last pixel in that slice.
   The pipeline comptes this pixel to be outside any slices and the output is NaN.
   To get around this, there's not test for the entire pipeline. Rather the tests are written to
   compute intermediate transforms using the table with verification data.

Both notes have been communicated to the INS team.

"""
from __future__ import absolute_import, division, unicode_literals, print_function

import numpy as np
from astropy.io import fits
from gwcs import wcs
from numpy.testing import utils

from ...datamodels.image import ImageModel
from .. import miri
from ..assign_wcs_step import AssignWcsStep


wcs_kw = {'wcsaxes': 2, 'ra_ref': 165, 'dec_ref': 54,
          'v2_ref': -8.3942412, 'v3_ref': -5.3123744, 'roll_ref': 37,
          'crpix1': 1024, 'crpix2': 1024,
          'cdelt1': .08, 'cdelt2': .08,
          'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
          'pc1_1': 1, 'pc1_2': 0, 'pc2_1': 0, 'pc2_2': 1,
          'pc3_1': 1, 'pc3_2': 0
          }

band_mapping = {'SHORT': 'A', 'MEDIUM': 'B', 'LONG': 'C'}


def create_hdul(detector, channel, band):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['instrume'] = 'MIRI'
    phdu.header['detector'] = detector
    phdu.header['CHANNEL'] = channel
    phdu.header['BAND'] = band
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2014-09-05'
    phdu.header['exp_type'] = 'MIR_MRS'
    phdu.header.update(wcs_kw)
    hdul.append(phdu)
    return hdul


def create_datamodel(hdul):
    im = ImageModel(hdul)
    ref = create_reference_files(im)

    pipeline = miri.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    im.meta.wcs = wcsobj
    return im


def create_reference_files(datamodel):
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)

    return refs


def run_test(model):
    wcsobj = model.meta.wcs
    for ch in model.meta.instrument.channel:
        ref_data = mrs_ref_data[ch + band_mapping[model.meta.instrument.band]]
        detector_to_alpha_beta = wcsobj.get_transform('detector', 'alpha_beta')
        ab_to_xan_yan = wcsobj.get_transform('alpha_beta', 'Xan_Yan').set_input(int(ch))
        ref_alpha = ref_data['alpha']
        ref_beta = ref_data['beta']
        ref_lam = ref_data['lam']
        #xyan_to_detector = wcsobj.get_transform('Xan_Yan', 'detector')
        #xin, yin = xyan_to_detector(ref_data['v2'], ref_data['v3'],
        #                                     ref_data['lam'])
        x, y = ref_data['x'], ref_data['y']
        for i, s in enumerate(ref_data['s']):
            sl = int(ch) * 100 + s
            alpha, beta, lam = detector_to_alpha_beta.set_input(sl)(x[i], y[i])
            utils.assert_allclose(alpha, ref_alpha[i], atol=10**-4)
            utils.assert_allclose(beta, ref_beta[i], atol=10**-4)
            utils.assert_allclose(lam, ref_lam[i], atol=10**-4)

        xan, yan, lam = ab_to_xan_yan(ref_alpha, ref_beta, ref_lam)
        utils.assert_allclose(xan, ref_data['v2'], atol=10**-4)
        utils.assert_allclose(yan, ref_data['v3'], atol=10**-4)
        utils.assert_allclose(lam, ref_data['lam'], atol=10**-4)

        #utils.assert_allclose(xin, x, atol=10**-5)
        #utils.assert_allclose(yin, y, atol=10**-5)


def test_miri_mrs_12A():
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="SHORT")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_12B():
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="MEDIUM")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_12C():
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="LONG")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_34A():
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="SHORT")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_34B():
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="MEDIUM")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_34C():
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="LONG")
    im = create_datamodel(hdul)
    run_test(im)


mrs_ref_data = {
    '1A': {'x': np.array([28.310396, 475.02154, 493.9777, 41.282537, 58.998266]),
           #'x': np.array([28.310396, 475.02154, 493.4777, 41.282537, 58.498266]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0, -1.66946, 1.65180, -1.70573, 1.70244]),
           'beta': np.array([0, -1.77210, -1.77210, 1.77210, 1.77210]),
           'lam': np.array([5.34437, 4.86642, 4.95325, 5.65296, 5.74349]),
           'v2': np.array([-8.39424, -8.41746, -8.36306, -8.42653, -8.37026]),
           'v3': np.array([-2.48763, -2.52081, -2.51311, -2.46269, -2.45395]),
           },
    '1B': {'x': np.array([28.648221, 475.07259, 493.98157, 41.559386, 59.738296]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0., -1.70796, 1.60161, -1.70854, 1.78261]),
           'beta': np.array([0., -1.77204, -1.77204, 1.77204, 1.77204]),
           'lam': np.array([6.17572, 5.62345, 5.72380, 6.53231, 6.63698]),
           'v2': np.array([-8.39426, -8.41808, -8.36368, -8.42682, -8.36899]),
           'v3': np.array([-2.48492, -2.51808, -2.51040, -2.46001, -2.45126])
           },
    '1C': {'x': np.array([30.461871, 477.23742, 495.96228, 43.905314, 60.995683]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0., -1.60587, 1.67276, -1.60766, 1.68720]),
           'beta': np.array([0., -1.77202, -1.77202, 1.77202, 1.77202]),
           'lam': np.array([7.04951, 6.42424, 6.53753, 7.45360, 7.57167]),
           'v2': np.array([-8.39357, -8.41570, -8.36165, -8.42457, -8.36996]),
           'v3': np.array([-2.48987, -2.52271, -2.51525, -2.46467, -2.45649])
           },
    '2A': {'x': np.array([992.158, 545.38386, 525.76143, 969.29711, 944.19303]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.11250, 2.10676, -2.17239, 2.10447]),
           'beta': np.array([0., -2.23775, -2.23775, 2.23775, 2.23775]),
           'lam': np.array([8.20797, 7.52144, 7.64907, 8.68677, 8.83051]),
           'v2': np.array([-8.39393, -8.42259, -8.35355, -8.43583, -8.36499]),
           'v3': np.array([-2.48181, -2.52375, -2.51357, -2.44987, -2.44022])
           },
    '2B': {'x': np.array([988.39977, 541.23447, 521.60207, 964.91753, 940.10325]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.10593, 2.10015, -2.08817, 2.10422]),
           'beta': np.array([0., -2.23781, -2.23781, 2.23781, 2.23781]),
           'lam': np.array([9.44205, 8.65341, 8.79991, 9.99257, 10.15795]),
           'v2': np.array([-8.39645, -8.42502, -8.35603, -8.43716, -8.36742]),
           'v3': np.array([-2.47773, -2.51972, -2.50938, -2.44554, -2.43626])
           },
    '2C': {'x': np.array([990.89693, 543.82344, 524.34514, 967.98318, 942.77564]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.07490, 2.11234, -2.14704, 2.14196]),
           'beta': np.array([0., -2.23778, -2.23778, 2.23778, 2.23778]),
           'lam': np.array([10.90225, 9.99162, 10.16079, 11.53780, 11.72887]),
           'v2': np.array([-8.39303, -8.42129, -8.35221, -8.43454, -8.36352]),
           'v3': np.array([-2.47869, -2.52052, -2.51036, -2.44668, -2.43712])
           },
    '3A': {'x': np.array([574.80828, 1001.0602, 984.6387, 547.27479, 518.89992]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0., -2.86745, 3.20982, -3.01230, 2.96643]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([12.5335, 13.49968, 13.33846, 11.77148, 11.52350]),
           'v2': np.array([-8.40590, -8.44849, -8.34906, -8.46070, -8.36174]),
           'v3': np.array([-2.48992, -2.54104, -2.52854, -2.44547, -2.43112])
           },
    '3B': {'x': np.array([574.26012, 1001.7349, 985.30166, 548.016, 519.98]),
           'y': np.array([512., 10., 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0, -3.17728, 2.92434, -3.29402, 2.60797]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([14.53997, 15.66039, 15.47355, 13.65622, 13.36833]),
           'v2': np.array([-8.40044, -8.44785, -8.34786, -8.46088, -8.36211]),
           'v3': np.array([-2.48588, -2.53771, -2.52512, -2.44219, -2.42776])
           },
    '3C': {'x': np.array([573.25446, 1000.21721, 983.92918, 546.00285, 518.2782]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0., -2.94573, 3.09057, -3.07810, 2.73161]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([16.79017, 18.08441, 17.86845, 15.76948, 15.43724]),
           'v2': np.array([-8.40205, -8.44574, -8.34664, -8.45859, -8.36196]),
           'v3': np.array([-2.48627, -2.53761, -2.52502, -2.44221, -2.42787]),
           },
    '4A': {'x': np.array([80.987181, 434.34987, 461.90855, 26.322503, 53.674656]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.74625, 3.72621, -3.94261, 3.62762]),
           'beta': np.array([-0.32802, -3.60821, -3.60821, 3.60821, 3.60821]),
           'lam': np.array([19.34914, 20.93078, 20.6464, 18.07975, 17.67221]),
           'v2': np.array([-8.38446, -8.43506, -8.31378, -8.46256, -8.33609]),
           'v3': np.array([-2.48058, -2.5444, -2.52426, -2.42449, -2.40839])
           },
    '4B': {'x': np.array([77.625553, 431.57061, 458.86869, 23.559111, 50.632416]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.64817, 3.73313, -3.73558, 3.74096]),
           'beta': np.array([-0.32802, -3.60821, -3.60821, 3.60821, 3.60821]),
           'lam': np.array([22.38267, 24.21212, 23.88327, 20.91426, 20.44279]),
           'v2': np.array([-8.38581, -8.43443, -8.3141, -8.46152, -8.33604]),
           'v3': np.array([-2.48185, -2.54526, -2.52568, -2.42513, -2.40959])
           },
    '4C': {'x': np.array([79.662059, 433.73384, 460.75026, 25.820431, 52.412219]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.61682, 3.69713, -3.66259, 3.69888]),
           'beta': np.array([-0.32802, -3.60819, -3.60819, 3.60819, 3.60819]),
           'lam': np.array([26.18343, 28.32354, 27.93894, 24.46574, 23.91417]),
           'v2': np.array([-8.38603, -8.43509, -8.31524, -8.45888, -8.33707]),
           'v3': np.array([-2.48315, -2.54647, -2.52661, -2.42721, -2.41060])
           }
}
