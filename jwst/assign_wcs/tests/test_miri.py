"""
Test MIRI MRS WCS transformation against IDT team data.

Notes:

1. Test data use CDP-8b placeholder values computed by D. Law for the time being.

"""
import numpy as np
from astropy.io import fits
from astropy.modeling.models import Scale, Shift, Identity
from gwcs import wcs
from numpy.testing import utils

from ...datamodels.image import ImageModel
from .. import miri
from ..assign_wcs_step import AssignWcsStep


wcs_kw = {'wcsaxes': 3, 'ra_ref': 165, 'dec_ref': 54,
          'v2_ref': -8.3942412, 'v3_ref': -5.3123744, 'roll_ref': 37,
          'crpix1': 1024, 'crpix2': 1024, 'crpix3': 0,
          'cdelt1': .08, 'cdelt2': .08, 'cdelt3': 1,
          'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN', 'ctype3': 'WAVE',
          'pc1_1': 1, 'pc1_2': 0, 'pc1_3': 0,
          'pc2_1': 0, 'pc2_2': 1, 'pc2_3': 0,
          'pc3_1': 0, 'pc3_2': 0, 'pc3_3': 1,
          'cunit1': 'deg', 'cunit2': 'deg', 'cunit3': 'um',
          }

band_mapping = {'SHORT': 'A', 'MEDIUM': 'B', 'LONG': 'C'}


def create_hdul(detector, channel, band):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['telescop'] = "JWST"
    phdu.header['filename'] = "test" + channel + band
    phdu.header['instrume'] = 'MIRI'
    phdu.header['detector'] = detector
    phdu.header['CHANNEL'] = channel
    phdu.header['BAND'] = band
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2017-09-05'
    phdu.header['exp_type'] = 'MIR_MRS'
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header.update(wcs_kw)
    hdul.append(phdu)
    hdul.append(scihdu)
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
        ab_to_v2v3 = wcsobj.get_transform('alpha_beta', 'v2v3').set_input(int(ch))
        v2v3_to_ab = wcsobj.get_transform('v2v3', 'alpha_beta').set_input(int(ch))
        ab_to_detector = wcsobj.get_transform('alpha_beta','detector')
        
        ref_alpha = ref_data['alpha']
        ref_beta = ref_data['beta']
        ref_lam = ref_data['lam']
        ref_v2 = ref_data['v2']
        ref_v3 = ref_data['v3']

        x, y = ref_data['x'], ref_data['y']
        for i, s in enumerate(ref_data['s']):
            sl = int(ch) * 100 + s
            alpha, beta, lam = detector_to_alpha_beta.set_input(sl)(x[i], y[i])
            utils.assert_allclose(alpha, ref_alpha[i], atol=0.05)
            utils.assert_allclose(beta, ref_beta[i], atol=0.05)
            utils.assert_allclose(lam, ref_lam[i], atol=0.05)

        v2, v3, lam = ab_to_v2v3(ref_alpha, ref_beta, ref_lam)
        utils.assert_allclose(v2, ref_v2, atol=0.05)
        utils.assert_allclose(v3, ref_v3, atol=0.05)
        utils.assert_allclose(lam, ref_lam, atol=0.05)

        # Test the reverse transform
        alpha_back, beta_back, lam_back = v2v3_to_ab(v2,v3,lam)
        utils.assert_allclose(alpha_back, ref_alpha, atol=0.05)
        utils.assert_allclose(beta_back, ref_beta, atol=0.05)
        utils.assert_allclose(lam_back, ref_lam, atol=0.05)
        
        for i, s in enumerate(ref_data['s']):
            sl = int(ch) * 100 + s
            x_back, y_back = ab_to_detector.set_input(sl)(alpha_back[i],beta_back[i],lam_back[i])
            utils.assert_allclose(x_back, x[i], atol=0.08)
            utils.assert_allclose(y_back, y[i], atol=0.08)


def test_miri_mrs_12A():
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="SHORT")
    im = create_datamodel(hdul)
    run_test(im)


def test_miri_mrs_12B():
    hdul = create_hdul(detector="MIRIFUSHORT", channel="12", band="MEDIUM")
    print(hdul[0].header)
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

# MRS test reference data
# These values are all CDP-8b placeholders mocked up by D. Law

mrs_ref_data = {
    '1A': {'x': np.array([123, 468]),
           'y': np.array([245, 1000]),
           's': np.array([9, 12]),
           'alpha': np.array([0.10731135, 1.12529977]),
           'beta': np.array([-0.35442029,  0.17721014]),
           'lam': np.array([5.11499695, 5.74016832]),
           'v2': np.array([-503.49916454, -502.56838429]),
           'v3': np.array([-318.40628359, -319.08393052]),
           },
    '1B': {'x': np.array([51, 244]),
           'y': np.array([1016,  476]),
           's': np.array([21, 17]),
           'alpha': np.array([0.38199789, 0.63723143]),
           'beta': np.array([1.77204177, 1.06322506]),
           'lam': np.array([6.62421656, 6.12990972]),
           'v2': np.array([-503.5315322 , -503.17667842]),
           'v3': np.array([-320.71381511, -320.05014591]),
           },
    '1C': {'x': np.array([127, 394]),
           'y': np.array([747, 111]),
           's': np.array([9, 3]),
           'alpha': np.array([0.37487296, -0.87620923]),
           'beta': np.array([-0.35440438, -1.4176175]),
           'lam': np.array([7.26440645, 6.52961571]),
           'v2': np.array([-503.19062471, -504.27276937]),
           'v3': np.array([-318.31059564, -317.08599443]),
           },
    '2A': {'x': np.array([574, 913]),
           'y': np.array([578, 163]),
           's': np.array([10, 16]),
           'alpha': np.array([0.02652122, -1.44523112]),
           'beta': np.array([0.27971819, 1.9580273]),
           'lam': np.array([8.22398597, 7.66495464]),
           'v2': np.array([-503.65420691, -505.38172957]),
           'v3': np.array([-319.37148692, -320.82933868]),
           },
    '2B': {'x': np.array([634, 955]),
           'y': np.array([749,  12]),
           's': np.array([11, 17]),
           'alpha': np.array([-1.31986085, -1.66029886]),
           'beta': np.array([0.5594521 , 2.23780842]),
           'lam': np.array([9.85535403, 8.65341739]),
           'v2': np.array([-505.18703764, -505.80250684]),
           'v3': np.array([-319.7057936 , -321.32425399]),
           },
    '2C': {'x': np.array([530, 884]),
           'y': np.array([965, 346]),
           's': np.array([1, 7]),
           'alpha': np.array([1.17219936, -0.13199122]),
           'beta': np.array([-2.23777695, -0.55944424]),
           'lam': np.array([11.68798183, 10.65732315]),
           'v2': np.array([-502.0634552 , -503.62291245]),
           'v3': np.array([-317.2417194 , -318.70820411]),
           },
    '3A': {'x': np.array([573, 913]),
           'y': np.array([851, 323]),
           's': np.array([8, 10]),
           'alpha': np.array([-1.0181757 ,  0.65295329]),
           'beta': np.array([-0.19490689,  0.58472067]),
           'lam': np.array([11.84245153, 12.96396074]),
           'v2': np.array([-505.35888594, -503.7824966]),
           'v3': np.array([-318.46913272, -319.4685406]),
           },
    '3B': {'x': np.array([606, 861]),
           'y': np.array([926, 366]),
           's': np.array([15, 11]),
           'alpha': np.array([-1.5124193 , -0.79361415]),
           'beta': np.array([2.53378956, 0.97453445]),
           'lam': np.array([13.60306079, 14.94878428]),
           'v2': np.array([-505.82191056, -504.9372123]),
           'v3': np.array([-321.34413558, -319.90108102]),
           },
    '3C': {'x': np.array([663, 852]),
           'y': np.array([822,  86]),
           's': np.array([14, 11]),
           'alpha': np.array([0.83845626, -1.00005387]),
           'beta': np.array([2.14397578, 0.97453445]),
           'lam': np.array([16.01468948, 17.97678143]),
           'v2': np.array([-503.52817761, -505.23700039]),
           'v3': np.array([-321.27004219, -319.84577337]),
           },
    '4A': {'x': np.array([448, 409]),
           'y': np.array([873,  49]),
           's': np.array([1, 7]),
           'alpha': np.array([-0.45466621, -1.07614592]),
           'beta': np.array([-3.60820915,  0.32801901]),
           'lam': np.array([18.05366191, 20.88016154]),
           'v2': np.array([-502.89806847, -504.25439193]),
           'v3': np.array([-315.86847223, -319.65622713]),
           },
    '4B': {'x': np.array([380, 260]),
           'y': np.array([926, 325]),
           's': np.array([2, 9]),
           'alpha': np.array([1.64217386, -1.70062938]),
           'beta': np.array([-2.95217555,  1.64009753]),
           'lam': np.array([20.69573674, 23.17990504]),
           'v2': np.array([-501.01720495, -505.23791555]),
           'v3': np.array([-316.76598039, -320.79546159]),
           },
    '4C': {'x': np.array([309, 114]),
           'y': np.array([941, 196]),
           's': np.array([3, 11]),
           'alpha': np.array([1.65440228, -0.87408042]),
           'beta': np.array([-2.29611932,  2.95215341]),
           'lam': np.array([24.17180582, 27.63402178]),
           'v2': np.array([-501.1647203 , -504.64107203]),
           'v3': np.array([-317.34628   , -322.10088837]),
           }
    
}
