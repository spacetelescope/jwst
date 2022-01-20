"""
Test MIRI MRS WCS transformation against IDT team data.

Notes:

1. Test data use CDP-8b placeholder values computed by D. Law for the time being.

"""
import numpy as np
from astropy.io import fits
from gwcs import wcs
from numpy.testing import assert_allclose

from jwst.datamodels import ImageModel
from jwst.assign_wcs import miri
from jwst.assign_wcs import AssignWcsStep


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
        ab_to_detector = wcsobj.get_transform('alpha_beta', 'detector')

        ref_alpha = ref_data['alpha']
        ref_beta = ref_data['beta']
        ref_lam = ref_data['lam']
        ref_v2 = ref_data['v2']
        ref_v3 = ref_data['v3']

        x, y = ref_data['x'], ref_data['y']
        for i, s in enumerate(ref_data['s']):
            sl = int(ch) * 100 + s
            alpha, beta, lam = detector_to_alpha_beta.set_input(sl)(x[i], y[i])
            assert_allclose(alpha, ref_alpha[i], atol=0.05)
            assert_allclose(beta, ref_beta[i], atol=0.05)
            assert_allclose(lam, ref_lam[i], atol=0.05)

        v2, v3, lam = ab_to_v2v3(ref_alpha, ref_beta, ref_lam)
        assert_allclose(v2, ref_v2, atol=0.05)
        assert_allclose(v3, ref_v3, atol=0.05)
        assert_allclose(lam, ref_lam, atol=0.05)

        # Test the reverse transform
        alpha_back, beta_back, lam_back = v2v3_to_ab(v2, v3, lam)
        assert_allclose(alpha_back, ref_alpha, atol=0.05)
        assert_allclose(beta_back, ref_beta, atol=0.05)
        assert_allclose(lam_back, ref_lam, atol=0.05)

        for i, s in enumerate(ref_data['s']):
            sl = int(ch) * 100 + s
            x_back, y_back = ab_to_detector.set_input(sl)(alpha_back[i], beta_back[i], lam_back[i])
            assert_allclose(x_back, x[i], atol=0.08)
            assert_allclose(y_back, y[i], atol=0.08)


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


# MRS test reference data
mrs_ref_data = {
    '1A': {'x': np.array([76.0, 354.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([0.05765538365149925, -0.017032619150995743]),
           'beta': np.array([-0.17721014379699995, -1.240471006579]),
           'lam': np.array([5.348546577257886, 5.5136420569934925]),
           'v2': np.array([-503.57285226785064, -503.4979806620663]),
           'v3': np.array([-318.5749892859028, -317.5090073056335]),
           },
    '1B': {'x': np.array([76.0, 355.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.012990737471741731, 0.10766447914943456]),
           'beta': np.array([-0.17720417669099997, -1.240429236837]),
           'lam': np.array([6.168310398808807, 6.358007642348213]),
           'v2': np.array([-503.643100332753, -503.37069816112813]),
           'v3': np.array([-318.72773306477103, -317.6938248759762]),
           },
    '1C': {'x': np.array([78.0, 356.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([0.02871804339196271, -0.028315822861031847]),
           'beta': np.array([-0.17720218765499984, -1.240415313585]),
           'lam': np.array([7.006608159574103, 7.218455147089075]),
           'v2': np.array([-503.5598371896608, -503.45975848303885]),
           'v3': np.array([-318.4367657801553, -317.3779485524358]),
           },
    '2A': {'x': np.array([574.0, 719.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([0.022862344416012093, 0.024104763006107532]),
           'beta': np.array([0.27971818633699996, -1.3985909316610001]),
           'lam': np.array([8.139463800053713, 8.423879719165456]),
           'v2': np.array([-503.65782416704644, -503.3907046961389]),
           'v3': np.array([-319.3709764579651, -317.71318662530217]),
           },
    '2B': {'x': np.array([570.0, 715.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.04101483043351095, -0.021964438108625473]),
           'beta': np.array([0.27972605223, -1.39863026115]),
           'lam': np.array([9.49091778668766, 9.826112199836349]),
           'v2': np.array([-503.872441161987, -503.58468453126545]),
           'v3': np.array([-319.6066193816802, -317.9526192173689]),
           },
    '2C': {'x': np.array([573.0, 718.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.08065540123411097, -0.07196315905207484]),
           'beta': np.array([0.2797221192789996, -1.3986105964070001]),
           'lam': np.array([10.909558387414732, 11.292658213110698]),
           'v2': np.array([-503.7062367371822, -503.4292038385116]),
           'v3': np.array([-319.54349206004053, -317.8886490566051]),
           },
    '3A': {'x': np.array([918.0, 827.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.14902640584477922, -0.1394111481404252]),
           'beta': np.array([0.5847206674920002, -1.7541620024759998]),
           'lam': np.array([12.586085291551054, 12.171803779467552]),
           'v2': np.array([-504.57532179184557, -504.3428404141017]),
           'v3': np.array([-319.3596209726561, -317.0363338552647]),
           },
    '3B': {'x': np.array([919.0, 827.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.18610616903060873, 0.05223448620927229]),
           'beta': np.array([0.5847206674920002, -1.7541620024759998]),
           'lam': np.array([14.60074101845329, 14.120353260795175]),
           'v2': np.array([-504.29128783278026, -503.81513623681207]),
           'v3': np.array([-319.5977726217362, -317.30169796071453]),
           },
    '3C': {'x': np.array([917.0, 826.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.08917305254544772, -0.09924683542340063]),
           'beta': np.array([0.5847206674920002, -1.7541620024759998]),
           'lam': np.array([16.860616228418674, 16.305648049347006]),
           'v2': np.array([-504.29179372150304, -504.06099473540036]),
           'v3': np.array([-319.5864222556306, -317.26146053061063]),
           },
    '4A': {'x': np.array([195.0, 232.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.18281231856817595, -0.10820926727846612]),
           'beta': np.array([2.2961330928359995, -1.640095066308]),
           'lam': np.array([19.42967253041467, 18.733785802367724]),
           'v2': np.array([-503.73916258138155, -502.9287085654886]),
           'v3': np.array([-321.7198475574414, -317.8596067111157]),
           },
    '4B': {'x': np.array([192.0, 229.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.03596952007447607, -0.10259402857181654]),
           'beta': np.array([2.2961365363689996, -1.640097525977]),
           'lam': np.array([22.47574268879503, 21.67074830984225]),
           'v2': np.array([-503.7051048327475, -502.9891450100565]),
           'v3': np.array([-321.6637327196876, -317.78403487305536]),
           },
    '4C': {'x': np.array([194.0, 231.0]),
           'y': np.array([512.0, 700.0]),
           's': np.array([10, 4]),
           'alpha': np.array([-0.0661930805678849, -0.01176625661012924]),
           'beta': np.array([2.296119318687, -1.640085227631]),
           'lam': np.array([26.292379242285914, 25.350694577065074]),
           'v2': np.array([-503.7171854824459, -502.9282547181127]),
           'v3': np.array([-321.57006077329663, -317.7252303132135]),
           }

}
