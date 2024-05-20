import pytest
import numpy as np
from jwst.badpix_selfcal.badpix_selfcal_step import BadpixSelfcalStep
from jwst.badpix_selfcal.badpix_selfcal import MRSPolyfitBackgroundFlagger, MRSMedfiltBackgroundFlagger
from jwst import datamodels as dm
from jwst.assign_wcs import AssignWcsStep, miri
from gwcs import wcs
from astropy.io import fits

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

hotpixel_intensity = 50
outlier_indices = [(100, 100), (300, 300), (500, 600), (1000, 900)]

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


def create_reference_files(datamodel):
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)

    return refs


@pytest.fixture(scope="module")
def background():
    """
    Create background IFUImageModel for testing. This is a mockup of the expected
    background data in a .rate file.
    Three components: random noise, low-order variability, and outliers
    """
    
    # random noise
    rng = np.random.default_rng(seed=77)
    shp = (1024, 1032)
    noise = rng.standard_normal(shp)

    # make 2-d polynomial representing background level
    c = np.array([[1, 3, 5], [2, 4, 6]])
    x = np.linspace(-1, 1, shp[0])
    y = np.linspace(-1, 1, shp[1])
    low_order_variability = np.polynomial.polynomial.polygrid2d(x, y, c)

    # add some outliers
    outliers = np.zeros(shp)
    for idx in outlier_indices:
        outliers[idx] += hotpixel_intensity
    # one negative one just for fun
    outliers[100, 100] = -hotpixel_intensity

    mock_data = low_order_variability + noise + outliers

    # build an IFUImageModel from these data and give it a wcs
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="LONG")
    hdul[1].data = mock_data

    im = dm.IFUImageModel(hdul)
    ref = create_reference_files(im)
    pipeline = miri.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    im.meta.wcs = wcsobj

    return im


@pytest.fixture(scope="module")
def sci():
    """Create science .rate file for testing"""
    return


@pytest.mark.parametrize("algorithm", [MRSPolyfitBackgroundFlagger, MRSMedfiltBackgroundFlagger])
def test_background_flagger_mrs(algorithm, background):
    """"""

    # first get wavelength array from WCS
    bg = background.data
    shp = bg.shape
    basex,basey = np.meshgrid(np.arange(shp[1]),np.arange(shp[0]))
    _,_,lam=background.meta.wcs.transform('detector','world',basex,basey)

    # pass into the MRSBackgroundFlagger and check it found the right pixels
    flagfrac = 0.001
    result = algorithm().flag_background(bg, lam, flagfrac)
    result_tuples = [(i,j) for i,j in zip(*result)]

    # check that the hot pixels were among those flagged
    for idx in outlier_indices:
        assert idx in result_tuples

    # check that the number of flagged pixels is as expected
    assert np.isclose(len(result_tuples)/bg.size, flagfrac*2, atol=0.0001)



# def test_badpix_selfcal_step(tmp_cwd, background, sci):

#     # Run the badpix_selfcal step
#     result = BadpixSelfcalStep.call(sci, background)

#     # Compare the calculated result with the expected result
#     assert np.allclose(result.data, expected_result.data)