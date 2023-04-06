"""
Unit test for Cube Build testing for various wcs functions
"""

import numpy as np
import math

from stdatamodels.jwst import datamodels

from jwst.cube_build import ifu_cube
from jwst.cube_build import coord
from jwst.cube_build import cube_build_wcs_util
from jwst.cube_build import instrument_defaults


shape = (101, 101)
xcenter = 50
ycenter = 50

slice_gap = np.zeros(shape)
slice_gap[:, 5:25] = 1
slice_gap[:, 30:50] = 2
slice_gap[:, 55:75] = 3
slice_gap[:, 80:] = 4


def dummy_wcs(x, y):
    """ Simple WCS for testing """

    global xcenter, ycenter, shape, slice_gap

    # for given shape and wcs this will result in
    # ra from 40.6 to 49.9 [x = -49, 49:  (x +1 -crpix1) * cdelt1 + crval1]
    # dec from 45.1 to 45.4
    # wave from 7.5 to 8.5

    crpix1 = xcenter
    crpix3 = 1.0
    cdelt1 = 0.1
    cdelt2 = 0.1
    cdelt3 = 0.01

    crval1 = 45.0
    crval2 = 45.0
    crval3 = 7.5

    dec = np.zeros(shape, dtype=float)
    ra = np.zeros(shape, dtype=float)
    wave = np.zeros(shape, dtype=float)

    wave = (y + 1 - crpix3) * cdelt3 + crval3
    index_x1 = np.where(slice_gap == 1)  # slice 1
    dec[index_x1] = crval2 + 1 * cdelt2

    index_x2 = np.where(slice_gap == 2)  # slice 2
    dec[index_x2] = crval2 + 2 * cdelt2

    index_x3 = np.where(slice_gap == 3)  # slice 3
    dec[index_x3] = crval2 + 3 * cdelt2

    index_x4 = np.where(slice_gap == 4)  # slice 4
    dec[index_x4] = crval2 + 4 * cdelt2
    ra = (x + 1 - crpix1) * cdelt1 + crval1

    index_nan = np.where(slice_gap == 0)
    dec[index_nan] = np.nan
    ra[index_nan] = np.nan
    wave[index_nan] = np.nan

    return ra, dec, wave


def test_coord_trans1():
    """ Test finding xi,eta and cos 90, ra 45 """

    crval1 = 45.0
    crval2 = 90.0
    diff_ra = 3.0  # in arc seconds
    diff_dec = 3.0  # in arc seconds
    ra = crval1 + diff_ra / 3600.0
    dec = crval2 + diff_dec / 3600.0

    # declination near 90 yields xi values close to 0
    # and an eta close to diff_dec

    xi, eta = coord.radec2std(crval1, crval2, ra, dec)
    assert math.isclose(xi, 0.0, abs_tol=0.001)
    assert math.isclose(eta, diff_ra, abs_tol=0.001)


def test_coord_trans2():
    """ Test finding ci,eta and cos 45, ra 45 """
    crval1 = 45.0
    crval2 = 45.0
    diff_ra = 5.0  # in arc seconds
    diff_dec = 5.0  # in arc seconds
    ra = crval1 + diff_ra / 3600.0
    dec = crval2 + diff_dec / 3600.0

    # both crval1 and crval2 = 45, gives h in equation = 1
    # and an eta close to diff_dec

    xi, eta = coord.radec2std(crval1, crval2, ra, dec)
    assert math.isclose(xi, -3.535, abs_tol=0.001)
    assert math.isclose(eta, diff_ra, abs_tol=0.001)


def test_coord_trans3():
    """ Test going from ra,dec -> xi,eta -> ra,dec """

    crval1 = 27.89
    crval2 = 56.08
    diff_ra = 5.0  # in arc seconds
    diff_dec = 5.0  # in arc seconds
    ra = crval1 + diff_ra / 3600.0
    dec = crval2 + diff_dec / 3600.0

    # both crval1 and crval2 = 45, gives h in equation = 1
    # and an eta close to diff_dec

    xi, eta = coord.radec2std(crval1, crval2, ra, dec)
    ra_test, dec_test = coord.std2radec(crval1, crval2, xi, eta)
    assert math.isclose(ra, ra_test, abs_tol=0.00001)
    assert math.isclose(dec, dec_test, abs_tol=0.00001)


def test_wrap_ra():
    """ Test function wrap_ra but all ra on same side of 0/360 border """

    # test 1  wrap ra should do nothing
    ra = np.zeros(5, dtype=float)
    ra[0] = 0.26
    ra[1] = 0.12
    ra[2] = 0.35
    ra[3] = 0.23
    ra[4] = 0.45

    ra_test = cube_build_wcs_util.wrap_ra(ra)
    assert np.all(ra == ra_test)

    # test 2 the 4th value should be changed to a negative
    ra[4] = 359.8
    ra_compare = ra.copy()
    ra_compare[4] = ra[4] - 360.0
    ra_test = cube_build_wcs_util.wrap_ra(ra)
    assert np.all(ra_compare == ra_test)

    # test 3 the 4th value will be changed to 362.0
    ra[0] = 356.8
    ra[1] = 350.0
    ra[2] = 358.9
    ra[3] = 359.6
    ra[4] = 2.0

    ra_compare = ra.copy()
    ra_compare[4] = 360.0 + ra[4]
    ra_test = cube_build_wcs_util.wrap_ra(ra)
    assert np.all(ra_compare == ra_test)


def test_setup_wcs():
    """ setting size of IFU given input min,max and cdelts """
    ra1 = 98.83006930071556
    dec1 = -66.8274397956464
    ra2 = 98.8334511693978
    dec2 = -66.82720255674548
    ra3 = 98.83217303911556
    dec3 = -66.82798122840644
    ra4 = 98.83139113841862
    dec4 = -66.82665445039441
    lambda_min = 6.420
    lambda_max = 7.511

    corner_a = []
    corner_b = []
    corner_a.append(ra1)
    corner_a.append(ra2)
    corner_a.append(ra3)
    corner_a.append(ra4)

    corner_b.append(dec1)
    corner_b.append(dec2)
    corner_b.append(dec3)
    corner_b.append(dec4)

    pars_cube = {
        'scale1': 0.0,
        'scale2': 0.0,
        'scalew': 0.0,
        'interpolation': 'pointcloud',
        'weighting': 'emsm',
        'weight_power': 2,
        'coord_system': 'skyalign',
        'rois': 0.0,
        'roiw': 0.0,
        'wavemin': lambda_min,
        'wavemax': lambda_max,
        'skip_dqflagging': False,
        'xdebug': None,
        'ydebug': None,
        'zdebug': None,
        'spaxel_debug': None}

    pipeline = 3
    input_model = None
    output_name_base = None
    output_type = None
    instrument = None
    list_par1 = None
    list_par2 = None
    master_table = None
    instrument_info = None
    thiscube = ifu_cube.IFUCubeData(
        pipeline,
        input_model,
        output_name_base,
        output_type,
        instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube)

    thiscube.cdelt1 = 0.13
    thiscube.cdelt2 = 0.13
    thiscube.cdelt3 = 0.001
    thiscube.linear_wavelength = True
    thiscube.set_geometry(corner_a, corner_b, lambda_min, lambda_max)

    assert thiscube.naxis1 == 41
    assert thiscube.naxis2 == 41
    assert thiscube.naxis3 == 1092


def test_footprint_miri():

    global shape

    input_model = datamodels.IFUImageModel()
    input_model.meta.instrument.name = 'MIRI'
    input_model.meta.instrument.detector = 'MIRIFULONG'
    input_model.meta.instrument.channel = '34'
    input_model.meta.instrument.band = 'SHORT'

    input_model.data = np.zeros(shape)
    input_model.meta.wcs = dummy_wcs

    this_channel = '3'
    coord_system = 'skyalign'
    instrument_info = instrument_defaults.InstrumentInfo()
    instrument_info.SetXSliceLimits(0, 101, this_channel)
    x1, x2 = instrument_info.GetMIRISliceEndPts(this_channel)

    corners = cube_build_wcs_util.find_corners_MIRI(input_model,
                                                    this_channel,
                                                    instrument_info,
                                                    coord_system)

    (ra_min, b1, ra_max, b2, a1, dec_min, a2, dec_max,
     lambda_min, lambda_max) = corners
    assert ra_min == 40.6
    assert ra_max == 50.1
    assert dec_min == 45.1
    assert dec_max == 45.4
    assert lambda_min == 7.5
    assert lambda_max == 8.5
