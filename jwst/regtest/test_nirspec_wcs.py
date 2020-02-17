import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel


@pytest.mark.bigdata
def test_nirspec_fixedslit_wcs(rtdata):

    input_file = 'jw00023001001_01101_00001_nrs1_rate.fits'
    rtdata.get_data('nirspec/test_wcs/' + input_file)
    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('rate', 'assign_wcs')
    rtdata.output = output

    rtdata.get_truth("truth/test_nirspec_wcs/" + output)

    im = ImageModel(rtdata.output)
    im_truth = ImageModel(rtdata.truth)

    for slit in ['S200A1']:
        wcs = nirspec.nrs_wcs_set_input(im, slit)
        wcs_truth = nirspec.nrs_wcs_set_input(im_truth, slit)
        assert_specwcs_grid_allclose(wcs, wcs_truth)


@pytest.mark.bigdata
def test_nirspec_mos_wcs(rtdata):

    input_file = 'msa_patt_num.fits'
    rtdata.get_data('nirspec/mos/V9621500100101_short_msa.fits')
    rtdata.get_data('nirspec/test_wcs/' + input_file)
    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('.fits', '_assign_wcs.fits')
    rtdata.output = output

    truth = input_file.replace('.fits', '_truth_assign_wcs.fits')
    rtdata.get_truth("truth/test_nirspec_wcs/" + truth)

    im = ImageModel(rtdata.output)
    im_truth = ImageModel(rtdata.truth)

    slits = nirspec.get_open_slits(im)
    names = [slit.name for slit in slits]
    for name in names:
        wcs = nirspec.nrs_wcs_set_input(im, name)
        wcs_truth = nirspec.nrs_wcs_set_input(im_truth, name)
        assert_specwcs_grid_allclose(wcs, wcs_truth)


params = [
    'jw00011001001_01120_00001_nrs1_rate.fits',
    'jw00011001001_01120_00002_nrs1_rate.fits',
    'jw00011001001_01120_00003_nrs2_rate.fits',
]
ids = ['ifu_nrs1', 'ifu_nrs1_opaque', 'ifu_nrs2']

@pytest.mark.parametrize("input_file", params, ids=ids)
@pytest.mark.bigdata
def test_nirspec_ifu_wcs(input_file, rtdata):

    rtdata.get_data('nirspec/test_wcs/' + input_file)
    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('rate.fits', 'assign_wcs.fits')
    rtdata.output = output

    rtdata.get_truth("truth/test_nirspec_wcs/" + output)

    im = ImageModel(rtdata.output)
    im_truth = ImageModel(rtdata.truth)

    for slice_ in range(29):
        wcs = nirspec.nrs_wcs_set_input(im, slice_)
        wcs_truth = nirspec.nrs_wcs_set_input(im_truth, slice_)
        assert_specwcs_grid_allclose(wcs, wcs_truth)


def assert_specwcs_grid_allclose(wcs, wcs_truth):
    __traceback__ = False

    # Compute RA, Dec, lambda values for each image array
    grid = grid_from_bounding_box(wcs.bounding_box)
    ra, dec, lam = wcs(*grid)
    grid_truth = grid_from_bounding_box(wcs_truth.bounding_box)
    ra_truth, dec_truth, lam_truth = wcs_truth(*grid_truth)

    # Compare the sky coordinates
    # equal_nan is used, because many of the entries are NaN,
    # due to the bounding_box being rectilinear while the
    # defined spectral traces are curved
    assert_allclose(ra, ra_truth, equal_nan=True)
    assert_allclose(dec, dec_truth, equal_nan=True)
    assert_allclose(lam, lam_truth, equal_nan=True)
