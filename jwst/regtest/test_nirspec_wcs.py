import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst import datamodels


@pytest.mark.bigdata
def test_nirspec_fixedslit_wcs(rtdata):
    """Test NIRSpec fixed slit wcs"""
    input_file = 'jw00023001001_01101_00001_nrs1_rate.fits'
    rtdata.get_data(f"nirspec/fs/{input_file}")
    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('rate', 'assign_wcs')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        # Check the 4 science slits
        for slit in ['S200A1', 'S200A2', 'S400A1', 'S1600A1']:
            wcs = nirspec.nrs_wcs_set_input(im, slit)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, slit)

            assert_wcs_grid_allclose(wcs, wcs_truth)


@pytest.mark.bigdata
def test_nirspec_mos_wcs(rtdata):
    """Test NIRSpec MOS wcs"""
    input_file = 'msa_patt_num.fits'
    # Get MSA meta file
    rtdata.get_data('nirspec/mos/V9621500100101_short_msa.fits')
    rtdata.get_data(f"nirspec/mos/{input_file}")
    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('.fits', '_assign_wcs.fits')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        names = [slit.name for slit in nirspec.get_open_slits(im)]
        for name in names:
            wcs = nirspec.nrs_wcs_set_input(im, name)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, name)

            assert_wcs_grid_allclose(wcs, wcs_truth)


@pytest.mark.parametrize("input_file",
    [
    'jw00011001001_01120_00001_nrs1_rate.fits',
    'jw00011001001_01120_00002_nrs1_rate.fits',
    'jw00011001001_01120_00003_nrs2_rate.fits',
    ],
    ids=['nrs1_f170lp', 'nrs1_opaque', 'nrs2_f170lp'])
@pytest.mark.bigdata
def test_nirspec_ifu_wcs(input_file, rtdata):
    """Test NIRSpec IFU wcs"""
    rtdata.get_data(f"nirspec/ifu/{input_file}")

    AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')

    output = input_file.replace('rate.fits', 'assign_wcs.fits')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        # Test several slices in the IFU, range(30)
        for slice_ in [0, 9, 16, 23, 29]:
            wcs = nirspec.nrs_wcs_set_input(im, slice_)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, slice_)

            assert_wcs_grid_allclose(wcs, wcs_truth)


def assert_wcs_grid_allclose(wcs, wcs_truth):
    """Assertion helper verifying the RA/DEC/(lam) are the same for 2 WCSs"""
    __tracebackhide__ = True
    # Compute RA, Dec[, lambda] values at each pixel in bounding box.
    grid = grid_from_bounding_box(wcs.bounding_box)
    grid_truth = grid_from_bounding_box(wcs_truth.bounding_box)
    # Store in tuple (RA, Dec, lambda)
    skycoords = wcs(*grid)
    skycoords_truth = wcs_truth(*grid_truth)

    # Compare each RA, Dec[, lambda] grid
    for n, (coord, coord_truth) in enumerate(zip(skycoords, skycoords_truth)):
        assert_allclose(coord, coord_truth,
            err_msg=f"for coordinate axis '{wcs.output_frame.axes_names[n]}'")
