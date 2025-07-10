import warnings

import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep, nirspec


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "input_file",
    ["jw02072002001_05101_00001_nrs1_rate.fits", "jw01309022001_04102_00001_nrs2_rate.fits"],
    ids=["full", "allslits"],
)
def test_nirspec_fixedslit_wcs(rtdata, input_file):
    """Test NIRSpec fixed slit wcs"""
    rtdata.get_data(f"nirspec/fs/{input_file}")
    AssignWcsStep.call(input_file, save_results=True, suffix="assign_wcs")

    output = input_file.replace("rate", "assign_wcs")
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        # Check the 4 science slits
        for slit in ["S200A1", "S200A2", "S400A1", "S1600A1"]:
            wcs = nirspec.nrs_wcs_set_input(im, slit)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, slit)

            assert_wcs_grid_allclose(wcs, wcs_truth)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "input_file,msa_file",
    [
        ("jw01345066001_05101_00001_nrs1_rate.fits", "jw01345066001_01_msa.fits"),
        ("jw02674004001_03101_00001_nrs1_rate.fits", "jw02674004001_01_msa.fits"),
    ],
    ids=["mos", "mos_fs"],
)
def test_nirspec_mos_wcs(rtdata, input_file, msa_file):
    """Test NIRSpec MOS wcs"""
    # Get MSA meta file
    rtdata.get_data(f"nirspec/mos/{msa_file}")
    rtdata.get_data(f"nirspec/mos/{input_file}")
    AssignWcsStep.call(input_file, save_results=True, suffix="assign_wcs")

    output = input_file.replace("rate", "assign_wcs")
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        # Get validated open slits from WCS transform
        transform = im_truth.meta.wcs.get_transform("gwa", "slit_frame")
        open_slits = transform.slits[:]
        names = [slit.name for slit in open_slits]

        for name in names:
            wcs = nirspec.nrs_wcs_set_input(im, name)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, name)

            assert_wcs_grid_allclose(wcs, wcs_truth)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "input_file",
    ["jw01251004001_03107_00001_nrs1_rate.fits", "jw01251004001_03107_00001_nrs2_rate.fits"],
    ids=["nrs1", "nrs2"],
)
def test_nirspec_ifu_wcs(rtdata, input_file):
    """Test NIRSpec IFU wcs"""
    rtdata.get_data(f"nirspec/ifu/{input_file}")

    AssignWcsStep.call(input_file, save_results=True, suffix="assign_wcs", nrs_ifu_slice_wcs=True)

    output = input_file.replace("rate", "assign_wcs")
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_wcs/{output}")

    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        # Test all the IFU slices
        for k in range(30):
            wcs = nirspec.nrs_wcs_set_input(im, k)
            wcs_truth = nirspec.nrs_wcs_set_input(im_truth, k)
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
        assert_allclose(
            coord, coord_truth, err_msg=f"for coordinate axis '{wcs.output_frame.axes_names[n]}'"
        )
