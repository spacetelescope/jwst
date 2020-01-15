import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel

test_data = [
    ('ifu_nrs1', 'jw00011001001_01120_00001_nrs1_rate.fits', 'jw00011001001_01120_00001_nrs1_assign_wcs.fits'),
    ('ifu_nrs1_opaque', 'jw00011001001_01120_00002_nrs1_rate.fits', 'jw00011001001_01120_00002_nrs1_assign_wcs.fits'),
    ('ifu_nrs2', 'jw00011001001_01120_00003_nrs2_rate.fits', 'jw00011001001_01120_00003_nrs2_assign_wcs.fits'),
    ('fs_nrs1', 'jw00023001001_01101_00001_nrs1_rate.fits', 'jw00023001001_01101_00001_nrs1_assign_wcs.fits')]

@pytest.mark.bigdata
@pytest.mark.parametrize("test_id, input_file, truth_file", test_data)
def test_nirspec_wcs(_jail, rtdata, test_id, input_file, truth_file):
    """
        Test of the AssignWcs step on 4 different NIRSpec exposures:
        1) IFU NRS1 exposure,
        2) IFU NRS1 exposure with FILTER=OPAQUE,
        3) IFU NRS2 exposure, and
        4) FS NRS1 exposure with 4 slits.
    """

    # Get the input and truth files
    rtdata.get_data('nirspec/test_wcs/' + input_file)
    rtdata.get_truth('truth/test_nirspec_wcs/' + truth_file)

    # Run the AssignWcs step
    result = AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')
    result.close()

    # Open the output and truth files
    im = ImageModel(result.meta.filename)
    im_ref = ImageModel(truth_file)

    if result.meta.exposure.type == 'NRS_FIXEDSLIT':

        # Loop over the 4 slit instances
        for slit in ['S200A1', 'S200A2', 'S400A1', 'S1600A1']:

            # Create WCS objects for each image
            wcs = nirspec.nrs_wcs_set_input(im, slit)
            wcs_ref = nirspec.nrs_wcs_set_input(im_ref, slit)

            # Compute RA, Dec, lambda values for each image array
            grid = grid_from_bounding_box(wcs.bounding_box)
            ra, dec, lam = wcs(*grid)
            ra_ref, dec_ref, lam_ref = wcs_ref(*grid)

            # Compare the sky coordinates
            assert_allclose(ra, ra_ref, equal_nan=True)
            assert_allclose(dec, dec_ref, equal_nan=True)
            assert_allclose(lam, lam_ref, equal_nan=True)

    else:

        # Create WCS objects for each image
        wcs = nirspec.nrs_wcs_set_input(im, 0)
        wcs_ref = nirspec.nrs_wcs_set_input(im_ref, 0)

        # Compute RA, Dec, lambda values for each image array
        grid = grid_from_bounding_box(wcs.bounding_box)
        ra, dec, lam = wcs(*grid)
        ra_ref, dec_ref, lam_ref = wcs_ref(*grid)

        # Compare the sky coordinates
        # equal_nan is used, because many of the entries are NaN,
        # due to the bounding_box being rectilinear while the
        # defined spectral traces are curved
        assert_allclose(ra, ra_ref, equal_nan=True)
        assert_allclose(dec, dec_ref, equal_nan=True)
        assert_allclose(lam, lam_ref, equal_nan=True)
