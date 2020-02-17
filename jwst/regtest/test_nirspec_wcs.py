import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel

test_data = {'ifu_nrs1':'jw00011001001_01120_00001_nrs1_rate.fits',
             'ifu_nrs1_opaque':'jw00011001001_01120_00002_nrs1_rate.fits',
             'ifu_nrs2': 'jw00011001001_01120_00003_nrs2_rate.fits',
             'fs_nrs1': 'jw00023001001_01101_00001_nrs1_rate.fits',
             'mos_nrs1':'msa_patt_num.fits'}

@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run assign_wcs"""
    def _run_pipe_with_file(input_file):
        rtdata = rtdata_module
        
        if input_file == 'msa_patt_num.fits':
            rtdata.get_data('nirspec/mos/V9621500100101_short_msa.fits')

        rtdata.get_data('nirspec/test_wcs/' + input_file)

        AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')
        
        return rtdata
    
    return _run_pipe_with_file

@pytest.mark.bigdata
def test_nirspec_fixedslit_wcs(run_pipeline):
    
    input_file = test_data['fs_nrs1']
    rtdata = run_pipeline(input_file)

    output = input_file.replace('rate', 'assign_wcs')
    rtdata.output = output
        
    rtdata.get_truth("truth/test_nirspec_wcs/" + output)
    
    im = ImageModel(rtdata.output)
    im_ref = ImageModel(rtdata.truth)
    
    print(rtdata.output)
    print(rtdata.truth)
    print(im.meta.exposure.type)

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

@pytest.mark.bigdata
def test_nirspec_mos_wcs(run_pipeline):
    
    input_file = test_data['mos_nrs1']
    rtdata = run_pipeline(input_file)

    output = input_file.replace('.fits', '_assign_wcs.fits')
    rtdata.output = output
    
    truth = input_file.replace('.fits', '_truth_assign_wcs.fits')
    rtdata.get_truth("truth/test_nirspec_wcs/" + truth)
    
    im = ImageModel(rtdata.output)
    im_ref = ImageModel(rtdata.truth)
    
    # Get the WCS test data
    slits = nirspec.get_open_slits(im)
    name = slits[0].name
    wcs = nirspec.nrs_wcs_set_input(im, name)

    # Get the WCS for truth data
    slits_ref = nirspec.get_open_slits(im)
    name_ref = slits_ref[0].name
    wcs_ref = nirspec.nrs_wcs_set_input(im_ref, name_ref)

    # Compute RA, Dec, lambda values for each image array
    grid = grid_from_bounding_box(wcs.bounding_box)
    ra, dec, lam = wcs(*grid)
    ra_ref, dec_ref, lam_ref = wcs_ref(*grid)

    # Compare the sky coordinates
    assert_allclose(ra, ra_ref, equal_nan=True)
    assert_allclose(dec, dec_ref, equal_nan=True)
    assert_allclose(lam, lam_ref, equal_nan=True)

@pytest.mark.bigdata
def test_nirspec_ifu_wcs(run_pipeline):
    
    input_file = test_data['ifu_nrs1']
    rtdata = run_pipeline(input_file)

    output = input_file.replace('rate.fits', 'assign_wcs.fits')
    rtdata.output = output
    
    rtdata.get_truth("truth/test_nirspec_wcs/" + output)
    
    im = ImageModel(rtdata.output)
    im_ref = ImageModel(rtdata.truth)
    
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

@pytest.mark.bigdata
def test_nirspec_ifu_nrs2_wcs(run_pipeline):
    
    input_file = test_data['ifu_nrs2']
    rtdata = run_pipeline(input_file)

    output = input_file.replace('rate.fits', 'assign_wcs.fits')
    rtdata.output = output
    
    rtdata.get_truth("truth/test_nirspec_wcs/" + output)
    
    im = ImageModel(rtdata.output)
    im_ref = ImageModel(rtdata.truth)
    
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

@pytest.mark.bigdata
def test_nirspec_ifu_opaque_wcs(run_pipeline):
    
    input_file = test_data['ifu_nrs1_opaque']
    rtdata = run_pipeline(input_file)

    output = input_file.replace('rate.fits', 'assign_wcs.fits')
    rtdata.output = output
    
    rtdata.get_truth("truth/test_nirspec_wcs/" + output)
    
    im = ImageModel(rtdata.output)
    im_ref = ImageModel(rtdata.truth)
    
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
