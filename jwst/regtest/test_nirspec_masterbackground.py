import os
import pytest
import numpy as np

from astropy.io.fits.diff import FITSDiff
from numpy.testing import assert_allclose
from typing import List

from gwcs.wcstools import grid_from_bounding_box
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.assign_wcs import nirspec
from jwst.stpipe import Step
from jwst import datamodels

pytestmark = pytest.mark.bigdata


@pytest.fixture
def run_nirspec_mbkg(jail, rtdata_module):
    """Fixture-factory that runs the MasterBackground step for NIRSpec MasterBackground tests."""
    rtdata = rtdata_module
    collect_pipeline_cfgs('config')

    def _set_pipeline_run(input_file: str, user_files: List[tuple] = None, asn: bool = False,
                          options: List[str] = None):
        """Set additional pipeline options for the different MasterBackground tests.
        Additional user_files and their inputs (inputs are assumed to be artifactory paths) should be given as tuples
        for the user_files argument.
        Set asn to True if getting an association.
        Add additional command line options with the options argument.
        """
        # Create optional arguments list using rtdata to find artifact paths
        user_files = [
            f'{arg}={rtdata.get_data(input_path)}' for arg, input_path in user_files
        ] if user_files is not None else []

        if options is None:
            options = []

        # Get the input data
        if asn:
            rtdata.get_asn(input_file)

        else:
            rtdata.get_data(input_file)

        # Set arguments and execute
        args = ['config/master_background.cfg', rtdata.input, '--save_results=True'] + user_files + options
        step = Step.from_cmdline(args)

        return rtdata, step

    return _set_pipeline_run


def test_nirspec_fs_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Run a test for NIRSpec FS data with a user-supplied background file."""
    rtdata, _ = run_nirspec_mbkg(
        'nirspec/fs/nrs_sci+bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/fs/nrs_bkg_user_clean_x1d.fits')]
    )

    rtdata.output = "nrs_sci+bkg_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_mbkg_user", "nrs_sci+bkg_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_ifu_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Test NIRSpec IFU data with a user-supplied background file."""
    rtdata, step = run_nirspec_mbkg(
        'nirspec/ifu/prism_sci_bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/ifu/prism_bkg_x1d.fits')]
    )

    rtdata.output = 'prism_sci_bkg_master_background.fits'

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_nirspec_ifu_mbkg_user", 'prism_sci_bkg_master_background.fits'))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # Test 2  compare the science  data with no background to the output from the masterBackground Subtraction step
    #  background subtracted science image.
    input_sci_cal_file = rtdata.get_data('nirspec/ifu/prism_sci_cal.fits')
    input_sci_model = datamodels.open(input_sci_cal_file)

    # We don't want the slices gaps to impact the statisitic loop over the 30 Slices
    for i in range(30):
        slice_wcs = nirspec.nrs_wcs_set_input(input_sci_model, i)
        x, y = grid_from_bounding_box(slice_wcs.bounding_box)

        ra, dec, lam = slice_wcs(x, y)

        valid = np.isfinite(lam)

        result_slice_region = step.data[y.astype(int), x.astype(int)]
        sci_slice_region = input_sci_model.data[y.astype(int), x.astype(int)]

        sci_slice = sci_slice_region[valid]
        result_slice = result_slice_region[valid]
        sub = result_slice - sci_slice

        # check for outliers in the science image
        sci_mean = np.nanmean(sci_slice)
        sci_std = np.nanstd(sci_slice)

        upper = sci_mean + sci_std * 5.0
        lower = sci_mean - sci_std * 5.0

        mask_clean = np.logical_and(sci_slice < upper, sci_slice > lower)

        sub_mean = np.absolute(np.nanmean(sub[mask_clean]))
        atol = 2.0
        assert_allclose(sub_mean, 0, atol=atol)


def test_nirspec_mos_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Test NIRSpec MOS data with a user-supplied background file."""
    rtdata, _ = run_nirspec_mbkg(
        'nirspec/mos/nrs_mos_sci+bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/mos/nrs_mos_bkg_x1d.fits')]
    )

    rtdata.output = "nrs_mos_sci+bkg_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_nirspec_mos_mbkg_user", "nrs_mos_sci+master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_ifu_nod_mbkg(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Test NIRSpec IFU nodded data."""
    rtdata, step = run_nirspec_mbkg(
        'nirspec/ifu/nirspec_spec3_asn.json',
        asn=True,
        options=['--save_background=True']
    )

    rtdata.output = 'ifu_prism_source_off_NRS1_o001_masterbg.fits'

    # Get the truth file
    rtdata.get_truth(
        os.path.join("truth/test_nirspec_ifu_nod_mbkg", "ifu_prism_source_off_fix_NRS1_o001_masterbg.fits")
    )

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # test 2
    #  compare  background subtracted data  to truth files and check that the  cal_step master_background ran to
    #  completion
    outputs = []
    for model in step:
        assert model.meta.cal_step.master_background == 'COMPLETE'

        result_file = model.meta.filename.replace('cal', 'master_background')
        truth_file = rtdata.get_data(os.path.join('nirspec/ifu/nod/', result_file))

        outputs.append((result_file, truth_file))

    for output in outputs:
        diff = FITSDiff(*output, **fitsdiff_default_kwargs)
        assert diff.identical, diff.report()
