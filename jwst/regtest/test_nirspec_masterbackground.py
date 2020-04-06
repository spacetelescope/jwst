import os
from typing import List

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

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
        Step.from_cmdline(args)

        return rtdata

    return _set_pipeline_run


def test_nirspec_fs_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Run a test for NIRSpec FS data with a user-supplied background file."""
    rtdata = run_nirspec_mbkg(
        'nirspec/fs/nrs_sci+bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/fs/v2_nrs_bkg_user_clean_x1d.fits')]
    )

    output = "nrs_sci+bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_fs_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_ifu_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Test NIRSpec IFU data with a user-supplied background file."""
    rtdata = run_nirspec_mbkg(
        'nirspec/ifu/prism_sci_bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/ifu/prism_bkg_x1d.fits')]
    )

    output = "prism_sci_bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_mos_mbkg_user(run_nirspec_mbkg, fitsdiff_default_kwargs):
    """Test NIRSpec MOS data with a user-supplied background file."""
    rtdata = run_nirspec_mbkg(
        'nirspec/mos/nrs_mos_sci+bkg_cal.fits',
        user_files=[('--user_background', 'nirspec/mos/v2_nrs_mos_bkg_x1d.fits')]
    )

    output = "nrs_mos_sci+bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_mos_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    'output_file',
    ['ifu_prism_source_on_NRS1_master_background.fits',
    'ifu_prism_source_off_NRS1_o001_masterbg.fits'],
    ids=["on-source", "off-source"]
)
def test_nirspec_ifu_mbkg_nod(run_nirspec_mbkg, fitsdiff_default_kwargs, output_file):
    """Test NIRSpec IFU nodded data."""
    rtdata = run_nirspec_mbkg(
        'nirspec/ifu/nirspec_spec3_asn.json',
        asn=True,
        options=['--save_background=True']
    )

    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_nod/{output_file}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
