import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

pytestmark = pytest.mark.bigdata


def test_nirspec_spec2_mbkg(rtdata, fitsdiff_default_kwargs):
    """Run spec2 with master background"""

    # Get data
    rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_msa.fits')
    rtdata.get_data('nirspec/mos/nrs_mos_with_bkgslits_rate.fits')
    collect_pipeline_cfgs('config')

    # Run the pipeline
    args = [
        'config/calwebb_spec2.cfg',
        rtdata.input,
        '--steps.master_background.skip=false'
    ]
    Step.from_cmdline(args)

    # Compare results
    rtdata.output = "nrs_mos_with_bkgslits_cal.fits"
    rtdata.get_truth('truth/test_nirspec_mos_mbkg_user/nrs_mos_with_bkgslits_cal.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_fs_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Run a test for NIRSpec FS data with a user-supplied background file."""

    # Get user-supplied background
    user_background = "v2_nrs_bkg_user_clean_x1d.fits"
    rtdata.get_data(f"nirspec/fs/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/fs/nrs_sci+bkg_cal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

    output = "nrs_sci+bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_fs_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_ifu_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Test NIRSpec IFU data with a user-supplied background file."""
    # Get user-supplied background
    user_background = "prism_bkg_x1d.fits"
    rtdata.get_data(f"nirspec/ifu/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/ifu/prism_sci_bkg_cal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

    output = "prism_sci_bkg_master_background.fits"
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_user/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nirspec_mos_mbkg_user(rtdata, fitsdiff_default_kwargs):
    """Test NIRSpec MOS data with a user-supplied background file."""
    # Get user-supplied background
    user_background = "v2_nrs_mos_bkg_x1d.fits"
    rtdata.get_data(f"nirspec/mos/{user_background}")

    # Get input data
    rtdata.get_data("nirspec/mos/nrs_mos_sci+bkg_cal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--user_background", user_background]
    Step.from_cmdline(args)

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
def test_nirspec_ifu_mbkg_nod(rtdata, fitsdiff_default_kwargs, output_file):
    """Test NIRSpec IFU prism nodded data."""
    # Get input data
    rtdata.get_asn("nirspec/ifu/nirspec_spec3_asn.json")

    collect_pipeline_cfgs("config")
    args = ["config/master_background.cfg", rtdata.input,
            "--save_background=True"]
    Step.from_cmdline(args)

    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(f"truth/test_nirspec_ifu_mbkg_nod/{output_file}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
