"""Regression tests for MIRI MRS straylight step"""
import pytest

from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth file paths
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs_straylight'


@pytest.mark.bigdata
def test_miri_mrs_straylight_clean_showers(rtdata, fitsdiff_default_kwargs):
    """Test running straylight with clean shower= True on an rate file."""

    rtdata.get_data("miri/mrs/jw01024001001_04101_00001_mirifulong_rate.fits")

    args = ['jwst.straylight.StraylightStep', rtdata.input,
            '--clean_showers=True',
            '--save_shower_model=True']
    Step.from_cmdline(args)
    rtdata.output = "jw01024001001_04101_00001_mirifulong_straylightstep.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_straylight/jw01024001001_04101_00001_mirifulong_straylightstep.fits')

    # Compare the results for straylight output with clean_showers turned on.
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # Compare the results for shower_model
    rtdata.output = "jw01024001001_04101_00001_mirifulong_shower_model.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_straylight/jw01024001001_04101_00001_mirifulong_shower_model.fits')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


