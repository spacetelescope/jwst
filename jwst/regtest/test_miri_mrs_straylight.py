"""Regression tests for MIRI MRS straylight step"""
import pytest

from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth file paths
INPUT_PATH = 'miri/mrs'
TRUTH_PATH = 'truth/test_miri_mrs_straylight'


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'suffix',
    ["straylightstep", "shower_model"],
)
def test_miri_mrs_straylight_clean_showers(rtdata, fitsdiff_default_kwargs,suffix):
    """Test running straylight with clean shower= True on an rate file."""

    filename = "jw01024001001_04101_00001_mirifulong_rate.fits"
    rtdata.get_data(INPUT_PATH + '/' + filename) 

    args = ['jwst.straylight.StraylightStep', rtdata.input,
            '--clean_showers=True',
            '--save_shower_model=True']
    Step.from_cmdline(args)

    output = "jw01024001001_04101_00001_mirifulong_" + suffix + ".fits"
    rtdata.output = output
    # Get the truth file
    rtdata.get_truth(f"{TRUTH_PATH}/{output}")
    
    # Compare the results for straylight output with clean_showers turned on.
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


