"""Test Extract1dStep on MIRI MRS point source"""
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.mark.bigdata
def test_miri_mrs_extract1d(rtdata, fitsdiff_default_kwargs):
    """Test running extract_1d on an s3d cube containing a point source"""
    rtdata.get_data("miri/mrs/miri_003_det_image_seq1_MIRIFUSHORT_12SHORTexp1_s3d.fits")

    args = ["jwst.extract_1d.Extract1dStep", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "miri_003_det_image_seq1_MIRIFUSHORT_12SHORTexp1_extract1dstep.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_extract1d/miri_003_det_image_seq1_MIRIFUSHORT_12SHORTexp1_extract1dstep.fits')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
