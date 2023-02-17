"""Test Extract1dStep on MIRI point source with user supplied center and extended source """
# test_miri_mrs_spec3 already tests point source extract1d source
# In this test we will test providing the extraction center  and extended source extraction.
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.mark.bigdata
def test_miri_mrs_extract1d_nominal(rtdata, fitsdiff_default_kwargs):
    """Test running extract_1d on an s3d cube containing a point source"""
    # input s3d are created using the same data that was used in
    # test_miri_mrs_spec3_ifushort: run calwebb_spec3 on
    # miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json

    rtdata.get_data("miri/mrs/jw01024-c1000_t002_miri_ch2-mediumlong_s3d.fits")

    args = ["jwst.extract_1d.Extract1dStep", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw01024-c1000_t002_miri_ch2-mediumlong_extract1dstep.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_extract1d/jw01024-c1000_t002_miri_ch2-mediumlong_extract1dstep.fits')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_mrs_extract1d_center(rtdata, fitsdiff_default_kwargs):
    """Test running extract_1d on an s3d cube containing a point source with user-supplied center"""
    # input s3d are created using the same data that was used in
    # input s3d are created using the same data that was used in
    # test_miri_mrs_spec3_ifushort: run calwebb_spec3 on
    # using  ch2 data

    rtdata.get_data("miri/mrs/jw01024-c1000_t002_miri_ch2-mediumlong_s3d.fits")

    args = ['jwst.extract_1d.Extract1dStep', rtdata.input,
            '--output_file=jw01024-c1000_t002_miri_ch2-mediumlong_center',
            '--center_xy=33,27']
    Step.from_cmdline(args)
    rtdata.output = "jw01024-c1000_t002_miri_ch2-mediumlong_center_extract1dstep.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_extract1d/jw01024-c1000_t002_miri_ch2-mediumlong_center_extract1dstep.fits')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_mrs_extract1d_extended(rtdata, fitsdiff_default_kwargs):
    """Test running extract_1d on an s3d cube for extended source"""
    # input s3d are created using the same data that was used in
    # s3d IFULONG data created by running calspec3 on  IFULONG data in
    # jw01355-o005_20230109t002554_spec3_00001_asn.json

    rtdata.get_data("miri/mrs/jw01355-o005_t010_miri_ch3-long_s3d.fits")

    args = ["jwst.extract_1d.Extract1dStep", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw01355-o005_t010_miri_ch3-long_extract1dstep.fits"

    # Get the truth file
    rtdata.get_truth('truth/test_miri_mrs_extract1d/jw01355-o005_t010_miri_ch3-long_extract1dstep.fits')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
