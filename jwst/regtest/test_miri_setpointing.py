import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.lib.set_telescope_pointing import add_wcs


@pytest.mark.bigdata
def test_miri_setpointing(rtdata, fitsdiff_default_kwargs):
    """
    Regression test of the set_telescope_pointing script on a level-1b MIRI image.
    """

    # Get the input level-1b file
    rtdata.get_data("miri/mrs/jw01282004001_02101_00001_mirifulong_uncal.fits")

    # The add_wcs function overwrites its input, so output = input
    rtdata.output = rtdata.input

    # Call the WCS routine, using the ENGDB_Service
    add_wcs(rtdata.input)

    # Compare the results
    rtdata.get_truth("truth/test_miri_setpointing/jw01282004001_02101_00001_mirifulong_uncal.fits")
    fitsdiff_default_kwargs["rtol"] = 1e-6
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
