import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.set_telescope_pointing import add_wcs
from jwst.lib import engdb_tools


@pytest.mark.bigdata
def test_miri_setpointing(_jail, rtdata, fitsdiff_default_kwargs):
    """
    Regression test of the set_telescope_pointing script on a level-1b MIRI image.
    """

    # Get SIAF PRD database file
    siaf_path = rtdata.get_data("common/prd.db")

    # Get the input level-1b file
    rtdata.get_data("miri/image/jw80600010001_02101_00001_mirimage_uncal.fits")

    # The add_wcs function overwrites its input, so output = input
    rtdata.output = rtdata.input

    # Call the WCS routine, using the ENGDB_Service
    # Note that there aren't any quaternion mnemonics in the ENGDB for the time
    # range of this exposure, so we set "allow_default" to tell add_wcs to use
    # default values for the pointing-related keywords. Keyword values retrieved
    # from the SIAF will be good.
    add_wcs(rtdata.input, allow_default=True,
            siaf_path=siaf_path,
            engdb_url=engdb_tools.ENGDB_BASE_URL)

    # Compare the results
    rtdata.get_truth("truth/test_miri_setpointing/jw80600010001_02101_00001_mirimage_uncal.fits")
    fitsdiff_default_kwargs['rtol'] = 1e-6
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
