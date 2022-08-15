import pytest
from astropy.io.fits.diff import FITSDiff
from pathlib import Path

from jwst.lib import engdb_tools
from jwst.lib.set_telescope_pointing import add_wcs
from jwst.lib.tests.engdb_mock import EngDB_Mocker

# Get the mock databases
DATA_PATH = Path(__file__).parents[1] / 'lib' / 'tests' / 'data'
ENGDB_PATH = DATA_PATH / 'engdb'


@pytest.mark.bigdata
def test_miri_setpointing(_jail, rtdata, engdb, fitsdiff_default_kwargs):
    """
    Regression test of the set_telescope_pointing script on a level-1b MIRI image.
    """

    # Get the input level-1b file
    rtdata.get_data("miri/image/jw80600010001_02101_00001_mirimage_uncal.fits")

    # The add_wcs function overwrites its input, so output = input
    rtdata.output = rtdata.input

    # Call the WCS routine, using the ENGDB_Service
    add_wcs(rtdata.input, allow_default=True, engdb_url='http://localhost', prd='PRDOPSSOC-055')

    # Compare the results
    rtdata.get_truth("truth/test_miri_setpointing/jw80600010001_02101_00001_mirimage_uncal.fits")
    fitsdiff_default_kwargs['rtol'] = 1e-6
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


# ########
# Fixtures
# ########
@pytest.fixture
def engdb():
    """Setup the mock engineering database"""
    with EngDB_Mocker(db_path=ENGDB_PATH):
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb
