"""Test suite for engdb_tools that require DB connection."""

import pytest
from astropy.time import Time

from jwst.lib.engdb_tools import ENGDB_Service

# Mark all tests in this module as slow due to remote DB connection
pytestmark = pytest.mark.slow

GOOD_MNEMONIC = "INRSI_GWA_Y_TILT_AVGED"
GOOD_STARTTIME = "2022-01-25 23:29:02.188"
GOOD_ENDTIME = "2022-01-26"

SHORT_STARTTIME = "2022-01-26 02:29:02.188"

BAD_MNEMONIC = "No_Such_MNEMONIC"
NODATA_STARTIME = "2014-01-01"
NODATA_ENDTIME = "2014-01-02"


class TestEngdbTools:
    """
    Class to test engdb_tools with DB connection.

    Because the DB constructor actually pings the service provider,
    we only want to do it once to prevent unnecessary server spam.
    """

    def setup_class(self):
        self.engdb = ENGDB_Service()

    def test_basic(self):
        assert self.engdb._get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)

    def test_values(self):
        engdb = self.engdb

        records = engdb._get_records(GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME)
        assert len(records) == 2

        values = engdb.get_values(GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME)
        assert len(values) == 10547
        assert values[0] == 0

        # test_values_with_bracket
        values = engdb.get_values(
            GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME, include_bracket_values=True
        )
        assert len(values) == 2
        assert values[1] == 0

    def test_values_with_time(self):
        values = self.engdb.get_values(
            GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME, include_obstime=True
        )
        assert len(values) >= 1
        assert isinstance(values[0], tuple)
        assert isinstance(values[0].obstime, Time)

    def test_novalues(self):
        values = self.engdb.get_values(GOOD_MNEMONIC, NODATA_STARTIME, NODATA_ENDTIME)
        assert len(values) == 0

    def test_meta(self):
        try:
            response = self.engdb.get_meta(GOOD_MNEMONIC)
        except NotImplementedError:
            pytest.skip("Test only valid with Direct EngDB connection.")
        assert response["Count"] == 1
        assert response["TlmMnemonics"][0]["TlmMnemonic"] == GOOD_MNEMONIC

    def test_unzip(self):
        """Test forunzipped versions of content."""
        values = self.engdb.get_values(
            GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME, include_obstime=True, zip_results=False
        )
        assert isinstance(values, tuple)
        assert len(values.obstime) == len(values.value)
