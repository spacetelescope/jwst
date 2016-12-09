"""
Set of various utilities to access the JWST
Engineering Database
"""

from astropy.time import Time
from collections import namedtuple
import logging
import re
import requests

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# #############################################
# Where is the engineering service? Its HERE!!!
# #############################################
ENGDB_HOST = 'http://iwjwdmsdemwebv.stsci.edu/'
ENGDB_BASE_URL = ''.join([
    ENGDB_HOST,
    'JWDMSEngFqAccB7_testFITSw/',
    'TlmMnemonicDataSrv.svc/',
])

# URI paths necessary to access the db
ENGDB_DATA = 'Data/'
ENGDB_DATA_XML = 'xml/Data/'
ENGDB_METADATA = 'MetaData/TlmMnemonics/'
ENGDB_METADATA_XML = 'xml/MetaData/TlmMnemonics/'

__all__ = [
    'ENGDB_Service',
    'EngDB_Value'
]

# Define the returned value tuple.
EngDB_Value = namedtuple('EngDB_Value', ['obstime', 'value'])


class ENGDB_Service(object):
    """
    Set of various utilities to access the JWST
    Engineering Database

    Parameters
    ----------
    base_url: str
        The base url for the engineering RESTful service

    default_format: str
        The format the results of the data should be returned.
        If 'dict', the result will be in Python dict format.

    Attributes
    ----------
    response: `requests.response`
        The results of the last query.

    starttime: `astropy.time.Time`
        The start time of the last query.

    endtime: `astropy.time.Time`

    base_url: str
        The base URL for the engineering service.

    default_format: str
        The format to retrieve from the service.
        This is not the format of the returned data.
    """

    def __init__(self, base_url=ENGDB_BASE_URL, default_format='dict'):
        self.base_url = base_url
        self.default_format = default_format

        # Check for aliveness
        response = requests.get(''.join([
            self.base_url,
            self.default_format,
            ENGDB_METADATA
        ]))
        response.raise_for_status()

    @property
    def default_format(self):
        return self._default_format

    @default_format.setter
    def default_format(self, result_format):
        result_format += '/'
        if result_format == 'dict/':
            result_format = ''
        self._default_format = result_format

    def get_records(
            self,
            mnemonic,
            starttime,
            endtime,
            result_format=None,
            time_format=None
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic: str
            The engineering mnemonic to retrieve

        starttime: str or astropy.time.Time
            The, inclusive, start time to retireve from.

        endttime: str or astropy.time.Time
            The, inclusive, end time to retireve from.

        result_format: str
            The format to request from the service.
            If None, the `default_format` is used.

        time_format: str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        Returns
        -------
        records: dict
            Returns the dict of the request. This includes all
            the data returned form the DB concerning the requested
            mnemonic.

        Raises
        ------
        requests.exceptions.HTTPError
            Either a bad URL or non-existant mnemonic.
        """
        if result_format is None:
            result_format = self.default_format

        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)

        # Build the URL
        query = ''.join([
            self.base_url,
            result_format,
            'Data/',
            mnemonic,
            '?sTime=',
            starttime.iso,
            '&eTime=',
            endtime.iso,
        ])
        logger.debug('Query URL="{}"'.format(query))

        # Make our request
        response = requests.get(query)
        logger.debug('Response="{}"'.format(response))
        response.raise_for_status()

        # That's all folks
        self.response = response
        self.starttime = starttime
        self.endtime = endtime
        return response.json()

    def get_values(
            self,
            mnemonic,
            starttime,
            endtime,
            time_format=None,
            include_obstime=False
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic: str
            The engineering mnemonic to retrieve

        starttime: str or astropy.time.Time
            The, inclusive, start time to retireve from.

        endttime: str or astropy.time.Time
            The, inclusive, end time to retireve from.

        time_format: str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        include_obstime: bool
            If True, the return values will be a 2-tuple of
            (astropy.time.Time, value)

        Returns
        -------
        values: list
            Returns the list of values. If `include_obstime` is True
            the values will be a 2-tuple of (astropy.time.Time, value)
            instead of just value.

        Raises
        ------
        requests.exceptions.HTTPError
            Either a bad URL or non-existant mnemonic.

        ValueError
            Mnemonic is found but has no data.
        """
        records = self.get_records(
            mnemonic=mnemonic,
            starttime=starttime,
            endtime=endtime,
            time_format=time_format
        )

        # Records returned are apparent not strictly correlated with
        # observation time. So, need to filter further.
        db_starttime = extract_db_time(records['ReqSTime'])
        db_endttime = extract_db_time(records['ReqETime'])
        results = []
        if records['Data'] is None:
            raise ValueError('Mnemonic {} has no data'.format(mnemonic))
        for record in records['Data']:
            obstime = extract_db_time(record['ObsTime'])
            if obstime >= db_starttime and obstime <= db_endttime:
                value = record['EUValue']
                if include_obstime:
                    result = EngDB_Value(
                        obstime=Time(obstime / 1000., format='unix'),
                        value=value
                    )
                else:
                    result = value
                results.append(result)

        return results

    def get_meta(self, mnemonic, result_format=None):
        """Get the menonics meta info

        Parameters
        ----------
        mnemonic: str
            The engineering mnemonic to retrieve

        result_format: str
            The format to request from the service.
            If None, the `default_format` is used.
        """
        if result_format is None:
            result_format = self.default_format

        query = ''.join([
            self.base_url,
            result_format,
            ENGDB_METADATA,
            mnemonic
        ])
        logger.debug('Query URL="{}"'.format(query))

        # Make our request
        response = requests.get(query)
        logger.debug('Response="{}"'.format(response))
        response.raise_for_status()

        # That's all folks
        self.response = response
        return response.json()


# #########
# Utilities
# #########
def extract_db_time(db_date):
    """Extract date from date string in the Database

    Parameters
    ----------
    db_date: str
        The string from a date field in the database.

    Returns
    -------
    time: int
        The UNIX time in milliseconds

    Notes
    -----
    The format of time in the database is

        /Date(1234567890123+1234)/

    where the plus could be a minus. What is returned is
    the actual 13 digit number before the plus/minus.
    """
    match = re.match(r'\/Date\((\d{13})(.*)\)\/', db_date)
    milliseconds = int(match.group(1))

    return milliseconds
