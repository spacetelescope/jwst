"""
Set of various utilities to access the JWST
Engineering Database
"""

from astropy.time import Time
from collections import namedtuple
import logging
from os import getenv
import re
import requests

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# #############################################
# Where is the engineering service? Its HERE!!!
# #############################################

# This is currently set to the D-string hostname:
ENGDB_HOST = 'http://twjwdmsemwebag.stsci.edu/'
ENGDB_BASE_URL = ''.join([
    ENGDB_HOST,
    'JWDMSEngFqAccSide2/',
    'TlmMnemonicDataSrv.svc/',
])

# URI paths necessary to access the db
ENGDB_DATA = 'Data/'
ENGDB_DATA_XML = 'xml/Data/'
ENGDB_METADATA = 'MetaData/TlmMnemonics/'
ENGDB_METADATA_XML = 'xml/MetaData/TlmMnemonics/'

__all__ = [
    'ENGDB_Service'
]

# Define the returned value tuple.
_EngDB_Value = namedtuple('EngDB_Value', ['obstime', 'value'])


class _Value_Collection():
    """Engineering Value Collection

    Parameters
    ----------
    include_obstime: bool
        If `True`, the return values will include observation
        time as `astropy.time.Time`. See `zip` for further details.

    zip: bool
        If `True` and `include_obstime` is `True`, the return values
        will be a list of 2-tuples. If false, the return will
        be a single 2-tuple, where each element is a list.


    Attributes
    ----------
    collection: [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
        Returns the list of values.
        See `include_obstime` and `zip` for modifications.
    """
    def __init__(self, include_obstime=False, zip=True):
        self._include_obstime = include_obstime
        self._zip = zip
        if zip:
            self.collection = []
        else:
            self.collection = _EngDB_Value([], [])

    def append(self, obstime, value):
        """Append value to collection

        Parameters
        ----------
        obstime: int(milliseconds)
            Observation time as returned from the engineering
            db, in milliseconds

        value: numeric
            Value from db.

        Notes
        -----
        The `obstime` is converted to an `astropy.time.Time`
        """
        if self._include_obstime:
            obstime = Time(obstime / 1000., format='unix')
            if self._zip:
                self.collection.append(
                    _EngDB_Value(obstime, value)
                )
            else:
                self.collection.obstime.append(obstime)
                self.collection.value.append(value)
        else:
            self.collection.append(value)


class ENGDB_Service():
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

    def __init__(self, base_url=None, default_format='dict'):
        if base_url is None:
            base_url = getenv('ENG_BASE_URL', ENGDB_BASE_URL)
        if base_url[-1] !='/':
            base_url += '/'
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

        Notes
        -----
        The engineering service always returns the bracketing entries
        before and after the requested time range.
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
            include_obstime=False,
            include_bracket_values=False,
            zip=True
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic: str
            The engineering mnemonic to retrieve

        starttime: str or `astropy.time.Time`
            The, inclusive, start time to retireve from.

        endttime: str or `astropy.time.Time`
            The, inclusive, end time to retireve from.

        time_format: str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        include_obstime: bool
            If `True`, the return values will include observation
            time as `astropy.time.Time`. See `zip` for further details.

        include_bracket_values: bool
            The DB service, by default, returns the bracketing
            values outside of the requested time. If `True`, include
            these values.

        zip: bool
            If `True` and `include_obstime` is `True`, the return values
            will be a list of 2-tuples. If false, the return will
            be a single 2-tuple, where each element is a list.

        Returns
        -------
        values: [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
            Returns the list of values. See `include_obstime` and `zip` for modifications.

        Raises
        ------
        requests.exceptions.HTTPError
            Either a bad URL or non-existant mnemonic.
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
        results = _Value_Collection(
            include_obstime=include_obstime,
            zip=zip
        )
        if records['Data'] is not None:
            for record in records['Data']:
                obstime = extract_db_time(record['ObsTime'])
                if not include_bracket_values:
                    if obstime < db_starttime or obstime > db_endttime:
                        continue
                value = record['EUValue']
                results.append(obstime, value)

        return results.collection

    def get_meta(self, mnemonic='', result_format=None):
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
