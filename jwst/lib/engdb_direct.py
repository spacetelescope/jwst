"""
Access the JWST Engineering Mnemonic Database through direct connection
"""

from astropy.time import Time
import logging
from os import getenv
import re
import requests
from requests.adapters import HTTPAdapter, Retry

from .engdb_lib import EngDB_Value, EngdbABC, FORCE_STATUSES, RETRIES, TIMEOUT

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# #############################################
# Where is the engineering service? Its HERE!!!
# #############################################

# URI paths necessary to access the db
ENGDB_DATA = 'Data/'
ENGDB_DATA_XML = 'xml/Data/'
ENGDB_METADATA = 'MetaData/TlmMnemonics/'
ENGDB_METADATA_XML = 'xml/MetaData/TlmMnemonics/'

__all__ = [
    'EngdbDirect'
]


class EngdbDirect(EngdbABC):
    """
    Access the JWST Engineering Database through direct connection

    Parameters
    ----------
    base_url : str
        The base url for the engineering RESTful service

    default_format : str
        The format the results of the data should be returned from the service.
        If 'dict', the result will be in Python dict format.

    service_kwargs : dict
        Service-specific keyword arguments that are not relevant to this implementation
        of EngdbABC.
    """

    #: The base URL for the engineering service.
    base_url = None

    #: The format the results of the data should be returned from the service.
    default_format = None

    #: The end time of the last query.
    endtime = None

    #: The results of the last query.
    response = None

    #: The start time of the last query.
    starttime = None

    def __init__(self, base_url=None, default_format='dict', **service_kwargs):
        logger.debug('kwargs not used by this service: %s', service_kwargs)

        self.configure(base_url=base_url)

        self.default_format = default_format

        self.set_session()

        # Check for aliveness
        response = self._session.get(''.join([
            self.base_url,
            self.default_format,
            ENGDB_METADATA
        ]), timeout=self.timeout)
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

    def configure(self, base_url=None):
        """Configure from parameters and environment

        Parameters
        ----------
        base_url : str
            The base url for the engineering RESTful service
        """
        # Determine the database to use.
        if base_url is None:
            base_url = getenv('ENG_BASE_URL')
        if not base_url:
            raise RuntimeError('No engineering database URL given.')
        if base_url[-1] != '/':
            base_url += '/'
        self.base_url = base_url

        # Get various timeout parameters
        self.retries = getenv('ENG_RETRIES', RETRIES)
        self.timeout = getenv('ENG_TIMEOUT', TIMEOUT)

    def get_meta(self, mnemonic='', result_format=None):
        """Get the mnemonics meta info

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve

        result_format : str
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
        response = self._session.get(query, timeout=self.timeout)
        logger.debug('Response="{}"'.format(response))
        response.raise_for_status()

        # That's all folks
        self.response = response
        return response.json()

    def get_values(
            self,
            mnemonic,
            starttime,
            endtime,
            time_format=None,
            include_obstime=False,
            include_bracket_values=False,
            zip_results=True
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve

        starttime : str or `astropy.time.Time`
            The, inclusive, start time to retrieve from.

        endtime : str or `astropy.time.Time`
            The, inclusive, end time to retrieve from.

        time_format : str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        include_obstime : bool
            If `True`, the return values will include observation
            time as `astropy.time.Time`. See `zip_results` for further details.

        include_bracket_values : bool
            The DB service, by default, returns the bracketing
            values outside of the requested time. If `True`, include
            these values.

        zip_results : bool
            If `True` and `include_obstime` is `True`, the return values
            will be a list of 2-tuples. If false, the return will
            be a single 2-tuple, where each element is a list.

        Returns
        -------
        values : [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
            Returns the list of values. See `include_obstime` and `zip_results` for modifications.

        Raises
        ------
        requests.exceptions.HTTPError
            Either a bad URL or non-existant mnemonic.
        """
        records = self._get_records(
            mnemonic=mnemonic,
            starttime=starttime,
            endtime=endtime,
            time_format=time_format
        )

        # Records returned are apparent not strictly correlated with
        # observation time. So, need to filter further.
        db_starttime = extract_db_time(records['ReqSTime'])
        db_endtime = extract_db_time(records['ReqETime'])
        results = _ValueCollection(
            include_obstime=include_obstime,
            zip_results=zip_results
        )
        if records['Data'] is not None:
            for record in records['Data']:
                obstime = extract_db_time(record['ObsTime'])
                if not include_bracket_values:
                    if obstime < db_starttime or obstime > db_endtime:
                        continue
                value = record['EUValue']
                results.append(obstime, value)

        return results.collection

    def set_session(self):
        """Setup HTTP session"""
        s = requests.Session()
        retries = Retry(total=10, backoff_factor=1.0, status_forcelist=FORCE_STATUSES, raise_on_status=True)
        s.mount('https://', HTTPAdapter(max_retries=retries))

        self._session = s

    def _get_records(
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
        mnemonic : str
            The engineering mnemonic to retrieve

        starttime : str or astropy.time.Time
            The, inclusive, start time to retrieve from.

        endtime : str or astropy.time.Time
            The, inclusive, end time to retrieve from.

        result_format : str
            The format to request from the service.
            If None, the `default_format` is used.

        time_format : str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        Returns
        -------
        records : dict
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
        self.starttime = starttime
        self.endtime = endtime

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
        response = self._session.get(query, timeout=self.timeout)
        logger.debug('Response: %s', response)
        logger.debug('Response: %s', response.json())
        response.raise_for_status()

        # That's all folks
        self.response = response
        return response.json()


class _ValueCollection():
    """Engineering Value Collection

    Parameters
    ----------
    include_obstime : bool
        If `True`, the return values will include observation
        time as `astropy.time.Time`. See `zip_results` for further details.

    zip_results : bool
        If `True` and `include_obstime` is `True`, the return values
        will be a list of 2-tuples. If false, the return will
        be a single 2-tuple, where each element is a list.


    Attributes
    ----------
    collection : [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
        Returns the list of values.
        See `include_obstime` and `zip_results` for modifications.
    """
    def __init__(self, include_obstime=False, zip_results=True):
        self._include_obstime = include_obstime
        self._zip_results = zip_results

        self.obstimes = []
        self.values = []

    @property
    def collection(self):
        if self._include_obstime:
            obstimes = Time(self.obstimes, format='unix')
            obstimes.format = 'isot'
            if self._zip_results:
                collection = [EngDB_Value(t, v) for t, v in zip(obstimes, self.values)]
            else:
                collection = EngDB_Value(obstimes, self.values)
        else:
            collection = self.values
        return collection

    def append(self, obstime, value):
        """Append value to collection

        Parameters
        ----------
        obstime : int(milliseconds)
            Observation time as returned from the engineering
            db, in milliseconds

        value : numeric
            Value from db.

        Notes
        -----
        The `obstime` is converted to an `astropy.time.Time`
        """
        # Convert from milliseconds to seconds before appending.
        self.obstimes.append(obstime / 1000.)
        self.values.append(value)


# #########
# Utilities
# #########
def extract_db_time(db_date):
    """Extract date from date string in the Database

    Parameters
    ----------
    db_date : str
        The string from a date field in the database.

    Returns
    -------
    time : int
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
