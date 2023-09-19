"""
Access the JWST Engineering Mnemonic Database through MAST
"""
import json
import logging
from os import getenv
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter, Retry
from shutil import copy2

from astropy.table import Table
from astropy.time import Time
import numpy as np

from .engdb_lib import EngDB_Value, EngdbABC, FORCE_STATUSES, RETRIES, TIMEOUT, mnemonic_data_fname

__all__ = ['EngdbMast']

# Default MAST info.
MAST_BASE_URL = 'https://mast.stsci.edu'
API_URI = 'api/v0.1/Download/file'
SERVICE_URI = 'mast:jwstedb/'

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class EngdbMast(EngdbABC):
    """
    Access the JWST Engineering Database through MAST

    Parameters
    ----------
    base_url : str
        The base url for the engineering RESTful service. If not defined,
        the environmental variable ENG_BASE_URL is queried. Otherwise
        the default MAST website is used.

    token : str or None
        The MAST access token. If not defined, the environmental variable
        MAST_API_TOKEN is queried. A token is required.
        For more information, see 'https://auth.mast.stsci.edu/'

    service_kwargs : dict
        Service-specific keyword arguments that are not relevant to this implementation
        of EngdbABC.

    Raises
    ------
    RuntimeError
        Any and all failures with connecting with the MAST server.
    """

    #: The base URL for the engineering service.
    base_url = None

    #: The end time of the last query.
    endtime = None

    #: The results of the last query.
    response = None

    #: Number of retries to attempt to contact the service
    retries = RETRIES

    #: The start time of the last query.
    starttime = None

    #: Network timeout when communicating with the service
    timeout = TIMEOUT

    #: MAST Token
    token = None

    def __init__(self, base_url=None, token=None, **service_kwargs):
        logger.debug('kwargs not used by this service: %s', service_kwargs)

        self.configure(base_url=base_url, token=token)

        # Check for basic aliveness.
        try:
            resp = requests.get(self.base_url + 'api/', timeout=self.timeout)
        except requests.exceptions.ConnectionError as exception:
            raise RuntimeError(f'MAST url: {self.base_url} is unreachable.') from exception
        if resp.status_code != 200:
            raise RuntimeError(f'MAST url: {self.base_url} is not available. Returned HTTPS status {resp.status_code}')

        # Basics are covered. Finalize initialization.
        self.set_session()

    def cache(self, mnemonics, starttime, endtime, cache_path):
        """Cache results for the list of mnemonics

        Parameters
        ----------
        mnemonics : iterable
            List of mnemonics to retrieve

        starttime : str or astropy.time.Time
            The, inclusive, start time to retrieve from.

        endtime : str or astropy.time.Time
            The, inclusive, end time to retrieve from.

        cache_path : str or Path-like
            Path of the cache directory.
        """
        cache_path = Path(cache_path)
        cache_path.mkdir(parents=True, exist_ok=True)

        for mnemonic in mnemonics:
            records = self._get_records(mnemonic, starttime, endtime)
            records.write(cache_path / f'{mnemonic}.ecsv', format='ascii.ecsv')

    def cache_as_local(self, mnemonics, starttime, endtime, cache_path):
        """Cache results for the list of mnemonics, but in the EngdbLocal format

        The target format is native to what the EngdbDirect service provides.

        Parameters
        ----------
        mnemonics : iterable
            List of mnemonics to retrieve

        starttime : str or astropy.time.Time
            The, inclusive, start time to retrieve from.

        endtime : str or astropy.time.Time
            The, inclusive, end time to retrieve from.

        cache_path : str or Path-like
            Path of the cache directory.
        """
        cache_path = Path(cache_path)
        cache_path.mkdir(parents=True, exist_ok=True)

        # Get mnemonic data.
        for mnemonic in mnemonics:
            records = self._get_records(mnemonic, starttime, endtime)

            target = dict()
            target['TlmMnemonic'] = mnemonic.upper()
            target['AllPoints'] = 1
            target['Count'] = len(records)
            target['Data'] = list()
            for record in records:
                t = Time(record['MJD'], format='mjd')
                t = int(t.unix * 1000.)
                v = record['euvalue']
                if record['sqldataType'] in ['int', 'tinyint']:
                    v = int(v)
                entry = {'ObsTime': f'/Date({t}+0000)/', 'EUValue': v}
                target['Data'].append(entry)

            with open(cache_path / mnemonic_data_fname(mnemonic), 'w') as fp:
                json.dump(target, fp)

        # When used with EngDB_Mocker, a `meta.json` needs to exist.
        # A default is saved within the package.
        metas_path = Path(__file__).parent / 'tests/data/meta_for_mock.json'
        copy2(metas_path, cache_path / 'meta.json')

    def configure(self, base_url=None, token=None):
        """Configure from parameters and environment

        Parameters
        ----------
        base_url : str
            The base url for the engineering RESTful service. If not defined,
            the environmental variable ENG_BASE_URL is queried. Otherwise
            the default MAST website is used.

        token : str or None
            The MAST access token. If not defined, the environmental variable
            MAST_API_TOKEN is queried. A token is required.
            For more information, see 'https://auth.mast.stsci.edu/'
        """
        # Determine the database to use
        if base_url is None:
            base_url = getenv('ENG_BASE_URL', MAST_BASE_URL)
        if base_url[-1] != '/':
            base_url += '/'
        self.base_url = base_url

        # Get the token
        if token is None:
            token = getenv('MAST_API_TOKEN', None)
        if token is None:
            raise RuntimeError('No MAST token provided but is required. See https://auth.mast.stsci.edu/ for more information.')
        self.token = token

        # Get various timeout parameters
        self.retries = getenv('ENG_RETRIES', RETRIES)
        self.timeout = getenv('ENG_TIMEOUT', TIMEOUT)

    def get_meta(self, *kwargs):
        """Get the mnemonics meta info

        The MAST interface does not provide any meta.
        """
        raise NotImplementedError('MAST Engineering AUI does not provide a meta service')

    def get_values(self, mnemonic, starttime, endtime,
                   time_format=None, include_obstime=False, include_bracket_values=False, zip_results=True):
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
            Returns the list of values. See `include_obstime` and `zip` for modifications.
        """
        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)

        records = self._get_records(mnemonic=mnemonic, starttime=starttime,
                                    endtime=endtime, time_format=time_format)

        # If desired, remove bracket or outside of timeframe entries.
        if not include_bracket_values:
            selection = np.logical_and(records['MJD'] >= starttime.mjd,
                                       records['MJD'] <= endtime.mjd)
            records = records[selection]

        # Reformat to the desired list formatting.
        results = _ValueCollection(
            include_obstime=include_obstime,
            zip_results=zip_results
        )
        values = records['euvalue']
        obstimes = Time(records['MJD'], format='mjd')
        for obstime, value in zip(obstimes, values):
            results.append(obstime, value)

        return results.collection

    def set_session(self):
        """Setup HTTP session"""
        self._req = requests.Request(method='GET',
                                     url=self.base_url + API_URI,
                                     headers={'Authorization': f'token {self.token}'})

        s = requests.Session()
        retries = Retry(total=self.retries, backoff_factor=1.0, status_forcelist=FORCE_STATUSES, raise_on_status=True)
        s.mount('https://', HTTPAdapter(max_retries=retries))

        self._session = s

    def _get_records(
            self,
            mnemonic,
            starttime,
            endtime,
            time_format=None,
            **other_kwargs
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

        other_kwargs : dict
            Keyword arguments not relevant to this implementation.

        Returns
        -------
        records : `astropy.Table`
            Returns the resulting table.

        Notes
        -----
        The engineering service always returns the bracketing entries
        before and after the requested time range.
        """
        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)
        self.starttime = starttime
        self.endtime = endtime

        # Make the request
        mnemonic = mnemonic.strip()
        mnemonic = mnemonic.upper()
        starttime_fmt = starttime.strftime('%Y%m%dT%H%M%S')
        endtime_fmt = endtime.strftime('%Y%m%dT%H%M%S')
        uri = f'{mnemonic}-{starttime_fmt}-{endtime_fmt}.csv'
        self._req.params = {'uri': SERVICE_URI + uri}
        prepped = self._session.prepare_request(self._req)
        settings = self._session.merge_environment_settings(prepped.url, {}, None, None, None)
        logger.debug('Query: %s', prepped.url)
        self.response = self._session.send(prepped, timeout=self.timeout, **settings)
        self.response.raise_for_status()
        logger.debug('Response: %s', self.response)
        logger.debug('Response test: %s', self.response.text)

        # Convert to table.
        r_list = self.response.text.split('\r\n')
        table = Table.read(r_list, format='ascii.csv')

        return table


class _ValueCollection:
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
        if zip_results:
            self.collection = []
        else:
            self.collection = EngDB_Value([], [])

    def append(self, obstime, value):
        """Append value to collection

        Parameters
        ----------
        obstime : `astropy.time.Time`
            Observation time as returned from the engineering

        value : numeric
            Value from db.
        """
        # Make all the times readable
        obstime.format = 'isot'

        # Append
        if self._include_obstime:
            if self._zip_results:
                self.collection.append(
                    EngDB_Value(obstime, value)
                )
            else:
                self.collection.obstime.append(obstime)
                self.collection.value.append(value)
        else:
            self.collection.append(value)
