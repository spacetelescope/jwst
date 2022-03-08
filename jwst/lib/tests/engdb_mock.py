"""
Tools to mock the JWST Engineering Database
"""
from copy import copy
from functools import lru_cache
import json
import os
import re
import requests_mock

from astropy.time import Time

from jwst.lib import engdb_direct
from jwst.lib.engdb_lib import mnemonic_data_fname

__all__ = [
    'ENGDB_PATH',
    'cache_engdb'
]

# Setup for local testing cache
ENGDB_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'data',
    'engdb'
)
MNEMONICS_TO_CACHE = [
    'SA_ZATTEST1',
    'SA_ZATTEST2',
    'SA_ZATTEST3',
    'SA_ZATTEST4',
    'SA_ZRFGS2J11',
    'SA_ZRFGS2J21',
    'SA_ZRFGS2J31',
    'SA_ZRFGS2J12',
    'SA_ZRFGS2J22',
    'SA_ZRFGS2J32',
    'SA_ZRFGS2J13',
    'SA_ZRFGS2J23',
    'SA_ZRFGS2J33',
    'SA_ZADUCMDX',
    'SA_ZADUCMDY',
    'SA_ZFGGSCMDX',
    'SA_ZFGGSCMDY',
    'SA_ZFGGSPOSX',
    'SA_ZFGGSPOSY',
    'SA_ZFGDETID',
    'INRSI_GWA_Y_TILT_AVGED',
]

# Path templates
META = 'meta.json'
DATA = '_data.json'


class EngDB_Mocker(requests_mock.Mocker):
    """
    Setup a mock of the JWST Engineering DB.
    """

    def __init__(self, *args, real_http=True, **kwargs):
        db_path = kwargs.pop('db_path', ENGDB_PATH)
        super(EngDB_Mocker, self).__init__(*args, real_http=real_http, **kwargs)

        # Setup the local engineering cache
        self.cache = EngDB_Local(db_path)

        # Setup from meta query
        meta_query = re.compile(''.join([
            engdb_direct.ENGDB_METADATA,
            '.*'
        ]))
        self.get(meta_query, json=self.response_meta)

        # Setup to return a general data query
        data_query = re.compile(''.join([
            engdb_direct.ENGDB_DATA,
            r'.+\?sTime=.+\&eTime=.+'
        ]))
        self.get(data_query, json=self.response_data)

    def response_meta(self, request, context):
        """
        Respond with the meta data

        Parameters
        ----------
        request, context:
            request-mock parameters

        Returns
        -------
        response: json
            The expected JSON response

        Modifies
        --------
        context
            `Update status_code` with the HTTP expected status
        """

        # Get the mnemonic substring.
        mnemonic = os.path.basename(request.path)
        results = self.cache.fetch_meta(mnemonic)
        context.status_code = 200
        return results

    def response_data(self, request, context):
        """
        Respond with the mnemonic data

        Parameters
        ----------
        request, context:
            request-mock parameters

        Returns
        -------
        response: json
            The expected JSON response

        Modifies
        --------
        context
            `Update status_code` with the HTTP expected status
        """
        mnemonic = os.path.basename(request.path)
        data = self.cache.fetch_data(
            mnemonic,
            request.qs['stime'][0],
            request.qs['etime'][0]
        )
        context.status_code = 200
        return data

    def __enter__(self):
        """Setup environment for the context

        Remove MAST_API_TOKEN to ensure the EngdbMast service is not used.
        """
        self._environ = os.environ.copy()
        try:
            del os.environ['MAST_API_TOKEN']
        except KeyError:
            pass
        super().__enter__()

    def __exit__(self, type, value, traceback):
        """Restore the environment"""
        super().__exit__(type, value, traceback)
        os.environ = self._environ


class EngDB_Local():
    """
    Fetch engineering data from the local cache
    """

    def __init__(self, db_path=''):
        self.db_path = db_path
        self.json_load.cache_clear()

    def fetch_data(self, mnemonic, starttime, endtime):
        """
        Get data for a mnemonic.

        Parameters
        mnemonic:z str
            The engineering mnemonic to retrieve

        starttime: str or astropy.time.Time
            The, inclusive, start time to retrieve from.

        endtime: str or astropy.time.Time
            The, inclusive, end time to retrieve from.

        Returns
        -------
        mnemonic_data: dict
            The returned structure
        """
        db_data = self.json_load(mnemonic)

        # Find the range of data
        stime = Time(starttime, format='iso')
        etime = Time(endtime, format='iso')
        stime_mil = int(stime.unix * 1000)
        etime_mil = int(etime.unix * 1000)

        data = db_data['Data']
        start_idx = 0
        try:
            while engdb_direct.extract_db_time(data[start_idx]['ObsTime']) < stime_mil:
                start_idx += 1
        except IndexError:
            pass
        end_idx = start_idx
        try:
            while engdb_direct.extract_db_time(data[end_idx]['ObsTime']) <= etime_mil:
                end_idx += 1
        except IndexError:
            end_idx -= 1

        # Real database always returns a bracketed result.
        start_idx = max(0, start_idx - 1)
        end_idx = min(len(data), end_idx + 1)

        # Construct the resultant data
        data_return = [
            data[idx]
            for idx in range(start_idx, end_idx)
        ]

        # Construct the return structure
        result = copy(db_data)
        result['ReqSTime'] = '/Date({:013d}+0000)/'.format(stime_mil)
        result['ReqETime'] = '/Date({:013d}+0000)/'.format(etime_mil)
        result['Count'] = len(data_return)
        result['Data'] = data_return

        return result

    def fetch_meta(self, mnemonic_substr=''):
        """
        Get the meta for the match to the mnemonic

        Parameters
        ----------
        mnemonic_substr: str
            The substring to match in all available mnemonics

        Returns
        -------
        The Meta structure which matches the live engineering
        RESTful interface
        """
        with open(
                os.path.join(
                    self.db_path,
                    META
                ),
                'r'
        ) as fp:
            meta = json.load(fp)

        tlmmnemonics = meta['TlmMnemonics']

        # Now look for only the requested mnemonic substring
        mnemonic_substr = mnemonic_substr.lower()
        mnemonics = [
            tlmmnemonics[idx]
            for idx, mnemonic in enumerate(tlmmnemonics)
            if mnemonic_substr in mnemonic['TlmMnemonic'].lower()
        ]

        # Construct the return structure
        result = copy(meta)
        result['Count'] = len(mnemonics)
        result['TlmMnemonics'] = mnemonics

        return result

    @lru_cache(maxsize=128)
    def json_load(self, mnemonic):
        with open(
                os.path.join(
                    self.db_path,
                    mnemonic_data_fname(mnemonic)),
                'r') as fp:
            db_data = json.load(fp)
        return db_data


# #########
# Utilities
# #########
def cache_engdb(
        mnemonics=MNEMONICS_TO_CACHE,
        starttime='2016-01-18T15:40:00',
        endtime='2016-01-18T15:40:05',
        db_path=ENGDB_PATH
):
    """Create a local cache of the database

    Parameters
    ----------
    mnemonics: iterable
        List of mnemonics to retrieve

    starttime: str or astropy.time.Time
        The, inclusive, start time to retrieve from.

    endtime: str or astropy.time.Time
        The, inclusive, end time to retrieve from.

    path: str
        Path of the cache directory.

    Notes
    -----
    Suggested shell command is
    `python -c 'from jwst.lib.tests.engdb_mock import cache_engdb; cache_engdb()'`
    """
    if not os.path.exists(db_path):
        os.mkdir(db_path)

    edb = engdb_direct.EngdbDirect()

    # Get the meta info for all mnemonics regardless.
    meta = edb.get_meta()
    with open(os.path.join(db_path, META), 'w') as fp:
        json.dump(meta, fp)

    for mnemonic in mnemonics:
        records = edb._get_records(mnemonic, starttime, endtime)

        # Remove the request times. These are filled back in
        # during retrieval.
        del records['ReqSTime']
        del records['ReqETime']

        with open(
                os.path.join(db_path, mnemonic_data_fname(mnemonic)),
                'w'
        ) as fp:
            json.dump(records, fp)
