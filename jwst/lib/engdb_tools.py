"""Access the JWST Engineering Database

Access can be either through the public MAST API
or by direct connection to the database server.
"""
import logging
from os import getenv

from .engdb_direct import EngdbDirect
from .engdb_mast import EngdbMast

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# This is currently set to the OPS database.
ENGDB_HOST = 'http://pwjwdmsemwebag.stsci.edu/'
ENGDB_BASE_URL = ''.join([
    ENGDB_HOST,
    'JWDMSEngFqAcc/',
    'TlmMnemonicDataSrv.svc/',
])


class ENGDB_Service:
    """Access the JWST Engineering Database

    Access can be either through the public MAST API
    or by direct connection to the database server.

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

    base_url: str or None
        The base URL for the engineering service.
        If None, the MAST AUI will be tried first, then
        the direct service.

    default_format: str
        The format to retrieve from the service.
        This is not the format of the returned data.
    """
    def __init__(self, base_url=None, default_format='dict'):

        # Determine the database to use
        for db_attempt in [EngdbMast, EngdbDirect]:
            try:
                db = db_attempt(base_url=base_url, default_format=default_format)
            except RuntimeError as excp:
                logger.debug('Service %s cannot use base_url %s.', db_attempt, base_url)
                logger.debug('Exception: %s', excp)
            else:
                # Found the working service. Continue on.
                break
        else:
            raise RuntimeError('Base URL of %s cannot be accessed.')

        # Service is in hand, continue initialization.
        self.db = db

    def __getattr__(self, name):
        """Attempt access through the underlying database object

        Parameters
        ----------
        name : str
            The attribute to access.
        """
        return getattr(self.db, name)
