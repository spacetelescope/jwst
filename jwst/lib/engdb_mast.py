"""
Access the JWST Engineering Mnemonic Database through MAST
"""
import logging
from os import getenv
import requests

# Default MAST base url
MAST_BASE_URL = 'https://mast.stsci.edu'

# Aliveness query
ALIVE_QUERY = 'SA_ZATTEST2-20210522T000000-20210522T000001.csv'

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class EngdbMast():
    """
    Access the JWST Engineering Database through MAST

    Parameters
    ----------
    base_url : str
        The base url for the engineering RESTful service. If not defined,
        the environmental variable ENG_BASE_URL is queried. Otherwise
        the default MAST website is used.

    token : str
        The MAST access token. If not defined, the environmental variable
        MAST_API_TOKEN is queried. A token is required.
        For more information, see 'https://auth.mast.stsci.edu/'

    db_kwargs : dict
        Other keyword arguments that may be needed by other versions of the service
        but are not used by this service.

    Attributes
    ----------
    response: `requests.response`
        The results of the last query.

    starttime: `astropy.time.Time`
        The start time of the last query.

    endtime: `astropy.time.Time`

    base_url: str
        The base URL for the engineering service.
    """
    def __init__(self, base_url=None, token=None, **db_kwargs):

        # Determine the database to use
        if base_url is None:
            base_url = getenv('ENG_BASE_URL', MAST_BASE_URL)
        if base_url[-1] != '/':
            base_url += '/'

        # Get the token
        if token is None:
            token = getenv('MAST_API_TOKEN', None)
        if token is None:
            raise RuntimeError('No MAST token provided but is required. See https://auth.mast.stsci.edu/ for more information.')

        # Check for basic aliveness.
        try:
            resp = requests.get(base_url + 'api/')
        except requests.exceptions.ConnectionError as exception:
            raise RuntimeError(f'MAST url: {base_url} is unreachable.') from exception
        if resp.status_code != 200:
            raise RuntimeError(f'MAST url: {base_url} is not available. Returned HTTPS status {resp.status_code}')
