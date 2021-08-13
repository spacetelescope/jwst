"""Access the JWST Engineering Database

Access can be either through the public MAST API
or by direct connection to the database server.
"""
import logging

from .engdb_direct import EngdbDirect
from .engdb_mast import EngdbMast

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def ENGDB_Service(base_url=None, default_format='dict'):
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

    Returns
    -------
    service : EngdbABC
        The engineering database service to use.
    """

    # Determine the database to use
    for db_attempt in [EngdbMast, EngdbDirect]:
        try:
            service = db_attempt(base_url=base_url, default_format=default_format)
        except RuntimeError as excp:
            logger.debug('Service %s cannot use base_url %s.', db_attempt, base_url)
            logger.debug('Exception: %s', excp)
        else:
            # Found the working service. Continue on.
            break
    else:
        raise RuntimeError('Base URL of %s cannot be accessed.')

    # Service is in hand.
    return service
