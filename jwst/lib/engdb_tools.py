"""Access the JWST Engineering Mnemonic Database

The engineering mnemonics are provided by multiple services, all of which require a level of authentication.

For non-operational use, the providing service is through the MAST AUI website

    https://mast.stsci.edu

Authorization can be requested through

    https://auth.mast.stsci.edu/

Interface
=========

The primary entry point is the function `jwst.lib.engdb_tools.ENGDB_Service`.
This function returns a `jwst.lib.engdb_lib.EngdbABC` connection object. Using
this object, values for a mnemonic covering a specified time range can be
retrieved using the `~jwst.lib.engdb_lib.EngdbABC.get_values` method.

By default, only values inclusively between the time end points are returned.
Depending on the frequency a mnemonic is updated, there can be no values. If
values are always desired, the nearest, bracketing values outside the time
range can be requested.

Warning
=======

Many mnemonics are updated very quickly, up to 16Hz. When in doubt, specify a
very short time frame, and request bracketing values. Otherwise, the request
can return a very large amount of data, risking timeout, unnecessary memory
consumption, or access restrictions.

Examples
========

The typical workflow is as follows:

.. code-block:: python

    from jwst.lib.engdb_tools import ENGDB_Service

    service = ENGDB_Service()  # By default, will use the public MAST service.

    values = service.get_values('sa_zattest2', '2021-05-22T00:00:00', '2021-05-22T00:00:01')

Environmental Variables
=======================

ENG_BASE_URL
    If no URL is specified in code or by command line parameters, this value is used.
    If not defined, a default, as defined by the individual services, will be attempted.

MAST_API_TOKEN
    If no token is provided in code or by command line parameters, this value will be used.
    `~jwst.lib.engdb_mast.EngdbMast` service requires a token to be provided.
    See https://auth.mast.stsci.edu/ for more information.

ENG_RETRIES
    Number of attempts to make when connecting to the service. Default is 10.

ENG_TIMEOUT
    Number of seconds before timing out a network connection. Default is 600 seconds (10 minutes)
"""
import logging

from .engdb_direct import EngdbDirect
from .engdb_mast import EngdbMast

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def ENGDB_Service(base_url=None, **service_kwargs):
    """Access the JWST Engineering Database

    Access can be either through the public MAST API
    or by direct connection to the database server.

    Parameters
    ----------
    base_url : str or None.
        The base url for the engineering RESTful service

    service_kwargs : dict
        Service-specific keyword arguments. Refer to the concrete implementations
        of `~jwst.lib.engdb_lib.EngdbABC`.

    Returns
    -------
    service : `~jwst.lib.engdb_lib.EngdbABC`
        The engineering database service to use.
    """

    # Determine the database to use
    for db_attempt in [EngdbMast, EngdbDirect]:
        try:
            service = db_attempt(base_url=base_url, **service_kwargs)
        except RuntimeError as excp:
            logger.debug('Service %s cannot use base_url %s.', db_attempt, base_url)
            logger.debug('Exception: %s', excp)
        else:
            # Found the working service. Continue on.
            break
    else:
        raise RuntimeError(f'Base URL of {base_url} cannot be accessed.')

    # Service is in hand.
    return service
