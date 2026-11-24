"""Access the JWST Engineering Mnemonic Database through direct connection."""

import logging
from os import getenv

import requests
from astropy.time import Time
from requests.adapters import HTTPAdapter, Retry

from jwst.lib.engdb_lib import FORCE_STATUSES, RETRIES, TIMEOUT, EngDB_Value, EngdbABC

# Configure logging
logger = logging.getLogger(__name__)

# #############################################
# Where is the engineering service? Its HERE!!!
# #############################################


# The big difference between old and new query service is:
#  old: ObsTime in the JSON-formatted response string was milliseconds from
#      UNIX epoch, e.g., Date(1453041183719+0000)
#  new: ObsTime in the JSON-formatted response string is an ISO-formatted string:
#      YYYY-MM-DDThh:mm:ss.nn
# The response is always a string in JSON format.
#  I believe that was always the default format, with nothing specified in the call.
#  The only other option was 'xml', which instead returned the results in XML format.
#  I do not believe that XML response is still supported, but it looks like it was not used anyway.

# Engineering DB Query Service Call (currently using the Data Archive API (DAAPI))
# base URL ENV VAR
#  ENGDB_DAAPI_QUERY = "https://{server}/api/query/jwst/edp/mnemonics/"
# EngDb access type: fqa=frequent access OR spa=sporadic access (CAL and SDP query fqa)
# Use an env var for this? --> ENGDB_ACCESS
# A full EngB Query Service Call endpoint for a mnemonic retrieval is:
#  ENGDB_DAAPI_QUERY + ENGDB_ACCESS/ + call,
#  where call is:
#   data?mnemonic={mnemonic}&s_time={startTimeIso}&e_time={endTimeIso}
#   metadata?mnemonic={mnemonic}
#   metadata    # (retrieve metadata for all eng mnemonics)
ENGDB_DATA = "/data"
ENGDB_METADATA = "/metadata"


__all__ = ["EngdbDirect"]


class EngdbDirect(EngdbABC):
    """
    Access the JWST Engineering Database through direct connection.

    Retrieve mnemonic metadata from the JWST Engineering Database using the DAAPI Query Service.

    Parameters
    ----------
    base_url : str
        The base url for the engineering database query service.

    **service_kwargs : dict
        Service-specific keyword arguments that are not relevant to this implementation
        of EngdbABC.
    """

    #: The base URL for the engineering service.
    base_url = None

    #: The end time of the last query.
    endtime = None

    #: The results of the last query.
    response = None

    #: The start time of the last query.
    starttime = None

    def __init__(self, base_url=None, **service_kwargs):
        logger.debug("kwargs not used by this service: %s", service_kwargs)

        self.configure(base_url=base_url)

        self.set_session()

        response = self._session.get("".join([self.base_url, ENGDB_METADATA]), timeout=self.timeout)
        response.raise_for_status()

    def configure(self, base_url=None):
        """
        Configure from parameters and environment.

        Parameters
        ----------
        base_url : str
            The base url for the engineering database query service.
        """
        # Setup for a retrieval using the engineering database query service
        if base_url is None:
            base_url = getenv("ENGDB_DAAPI_QUERY")
        if not base_url:
            raise RuntimeError("No engineering database URL given.")
        if base_url[-1] != "/":
            base_url += "/"
        # Accommodate differing access methods
        # Default to fqa for operations
        access = getenv("ENGDB_ACCESS")
        if access is None:
            base_url = base_url + "fqa"
        else:
            base_url = base_url + access
        self.base_url = base_url + "/" if base_url[-1] != "/" else base_url

        # Get various timeout parameters
        self.retries = getenv("ENG_RETRIES", RETRIES)
        self.timeout = getenv("ENG_TIMEOUT", TIMEOUT)

    def get_meta(self, mnemonic=""):
        """
        Get the mnemonics meta info.

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve.

        Returns
        -------
        response_json : dict
            Meta info as JSON.
        """
        query = "".join([self.base_url, ENGDB_METADATA, "?mnemonic=", mnemonic])
        logger.debug('Query URL="%s"', query)

        # Make our request
        response = self._session.get(query, timeout=self.timeout)
        logger.debug('Response="%s"', response)
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
        zip_results=True,
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve.

        starttime : str or `astropy.time.Time`
            The, inclusive, start time to retrieve from.

        endtime : str or `astropy.time.Time`
            The, inclusive, end time to retrieve from.

        time_format : str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        include_obstime : bool
            If `True`, the return values will include observation
            time as `astropy.time.Time`. See ``zip_results`` for further details.

        include_bracket_values : bool
            The DB service, by default, returns the bracketing
            values outside of the requested time. If `True`, include
            these values.

        zip_results : bool
            If `True` and ``include_obstime`` is `True`, the return values
            will be a list of 2-tuples. If false, the return will
            be a single 2-tuple, where each element is a list.

        Returns
        -------
        values : [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
            Returns the list of values. See ``include_obstime``
            and ``zip_results`` for modifications.

        Raises
        ------
        requests.exceptions.HTTPError
            Either a bad URL or non-existent mnemonic.
        """
        records = self._get_records(
            mnemonic=mnemonic, starttime=starttime, endtime=endtime, time_format=time_format
        )

        # The retrieved rows include:
        #   all rows where the observed time is within starttime to endtime, inclusive
        #   AND the latest row prior to the starttime
        #   AND the earliest row beyond the endtime
        # All times are retrieved as ISOT-formatted strings, converted to astropy.time.
        # Records returned are apparent not strictly correlated with
        # observation time. So, need to filter further.
        db_starttime = Time(records["ReqSTime"], format="isot")
        db_endtime = Time(records["ReqETime"], format="isot")
        results = _ValueCollection(include_obstime=include_obstime, zip_results=zip_results)
        if records["Data"] is not None:
            for record in records["Data"]:
                obstime = Time(record["ObsTime"], format="isot")
                if not include_bracket_values:
                    if obstime < db_starttime or obstime > db_endtime:
                        continue
                value = record["EUValue"]
                results.append(obstime, value)

        return results.collection

    def set_session(self):
        """Set up HTTP session."""
        s = requests.Session()
        retries = Retry(
            total=10, backoff_factor=1.0, status_forcelist=FORCE_STATUSES, raise_on_status=True
        )
        s.mount("https://", HTTPAdapter(max_retries=retries))

        self._session = s

    def _get_records(self, mnemonic, starttime, endtime, time_format=None):
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
            Either a bad URL or non-existent mnemonic.

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

        query = "".join(
            [
                self.base_url,
                ENGDB_DATA,
                "?mnemonic=",
                mnemonic,
                "&s_time=",
                starttime.iso,
                "&e_time=",
                endtime.iso,
            ]
        )
        logger.debug('Query URL="%s"', query)

        # Make our request
        response = self._session.get(query, timeout=self.timeout)
        logger.debug("Response: %s", response)
        logger.debug("Response: %s", response.json())
        response.raise_for_status()

        # That's all folks
        self.response = response
        return response.json()


class _ValueCollection:
    """
    Engineering Value Collection.

    Parameters
    ----------
    include_obstime : bool
        If `True`, the return values include observed time
        as `astropy.time.Time`. See `zip_results` for further details.

    zip_results : bool
        If `True` and `include_obstime` is `True`, a single list of 2-tuples is returned.
        If false, a single 2-tuple is returned, where each element is a list.

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
        """
        Define a collection to hold Data from the query response.

        Returns
        -------
        list
            Collection of values from EngDB.
        """
        if self._include_obstime:
            if self._zip_results:
                collection = [
                    EngDB_Value(t, v) for t, v in zip(self.obstimes, self.values, strict=False)
                ]
            else:
                collection = EngDB_Value(self.obstimes, self.values)
        else:
            collection = self.values
        return collection

    def append(self, obstime, value):
        """
        Append an astropy.time observed time and EU value to collection.

        Parameters
        ----------
        obstime : astropy.time
            Observed time of the EU value

        value : numeric
            Mnemonic EU value at the observed time (may be string or numeric)
        """
        self.obstimes.append(obstime)
        self.values.append(value)
