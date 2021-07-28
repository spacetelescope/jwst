"""SIAF Database Access

Provide a common interface to different versions of the SIAF.

Under operations, the SIAF is found in a sqlite database.
Otherwise, use the standard interface defined by the `pysiaf` package
"""
import logging
import os
from pathlib import Path

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class SiafDb:
    """SIAF Database Access

    Provide a common interface to different versions of the SIAF.

    Under operations, the SIAF is found in a sqlite database.
    Otherwise, use the standard interface defined by the `pysiaf` package

    Parameters
    ----------
    source : None, str, or a file-like object
        The SIAF database source. See notes for more details.

    useafter : None or str
        Only use entries on or after this date. Date should be in ISO format.
        If None, the useafter date is set to '1900-01-01'.

    Notes
    -----
    The interpretation of `source` is as follows:

    If None, the environmental 'XML_DATA' is queried for a value.
    If None, then the `pysiaf` package is used.
    If a string, the string is treated as a path.
    If that path is to a folder, the `pysiaf` package is used with the folder
        as the XML source folder. See the `pysiaf` package for more information.
    Finally, an attempt is made to open the path as a sqlite database.
    Otherwise, fail.
    """
    def __init__(self, source=None, useafter='1900-01-01'):

        # If no source, retrieve the environmental XML_DATA
        if source is None:
            source = os.environ.get('XML_DATA', None)
            if source is not None:
                source = Path(source) / 'prd.db'

        # Attempt to access source as a pysiaf source
        try:
            source = SiafDbPySiaf(source, useafter)
        except ValueError:
            # Source is incompatible.
            logger.debug('Could not open as a pysiaf object: %s', source)
        else:
            self._source = source
            return

        # Attempt to access source as an sqlite database.
        self._source = SiafDbSqlite(source, useafter)


class SiafDbPySiaf:
    """Use pysiaf as the source of siaf information

    Parameters
    ----------
    source : None, str, or a file-like object
        The SIAF database source. See notes for more details.

    useafter : None or str
        Only use entries on or after this date. Date should be in ISO format.
        If None, the useafter date is set to '1900-01-01'.

    Notes
    -----
    The interpretation of `source` is as follows:

    If None, then the `pysiaf` package is used.
    If a string, the string is treated as a path.
    If that path is to a folder, the `pysiaf` package is used with the folder
        as the XML source folder. See the `pysiaf` package for more information.
    Otherwise, fail.
    """
    def __init__(self, source=None, useafter='1900-01-01'):
        import pysiaf  # noqa: Ensure that the pysiaf package is available.

        if source is not None:
            source = Path(source)
            if not source.is_dir():
                raise ValueError('Source %s: Needs to be a folder for use with pysiaf')
        self._source = source


class SiafDbSqlite:
    """Use a sqlite db as the source of siaf information

    Parameters
    ----------
    source : str, or a file-like object
        The SIAF database source.

    useafter : None or str
        Only use entries on or after this date. Date should be in ISO format.
        If None, the useafter date is set to '1900-01-01'.
    """
    def __init__(self, source, useafter='1900-01-01'):
        import sqlite3

        source = Path(source)
        if not source.exists():
            raise ValueError('Source: %s does not exist.', source)
        source = f'file:{str(source)}?mode=ro'

        self._source = sqlite3.connect(source, uri=True)
        self._cursor = self._source.cursor()

    def __del__(self):
        """Close out any connections"""
        if self._source:
            self._source.close()
