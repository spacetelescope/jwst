"""Module to handle association pools."""

from collections import UserDict

from astropy.io.ascii import convert_numpy
from astropy.table import Table

__all__ = ["AssociationPool", "PoolRow"]

DEFAULT_DELIMITER = "|"
DEFAULT_FORMAT = "ascii"


class AssociationPool(Table):
    """
    Association pool that is built off an astropy table.

    This class is essentially a `~astropy.table.Table` with the
    following modified behaviors:

    - ASCII tables with a default delimiter of ``|``.
    - All values are read in as strings.
    """

    def __init__(self, *args, **kwargs):
        super(AssociationPool, self).__init__(*args, **kwargs)
        self.meta["pool_file"] = self.meta.get("pool_file", "in-memory")

    @classmethod
    def read(cls, filename, delimiter=DEFAULT_DELIMITER, fmt=DEFAULT_FORMAT, **kwargs):
        """
        Read in an association pool file.

        Parameters
        ----------
        filename : str
            File path to read in as a table.
        delimiter : str
            Character used to delineate columns.
        fmt : str
            The format of the input file.

        Returns
        -------
        `AssociationPool`
            The association pool representation of the file.
        """
        table = super(AssociationPool, cls).read(
            filename, delimiter=delimiter, format=fmt, converters=convert_to_str, **kwargs
        )

        # If anything has been masked, just fill
        table = table.filled("null")

        # Lowercase the column names
        # Note: Cannot do in-place because modifying the
        #       list while iterating.
        columns = [column for name, column in table.columns.items()]
        for c in columns:
            c.name = c.name.lower()

        table.meta["pool_file"] = filename
        return table

    def write(self, *args, **kwargs):
        """
        Write the association pool to a file.

        Parameters
        ----------
        *args, **kwargs
            Arguments and keywords that :func:`astropy.io.ascii.write` can accept.
        """
        delimiter = kwargs.pop("delimiter", DEFAULT_DELIMITER)
        fmt = kwargs.pop("fmt", None)
        if fmt is None:
            fmt = kwargs.pop("format", DEFAULT_FORMAT)
        try:
            super(AssociationPool, self).write(*args, delimiter=delimiter, format=fmt, **kwargs)
        except TypeError:
            # Most likely caused by the actual `write` called
            # does not handle `delimiter`. `jsviewer` is one
            # such format.
            # So, try again without a delimiter.
            super(AssociationPool, self).write(*args, format=fmt, **kwargs)


class PoolRow(UserDict):
    """
    A row from an association pool.

    Class to create a copy of an `AssociationPool` row without copying
    all of the `astropy.table.Row` private attributes.
    """

    def __init__(self, init=None):
        dict_init = dict(init)
        super().__init__(dict_init)
        try:
            self.meta = init.meta
        except AttributeError:
            self.meta = {}


def _convert_to_str():
    func, type_ = convert_numpy(str)

    def convert_func(vals):
        """
        Lowercase the conversion.

        Parameters
        ----------
        vals : list of str
            The list of strings to lowercase.

        Returns
        -------
        list of str
            The lowercase strings.
        """
        results = func(vals)
        results = [result.lower() for result in results]
        return results

    return [(convert_func, type_)]


convert_to_str = {"*": _convert_to_str()}
