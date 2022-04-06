"""
Association Pools
"""

from pkg_resources import parse_version, get_distribution

from collections import UserDict

from astropy.io.ascii import convert_numpy
from astropy.table import Row, Table

__all__ = ['AssociationPool']

DEFAULT_DELIMITER = '|'
DEFAULT_FORMAT = 'ascii'


class AssociationPool(Table):
    """Association Pool

    An ``AssociationPool`` is essentially an astropy Table with the
    following default behaviors:

    - ASCII tables with a default delimiter of `|`
    - All values are read in as strings
    """
    def __init__(self, *args, **kwargs):
        super(AssociationPool, self).__init__(*args, **kwargs)
        self.meta['pool_file'] = self.meta.get('pool_file', 'in-memory')

    @classmethod
    def read(
            cls,
            filename,
            delimiter=DEFAULT_DELIMITER,
            format=DEFAULT_FORMAT,
            **kwargs
    ):
        """Read in a Pool file

        Parameters
        ----------
        filename : str
            File path to read in as a table.

        delimiter : str
            Character used to delineate columns.

        format : str
            The format of the input file.

        Returns
        -------
        AssociationPool
            The ``AssociationPool`` representation of the file.
        """
        table = super(AssociationPool, cls).read(
            filename, delimiter=delimiter,
            format=format,
            converters=convert_to_str, **kwargs
        )

        # If anything has been masked, just fill
        table = table.filled('null')

        # Lowercase the column names
        # Note: Cannot do in-place because modifying the
        #       list while iterating.
        columns = [column for name, column in table.columns.items()]
        for c in columns:
            c.name = c.name.lower()

        table.meta['pool_file'] = filename
        return table

    def write(self, *args, **kwargs):
        """Write the pool to a file.

        Parameters
        ----------
        output : str, file-like
            The output file or file-like object.

        delimiter : str
            The string to use to delineate columns.
            Default is '|'.

        format : str
            The format the file should be written in.
            Default is 'ascii'.

        args, kwargs : obj
            Other parameters that ``astropy.io.ascii.write`` can accept.
        """
        delimiter = kwargs.pop('delimiter', DEFAULT_DELIMITER)
        format = kwargs.pop('format', DEFAULT_FORMAT)
        try:
            super(AssociationPool, self).write(
                *args, delimiter=delimiter, format=format, **kwargs
            )
        except TypeError:
            # Most likely caused by the actual `write` called
            # does not handle `delimiter`. `jsviewer` is one
            # such format.
            # So, try again without a delimiter.
            super(AssociationPool, self).write(
                *args, format=format, **kwargs
            )

class PoolRow(UserDict):
    """A row from an AssociationPool

    Class to create a copy of an AssociationPool row without copying
    all of the astropy.Table.Row private attributes.
    """
    def __init__(self, init=None):
        dict_init = dict(init)
        super().__init__(dict_init)
        try:
            self.meta = init.meta
        except AttributeError:
            self.meta = dict()


def _convert_to_str():
    func, type_ = convert_numpy(str)

    def convert_func(vals):
        """Lowercase the conversion"""
        results = func(vals)
        results = [result.lower() for result in results]
        return results

    return [(convert_func, type_)]


class _ConvertToStr(dict):
    def __getitem__(self, k):
        return _convert_to_str()

    def get(self, k, default=None):
        return self.__getitem__(k)


if parse_version(get_distribution('astropy').version) >= parse_version('5.0.dev'):
    convert_to_str = {'*': _convert_to_str()}
else:
    convert_to_str = _ConvertToStr()
