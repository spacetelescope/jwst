import numpy as np


class _MissingValueType:
    pass


_MISSING_VALUE = _MissingValueType()


def _convert_dtype(value):
    """
    Convert numarray column dtype into YAML-compatible format description.

    Parameters
    ----------
    value : str, bool, or numpy.dtype
        The datatype to convert.

    Returns
    -------
    new_dtype : str or list
        The converted datatype.
    """
    if "U" in value:
        # working with a string description
        str_len = int(value[value.find("U") + 1 :])
        new_dtype = ["ascii", str_len]
    elif value == "bool":
        # cast all bool to int8 to avoid issues on write
        new_dtype = "int8"
    else:
        new_dtype = str(value)

    return new_dtype


def table_to_schema(table):
    """
    Construct a schema for a table.

    Convert a "table" (a structured ndarray) to a stdatamodels
    sub-schema that will allow the "table" to be stored to a fits
    extension HDRTAB.

    Parameters
    ----------
    table : ndarray
        The structured array containing the data (and datatype).

    Returns
    -------
    subschema : dict
        Stdatamodels for the "table" datatype.
    """
    return {
        "title": "Combined header table",
        "fits_hdu": "HDRTAB",
        "datatype": [
            {
                "name": col_name,
                "datatype": _convert_dtype(str(table.dtype[col_name])),
            }
            for col_name in table.dtype.fields
        ],
    }


class TableBuilder:
    """
    Class to build a metadata table.

    Used to incrementally build a metadata "table" (a numpy
    structured array) containing metadata from several models.

    >>> tb = TableBuilder({"meta.filename": "FN"})
    >>> tb.header_to_row({"meta.filename": "foo.fits"})
    >>> tb.build_table()
    rec.array([('foo.fits',)],
         dtype=[('FN', '<U8')])
    """

    def __init__(self, attr_to_column):
        """
        Create a new `TableBuilder`.

        Parameters
        ----------
        attr_to_column : dict
            A one-to-one mapping of attribute names (as
            dotted paths like "meta.filename") to column
            names (strings).
        """
        self.attr_to_column = attr_to_column
        self.columns = {col: [] for col in self.attr_to_column.values()}
        if len(attr_to_column) != len(self.columns):
            raise ValueError(f"Invalid attr_to_column, mapping is not 1-to-1: {attr_to_column}")

    def header_to_row(self, header):
        """
        Add metadata in header to the table.

        This function will add a complete row for each
        header. If header is missing a required attribute
        (as defined in the mapping provided to
        `TableBuilder.__init__`) the column with this
        missing value for this row will contain a ``nan``.

        Parameters
        ----------
        header : dict
            Often produced from ``Datamodel.to_flat_dict``.
        """
        row = {}
        for attr in self.attr_to_column:
            if attr in header:
                row[attr] = header[attr]
        self._add_row(row)

    def _add_row(self, row):
        for attr, col in self.attr_to_column.items():
            self.columns[col].append(row[attr] if attr in row else _MISSING_VALUE)

    def build_table(self):
        """
        Build the final "table".

        Returns
        -------
        table : numpy.ndarray
            Structured array containing fields with datatypes
        """
        arrays = []
        table_dtype = []
        for col, items in self.columns.items():
            if all(i is _MISSING_VALUE for i in items):
                continue
            arrays.append(np.array([np.nan if i is _MISSING_VALUE else i for i in items]))
            table_dtype.append((col, arrays[-1].dtype))
        return np.rec.fromarrays(arrays, dtype=table_dtype)
