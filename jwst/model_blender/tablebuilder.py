import numpy as np


def _convert_dtype(value):
    """Convert numarray column dtype into YAML-compatible format description"""
    if 'U' in value:
        # working with a string description
        str_len = int(value[value.find('U') + 1:])
        new_dtype = ['ascii', str_len]
    elif value == 'bool':
        # cast all bool to int8 to avoid issues on write
        new_dtype = 'int8'
    else:
        new_dtype = str(value)

    return new_dtype


def table_to_schema(table):
    return {
        'title': 'Combined header table',
        'fits_hdu': 'HDRTAB',
        'datatype': [
            {
                'name': col_name,
                'datatype': _convert_dtype(str(table.dtype[col_name])),
            }
            for col_name in table.dtype.fields
        ],
    }


class TableBuilder:
    def __init__(self, attr_to_column):
        self.attr_to_column = attr_to_column
        self.columns = {col: [] for col in self.attr_to_column.values()}
        self.masks = {col: [] for col in self.columns}
        if len(attr_to_column) != len(self.columns):
            raise ValueError(f"Invalid attr_to_column, mapping is not 1-to-1: {attr_to_column}")

    def header_to_row(self, header):
        row = {}
        for attr in self.attr_to_column:
            if attr in header:
                row[attr] = header[attr]
        self._add_row(row)

    def _add_row(self, row):
        for attr, col in self.attr_to_column.items():
            if attr in row:
                self.columns[col].append(row[attr])
                self.masks[col].append(False)
            else:
                self.columns[col].append(np.nan)
                self.masks[col].append(True)

    def build_table(self):
        arrays = []
        table_dtype = []
        for col, items in self.columns.items():
            arrays.append(np.ma.array(items, mask=self.masks[col]))
            table_dtype.append((col, arrays[-1].dtype))
        # TODO loses masks... but so did the old code?
        return np.rec.fromarrays(arrays, dtype=table_dtype)
