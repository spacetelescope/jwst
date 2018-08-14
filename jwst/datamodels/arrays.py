"""
This module contains functions for supporting Numpy arrays and
structured arrays.
"""

# TODO: Add __all__ members to modules

import re

import numpy as np

def _schema_dtype_to_numpy_dtype_single(dtype):
    names = {
        'bool8': np.int8,
        'int8': np.int8,
        'int16': np.int16,
        'int32': np.int32,
        'int64': np.int64,
        'uint8': np.uint8,
        'uint16': np.uint16,
        'uint32': np.uint32,
        'uint64': np.uint64,
        'float32': np.float32,
        'float64': np.float64,
        'complex64': np.complex64,
        'complex128': np.complex128,
        }

    if isinstance(dtype, str):
        if dtype in names:
            return names[dtype]
        elif re.match('string[0-9]+', dtype):
            size = int(dtype[6:])
            return 'S{0}'.format(size)
        else:
            raise ValueError("invalid dtype {0!r}".format(dtype))
    else:
        if dtype not in names.values():
            raise ValueError("invalid dtype {0!r}".format(dtype))
        return dtype


def _parse_column(i, entry):
    if isinstance(entry, str):
        col_dtype = _schema_dtype_to_numpy_dtype_single(entry)
        col_name = b'f{0}'.format(i)
        col_shape = 1
    elif isinstance(entry, dict):
        if 'dtype' in entry:
            col_dtype = _schema_dtype_to_numpy_dtype_single(
                entry['dtype'])
        else:
            col_dtype = np.float_

        if 'name' in entry:
            col_name = entry['name']
            if isinstance(col_name, bytes):
                col_name = col_name.encode('ascii')
            if not isinstance(col_name, str):
                raise ValueError(
                    "Invalid name {0} in dtype".format(col_name))
        else:
            col_name = b'f{0}'.format(i)

        if 'title' in entry:
            col_title = entry['title']
            if not isinstance(col_title, str):
                raise TypeError(
                    "Invalid title {0} in dtype".format(col_title))
            col_name = (col_name, col_title)

        if 'shape' in entry:
            col_shape = entry['shape']
            try:
                col_shape = int(col_shape)
            except TypeError:
                try:
                    col_shape = tuple(col_shape)
                except TypeError:
                    raise ValueError(
                        "Invalid shape {0} in dtype".format(col_shape))
                else:
                    for x in col_shape:
                        try:
                            x = int(x)
                        except TypeError:
                            raise ValueError(
                                "Invalid shape {0} in dtype".format(col_shape))
        else:
            col_shape = 1

    return (bytes(col_name), col_dtype, col_shape)


def schema_dtype_to_numpy_dtype(dtype):
    """
    Given a schema dtype, returns a numpy dtype.
    """
    if isinstance(dtype, str):
        return _schema_dtype_to_numpy_dtype_single(dtype)
    elif isinstance(dtype, list):
        new_dtype = []
        for i, entry in enumerate(dtype):
            new_dtype.append(
                _parse_column(i, entry))
        return new_dtype
