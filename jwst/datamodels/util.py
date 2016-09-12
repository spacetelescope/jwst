"""
Various utility functions and data types
"""
from __future__ import absolute_import, unicode_literals, division, print_function

import sys

import numpy as np

from astropy.extern import six

def can_broadcast(a, b):
    """
    Given two shapes, returns True if they are broadcastable.
    """
    for i in range(1, min(len(a), len(b)) + 1):
        adim = a[-i]
        bdim = b[-i]

        if not (adim == 1 or bdim == 1 or adim == bdim):
            return False

    return True


def to_camelcase(token):
    return ''.join(x.capitalize() for x in token.split('_-'))


if six.PY3:
    def fits_header_name(name):
        """
        Returns a FITS header name in the correct form for the current
        version of Python.
        """
        if isinstance(name, bytes):
            return name.decode('ascii')
        return name
else:
    def fits_header_name(name):
        """
        Returns a FITS header name in the correct form for the current
        version of Python.
        """
        if isinstance(name, unicode):
            return name.encode('ascii')
        return name


def gentle_asarray(a, dtype):
    """
    Performs an asarray that doesn't cause a copy if the byteorder is
    different.  It also ignores column name differences -- the
    resulting array will have the column names from the given dtype.
    """
    out_dtype = np.dtype(dtype)
    if isinstance(a, np.ndarray):
        in_dtype = a.dtype
        # Non-table array
        if in_dtype.fields is None and out_dtype.fields is None:
            if np.can_cast(in_dtype, out_dtype, 'equiv'):
                return a
            else:
                return np.asarray(a, dtype=out_dtype)
        elif in_dtype.fields is not None and out_dtype.fields is not None:
            if in_dtype == out_dtype:
                return a
            if len(in_dtype) != len(out_dtype):
                raise ValueError(
                    "Wrong number of columns.  Expected {0}, got {1}".format(
                        len(out_dtype), len(in_dtype)))
            new_dtype = []
            # We cannot change the names in the dtype record because
            # they are separately stored in the table _coldefs attribute
            if hasattr(in_dtype, 'names') and hasattr(out_dtype, 'names'):
                out_dtype.names = in_dtype.names
            for i in range(len(out_dtype.fields)):
                in_type = in_dtype[i]
                out_type = out_dtype[i]
                if in_type.subdtype is None:
                    type_str = in_type.str
                else:
                    type_str = in_type.subdtype[0].str
                if np.can_cast(in_type, out_type, 'equiv'):
                    new_dtype.append(
                        (out_dtype.names[i],
                         type_str,
                         in_type.shape))
                else:
                    return np.asarray(a, dtype=out_dtype)
            return a.view(dtype=np.dtype(new_dtype))
    else:
        try:
            a = np.asarray(a, dtype=out_dtype)
        except:
            raise ValueError("Can't convert {0!s} to ndarray".fornat(type(a)))
        return a

def get_short_doc(schema):
    title = schema.get('title', None)
    description = schema.get('description', None)
    if description is None:
        description = title or ''
    else:
        if title is not None:
            description = title + '\n\n' + description
    return description.partition('\n')[0]


def ensure_ascii(s):
    if isinstance(s, six.text_type):
        s = s.encode('ascii', 'replace')
        if six.PY3:
            s = s.decode('ascii')
    return s
