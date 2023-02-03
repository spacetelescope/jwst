# Copyright (C) 2010 Association of Universities for Research in Astronomy(AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

import numpy as np
from numpy import ma

from stdatamodels import DataModel
from stdatamodels.jwst import datamodels


class _KeywordMapping:
    """
    A simple class to verify and store information about each mapping
    entry.
    """

    def __init__(self, src_kwd, dst_name, agg_func=None, error_type="ignore",
                 error_value=np.nan):
        if not isinstance(src_kwd, str):
            raise TypeError(
                "The source keyword name must be a string")

        if not isinstance(dst_name, str):
            raise TypeError(
                "The destination name must be a string")

        if agg_func is not None:
            try:
                for i in agg_func:
                    if not hasattr(i, '__call__'):
                        raise TypeError(
                            "The aggregating function must be a callable " +
                            "object, None or a sequence of callables")
                self.agg_func_is_sequence = True
            except TypeError:
                if not hasattr(agg_func, '__call__'):
                    raise TypeError(
                        "The aggregating function must be a callable object, "
                        "None or a sequence of callables")
                self.agg_func_is_sequence = False

        if error_type not in ('ignore', 'raise', 'constant'):
            raise ValueError(
                "The error type must be either 'ignore', 'raise' or 'constant'")

        self.src_kwd = src_kwd
        self.dst_name = dst_name
        self.agg_func = agg_func
        self.error_type = error_type
        self.error_value = error_value


def metablender(input_models, spec):
    """
    Given a list of datamodels, aggregate metadata attribute values and
    create a table made up of values from a number of metadata instances,
    according to the given specification.

    **Parameters:**

    - *input_models* is a sequence where each element is either:

      - a `datamodels.DataModel` instance or sub-class

      - a string giving the *filename* for the input_model

    - *spec* is a list defining which keyword arguments are to be
      aggregated and how.  Each element in the list should be a
      sequence with 2 to 5 elements of the form:

        (*src_keyword*, *dst_name*, *function*, *error_type*, *error_value*)

      - *src_keyword* is the keyword to pull values from.  It is
        case-insensitive.

      - *dst_name* is the name to use as a dictionary key or column
        name for the destination values.

      - *function* (optional).  If function is not None, the values
        from the source are aggregated and returned in the
        *aggregate_dict*.  If function is None (or the tuple contains
        only 2 elements), all values are stored as a column with the
        name *dst_name* in the result *table*.

        If not None, *function* should be a callable object that takes
        a sequence of values and returns an aggregate result.  If the
        function returns None, no values will be added to the
        aggregate dictionary.  There are many functions in Numpy that
        are directly useful as an aggregating function, for example:

          - mean: `numpy.mean`

          - median: `numpy.median`

          - maximum: `numpy.max`

          - minimum: `numpy.min`

          - sum: `numpy.sum`

          - standard deviation: `numpy.std`

        Lambda functions are also often useful:

          - first: ``lambda x: x[0]``

          - last: ``lambda x: x[-1]``

        Additionally, *function* may be a tuple, where each member is
        itself a callable object.  The result will be a tuple
        containing results from each of the given functions.  For
        instance, to aggregate a range of values, i.e. both the
        minimum and maximum values, use the following as *function*:
        ``(numpy.min, numpy.max)``.

      - *error_type* (optional) defines how missing or syntax-errored
        values are handled.  It may be one of the following:

        - 'ignore': missing or unparsable values are ignored.  They
          are not included in the list of values passed to the
          aggregating function.  In the result *table*, missing values
          are masked out.

        - 'raise': missing or unparsable values raise a `ValueError`
          exception.

        - 'constant': missing or unparsable values are replaced with a
          constant, given by the *error_value* field.

      - *error_value* (optional) is the constant value to be used for
        missing or unparsable values when *error_type* is set to
        'constant'.  When not provided, it defaults to `NaN`.

    **Returns:**

    A 2-tuple of the form (*aggregate_dict*, *table*) where:

    - *aggregate_dict* is a dictionary of where the keys come from
      *dst_name* and the values are the aggregated values as run_KeywordMapping
      through *function*.

    - *table* is a masked Numpy structured array where the column
      names come from *dst_name* and the column contains the values
      from *src_keyword* for all of the given headers.  Missing values
      are masked out.
    """
    mappings = [_KeywordMapping(*x) for x in spec]
    data = [[] for x in spec]
    data_masks = [[] for x in spec]

    # Read in data
    for model in input_models:
        if not isinstance(model, DataModel):
            if not isinstance(model, str):
                raise TypeError(
                    "Each entry in the headers list must be either a " +
                    "datamodels.DataModel instance or a filename (str)")
            model = datamodels.open(model)
        header = model.to_flat_dict()
        filename = header['meta.filename']
        for i, mapping in enumerate(mappings):
            if mapping.src_kwd in header:
                value = model[mapping.src_kwd]
            elif mapping.error_type == 'raise':
                raise ValueError(
                    "%s is missing keyword '%s'" %
                    (filename, mapping.src_kwd))
            elif mapping.error_type == 'constant':
                value = mapping.error_value
            else:
                value = None

            if mapping.agg_func is None:
                if value is None:
                    data[i].append(np.nan)
                    data_masks[i].append(True)
                else:
                    data[i].append(value)
                    data_masks[i].append(False)
            else:
                if value is not None:
                    data[i].append(value)

    # Aggregate data into dictionary
    results = {}
    for i, mapping in enumerate(mappings):
        if data[i] == []:
            result = None
            continue
        if mapping.agg_func is not None:
            if mapping.agg_func_is_sequence:
                result = []
                for func in mapping.agg_func:
                    result.append(func(data[i]))
                result = tuple(result)
            else:
                result = mapping.agg_func(data[i])
            if result is not None:
                results[mapping.dst_name] = result

    # Aggregate data into table
    dtype = []
    arrays = []

    # Use Numpy to "guess" a data type for each of the columns
    for i, mapping in enumerate(mappings):
        if mapping.agg_func is None:
            array = np.array(data[i])
            if np.issubdtype(np.int32, array.dtype):
                # see about recasting as int32
                if not np.any(array / (2**31 - 1) > 1.):
                    array = array.astype(np.int32)
            dtype.append((mapping.dst_name, array.dtype))
            arrays.append(array)

    if len(dtype):
        # Combine the columns into a structured array
        table = ma.empty((len(input_models),), dtype=dtype)
        j = 0
        for i, mapping in enumerate(mappings):
            if mapping.agg_func is None:
                table.data[mapping.dst_name] = arrays[j]
                table.mask[mapping.dst_name] = data_masks[i]
                j += 1
    else:
        table = np.empty((0,))

    return results, table
