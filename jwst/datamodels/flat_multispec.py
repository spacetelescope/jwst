"""Utilities for re-organizing spectral products into a flat structure."""

import logging
import warnings
import numpy as np
from asdf.tags.core.ndarray import asdf_datatype_to_numpy_dtype


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def determine_vector_and_meta_columns(input_datatype, output_datatype):
    """
    Figure out which columns are vector-like and which are metadata.

    The vector-like columns are the ones defined in the input schema,
    and the metadata columns are the ones defined only in the output schema.
    The input and output datatypes are typically read from the schema as e.g.
    datatype = schema["properties"]["spec_table"]["datatype"].

    Parameters
    ----------
    input_datatype : list[dict]
        The datatype of the input model as read from the schema.
        Each inner dict should have at least the keys "name" and "datatype".
    output_datatype : list[dict]
        The datatype of the output model as read from the schema.
        Each inner dict should have at least the keys "name" and "datatype".

    Returns
    -------
    vector_columns : list[tuple]
        List of tuples containing the vector-like column names and their dtypes.
    meta_columns : list[tuple]
        List of tuples containing the metadata column names and their dtypes.
    """
    # Extract just names and dtypes, convert to numpy dtypes
    vector_colnames = np.array([col["name"] for col in input_datatype])
    vector_dtypes = np.array(
        [asdf_datatype_to_numpy_dtype(col["datatype"]) for col in input_datatype]
    )
    all_colnames = np.array([col["name"] for col in output_datatype])
    all_dtypes = np.array(
        [asdf_datatype_to_numpy_dtype(col["datatype"]) for col in output_datatype]
    )

    # Determine which columns are metadata
    is_meta = ~np.array([col in vector_colnames for col in all_colnames])
    meta_colnames = all_colnames[is_meta]
    meta_dtypes = all_dtypes[is_meta]

    # Construct the vector and meta column tuples
    vector_cols = list(zip(vector_colnames, vector_dtypes, strict=True))
    meta_cols = list(zip(meta_colnames, meta_dtypes, strict=True))

    return vector_cols, meta_cols


def make_empty_recarray(n_rows, n_spec, columns, is_vector, defaults=None):
    """
    Create an empty output table with the specified number of rows.

    Parameters
    ----------
    n_rows : int
        The number of rows in the output table; this is the maximum number of
        data points for any spectrum in the exposure.
    n_spec : int
        The number of spectra in the output table.
    columns : list[tuple]
        List of tuples containing the column names and their dtypes.
    is_vector : list[bool]
        List of booleans indicating whether each column is vector-like.
        If `True`, the column will be a 1D array of length `n_rows`.
        Otherwise, the column will be a scalar.
    defaults : list, optional
        List of default values for each column. If a column is vector-like,
        the default value will be repeated to fill the array.
        If a column is scalar, the default value will be used directly.
        If None, the array will be empty.

    Returns
    -------
    output_table : `~numpy.recarray`
        The empty output table with the specified shape and dtypes.
    """
    # build the data type
    fltdtype = []
    for i, (col, dtype) in enumerate(columns):
        if is_vector[i]:
            fltdtype.append((col, dtype, n_rows))
        else:
            fltdtype.append((col, dtype))

    arr = np.empty(n_spec, dtype=fltdtype)
    if defaults is None:
        return arr

    # fill the array with the default values
    for i, (col, dtype) in enumerate(columns):
        if is_vector[i]:
            arr[col] = np.full((n_spec, n_rows), defaults[i], dtype=dtype)
        else:
            arr[col] = np.full(n_spec, defaults[i], dtype=dtype)
    return arr


def populate_recarray(
    output_table, input_spec, n_rows, vector_columns, meta_columns, ignore_columns=None
):
    """
    Populate the output table in-place with data from the input spectrum.

    The output table is padded with NaNs to match the maximum number of
    data points for any spectrum in the exposure.
    The metadata columns are copied from the input spectrum assuming
    they have the same names as in the output table.

    Parameters
    ----------
    output_table : `~numpy.recarray`
        The output table to be populated with the spectral data.
    input_spec : `~jwst.datamodels.SpecModel` or `~jwst.datamodels.CombinedSpecModel`
        The input data model containing the spectral data.
    n_rows : int
        The number of rows in the output table; this is the maximum number of
        data points for any spectrum in the exposure.
    vector_columns : list[tuple]
        List of tuples containing the vector-like column names and their dtypes.
    meta_columns : list[tuple]
        List of tuples containing the metadata column names and their dtypes.
    ignore_columns : list[str], optional
        List of column names to ignore when copying data or metadata from the input
        spectrum to the output table. This is useful for columns that are not
        present in the input spectrum but are required in the output table,
        and are handled separately in the calling code.
    """
    if ignore_columns is None:
        ignore_columns = []
    input_table = input_spec.spec_table

    # Copy the data into the new table with NaN padding
    for col, _ in vector_columns:
        if col in ignore_columns:
            continue
        padded_data = np.full(n_rows, np.nan)
        padded_data[: input_table.shape[0]] = input_table[col]
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=RuntimeWarning, message="invalid value encountered in cast"
            )
            output_table[col] = padded_data

    # Copy the metadata into the new table
    # Metadata columns must have identical names to spec_meta columns
    problems = []
    for col, _ in meta_columns:
        if col in ignore_columns:
            continue

        spec_meta = getattr(input_spec, col.lower(), None)
        if spec_meta is None:
            problems.append(col.lower())
        else:
            output_table[col] = spec_meta

    if len(problems) > 0:
        log.warning(f"Metadata could not be determined from input spec_table: {problems}")
