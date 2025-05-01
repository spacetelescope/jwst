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
    columns : np.ndarray[tuple]
        Array of tuples containing the column names and their dtypes.
    is_vector : np.ndarray[bool]
        Array of booleans indicating whether each column is vector-like,
        same length as `columns`.
    """
    # Extract just names and dtypes, convert to numpy dtypes
    vector_colnames = np.array([col["name"] for col in input_datatype])
    all_colnames = np.array([col["name"] for col in output_datatype])
    all_dtypes = np.array(
        [asdf_datatype_to_numpy_dtype(col["datatype"]) for col in output_datatype]
    )
    all_cols = np.array(list(zip(all_colnames, all_dtypes, strict=True)))

    # Determine which columns are metadata
    is_vector = np.array([col in vector_colnames for col in all_colnames])

    return all_cols, is_vector


def make_empty_recarray(n_rows, n_spec, columns, is_vector, defaults=0):
    """
    Create an empty output table with the specified number of rows.

    Parameters
    ----------
    n_rows : int
        The number of rows in the output table; this is the maximum number of
        data points for any spectrum in the exposure.
    n_spec : int
        The number of spectra in the output table.
    columns : np.ndarray[tuple]
        Array of tuples containing the column names and their dtypes.
    is_vector : np.ndarray[bool]
        Array of booleans indicating whether each column is vector-like.
        If `True`, the column will be a 1D array of length `n_rows`.
        Otherwise, the column will be a scalar.
    defaults : list, np.ndarray, int, or float, optional
        List of default values for each column. If a column is vector-like,
        the default value will be repeated to fill the array.
        If a column is scalar, the default value will be used directly.
        If `int` or `float`, the same value will be used for all columns;
        string-type columns will be filled with the string representation of the value.

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
    if isinstance(defaults, (int, float)):
        arr[...] = defaults
        return arr

    # fill the array with the default values
    for i, (col, dtype) in enumerate(columns):
        if is_vector[i]:
            arr[col] = np.full((n_spec, n_rows), defaults[i], dtype=dtype)
        else:
            arr[col] = np.full(n_spec, defaults[i], dtype=dtype)
    return arr


def populate_recarray(output_table, input_spec, n_rows, columns, is_vector, ignore_columns=None):
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
    columns : np.ndarray[tuple]
        Array of tuples containing the column names and their dtypes.
    is_vector : np.ndarray[bool]
        Array of booleans indicating whether each column is vector-like,
    ignore_columns : list[str], optional
        List of column names to ignore when copying data or metadata from the input
        spectrum to the output table. This is useful for columns that are not
        present in the input spectrum but are required in the output table,
        and are handled separately in the calling code.
    """
    if ignore_columns is None:
        ignore_columns = []
    input_table = input_spec.spec_table

    vector_columns = columns[is_vector]
    meta_columns = columns[~is_vector]

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


def copy_column_units(input_model, output_model):
    """
    Copy units from input columns to output columns.

    Spectral tables in both input and output models must be
    in FITS record format. The output model is updated in place.

    Parameters
    ----------
    input_model : SpecModel
        Input spectral model containing vector columns in the
        ``spec_table`` attribute.
    output_model : DataModel
        Output spectral model containing a mix of vector columns
        and metadata columns in the ``spec_table`` attribute.
    """
    input_columns = input_model.spec_table.columns
    output_columns = output_model.spec_table.columns
    for col_name in input_columns.names:
        if col_name in output_columns.names:
            output_columns[col_name].unit = input_columns[col_name].unit
