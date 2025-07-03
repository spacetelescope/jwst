"""Utilities for re-organizing spectral products into a flat structure."""

import logging
from copy import deepcopy

import numpy as np
from asdf.tags.core.ndarray import asdf_datatype_to_numpy_dtype
from stdatamodels.jwst import datamodels

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


def populate_recarray(output_table, input_spec, columns, is_vector, ignore_columns=None):
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

    # Copy the data into the new table
    for col, _ in vector_columns:
        if col in ignore_columns:
            continue

        output_table[col][: input_table.shape[0]] = input_table[col]

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


def set_schema_units(model):
    """
    Give all columns in the model the units defined in the model schema.

    This gets around a bug/bad behavior in stdatamodels that units are not
    automatically assigned to the spec_table.

    Model is modified in place.

    Parameters
    ----------
    model : DataModel
        Any model containing a spec_table attribute.
    """
    data_type = model.schema["properties"]["spec_table"]["datatype"]
    for col in data_type:
        if "unit" in col:
            model.spec_table.columns[col["name"]].unit = col["unit"]


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


def copy_spec_metadata(input_model, output_model):
    """
    Copy spectral metadata from the input to the output spectrum.

    Values to be copied are any attributes of the input model,
    other than "meta" or "spec_table", e.g. "source_id", "name", etc.

    Parameters
    ----------
    input_model : DataModel or ObjectNode
        A spectral model, such as SpecModel or TSOSpecModel. If read
        in from a list of spectra, as in MultiSpecModel, the input model may be
        an ObjectNode rather than a full DataModel.
    output_model : DataModel
        A spectral model, such as SpecModel or TSOSpecModel. Updated in place
        with metadata from the input model.  The output model must be a full
        DataModel, not an ObjectNode.
    """
    copy_attributes = []
    for prop in output_model.schema["properties"]:
        if prop not in ["meta", "spec_table"]:
            copy_attributes.append(prop)
    for key in copy_attributes:
        if hasattr(input_model, key) and getattr(input_model, key) is not None:
            setattr(output_model, key, getattr(input_model, key))


def expand_table(spec):
    """
    Expand a table of spectra into a list of SpecModel objects.

    Parameters
    ----------
    spec : WFSSSpecModel, TSOMultiSpecModel, ObjectNode
        Any model containing a spec_table to expand into multiple spectra

    Returns
    -------
    list[SpecModel]
        A list of SpecModel objects, one for each spectrum in the input spec_table.
    """
    all_columns = np.array([str(x) for x in spec.spec_table.dtype.names])
    new_spec_list = []
    n_spectra = len(spec.spec_table)
    for i in range(n_spectra):
        # initialize a new SpecModel
        spec_row = spec.spec_table[i]
        n_elements = int(spec_row["N_ALONGDISP"])
        new_spec = datamodels.SpecModel()
        data_type = new_spec.schema["properties"]["spec_table"]["datatype"]
        columns_to_copy = np.array([col["name"] for col in data_type])

        # Copy over the vector columns from input spec_table to output spec_table
        spec_table = np.empty(n_elements, dtype=new_spec.spec_table.dtype)
        for col_name in columns_to_copy:
            spec_table[col_name] = spec_row[col_name][:n_elements]
        new_spec.spec_table = spec_table

        # Copy over the metadata columns from input spec_table to the spectrum's metadata
        meta_columns = all_columns[~np.isin(all_columns, columns_to_copy)].tolist()
        meta_columns.remove("N_ALONGDISP")
        for meta_key in meta_columns:
            try:
                setattr(new_spec, meta_key.lower(), spec_row[meta_key])
            except KeyError:
                pass

        # Copy over relevant metadata from the input model to the output model
        if hasattr(spec.meta, "wcs"):
            new_spec.meta.wcs = deepcopy(spec.meta.wcs)
        new_spec.meta.group_id = getattr(spec, "group_id", "")
        new_spec.meta.filename = getattr(spec, "filename", "")
        copy_spec_metadata(spec, new_spec)
        copy_column_units(spec, new_spec)

        new_spec_list.append(new_spec)

    return new_spec_list


def expand_flat_spec(input_model):
    """
    Create simple spectra from a flat spectral table.

    Parameters
    ----------
    input_model : TSOMultiSpecModel
        Spectral model containing spectra with a mix of vector columns
        and metadata columns in the ``spec_table`` attribute.
        Metadata columns will be dropped.

    Returns
    -------
    MultiSpecModel
        A set of simple spectra, one per extension.
    """
    output_model = datamodels.MultiSpecModel()
    for old_spec in input_model.spec:
        new_spec_list = expand_table(old_spec)
        for new_spec in new_spec_list:
            # Add the new spec to the output model
            output_model.spec.append(new_spec)

    # Update meta
    output_model.update(input_model, only="PRIMARY")

    # Copy int_times if present
    if hasattr(input_model, "int_times"):
        output_model.int_times = input_model.int_times.copy()

    return output_model
