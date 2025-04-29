"""Pipeline utilities objects."""

import logging
import warnings

import numpy as np
from asdf.tags.core.ndarray import asdf_datatype_to_numpy_dtype
from stdatamodels.properties import ObjectNode
from stdatamodels.jwst.datamodels import dqflags, JwstDataModel

from jwst.associations.lib.dms_base import TSO_EXP_TYPES


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def is_tso(model):
    """
    Check if data is a Time Series Observation (TSO).

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        Data to check.

    Returns
    -------
    is_tso : bool
       `True` if the model represents TSO data.
    """
    is_tso = False

    # Check on JWST-specific TSOVISIT flag
    try:
        is_tso = model.meta.visit.tsovisit
    except AttributeError:
        pass

    # Check on exposure types
    try:
        is_tso = is_tso or model.meta.exposure.type.lower() in TSO_EXP_TYPES
    except AttributeError:
        pass

    # Check on number of integrations
    try:
        if model.meta.exposure.nints is not None and model.meta.exposure.nints < 2:
            is_tso = False
    except AttributeError:
        pass

    # We've checked everything.
    return is_tso


def is_irs2(model):
    """
    Check whether the data are in IRS2 format.

    This currently assumes that only full-frame, near-infrared data can be
    taken using the IRS2 readout pattern.

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel` or ndarray
        Data to check.

    Returns
    -------
    status : bool
       `True` if the data are in IRS2 format.
    """
    if isinstance(model, np.ndarray):
        shape = model.shape
    else:
        try:
            shape = model.data.shape
        except AttributeError:
            return False

    max_length = 2048

    irs2_axis_length = max(shape[-1], shape[-2])

    if irs2_axis_length > max_length:
        return True
    else:
        return False


def match_nans_and_flags(input_model):
    """
    Ensure data, error, variance, and DQ are marked consistently for invalid data.

    Invalid data is assumed to be any pixel set to NaN in any one of the
    data, error, or variance arrays, or else set to the DO_NOT_USE flag
    in the DQ array.

    The input model is updated in place with NaNs or DO_NOT_USE flags, as
    appropriate, at all invalid data locations.

    Parameters
    ----------
    input_model : DataModel
        Input model containing some combination of data, dq, err, var_rnoise,
        var_poisson, and var_flat extensions. These extensions must all have
        matching dimensions if present.
    """
    # Check for datamodel input or slit instance
    if not isinstance(input_model, JwstDataModel) and not isinstance(input_model, ObjectNode):
        raise TypeError(f"Input {type(input_model)} is not a datamodel.")

    # Build up the invalid data flags from each available data extension.
    is_invalid = None
    data_shape = None
    nan_extensions = ["data", "err", "var_rnoise", "var_poisson", "var_flat"]
    for extension in nan_extensions:
        if not hasattr(input_model, extension):
            continue
        data = getattr(input_model, extension)
        if is_invalid is None:
            is_invalid = np.isnan(data)
            data_shape = data.shape
        else:
            if data.shape != data_shape:
                log.warning(
                    "Mismatched data shapes; skipping invalid data updates for extension '%s'",
                    extension,
                )
                continue
            is_invalid |= np.isnan(data)

    # Nothing to do if no extensions were found to update
    if is_invalid is None:
        return

    # Add in invalid flags from the DQ extension if present
    if hasattr(input_model, "dq"):
        do_not_use = (input_model.dq & dqflags.pixel["DO_NOT_USE"]).astype(bool)
        if input_model.dq.shape != data_shape:
            log.warning("Mismatched data shapes; skipping invalid data updates for extension 'dq'")
        else:
            is_invalid |= do_not_use

    # Update all the data extensions
    for extension in nan_extensions:
        if not hasattr(input_model, extension):
            continue
        data = getattr(input_model, extension)
        if data.shape != data_shape:
            continue
        data[is_invalid] = np.nan

    # Update the DQ extension
    if input_model.dq.shape == data_shape:
        input_model.dq[is_invalid] |= dqflags.pixel["DO_NOT_USE"]


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


def make_empty_recarray(n_rows, n_spec, vector_columns, meta_columns):
    """
    Create an empty output table with the specified number of rows.

    Parameters
    ----------
    n_rows : int
        The number of rows in the output table; this is the maximum number of
        data points for any spectrum in the exposure.
    n_spec : int
        The number of spectra in the output table.
    vector_columns : list[tuple]
        List of tuples containing the vector-like column names and their dtypes.
    meta_columns : list[tuple]
        List of tuples containing the metadata column names and their dtypes.

    Returns
    -------
    output_table : `~numpy.recarray`
        The empty output table with the specified shape and dtypes.
    """
    fltdtype = []
    for col, dtype in vector_columns:
        fltdtype.append((col, dtype, n_rows))
    for col, dtype in meta_columns:
        fltdtype.append((col, dtype))
    return np.empty(n_spec, dtype=fltdtype)


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
