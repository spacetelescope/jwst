import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels.utils import flat_multispec


def make_tso_specmodel(spec_list, segment=None):
    """
    Make a TSOSpecModel from a list of SpecModel.

    Parameters
    ----------
    spec_list : list of SpecModel
       Individual spectra for each integration.
    segment : int or None, optional
       The segment number for the input model.

    Returns
    -------
    TSOSpecModel
        A model containing all spectra in a flat table with vector
        columns for spectral data.
    """
    # Make a blank model
    tso_spec = datamodels.TSOSpecModel()

    # Compare input and output schema to get vector columns and meta columns
    input_schema = datamodels.SpecModel().schema
    in_cols = input_schema["properties"]["spec_table"]["datatype"]
    out_cols = tso_spec.schema["properties"]["spec_table"]["datatype"]
    all_cols, is_vector = flat_multispec.determine_vector_and_meta_columns(in_cols, out_cols)

    # Make an empty table to populate
    n_elements = [spec.spec_table["WAVELENGTH"].size for spec in spec_list]
    n_rows = np.max(n_elements)
    n_spectra = len(spec_list)
    defaults = tso_spec.schema["properties"]["spec_table"]["default"]
    spec_table = flat_multispec.make_empty_recarray(
        n_rows, n_spectra, all_cols, is_vector, defaults=defaults
    )

    # Populate the table from the input spectra
    time_keys = [
        "MJD-BEG",
        "MJD-AVG",
        "MJD-END",
        "TDB-BEG",
        "TDB-MID",
        "TDB-END",
    ]
    ignore_columns = ["N_ALONGDISP", "SEGMENT", "INT_NUM"] + time_keys
    for i in range(n_spectra):
        input_spec = spec_list[i]
        this_output = spec_table[i]

        flat_multispec.populate_recarray(
            this_output,
            input_spec,
            all_cols,
            is_vector,
            ignore_columns=ignore_columns,
        )

        # Update the special metadata columns
        this_output["N_ALONGDISP"] = n_elements[i]
        if segment is not None:
            this_output["SEGMENT"] = segment
        else:
            this_output["SEGMENT"] = 1

        # Set a default value for time keys - they'll be updated later if possible
        this_output["INT_NUM"] = i + 1
        for key in time_keys:
            this_output[key] = np.nan

    # Update the model with the new spec_table
    tso_spec.spec_table = spec_table

    # Add some units to the new columns
    flat_multispec.copy_column_units(spec_list[0], tso_spec)
    for column_name in time_keys:
        tso_spec.spec_table.columns[column_name].unit = "d"

    # Copy metadata from the first input_spec
    tso_spec.update(spec_list[0])
    if hasattr(spec_list[0].meta, "wcs"):
        tso_spec.meta.wcs = spec_list[0].meta.wcs

    return tso_spec
