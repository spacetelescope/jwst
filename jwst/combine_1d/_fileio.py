"""Utility functions for saving WFSS c1d products."""

import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.flat_multispec import (
    copy_column_units,
    determine_vector_and_meta_columns,
    make_empty_recarray,
    populate_recarray,
)


def save_wfss_c1d(results_list, output_filename):
    """
    Compile exposure-averaged sources into a single table and save to a file.

    The output c1d product will contain a single spec_table,
    which will end up in extension index 1 of the FITS file.

    Parameters
    ----------
    results_list : list[MultiSlitModel]
        List of MultiSlitModel objects to be combined into a single c1d file.
    output_filename : str
        Name of the output c1d file.
    """
    # determine shape of output table
    # each input model should have just one spec table
    n_sources = len(results_list)
    n_rows = max(len(model.spec[0].spec_table) for model in results_list)

    # figure out column names and dtypes
    input_datatype = dm.CombinedSpecModel().schema["properties"]["spec_table"]["datatype"]
    output_datatype = dm.WFSSMultiCombinedSpecModel().schema["properties"]["spec_table"]["datatype"]
    all_columns, is_vector = determine_vector_and_meta_columns(input_datatype, output_datatype)
    defaults = dm.WFSSMultiCombinedSpecModel().schema["properties"]["spec_table"]["default"]

    # create empty table
    fltdata = make_empty_recarray(n_rows, n_sources, all_columns, is_vector, defaults=defaults)

    # loop over sources to populate the table with data from the input spectrum
    for j, model in enumerate(results_list):
        populate_recarray(
            fltdata[j],
            model.spec[0],
            n_rows,
            all_columns,
            is_vector,
            ignore_columns=["NELEMENTS"],
        )
        # special handling for NELEMENTS because not defined in specmeta schema
        fltdata[j]["NELEMENTS"] = model.spec[0].spec_table.shape[0]

    # Create a new SpecModel to hold the combined data
    # with one SpecModel per exposure
    output_c1d = dm.WFSSMultiCombinedSpecModel()
    fltdata.sort(order=["SOURCE_ID"])
    output_c1d.spec_table = fltdata
    example_model = results_list[0]
    output_c1d.update(example_model)

    # copy units from any of the SpecModels (they should all be the same)
    copy_column_units(model.spec[0], output_c1d)

    output_c1d.save(output_filename)
