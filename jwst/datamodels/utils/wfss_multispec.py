"""Utilities for manipulating WFSS multi-spectral data."""

import numpy as np
import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.flat_multispec import (
    set_schema_units,
    copy_column_units,
    copy_spec_metadata,
    determine_vector_and_meta_columns,
    make_empty_recarray,
    populate_recarray,
    expand_table,
)


def make_wfss_multiexposure(input_list):
    """
    Compile a list of extracted sources into a single binary table.

    The output model will contain one binary table per exposure,
    with each table containing all sources extracted from that exposure
    (one row per source). The number of elements in each table row
    will be the same across all exposures, with NaNs used to pad
    shorter rows to match the longest row in the exposure.

    Parameters
    ----------
    input_list : MultiSpecModel or list[MultiSpecModel]
        List of MultiSpecModel objects to be combined.

    Returns
    -------
    output_x1d : WFSSMultiSpecModel
        The extract_1d product for WFSS modes.
    """
    if isinstance(input_list, dm.JwstDataModel):
        # if input is a single MultiSpecModel, convert to list
        input_list = [input_list]

    results_list = []
    for model in input_list:
        if isinstance(model, dm.WFSSMultiSpecModel):
            results_list.extend(wfss_multiexposure_to_multispec(model))
        elif isinstance(model, dm.MultiSpecModel):
            results_list.append(model)
        else:
            raise TypeError(
                "Input must be MultiSpecModel or WFSSMultiSpecModel, or a list of these models."
            )

    # first loop over source and exposure to figure out
    # final n_rows, n_exposures, n_sources
    exposure_counter = {}
    all_source_ids = []
    # for calwebb_spec3 the outer loop is over sources and the inner loop is over exposures
    # for calwebb_spec2 it is the opposite, but outer loop (exposures) should have just one element
    for model in results_list:
        for spec in model.spec:
            fname = getattr(spec.meta, "filename", None)
            exp_number = getattr(spec.meta, "group_id", None)

            # if this is the first time this exposure has been encountered,
            # create a new dictionary entry for it
            if exp_number not in exposure_counter.keys():
                n_rows = spec.spec_table.shape[0]
                exposure_counter[exp_number] = {
                    "n_rows": n_rows,
                    "filename": fname,
                    "exposure_time": model.meta.exposure.exposure_time,  # need for combine_1d
                    "integration_time": model.meta.exposure.integration_time,  # need for combine_1d
                    "spectral_order": spec.spectral_order,
                }
            else:
                # if this exposure has already been encountered,
                # check if number of rows is larger than the previous one
                exposure_counter[exp_number]["n_rows"] = max(
                    exposure_counter[exp_number]["n_rows"], spec.spec_table.shape[0]
                )

            all_source_ids.append(spec.source_id)

    all_source_ids = np.sort(np.unique(all_source_ids))
    n_sources = len(all_source_ids)

    exposure_numbers = list(exposure_counter.keys())
    n_exposures = len(exposure_numbers)
    n_rows_by_exposure = [exposure_counter[n]["n_rows"] for n in exposure_numbers]

    # Set up output table column names and dtypes
    # Use SpecModel.spectable to determine the vector-like columns
    # The additional metadata columns are all those that are defined in WFSSMultiSpecModel
    # but not in SpecModel
    input_datatype = dm.SpecModel().schema["properties"]["spec_table"]["datatype"]
    output_datatype = dm.WFSSSpecModel().schema["properties"]["spec_table"]["datatype"]
    all_columns, is_vector = determine_vector_and_meta_columns(input_datatype, output_datatype)
    defaults = dm.WFSSSpecModel().schema["properties"]["spec_table"]["default"]

    # loop over exposures to make empty tables for each exposure, which are initially populated
    # with the schema defaults (NaNs for float columns).
    # All tables should have the same source_id column for ease of indexing,
    # even if some sources are not present in some exposures.
    fltdata_by_exposure = []
    for i in range(n_exposures):
        n_rows = n_rows_by_exposure[i]
        flt_empty = make_empty_recarray(
            n_rows, n_sources, all_columns, is_vector, defaults=defaults
        )
        flt_empty["SOURCE_ID"] = all_source_ids
        fltdata_by_exposure.append(flt_empty)

    # Now loop through the models and populate the tables
    for model in results_list:
        for spec in model.spec:
            # ensure data goes to correct exposure table based on group_id attribute
            exp_num = spec.meta.group_id
            exposure_idx = exposure_numbers.index(exp_num)
            fltdata = fltdata_by_exposure[exposure_idx]
            n_rows = n_rows_by_exposure[exposure_idx]

            # ensure data goes to the correct source
            spec_idx = np.where(fltdata["SOURCE_ID"] == spec.source_id)[0][0]

            # populate the table with data from the input spectrum
            populate_recarray(
                fltdata[spec_idx],
                spec,
                all_columns,
                is_vector,
                ignore_columns=["SOURCE_ID", "N_ALONGDISP"],
            )

            # special handling for N_ALONGDISP because not defined in specmeta schema
            fltdata[spec_idx]["N_ALONGDISP"] = spec.spec_table.shape[0]

    # Finally, create a new WFSSMultiSpecModel to hold the combined data
    # with one WFSSMultiSpecModel table per exposure
    output_x1d = dm.WFSSMultiSpecModel()
    example_spec = results_list[0].spec[0]
    for i, exposure_number in enumerate(exposure_numbers):
        spec_table = fltdata_by_exposure[i]
        ext = dm.WFSSSpecModel(spec_table)

        # Set default units from the model schema
        set_schema_units(ext)
        # copy units from the example specmodel, overriding the schema defaults where applicable
        copy_column_units(example_spec, ext)

        # copy metadata
        ext.filename = exposure_counter[exposure_number]["filename"]
        ext.group_id = exposure_number
        ext.dispersion_direction = example_spec.dispersion_direction
        ext.spectral_order = exposure_counter[exposure_number]["spectral_order"]
        ext.exposure_time = exposure_counter[exposure_number]["exposure_time"]
        ext.integration_time = exposure_counter[exposure_number]["integration_time"]

        output_x1d.spec.append(ext)

    output_x1d.update(input_list[0], only="PRIMARY")
    return output_x1d


def wfss_multiexposure_to_multispec(input_model):
    """
    Transform a WFSSMultiSpecModel into a list of MultiSpecModel objects.

    Parameters
    ----------
    input_model : WFSSMultiSpecModel
        Input model to be reorganized.

    Returns
    -------
    output_list : list[MultiSpecModel]
        List of MultiSpecModel objects, one for each exposure in the input model.
    """
    # first extract all spectra as SpecModels in a flat list
    spec_list = []
    source_ids = []
    exposure_times = []
    integration_times = []
    for exp in input_model.spec:
        this_exp_list = expand_table(exp)
        source_ids.extend([spec.source_id for spec in this_exp_list])
        exposure_times.extend([exp.exposure_time for _ in this_exp_list])
        integration_times.extend([exp.integration_time for _ in this_exp_list])
        spec_list.extend(this_exp_list)

    # organize by unique source id such that there is one MultiSpecModel per source
    # with all exposures for that source in that model's model.spec
    unique_source_ids = set(source_ids)
    source_ids = np.array(source_ids)
    spec_list = np.array(spec_list)
    output_list = []
    for source_id in unique_source_ids:
        multispec = dm.MultiSpecModel()
        # BUG: currently there is no infrastructure for handling per-exposure weights
        # This is also a problem on main
        multispec.meta.exposure.exposure_time = exposure_times[0]
        multispec.meta.exposure.integration_time = integration_times[0]
        spec_this_id = spec_list[source_ids == source_id]
        multispec.spec.extend(spec_this_id)
        multispec.update(input_model, only="PRIMARY")
        output_list.append(multispec)

    return output_list


def make_wfss_multicombined(results_list):
    """
    Compile a list of exposure-averaged extracted sources into a single binary table.

    The output model will contain a single spec_table,
    which will end up in extension index 1 of the FITS file on save.

    Parameters
    ----------
    results_list : list[MultiCombinedSpecModel]
        List of MultiSpecModel objects to be combined.

    Returns
    -------
    output_c1d : WFSSMultiCombinedSpecModel
        The combined c1d product for WFSS modes.
    """
    # determine shape of output table.
    # Each input model should have one spec table per spectral order
    n_sources = len(results_list)

    # figure out column names and dtypes
    input_datatype = dm.CombinedSpecModel().schema["properties"]["spec_table"]["datatype"]
    output_schema = dm.WFSSMultiCombinedSpecModel().schema
    output_table_schema = output_schema["properties"]["spec"]["items"]["properties"]["spec_table"]
    output_datatype = output_table_schema["datatype"]
    all_columns, is_vector = determine_vector_and_meta_columns(input_datatype, output_datatype)
    defaults = output_table_schema["default"]

    # first loop over models and spec to figure out orders present, and number of rows per order
    order_rows = {}
    for model in results_list:
        for spec in model.spec:
            order = spec.spectral_order
            if order not in order_rows:
                order_rows[order] = spec.spec_table.shape[0]
            else:
                order_rows[order] = max(order_rows[order], spec.spec_table.shape[0])

    order_data = {}
    for j, model in enumerate(results_list):
        for spec in model.spec:
            # ensure data goes to the correct order
            order = spec.spectral_order
            n_rows = order_rows[order]
            if order not in order_data:
                order_data[order] = make_empty_recarray(
                    n_rows, n_sources, all_columns, is_vector, defaults=defaults
                )
            fltdata = order_data[order]

            # populate the table with data from the input spectrum
            populate_recarray(
                fltdata[j],
                spec,
                all_columns,
                is_vector,
                ignore_columns=["N_ALONGDISP"],
            )
            # special handling for N_ALONGDISP because not defined in specmeta schema
            fltdata[j]["N_ALONGDISP"] = spec.spec_table.shape[0]

    # Create a new model to hold the combined data table
    output_c1d = dm.WFSSMultiCombinedSpecModel()

    example_spec = results_list[0].spec[0]
    for order in order_data.keys():
        spec = dm.WFSSCombinedSpecModel()
        fltdata = order_data[order]
        fltdata.sort(order=["SOURCE_ID"])
        spec.spec_table = fltdata

        # Set default units from the model schema
        set_schema_units(spec)
        # copy units from any of the SpecModels (they should all be the same)
        copy_column_units(example_spec, spec)
        copy_spec_metadata(example_spec, spec)
        spec.spectral_order = order

        output_c1d.spec.append(spec)

    output_c1d.update(results_list[0], only="PRIMARY")
    return output_c1d
