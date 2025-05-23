"""Utilities for manipulating WFSS multi-spectral data."""

import numpy as np
import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.flat_multispec import (
    set_schema_units,
    copy_column_units,
    determine_vector_and_meta_columns,
    make_empty_recarray,
    populate_recarray,
    expand_table,
    NEW_NAMES,
    OLD_NAMES,
)


def make_wfss_multiexposure(input_list):
    """
    Combine all sources into a single table and save to a file.

    The output x1d product will have one extension per exposure.
    Each extension will contain a single table with one row per source.
    The table size is set by the maximum number of data points for any
    source in the exposure; the other sources will be end-padded with NaNs
    so their shape equals the maximum.

    Parameters
    ----------
    input_list : MultiSpecModel or list[MultiSpecModel]
        List of MultiSpecModel objects to be combined into a single x1d file.

    Returns
    -------
    output_x1d : WFSSMultiExposureSpecModel
        The combined x1d product for WFSS modes.
    """
    if isinstance(input_list, dm.JwstDataModel):
        # if input is a single MultiSpecModel, convert to list
        input_list = [input_list]

    results_list = []
    for model in input_list:
        if isinstance(model, dm.WFSSMultiExposureSpecModel):
            results_list.extend(wfss_multiexposure_to_multispec(model))
        elif isinstance(model, dm.MultiSpecModel):
            results_list.append(model)
        else:
            raise TypeError(
                "Input must be MultiSpecModel or WFSSMultiExposureSpecModel, "
                "or a list of these models."
            )

    # first loop over both source and exposure to figure out final n_rows, n_exposures, n_sources
    n_rows_by_exposure = []
    exposure_counter = {}
    # for calwebb_spec3 the outer loop is over sources and the inner loop is over exposures
    # for calwebb_spec2 it is the opposite, but outer loop (exposures) should have just one element
    for model in results_list:
        for spec in model.spec:
            fname = getattr(spec.meta, "filename", None)
            if hasattr(spec.meta, "observation"):
                exp_number = getattr(spec.meta.observation, "exposure_number", None)
            else:
                exp_number = None

            # if this is the first time this exposure has been encountered,
            # create a new dictionary entry for it
            if exp_number not in exposure_counter.keys():
                n_rows = spec.spec_table.shape[0]
                wcs = spec.meta.wcs
                exposure_counter[exp_number] = {
                    "n_rows": n_rows,
                    "n_sources": 1,
                    "filename": fname,
                    "wcs": wcs,
                    "exposure_time": model.meta.exposure.exposure_time,  # need for combine_1d
                    "integration_time": model.meta.exposure.integration_time,  # need for combine_1d
                }
            else:
                exposure_counter[exp_number]["n_sources"] += 1
                # if this exposure has already been encountered,
                # check if number of rows is larger than the previous one
                exposure_counter[exp_number]["n_rows"] = max(
                    exposure_counter[exp_number]["n_rows"], spec.spec_table.shape[0]
                )

    exposure_numbers = list(exposure_counter.keys())
    n_exposures = len(exposure_numbers)
    n_rows_by_exposure = [exposure_counter[n]["n_rows"] for n in exposure_numbers]
    n_sources_by_exposure = [exposure_counter[n]["n_sources"] for n in exposure_numbers]

    # Set up output table column names and dtypes
    # Use SpecModel.spectable to determine the vector-like columns
    # The additional metadata columns are all those that are defined in WFSSMultiSpecModel
    # but not in SpecModel
    input_datatype = dm.SpecModel().schema["properties"]["spec_table"]["datatype"]
    output_datatype = dm.WFSSMultiSpecModel().schema["properties"]["spec_table"]["datatype"]
    all_columns, is_vector = determine_vector_and_meta_columns(input_datatype, output_datatype)
    defaults = dm.WFSSMultiSpecModel().schema["properties"]["spec_table"]["default"]

    # loop over exposures to make tables for each exposure
    fltdata_by_exposure = []
    for i in range(n_exposures):
        n_rows = n_rows_by_exposure[i]
        n_sources = n_sources_by_exposure[i]
        flt_empty = make_empty_recarray(
            n_rows, n_sources, all_columns, is_vector, defaults=defaults
        )
        fltdata_by_exposure.append(flt_empty)

    # Now loop through the models and populate the tables
    # Need to index each exposure separately because they may have a different number of sources
    loop_index_by_exposure = [0] * n_exposures
    for model in results_list:
        for spec in model.spec:
            # ensure data goes to table corresponding to correct exposure based on filename
            exp_num = spec.meta.observation.exposure_number
            exposure_idx = exposure_numbers.index(exp_num)
            fltdata = fltdata_by_exposure[exposure_idx]
            n_rows = n_rows_by_exposure[exposure_idx]
            spec_idx = loop_index_by_exposure[exposure_idx]
            loop_index_by_exposure[exposure_idx] += 1

            # populate the table with data from the input spectrum
            populate_recarray(
                fltdata[spec_idx],
                spec,
                n_rows,
                all_columns,
                is_vector,
                ignore_columns=["NELEMENTS"] + NEW_NAMES,
            )

            # special handling for NELEMENTS because not defined in specmeta schema
            fltdata[spec_idx]["NELEMENTS"] = spec.spec_table.shape[0]

            # special handling for the names that we want renamed
            # at the request of the WFSS teams
            for old_name, new_name in zip(OLD_NAMES, NEW_NAMES, strict=True):
                fltdata[spec_idx][new_name] = getattr(spec, old_name.lower(), None)

    # Finally, create a new MultiExposureModel to hold the combined data
    # with one MultiSpecModel table per exposure
    output_x1d = dm.WFSSMultiExposureSpecModel()
    example_spec = results_list[0].spec[0]
    for i, exposure_number in enumerate(exposure_numbers):
        # Create a new extension for each exposure
        spec_table = fltdata_by_exposure[i]
        spec_table.sort(order="SOURCE_ID")
        ext = dm.WFSSMultiSpecModel(spec_table)

        # Set default units from the model schema
        set_schema_units(ext)
        # copy units from the example specmodel, overriding the schema defaults where applicable
        copy_column_units(example_spec, ext)

        # copy metadata
        ext.filename = exposure_counter[exposure_number]["filename"]
        ext.meta.wcs = exposure_counter[exposure_number]["wcs"]
        ext.exposure_number = exposure_number
        ext.meta.exposure.exposure_time = exposure_counter[exposure_number]["exposure_time"]
        ext.meta.exposure.integration_time = exposure_counter[exposure_number]["integration_time"]

        output_x1d.exposures.append(ext)

    # Save the combined results to a file using first input model for metadata
    example_model = results_list[0]
    output_x1d.update(example_model, only="PRIMARY")
    return output_x1d


def wfss_multiexposure_to_multispec(input_model):
    """
    Transform a WFSSMultiExposureSpecModel into a list of MultiSpecModel objects.

    Parameters
    ----------
    input_model : WFSSMultiExposureSpecModel
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
    for exp in input_model.exposures:
        this_exp_list = expand_table(exp)
        source_ids.extend([spec.source_id for spec in this_exp_list])
        exposure_times.extend([exp.meta.exposure.exposure_time for _ in this_exp_list])
        integration_times.extend([exp.meta.exposure.integration_time for _ in this_exp_list])
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
        output_list.append(multispec)

    return output_list


def make_wfss_multicombined(results_list):
    """
    Compile exposure-averaged sources into a single table and save to a file.

    The output c1d product will contain a single spec_table,
    which will end up in extension index 1 of the FITS file.

    Parameters
    ----------
    results_list : list[MultiCombinedSpecModel]
        List of MultiSlitModel objects to be combined into a single c1d file.

    Returns
    -------
    output_c1d : WFSSMultiCombinedSpecModel
        The combined c1d product for WFSS modes.
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

    set_schema_units(output_c1d)
    # copy units from any of the SpecModels (they should all be the same)
    copy_column_units(model.spec[0], output_c1d)

    return output_c1d
