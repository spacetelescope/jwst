"""Utilities for manipulating WFSS multi-spectral data."""

from collections import defaultdict

import numpy as np
import stdatamodels.jwst.datamodels as dm

from jwst.datamodels.utils.flat_multispec import (
    copy_column_units,
    copy_spec_metadata,
    determine_vector_and_meta_columns,
    expand_table,
    expand_wfss_table,
    make_empty_recarray,
    populate_recarray,
    set_schema_units,
)

__all__ = [
    "make_wfss_multiexposure",
    "make_wfss_multiexposure_spec3",
    "wfss_multiexposure_to_multispec",
    "make_wfss_multicombined",
    "wfss_multispec_to_source",
]


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
    input_list : `~stdatamodels.jwst.datamodels.MultiSpecModel` or list[MultiSpecModel]
        List of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects to be combined.

    Returns
    -------
    output_x1d : `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
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
    all_source_ids = set()
    # for calwebb_spec3 the outer loop is over sources and the inner loop is over exposures
    # for calwebb_spec2 it is the opposite, but outer loop (exposures) should have just one element
    for model in results_list:
        for spec in model.spec:
            fname = getattr(spec.meta, "filename", None)
            exp_number = getattr(spec.meta, "group_id", None)

            # if this is the first time this exposure has been encountered,
            # create a new dictionary entry for it
            n_rows = spec.spec_table.shape[0]
            if exp_number not in exposure_counter:
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
                    exposure_counter[exp_number]["n_rows"], n_rows
                )

            all_source_ids.add(spec.source_id)

    all_source_ids = sorted(all_source_ids)
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


def make_wfss_multiexposure_spec3(input_list):
    """
    Compile a list of extracted sources into a single binary table for calwebb_spec3.

    The output model will contain one binary table per exposure,
    with each table containing all sources extracted from that exposure
    (one row per source). The number of elements in each table row
    will be the same across all exposures, with NaNs used to pad
    shorter rows to match the longest row in the exposure.

    Parameters
    ----------
    input_list : list[MultiSpecModel]
        List of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects to be combined.

    Returns
    -------
    output_x1d : `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
        The extract_1d product for WFSS modes.
    """
    # First loop over source and exposure to figure out final parameters
    specs_db = defaultdict()  # (group_id, source_id): spec
    all_source_ids = set()
    exposure_nrows = defaultdict(int)
    # For calwebb_spec3 the outer loop is over sources and the inner loop is over exposures
    for model in input_list:
        for spec in model.spec:
            all_source_ids.add(spec.source_id)
            exposure_nrows[spec.group_id] = max(
                exposure_nrows[spec.group_id], spec.spec_table.shape[0]
            )
            specs_db[(spec.group_id, spec.source_id)] = spec

    all_source_ids = sorted(all_source_ids)
    n_sources = len(all_source_ids)

    all_columns = [
        (x["name"], x["datatype"])
        for x in dm.WFSSSpecModel().schema["properties"]["spec_table"]["datatype"]
    ]
    all_columns[20] = ("SOURCE_TYPE", "U20")  # schema incompatible with numpy
    vec_cols = (
        "WAVELENGTH",
        "FLUX",
        "FLUX_ERROR",
        "FLUX_VAR_POISSON",
        "FLUX_VAR_RNOISE",
        "FLUX_VAR_FLAT",
        "SURF_BRIGHT",
        "SB_ERROR",
        "SB_VAR_POISSON",
        "SB_VAR_RNOISE",
        "SB_VAR_FLAT",
        "DQ",
        "BACKGROUND",
        "BKGD_ERROR",
        "BKGD_VAR_POISSON",
        "BKGD_VAR_RNOISE",
        "BKGD_VAR_FLAT",
        "NPIXELS",
    )
    is_vector = [True if col[0] in vec_cols else False for col in all_columns]
    defaults = dm.WFSSSpecModel().schema["properties"]["spec_table"]["default"]

    # Finally, create a new WFSSMultiSpecModel to hold the combined data
    # with one WFSSMultiSpecModel table per exposure
    output_x1d = dm.WFSSMultiSpecModel()
    # WCS is needed to combine S_REGION in calwebb_spec3
    if (not getattr(output_x1d.meta, "wcs", None)) and hasattr(input_list[0].meta, "wcs"):
        output_x1d.meta.wcs = input_list[0].meta.wcs
    for exposure_number in sorted(exposure_nrows.keys()):
        n_rows = exposure_nrows[exposure_number]
        spec_table = make_empty_recarray(
            n_rows, n_sources, all_columns, is_vector, defaults=defaults
        )
        spec_table["SOURCE_ID"] = all_source_ids
        first_loop = True

        for src_id in all_source_ids:
            key = (exposure_number, src_id)
            if key not in specs_db:
                continue
            spec = specs_db[key]

            # ensure data goes to the correct source
            spec_idx = np.where(spec_table["SOURCE_ID"] == src_id)[0][0]

            # populate the table with data from the input spectrum
            populate_recarray(
                spec_table[spec_idx],
                spec,
                all_columns,
                is_vector,
                ignore_columns=["SOURCE_ID", "N_ALONGDISP"],
            )

            # special handling for N_ALONGDISP because not defined in specmeta schema
            spec_table[spec_idx]["N_ALONGDISP"] = spec.spec_table.shape[0]

            if first_loop:
                example_spec = spec

        ext = dm.WFSSSpecModel(spec_table)
        # Set default units from the model schema
        set_schema_units(ext)
        # copy units from the example specmodel, overriding the schema defaults where applicable
        copy_column_units(example_spec, ext)
        # copy metadata
        ext.filename = example_spec.filename
        ext.group_id = exposure_number
        ext.dispersion_direction = example_spec.dispersion_direction
        ext.spectral_order = example_spec.spectral_order
        ext.exposure_time = example_spec.exposure_time
        ext.integration_time = example_spec.integration_time
        ext.s_region = example_spec.s_region
        output_x1d.spec.append(ext)

    output_x1d.update(input_list[0], only="PRIMARY")
    return output_x1d


def wfss_multiexposure_to_multispec(input_model):
    """
    Transform a `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel` into
    a list of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
        Input model to be reorganized.

    Returns
    -------
    output_list : list[MultiSpecModel]
        List of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects,
        one for each source ID in the input model.
    """  # noqa: D205  # numpydoc ignore=SS06
    # first extract all spectra as SpecModels in a flat list
    spec_list = []
    source_ids = []
    first_loop = True
    exposure_time = 0
    integration_time = 0
    for exp in input_model.spec:
        this_exp_list = expand_table(exp)
        source_ids.extend([spec.source_id for spec in this_exp_list])
        if first_loop:
            exposure_time = exp.exposure_time
            integration_time = exp.integration_time
            first_loop = False
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
        multispec.meta.exposure.exposure_time = exposure_time
        multispec.meta.exposure.integration_time = integration_time
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
        List of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects
        to be combined.

    Returns
    -------
    output_c1d : `~stdatamodels.jwst.datamodels.WFSSMultiCombinedSpecModel`
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


def wfss_multispec_to_source(inputs):
    """
    Reformat exposure-based WFSS data to source-based.

    Parameters
    ----------
    inputs : list of `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
        List of model instances to reformat.

    Returns
    -------
    output_list : list[MultiSpecModel]
        List of `~stdatamodels.jwst.datamodels.MultiSpecModel` objects,
        one for each source ID in the input model.
    """
    # first extract all spectra as SpecModels in a flat list
    spec_list = []
    source_ids = []
    first_loop = True
    exposure_time = 0
    integration_time = 0
    for input_model in inputs:
        for exp in input_model.spec:
            this_exp_list = expand_wfss_table(exp)
            for spec in this_exp_list:
                if first_loop:
                    exposure_time = spec.exposure_time
                    integration_time = spec.integration_time
                    first_loop = False
                spec.source_type = spec.spec_table["SOURCE_TYPE"][0]
                spec_list.append(spec)
                source_ids.append(spec.source_id)

    # organize by unique source id such that there is one MultiSpecModel per source
    # with all exposures for that source in that model's model.spec
    unique_source_ids = sorted(set(source_ids))
    source_ids = np.array(source_ids)
    spec_list = np.array(spec_list)
    output_list = []
    for source_id in unique_source_ids:
        multispec = dm.MultiSpecModel()
        # BUG: currently there is no infrastructure for handling per-exposure weights
        # This is also a problem on main
        multispec.meta.exposure.exposure_time = exposure_time
        multispec.meta.exposure.integration_time = integration_time
        spec_this_id = spec_list[source_ids == source_id]
        multispec.spec.extend(spec_this_id)
        multispec.update(input_model, only="PRIMARY")
        output_list.append(multispec)

    return output_list
