"""Utility functions for saving WFSS x1d products."""

import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.flat_multispec import (
    copy_column_units,
    determine_vector_and_meta_columns,
    make_empty_recarray,
    populate_recarray,
)


def reorder_wfss(input_model):
    """
    Make spec2 x1d product look like spec3 x1d product so save_wfss_x1d can be used.

    The spec2pipeline is called on a single science exposure at a time,
    and the extracted spectra from the extract_1d step are saved in a MultiSpecModel
    with the .spec attribute containing a list of spectra for each source.
    This is opposite to the spec3 pipeline, which returns a list of MultiSpecModels
    with one MultiSpecModel per source and each spectrum in the .spec attribute
    corresponding to an exposure.
    This function reorders the input model to match the spec3 pipeline output
    so that the save_wfss_x1d function can be used to save the x1d product.

    Parameters
    ----------
    input_model : MultiSpecModel
        Input model containing the extracted spectra from the extract_1d step
        in calwebb_spec2.

    Returns
    -------
    list[MultiSpecModel]
        A list of MultiSpecModels, each with a single spectrum corresponding to one
        extracted source in the exposure.
    """
    results_list = []
    for spec in input_model.spec:
        # Create a new MultiSpecModel for each source
        new_spec = dm.MultiSpecModel()
        new_spec.spec.append(spec)
        results_list.append(new_spec)

    return results_list


def save_wfss_x1d(results_list, output_filename):
    """
    Combine all sources into a single table and save to a file.

    The output x1d product will have one extension per exposure.
    Each extension will contain a single table with one row per source.
    The table size is set by the maximum number of data points for any
    source in the exposure; the other sources will be end-padded with NaNs
    so their shape equals the maximum.

    Parameters
    ----------
    results_list : list[MultiSlitModel]
        List of MultiSlitModel objects to be combined into a single x1d file.
    output_filename : str
        Name of the output x1d file.
    """
    # first loop over both source and exposure to figure out final n_rows, n_exposures, n_sources
    n_rows_by_exposure = []
    exposure_counter = {}
    for model in results_list:  # loop over sources
        for spec in model.spec:  # loop over exposures
            # When this is called from calwebb_spec2, there
            fname = getattr(spec.meta, "filename", None)
            if hasattr(spec.meta, "observation"):
                exp_number = getattr(spec.meta.observation, "exposure_number", None)
            else:
                exp_number = None

            # if this is the first time this exposure has been encountered,
            # create a new dictionary entry for it
            if exp_number not in exposure_counter.keys():
                n_rows = spec.spec_table.shape[0]
                exposure_counter[exp_number] = {"n_rows": n_rows, "n_sources": 1, "filename": fname}
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
        # inner loop over exposures
        for spec in model.spec:
            # ensure data goes to table corresponding to correct exposure based on filename
            exp_num = spec.meta.observation.exposure_number
            exposure_idx = exposure_numbers.index(exp_num)
            fltdata = fltdata_by_exposure[exposure_idx]
            n_rows = n_rows_by_exposure[exposure_idx]
            spec_idx = loop_index_by_exposure[exposure_idx]
            loop_index_by_exposure[exposure_idx] += 1

            # rename some spec metadata to clarify where they came from
            # at the request of the WFSS teams
            name_roots = ["XSTART", "YSTART", "XSTOP", "YSTOP"]
            new_names = ["EXTRACT1D_" + name for name in name_roots]
            old_names = ["EXTRACTION_" + name for name in name_roots]
            for old_name, new_name in zip(old_names, new_names, strict=True):
                setattr(spec, new_name.lower(), getattr(spec, old_name.lower(), None))

            # populate the table with data from the input spectrum
            populate_recarray(
                fltdata[spec_idx],
                spec,
                n_rows,
                all_columns,
                is_vector,
                ignore_columns=["NELEMENTS"],
            )

            # special handling for NELEMENTS because not defined in specmeta schema
            fltdata[spec_idx]["NELEMENTS"] = spec.spec_table.shape[0]

    # Finally, create a new MultiExposureModel to hold the combined data
    # with one MultiSpecModel table per exposure
    output_x1d = dm.WFSSMultiExposureSpecModel()
    example_spec = results_list[0].spec[0]
    for i, exposure_number in enumerate(exposure_numbers):
        # Create a new extension for each exposure
        spec_table = fltdata_by_exposure[i]
        spec_table.sort(order="SOURCE_ID")
        ext = dm.WFSSMultiSpecModel(spec_table)

        # copy units from any of the SpecModels (they should all be the same)
        copy_column_units(example_spec, ext)

        # copy metadata
        fname = exposure_counter[exposure_number]["filename"]
        # ext.meta.filename = fname
        # ext.meta.observation.exposure_number = exposure_number

        output_x1d.exposures.append(ext)
        output_x1d.exposures[-1].filename = fname
        output_x1d.exposures[-1].exposure_number = int(exposure_number)

    # Save the combined results to a file using first input model for metadata
    example_model = results_list[0]
    output_x1d.update(example_model, only="PRIMARY")
    output_x1d.save(output_filename)
    output_x1d.close()
