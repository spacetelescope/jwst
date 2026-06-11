"""Conversion of Stage 2b exposure-based data products to Stage 3 source-based data products."""

import logging
from collections import defaultdict

import numpy as np
import stdatamodels.jwst.datamodels as dm
from stdatamodels.jwst.datamodels import MultiExposureModel
from stdatamodels.properties import merge_tree

from jwst.datamodels import SourceModelContainer
from jwst.datamodels.utils.flat_multispec import copy_column_units, copy_spec_metadata

__all__ = ["exp_to_source", "multislit_to_container", "wfss_multispec_to_source"]

log = logging.getLogger(__name__)


def exp_to_source(inputs):
    """
    Reformat exposure-based MSA data to source-based.

    Parameters
    ----------
    inputs : list of `~stdatamodels.jwst.datamodels.MultiSlitModel`
        List of `~stdatamodels.jwst.datamodels.MultiSlitModel` instances to reformat.

    Returns
    -------
    multiexposures : dict
        Returns a dictionary of `~stdatamodels.jwst.datamodels.MultiExposureModel`
        instances wherein each
        instance contains slits belonging to the same source.
        The key is the ``source_id`` of each source.
    """
    result = defaultdict(MultiExposureModel)

    for exposure in inputs:
        log.info(f"Reorganizing data from exposure {exposure.meta.filename}")

        for slit in exposure.slits:
            if slit.source_name is None or str(slit.source_name).strip() == "":
                # All MultiSlit data other than NIRSpec MOS get sorted by
                # source_id (source_name is not populated)
                key = slit.source_id
            else:
                # NIRSpec MOS slits get sorted by source_name
                key = slit.source_name
            log.debug(f"Copying source {key}")
            result_slit = result[str(key)]
            result_slit.exposures.append(slit)

            # store values for later use (after merge_tree)
            # these values are incorrectly getting overwritten by
            # the top model.
            slit_bunit = slit.meta.bunit_data
            slit_bunit_err = slit.meta.bunit_err
            slit_model = slit.meta.model_type
            slit_wcsinfo = slit.meta.wcsinfo.instance
            slit_exptype = None
            if hasattr(slit.meta, "exposure"):
                if hasattr(slit.meta.exposure, "type"):
                    slit_exptype = slit.meta.exposure.type

            # Before merge_tree the slits have a model_type of SlitModel.
            # After merge_tree it is overwritten with MultiSlitModel.
            # store the model type to undo overwriting of modeltype.
            merge_tree(result_slit.exposures[-1].meta.instance, exposure.meta.instance)

            result_slit.exposures[-1].meta.bunit_data = slit_bunit
            result_slit.exposures[-1].meta.bunit_err = slit_bunit_err
            result_slit.exposures[-1].meta.model_type = slit_model
            result_slit.exposures[-1].meta.wcsinfo = slit_wcsinfo

            # make sure top-level exposure type matches slit exposure type
            # (necessary for NIRSpec fixed slits defined as part of an MSA file)
            if slit_exptype is not None:
                result_slit.exposures[-1].meta.exposure.type = slit_exptype
                result_slit.meta.exposure.type = slit_exptype
                log.debug(f"Input exposure type: {exposure.meta.exposure.type}")
                log.debug(f"Output exposure type: {result_slit.meta.exposure.type}")

            if result_slit.meta.instrument.name is None:
                result_slit.update(exposure)
            result_slit.meta.filename = None  # Resulting merged data doesn't come from one file

        exposure.close()

    # Turn off the default factory
    result.default_factory = None

    return result


def multislit_to_container(inputs):
    """
    Reformat exposure-based MSA data to source-based containers.

    Parameters
    ----------
    inputs : list or `~jwst.datamodels.container.ModelContainer` of \
             `~stdatamodels.jwst.datamodels.MultiSlitModel`
        List of `~stdatamodels.jwst.datamodels.MultiSlitModel`
        instances to reformat

    Returns
    -------
    containers : dict
        Returns a dictionary of `~jwst.datamodels.container.ModelContainer`
        instances wherein each instance contains
        `~stdatamodels.jwst.datamodels.ImageModel` of slits belonging
        to the same source.
        The key is the ID of each slit, i.e., ``source_id``.
    """
    containers = exp_to_source(inputs)
    for container_id in containers:
        containers[container_id] = SourceModelContainer(containers[container_id])

    return containers


def _expand_wfss_table(spec):
    """
    Expand a table of spectra into a list of WFSSSpecModel objects.

    This is like :func:`expand_table` but for ``calwebb_spec3``
    where we need to break it down per source.

    Parameters
    ----------
    spec : `~stdatamodels.jwst.datamodels.WFSSSpecModel`
        WFSS model containing a spec_table to expand into multiple spectra

    Returns
    -------
    list[WFSSSpecModel]
        A list of `~stdatamodels.jwst.datamodels.WFSSSpecModel` objects,
        one for each spectrum in the input spec_table.
    """
    new_spec_list = []
    n_spectra = len(spec.spec_table)
    for i in range(n_spectra):
        # initialize a new SpecModel
        spec_row = spec.spec_table[i]
        n_elements = int(spec_row["N_ALONGDISP"])
        new_spec = dm.WFSSSpecModel()
        data_type = new_spec.schema["properties"]["spec_table"]["datatype"]
        columns_to_copy = np.array([col["name"] for col in data_type])

        # Copy over the vector columns from input spec_table to output spec_table
        spec_table = np.empty(n_elements, dtype=new_spec.get_dtype("spec_table"))
        for col_name in columns_to_copy:
            if isinstance(spec_row[col_name], np.ndarray) and spec_row[col_name].ndim == 1:
                spec_table[col_name] = spec_row[col_name][:n_elements]
            else:
                spec_table[col_name] = spec_row[col_name]
        new_spec.spec_table = spec_table
        new_spec.filename = spec.filename
        new_spec.group_id = spec.group_id
        new_spec.dispersion_direction = spec.dispersion_direction
        new_spec.spectral_order = spec.spectral_order
        new_spec.exposure_time = spec.exposure_time
        new_spec.integration_time = spec.integration_time
        new_spec.s_region = spec.s_region
        copy_spec_metadata(spec, new_spec)
        copy_column_units(spec, new_spec)

        source_ids = set(spec_table["SOURCE_ID"])
        if len(source_ids) != 1:
            raise ValueError(f"SOURCE_ID not unique: {source_ids}")
        new_spec.source_id = spec_table["SOURCE_ID"][0]

        new_spec_list.append(new_spec)

    return new_spec_list


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
            this_exp_list = _expand_wfss_table(exp)
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
