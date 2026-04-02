"""Conversion of Stage 2b exposure-based data products to Stage 3 source-based data products."""

import logging
from collections import defaultdict

from stdatamodels.jwst.datamodels import MultiExposureModel
from stdatamodels.properties import merge_tree

from jwst.datamodels import SourceModelContainer

__all__ = ["exp_to_source", "multislit_to_container"]

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
            if slit.source_name is None:
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
