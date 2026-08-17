"""Flag pixels affected by open MSA shutters in NIRSpec exposures."""

import logging
from pathlib import Path

from crds import api as crds_api
from stpipe.crds_client import get_context_used, reference_uri_to_cache_path

from jwst.assign_wcs import AssignWcsStep
from jwst.msaflagopen import msaflag_open
from jwst.stpipe import Step

__all__ = ["MSAFlagOpenStep"]

log = logging.getLogger(__name__)


class MSAFlagOpenStep(Step):
    """Flag pixels affected by MSA failed open shutters."""

    class_alias = "msa_flagging"

    spec = """
    """  # noqa: E501

    reference_file_types = ["msaoper"]

    def process(self, input_data):
        """
        Flag data affected by open shutters.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.ImageModel`
            Science data to be corrected.

        Returns
        -------
        `~stdatamodels.jwst.datamodels.ImageModel`
            Science data with DQ array modified.
        """
        # Open the input data model
        output_model = self.prepare_output(input_data)

        self.reference_name = self.get_reference_file(output_model, "msaoper")
        log.info("Using reference file %s", self.reference_name)

        # Check for a valid reference file
        if self.reference_name == "N/A":
            log.warning("No reference file found")
            log.warning("Step will be skipped")
            output_model.meta.cal_step.msa_flagging = "SKIPPED"
            return output_model

        # Get the reference file names for constructing the WCS pipeline
        wcs_reffile_names = create_reference_filename_dictionary(output_model)

        # Do the DQ flagging
        output_model = msaflag_open.do_correction(
            output_model, self.reference_name, wcs_reffile_names
        )

        # set the step status to complete
        output_model.meta.cal_step.msa_flagging = "COMPLETE"

        return output_model


def create_reference_filename_dictionary(input_model):
    """
    Get all relevant WCS reference files.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Input data with WCS assigned.

    Returns
    -------
    ref_files : dict
        Dictionary of reference files.  Keys are CRDS reference types; values
        are file paths or 'N/A'.
    """
    ref_files = {}
    for ref_type in AssignWcsStep.reference_file_types:
        ref_file = getattr(input_model.meta.ref_file, ref_type)
        ref_files[ref_type] = ref_file.name

    # Convert CRDS URIs to local paths, fetching exact recorded references
    # when they are not already present in the local cache.
    for key, ref_file in ref_files.items():
        if isinstance(ref_file, str) and ref_file.startswith("crds://"):
            ref_files[key] = _reference_uri_to_local_path(ref_file, input_model)
    return ref_files


def _reference_uri_to_local_path(reference_uri, input_model):
    """Resolve a CRDS URI to a readable local reference-file path."""
    cache_path = reference_uri_to_cache_path(reference_uri, input_model.crds_observatory)
    if Path(cache_path).exists():
        return cache_path

    context = getattr(input_model.meta.ref_file.crds, "context_used", None)
    if not context:
        context = get_context_used(input_model.crds_observatory)

    basename = reference_uri.removeprefix("crds://")
    log.info("Fetching reference file %s using CRDS context %s", basename, context)
    return crds_api.dump_references(context, [basename])[basename]
