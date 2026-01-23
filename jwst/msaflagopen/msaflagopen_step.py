import logging

from stpipe.crds_client import reference_uri_to_cache_path

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
    input_model : DataModel
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

    # Convert from crds protocol to absolute filenames
    for key in ref_files.keys():
        if ref_files[key].startswith("crds://"):
            ref_files[key] = reference_uri_to_cache_path(
                ref_files[key], input_model.crds_observatory
            )
    return ref_files
