from stdatamodels.jwst import datamodels

from stpipe.crds_client import reference_uri_to_cache_path

from jwst.stpipe import Step
from jwst.assign_wcs import AssignWcsStep
from jwst.msaflagopen import msaflag_open

__all__ = ["MSAFlagOpenStep"]


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
        input_data : DataModel or str
            Science data to be corrected.

        Returns
        -------
        DataModel
            Science data with DQ array modified.
        """
        # Open the input data model
        with datamodels.open(input_data) as input_model:
            self.reference_name = self.get_reference_file(input_model, "msaoper")
            self.log.info("Using reference file %s", self.reference_name)

            # Check for a valid reference file
            if self.reference_name == "N/A":
                self.log.warning("No reference file found")
                self.log.warning("Step will be skipped")
                result = input_model.copy()
                result.meta.cal_step.msa_flagging = "SKIPPED"
                return result

            # Get the reference file names for constructing the WCS pipeline
            wcs_reffile_names = create_reference_filename_dictionary(input_model)

            # Do the DQ flagging
            result = msaflag_open.do_correction(input_model, self.reference_name, wcs_reffile_names)

            # set the step status to complete
            result.meta.cal_step.msa_flagging = "COMPLETE"

        return result


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
