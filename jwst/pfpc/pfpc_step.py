import logging

from crds.core.exceptions import CrdsLookupError
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.pfpc.pfpc import apply_correction, combine_dithers, find_correction, process_exposures
from jwst.stpipe import Step, record_step_status

__all__ = ["PFPCStep"]

log = logging.getLogger(__name__)


class PFPCStep(Step):
    """Apply a point fixed pattern correction (PFPC) to dithered point source spectra."""

    class_alias = "pfpc"

    spec = """
        skip = boolean(default=True) # By default, skip the step.
        suffix = string(default='pfpc') # Set the output suffix.
        output_use_model = boolean(default=True) # Use input filenames in the output models
    """  # noqa: E501

    reference_file_types = ["pfpc"]

    def process(self, input_data):
        """
        Extract spectra and apply a point fixed pattern correction (PFPC).

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel` \
                     or `~jwst.datamodels.container.ModelContainer`
            The input filename, datamodel, or container of 2D spectral data.
            If the input is a container, each contained model is separately
            processed.

        Returns
        -------
        output_model : `~jwst.datamodels.container.ModelContainer`
            Extracted spectra with fixed pattern correction applied. If
            no correction can be applied, an empty container is returned.
        """
        # Read in models. In pipeline context, these may be modified
        # in place to set a 'FAILED' status if the step does not succeed,
        # but they are not returned by the step.
        output_model = self.prepare_output(input_data)
        if isinstance(output_model, ModelContainer):
            models = output_model
        else:
            models = [output_model]

        # Get output file name from self or parent if available
        output_file = self.search_attr("output_file")
        if output_file is None:
            # Set up output path name to include the ASN ID if available
            self.add_asn_id_to_output_name(models)

        # Make an empty container for output: if no products can be made,
        # none are returned
        output_container = ModelContainer()

        # Get the reference file from the first model
        try:
            pfpc_file = self.get_reference_file(models[0], "pfpc")
        except CrdsLookupError:
            pfpc_file = "N/A"
        if pfpc_file == "N/A":
            log.warning("No PFPC reference file found.")
            record_step_status(output_model, "pfpc", success=False)
            return output_container

        with datamodels.open(pfpc_file) as pfpc_model:
            pfpc_table = pfpc_model.pfpc_table

        # Check the first input against the PFPC table and skip processing if no match
        ignore_keys = ["channel"]  # channel may not exactly match before extraction
        first_correction = find_correction(
            models[0], pfpc_table, ignore_keys=ignore_keys, require_one=False
        )
        if first_correction is None:
            log.warning("No matching PFPC correction found for input models.")
            record_step_status(output_model, "pfpc", success=False)
            return output_container

        # TODO - Also check that all input models are point sources and did TA

        # Extract spectra from each exposure
        spec_by_exposure = process_exposures(models, output_file)

        # Apply the PFPC correction at each dither position
        corrected_spec = apply_correction(spec_by_exposure, pfpc_table)

        # Combine spectra across all dithers for each band
        if len(corrected_spec) > 0:
            combined_spec = combine_dithers(corrected_spec)
        else:
            log.warning("No corrected spectra were created.")
            record_step_status(output_model, "pfpc", success=False)
            return output_container

        # Update step status and assemble output
        for spec in combined_spec:
            spec.meta.cal_step.pfpc = "COMPLETE"
            output_container.append(spec)

        # Close the input models if necessary: returned models are newly created
        if output_model is not input_data:
            output_model.close()

        return output_container
