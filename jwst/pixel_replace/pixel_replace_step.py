"""Estimate missing pixel values in spectral data."""

import logging

from jwst import datamodels
from jwst.adaptive_trace_model.adaptive_trace_model_step import AdaptiveTraceModelStep
from jwst.lib.basic_utils import disable_logging
from jwst.pixel_replace.pixel_replace import PixelReplacement
from jwst.stpipe import Step, query_step_status, record_step_status

__all__ = ["PixelReplaceStep"]

log = logging.getLogger(__name__)


class PixelReplaceStep(Step):
    """Replace flagged bad pixels prior to spectral extraction."""

    class_alias = "pixel_replace"

    spec = """
        algorithm = option("fit_profile", "mingrad", "trace_model", default="mingrad") # Replacement algorithm
        n_adjacent_cols = integer(default=3) # Number of adjacent columns to use in profile creation
        skip = boolean(default=True) # Step must be turned on by parameter reference or user
        output_use_model = boolean(default=True) # Use input filenames in the output models
    """  # noqa: E501

    def process(self, input_data):
        """
        Execute the step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel`
            The file name or input datamodel containing spectral data in need of
            pixel replacement.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
            This will be the same as the input model if the step was skipped;
            otherwise, it will be a model containing data arrays with estimated fluxes
            for any bad pixels, now flagged as FLUX_ESTIMATED.
        """
        output_model = self.prepare_output(input_data)

        # If more than one 2d spectrum exists in input, call replacement

        if isinstance(
            output_model,
            datamodels.MultiSlitModel
            | datamodels.SlitModel
            | datamodels.ImageModel
            | datamodels.IFUImageModel
            | datamodels.CubeModel,
        ):
            log.debug(f"Input is a {str(output_model)}.")
        elif isinstance(output_model, datamodels.ModelContainer):
            log.debug("Input is a ModelContainer.")
        else:
            log.error(f"Input is of type {str(output_model)} for which")
            log.error("pixel_replace does not have an algorithm.")
            log.error("Pixel replacement will be skipped.")
            record_step_status(output_model, "pixel_replace", success=False)
            return output_model

        # Set up output path name to include the ASN ID
        # if associations are involved
        self.add_asn_id_to_output_name(output_model)

        # Create a trace model if needed
        if self.algorithm == "trace_model":
            atm_status = query_step_status(output_model, "adaptive_trace_model")
            if atm_status in ("NOT SET", "SKIPPED", None):
                log.info(
                    "The algorithm is 'trace_model' but the adaptive_trace_model "
                    "step has not been completed."
                )
                log.info("Fitting a trace model to the input data")
                with disable_logging(level=logging.INFO):
                    try:
                        output_model = AdaptiveTraceModelStep.call(output_model, oversample=1)
                    except ValueError as err:
                        log.error("Processing failed with ValueError: %s", str(err))

                # Check status again
                atm_status = query_step_status(output_model, "adaptive_trace_model")

            if atm_status != "COMPLETE":
                log.warning(
                    "The algorithm is 'trace_model' but the adaptive_trace_model "
                    "step failed. Defaulting to the 'mingrad' method instead."
                )
                self.algorithm = "mingrad"

        # Parameters to pass
        pars = {
            "algorithm": self.algorithm,
            "n_adjacent_cols": self.n_adjacent_cols,
        }

        # calwebb_spec3 case / ModelContainer
        if isinstance(output_model, datamodels.ModelContainer):
            # Check models to confirm they are the correct type
            for i, model in enumerate(output_model):
                if isinstance(
                    model,
                    datamodels.MultiSlitModel
                    | datamodels.SlitModel
                    | datamodels.ImageModel
                    | datamodels.IFUImageModel
                    | datamodels.CubeModel,
                ):
                    log.debug(f"Input is a {str(model)}.")
                    replacement = PixelReplacement(model, **pars)
                    replacement.replace()
                    success = True
                else:
                    log.error(f"Input is of type {str(model)} for which")
                    log.error("pixel_replace does not have an algorithm.")
                    log.error("Pixel replacement will be skipped.")
                    success = False
                record_step_status(output_model[i], "pixel_replace", success=success)

            return output_model

        # calwebb_spec2 case / single input model
        else:
            # Make copy of input to prevent overwriting
            replacement = PixelReplacement(output_model, **pars)
            replacement.replace()
            record_step_status(output_model, "pixel_replace", success=True)
            return output_model
