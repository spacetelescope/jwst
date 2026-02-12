import logging
from functools import partial

from jwst import datamodels
from jwst.pixel_replace.pixel_replace import PixelReplacement
from jwst.stpipe import Step, record_step_status

__all__ = ["PixelReplaceStep"]

log = logging.getLogger(__name__)


class PixelReplaceStep(Step):
    """Replace flagged bad pixels prior to spectral extraction."""

    class_alias = "pixel_replace"

    spec = """
        algorithm = option("fit_profile", "mingrad", "N/A", default="fit_profile") # Replacement algorithm
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
            output_model.meta.cal_step.pixel_replace = "SKIPPED"
            return output_model

        pars = {
            "algorithm": self.algorithm,
            "n_adjacent_cols": self.n_adjacent_cols,
        }

        # calwebb_spec3 case / ModelContainer
        if isinstance(output_model, datamodels.ModelContainer):
            # Set up output path name to include the ASN ID
            # if associations are involved
            asn_id = None
            try:
                asn_id = output_model.asn_table["asn_id"]
            except (AttributeError, KeyError):
                pass
            if asn_id is None:
                asn_id = self.search_attr("asn_id")
            if asn_id is not None:
                _make_output_path = self.search_attr("_make_output_path", parent_first=True)
                self._make_output_path = partial(_make_output_path, asn_id=asn_id)

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
