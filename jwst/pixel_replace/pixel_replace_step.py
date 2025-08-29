from functools import partial

from jwst import datamodels
from jwst.pixel_replace.pixel_replace import PixelReplacement
from jwst.stpipe import Step, record_step_status

__all__ = ["PixelReplaceStep"]


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
        input_data : datamodel, str
            The input datamodel or filename containing
            spectral data in need of pixel replacement.

        Returns
        -------
        JWST DataModel
            This will be `input` if the step was skipped; otherwise,
            it will be a model containing data arrays with estimated fluxes
            for any bad pixels, now flagged as TO-BE-DETERMINED (DQ bit 7?).
        """
        with datamodels.open(input_data) as input_model:
            # If more than one 2d spectrum exists in input, call replacement

            if isinstance(
                input_model,
                datamodels.MultiSlitModel
                | datamodels.SlitModel
                | datamodels.ImageModel
                | datamodels.IFUImageModel
                | datamodels.CubeModel,
            ):
                self.log.debug(f"Input is a {str(input_model)}.")
            elif isinstance(input_model, datamodels.ModelContainer):
                self.log.debug("Input is a ModelContainer.")
            else:
                self.log.error(f"Input is of type {str(input_model)} for which")
                self.log.error("pixel_replace does not have an algorithm.")
                self.log.error("Pixel replacement will be skipped.")
                result = input_model.copy()
                result.meta.cal_step.pixel_replace = "SKIPPED"
                return result

            pars = {
                "algorithm": self.algorithm,
                "n_adjacent_cols": self.n_adjacent_cols,
            }

            # calwebb_spec3 case / ModelContainer
            if isinstance(input_model, datamodels.ModelContainer):
                output_model = input_model

                # Set up output path name to include the ASN ID
                # if associations are involved
                asn_id = None
                try:
                    asn_id = input_model.asn_table["asn_id"]
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
                        self.log.debug(f"Input is a {str(model)}.")
                        replacement = PixelReplacement(model, **pars)
                        replacement.replace()
                        output_model[i] = replacement.output
                        success = True
                    else:
                        self.log.error(f"Input is of type {str(model)} for which")
                        self.log.error("pixel_replace does not have an algorithm.")
                        self.log.error("Pixel replacement will be skipped.")
                        output_model[i] = model.copy()
                        success = False
                    record_step_status(output_model[i], "pixel_replace", success=success)

                return output_model

            # calwebb_spec2 case / single input model
            else:
                # Make copy of input to prevent overwriting
                result = input_model.copy()
                replacement = PixelReplacement(result, **pars)
                replacement.replace()
                record_step_status(replacement.output, "pixel_replace", success=True)
                result = replacement.output
                return result
