import logging

from stdatamodels.jwst import datamodels

from jwst.combine_1d import combine1d
from jwst.datamodels.utils.wfss_multispec import (
    make_wfss_multicombined,
    wfss_multiexposure_to_multispec,
)
from jwst.stpipe import Step

__all__ = ["Combine1dStep"]

log = logging.getLogger(__name__)


class Combine1dStep(Step):
    """Combine 1D spectra."""

    class_alias = "combine_1d"

    spec = """
    exptime_key = string(default="exposure_time") # Metadata key to use for weighting
    sigma_clip = float(default=None) # Factor for clipping outliers
    """  # noqa: E501

    def process(self, input_data):
        """
        Combine the input data.

        Parameters
        ----------
        input_data : str, `~jwst.datamodels.container.ModelContainer`, \
                     `~stdatamodels.jwst.datamodels.MultiSpecModel`, \
                     `~stdatamodels.jwst.datamodels.TSOMultiSpecModel`, or \
                     `~stdatamodels.jwst.datamodels.MRSMultiSpecModel`
            Input is expected to be an association file name, ModelContainer,
            or multi-spectrum model containing multiple spectra to be combined.
            Individual members of the association or container are expected
            to be multi-spectrum model instances.

        Returns
        -------
        output_spectrum : MultiCombinedSpecModel
            A single combined 1D spectrum.
        """
        output_model = self.prepare_output(input_data)

        if isinstance(output_model, datamodels.WFSSMultiSpecModel):
            input_list = wfss_multiexposure_to_multispec(output_model)
            if len(input_list) == 1:
                # Single input: will be combined below
                output_model = input_list[0]
            else:
                # Multiple inputs: combine in a loop, then reconstitute the final model
                results_list = []
                for model in input_list:
                    result = combine1d.combine_1d_spectra(
                        model, self.exptime_key, sigma_clip=self.sigma_clip
                    )
                    if not result.meta.cal_step.combine_1d == "SKIPPED":
                        results_list.append(result)
                if not results_list:
                    log.error("No valid input spectra found in WFSSMultiSpecModel. Skipping.")
                    output_model.meta.cal_step.combine_1d = "SKIPPED"
                    return output_model
                result = make_wfss_multicombined(results_list)
                result.meta.cal_step.combine_1d = "COMPLETE"

                # Close any input models opened here before returning
                if output_model is not input_data:
                    output_model.close()

                return result

        try:
            result = combine1d.combine_1d_spectra(
                output_model, self.exptime_key, sigma_clip=self.sigma_clip
            )
        except TypeError:
            log.error("Invalid input model for combine_1d; skipping.")
            output_model.meta.cal_step.combine_1d = "SKIPPED"
            return output_model

        # The result is a new model: close any input models opened here
        if output_model is not input_data:
            output_model.close()

        return result
