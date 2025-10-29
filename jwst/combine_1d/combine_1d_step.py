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
        input_data : str, ModelContainer, MultiSpecModel, TSOMultiSpecModel, MRSMultiSpecModel
            Input is expected to be an association file name, ModelContainer,
            or multi-spectrum model containing multiple spectra to be combined.
            Individual members of the association or container are expected
            to be multi-spectrum model instances.

        Returns
        -------
        output_spectrum : MultiCombinedSpecModel
            A single combined 1D spectrum.
        """
        with datamodels.open(input_data) as input_model:
            if isinstance(input_model, datamodels.WFSSMultiSpecModel):
                input_list = wfss_multiexposure_to_multispec(input_model)
                if len(input_list) == 1:
                    input_model = input_list[0]
                else:
                    results_list = []
                    for model in input_list:
                        result = combine1d.combine_1d_spectra(
                            model, self.exptime_key, sigma_clip=self.sigma_clip
                        )
                        if not result.meta.cal_step.combine_1d == "SKIPPED":
                            results_list.append(result)
                    if not results_list:
                        log.error("No valid input spectra found in WFSSMultiSpecModel. Skipping.")
                        result = input_model.copy()
                        result.meta.cal_step.combine_1d = "SKIPPED"
                        return result
                    result = make_wfss_multicombined(results_list)
                    result.meta.cal_step.combine_1d = "COMPLETE"
                    return result
            try:
                result = combine1d.combine_1d_spectra(
                    input_model, self.exptime_key, sigma_clip=self.sigma_clip
                )
            except TypeError:
                log.error("Invalid input model for combine_1d; skipping.")
                result = input_model.copy()
                result.meta.cal_step.combine_1d = "SKIPPED"

        return result
