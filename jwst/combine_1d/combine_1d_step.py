import logging

from stdatamodels.jwst import datamodels

from jwst.combine_1d import combine1d
from jwst.datamodels.utils.wfss_multispec import (
    make_wfss_multicombined,
)
from jwst.stpipe import Step, record_step_status

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
                     `~stdatamodels.jwst.datamodels.TSOMultiSpecModel`, \
                     `~stdatamodels.jwst.datamodels.MRSMultiSpecModel`, or \
                     `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
            Input is expected to be an association file name, ModelContainer,
            or multi-spectrum model containing multiple spectra to be combined.
            Individual members of the association or container are expected
            to be multi-spectrum model instances.

        Returns
        -------
        output_spectrum : `~stdatamodels.jwst.datamodels.MultiCombinedSpecModel`
            A single combined 1D spectrum.
        """
        output_model = self.prepare_output(input_data)

        try:
            result = combine1d.combine_1d_spectra(
                output_model, self.exptime_key, sigma_clip=self.sigma_clip
            )

            # FIXME: Absorb this into combine1d natively.
            if isinstance(output_model, datamodels.WFSSMultiSpecModel):
                result = make_wfss_multicombined([result])
                result.meta.cal_step.combine_1d = "COMPLETE"

        except TypeError as err:
            log.error("%s; skipping.", str(err))
            record_step_status(output_model, "combine_1d", status="SKIPPED")
            return output_model

        # The result is a new model: close any input models opened here
        if output_model is not input_data:
            output_model.close()

        return result
