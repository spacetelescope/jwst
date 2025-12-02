import logging

from stdatamodels.jwst import datamodels

from jwst.fit_profile.fit_profile import fit_and_oversample_ifu
from jwst.stpipe import Step

__all__ = ["FitProfileStep"]

log = logging.getLogger(__name__)


class FitProfileStep(Step):
    """Fit a spatial profile to a spectral image."""

    class_alias = "fit_profile"

    spec = """
    threshsig = float(default=10) # Limiting sigma for fitting splines
    slopelim = float(default=0.1) # Slope limit for using splines in oversample
    oversample = float(default=1.0) # Use the profile fit to oversample the data by this factor
    # psfoptimal = boolean(default=False)
    """  # noqa: E501

    def process(self, input_data):
        """
        Execute the step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel`
            The input datamodel or filename containing spectral data.

        Returns
        -------
        `~stdatamodels.jwst.datamodels.JwstDataModel`
            The input model, updated with the fit profile.
        """
        output_model = self.prepare_output(input_data)

        if not isinstance(output_model, datamodels.IFUImageModel):
            log.warning("The fit_profile step is only implemented for IFU data.")
            log.warning("Skipping processing for datamodel type %s.", str(output_model))
            output_model.meta.cal_step.fit_profile = "SKIPPED"
            return output_model

        output_model = fit_and_oversample_ifu(
            output_model,
            threshsig=self.threshsig,
            slopelim=self.slopelim,
            oversample_factor=self.oversample,
        )

        output_model.meta.cal_step.fit_profile = "COMPLETE"
        return output_model
