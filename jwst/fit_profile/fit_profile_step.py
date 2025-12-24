import logging

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.fit_profile.fit_profile import fit_and_oversample
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
    skip = boolean(default=True) # By default, skip the step.
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
        if isinstance(output_model, ModelContainer):
            models = output_model
        else:
            models = [output_model]

        for model in models:
            if not isinstance(model, datamodels.IFUImageModel):
                log.warning("The fit_profile step is only implemented for IFU data.")
                log.warning("Skipping processing for datamodel type %s.", str(output_model))
                model.meta.cal_step.fit_profile = "SKIPPED"
                continue

            # Update the model in place
            log.info("Fitting profile for %s", model.meta.filename)
            fit_and_oversample(
                model,
                threshsig=self.threshsig,
                slopelim=self.slopelim,
                oversample_factor=self.oversample,
            )

            model.meta.cal_step.fit_profile = "COMPLETE"

        return output_model
