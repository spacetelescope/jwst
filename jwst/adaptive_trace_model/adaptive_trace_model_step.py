import logging

from stdatamodels.jwst import datamodels

from jwst.adaptive_trace_model.trace_model import fit_and_oversample
from jwst.datamodels import ModelContainer
from jwst.stpipe import Step

__all__ = ["AdaptiveTraceModelStep"]

log = logging.getLogger(__name__)


class AdaptiveTraceModelStep(Step):
    """Fit an adaptive trace model to a spectral image."""

    class_alias = "adaptive_trace_model"

    spec = """
    threshsig = float(default=10) # Limiting sigma for fitting splines
    oversample = float(default=1.0) # Use the trace model to oversample the data by this factor
    slopelim = float(default=0.1) # Slope limit for using splines in oversample
    skip = boolean(default=True) # By default, skip the step.
    """  # noqa: E501

    def process(self, input_data):
        """
        Fit an adaptive trace model to a spectral image.

        For each spectral region in the input (i.e. a slice or a slit), the algorithm
        fits a univariate basis spline to the spatial data within a window of each
        dispersion element. A spectral region may be ignored if its signal-to-noise
        ratio is too low. This decision is controlled by the ``threshsig`` parameter:
        lower values will create spline models for more slices.

        After fitting, the set of spline models is then evaluated at each pixel to
        create a model of the spectral trace at all valid data points.

        If no oversample factor is specified, the only change to the datamodel
        is to attach the resulting spectral trace image, in the ``trace_model``
        attribute.

        If an oversample factor is specified, the input flux, error, and variance
        images are replaced with interpolated data scaled to a new pixel grid. The DQ
        and wavelength images are also oversampled and the trace model image is
        attached, but any additional images (e.g. pathloss corrections) are not
        propagated.

        To perform this oversampling, the trace model is used in combination with a
        linear interpolation at each dispersion element.  Which model is used at each
        pixel depends on some heuristic decisions derived from the data: bright, compact
        source regions generally use the spline model; faint, diffuse regions use the linear
        interpolation.  The ``slopelim`` parameter can be modified to impact the decision
        on whether sources are bright and compact enough to use the spline models. Lower
        values will use the spline model for fainter sources.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.IFUImageModel` \
                     or `~jwst.datamodels.container.ModelContainer`
            The input datamodel, container, or filename containing spectral data.
            If the input is a container, each contained model is separately
            processed.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.IFUImageModel` \
                       or `~jwst.datamodels.container.ModelContainer`
            The input model, updated with the fit trace model.
        """
        output_model = self.prepare_output(input_data)
        if isinstance(output_model, ModelContainer):
            models = output_model
        else:
            models = [output_model]

        for model in models:
            if not isinstance(model, datamodels.IFUImageModel):
                log.warning("The adaptive_trace_model step is only implemented for IFU data.")
                log.warning("Skipping processing for datamodel type %s.", str(output_model))
                model.meta.cal_step.adaptive_trace_model = "SKIPPED"
                continue

            # Update the model in place
            log.info("Fitting trace model for %s", model.meta.filename)
            fit_and_oversample(
                model,
                threshsig=self.threshsig,
                slopelim=self.slopelim,
                oversample_factor=self.oversample,
            )

            model.meta.cal_step.adaptive_trace_model = "COMPLETE"

        return output_model
