import logging
from functools import partial

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
    fit_threshold = float(default=10.0) # Limiting sigma for fitting splines
    oversample = float(default=1.0) # Use the trace model to oversample the data by this factor
    slope_limit = float(default=0.1) # Slope limit for using splines in oversample
    psf_optimal = boolean(default=False) # Model the target as a simple point source; ignore slope_limit and fit_threshold.
    save_intermediate_results = boolean(default=False)  # Save the full spline model and residuals.
    skip = boolean(default=True) # By default, skip the step.
    output_use_model = boolean(default=True) # Use input filenames in the output models
    """  # noqa: E501

    def process(self, input_data):
        """
        Fit an adaptive trace model to a spectral image.

        For each spectral region in the input (i.e. a slice or a slit), the algorithm
        fits a univariate basis spline to the spatial data within a window of each
        dispersion element. A spectral region may be ignored if its signal-to-noise
        ratio is too low. This decision is controlled by the ``fit_threshold`` parameter:
        lower values will create spline models for more slices.

        After fitting, a further check is performed, using the ``slope_limit`` parameter
        to identify bright, compact source regions for which the spline models are likely
        to be valid. The set of spline models is then evaluated at each pixel within a
        compact source region to create a model of the spectral trace.

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
        source regions use the spline model; faint, diffuse regions use the linear
        interpolation.  The ``slope_limit`` parameter can be modified to impact the decision
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

            # Set up output path name to include the ASN ID if associations are involved
            # TODO: This check is also performed in pixel_replace and outlier_detection.
            #       It should be moved to a shared location instead.
            asn_id = None
            try:
                asn_id = models.asn_table["asn_id"]
            except (AttributeError, KeyError):
                pass
            if asn_id is None:
                asn_id = self.search_attr("asn_id")
            if asn_id is not None:
                _make_output_path = self.search_attr("_make_output_path", parent_first=True)
                self._make_output_path = partial(_make_output_path, asn_id=asn_id)

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
            results = fit_and_oversample(
                model,
                fit_threshold=self.fit_threshold,
                slope_limit=self.slope_limit,
                oversample_factor=self.oversample,
                psf_optimal=self.psf_optimal,
                return_intermediate_models=self.save_intermediate_results,
            )

            model.meta.cal_step.adaptive_trace_model = "COMPLETE"
            if self.save_intermediate_results:
                _, full_spline, used_spline, linear, residual = results
                basepath = model.meta.filename

                outpath = self.make_output_path(basepath=basepath, suffix="spline_full")
                full_spline.save(outpath)
                full_spline.close()
                log.info(f"Saved full spline model in {outpath}")

                outpath = self.make_output_path(basepath=basepath, suffix="spline_used")
                used_spline.save(outpath)
                used_spline.close()
                log.info(f"Saved spline model for compact sources in {outpath}")

                if linear is not None:
                    outpath = self.make_output_path(basepath=basepath, suffix="linear_interp")
                    linear.save(outpath)
                    linear.close()
                    log.info(f"Saved linearly interpolated data in {outpath}")
                else:
                    log.info(
                        f"No linearly interpolated data to save for oversample={self.oversample}"
                    )
                if residual is not None:
                    outpath = self.make_output_path(basepath=basepath, suffix="spline_residual")
                    residual.save(outpath)
                    residual.close()
                    log.info(f"Saved spline residuals in {outpath}")
                else:
                    log.info(f"No spline residuals to save for oversample={self.oversample}")

        return output_model
