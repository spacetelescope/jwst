import logging

from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.wfss_contam import wfss_contam

__all__ = ["WfssContamStep"]

log = logging.getLogger(__name__)


class WfssContamStep(Step):
    """Perform contamination correction of WFSS spectra."""

    class_alias = "wfss_contam"

    spec = """
        save_simulated_image = boolean(default=False)  # Save full-frame simulated image
        save_contam_images = boolean(default=False)  # Save source contam estimates
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none')
        skip = boolean(default=True)
        brightest_n = integer(default=None)
        wl_oversample = integer(default=2) # oversampling factor for wavelength grid
        max_pixels_per_chunk = integer(default=50000) # max number of pixels to disperse at once
    """  # noqa: E501

    reference_file_types = ["photom", "wavelengthrange"]

    def process(self, input_model):
        """
        Run the WFSS contamination correction step.

        Parameters
        ----------
        input_model : `~jwst.datamodels.MultiSlitModel`
            The input data model containing 2-D cutouts for each identified source.

        Returns
        -------
        output_model : `~jwst.datamodels.MultiSlitModel`
            A copy of the input_model with contamination removed
        """
        with datamodels.open(input_model) as dm:
            # Get the wavelengthrange ref file
            waverange_ref = self.get_reference_file(dm, "wavelengthrange")
            log.info(f"Using WAVELENGTHRANGE reference file {waverange_ref}")
            waverange_model = datamodels.WavelengthrangeModel(waverange_ref)

            # Get the photom ref file
            photom_ref = self.get_reference_file(dm, "photom")
            log.info(f"Using PHOTOM reference file {photom_ref}")
            photom_model = datamodels.open(photom_ref)

            result, simul, contam, simul_slits = wfss_contam.contam_corr(
                dm,
                waverange_model,
                photom_model,
                self.maximum_cores,
                brightest_n=self.brightest_n,
                oversample_factor=self.wl_oversample,
                max_pixels_per_chunk=self.max_pixels_per_chunk,
            )

            # Save intermediate results, if requested
            if self.save_simulated_image:
                simul_path = self.save_model(simul, suffix="simul", force=True)
                log.info(f'Full-frame simulated grism image saved to "{simul_path}"')
                simul_slits_path = self.save_model(simul_slits, suffix="simul_slits", force=True)
                log.info(f'Simulated slits saved to "{simul_slits_path}"')
            if self.save_contam_images:
                contam_path = self.save_model(contam, suffix="contam", force=True)
                log.info(f'Contamination estimates saved to "{contam_path}"')

        return result
