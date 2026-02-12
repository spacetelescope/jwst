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
        maximum_cores = string(default='1')
        skip = boolean(default=True)
        orders = list(default=None)  # Spectral orders to process, e.g. 1, or 1,2,3
        magnitude_limit = float(default=None) # Isophotal AB magnitude limit for sources to be included in the contamination correction
        wl_oversample = integer(default=2) # oversampling factor for wavelength grid
        max_pixels_per_chunk = integer(default=50000) # max number of pixels to disperse at once
    """  # noqa: E501

    reference_file_types = ["photom", "wavelengthrange"]

    def process(self, input_data):
        """
        Run the WFSS contamination correction step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.MultiSlitModel`
            The input file name or data model containing 2-D cutouts for
            each identified source.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
            A copy of the input data model, with contamination removed.
        """
        output_model = self.prepare_output(input_data)

        # Get the wavelengthrange ref file
        waverange_ref = self.get_reference_file(output_model, "wavelengthrange")
        log.info(f"Using WAVELENGTHRANGE reference file {waverange_ref}")

        # Get the photom ref file
        photom_ref = self.get_reference_file(output_model, "photom")
        log.info(f"Using PHOTOM reference file {photom_ref}")

        # Get requested orders
        orders = [int(o) for o in self.orders] if self.orders else None

        with (
            datamodels.WavelengthrangeModel(waverange_ref) as waverange_model,
            datamodels.open(photom_ref) as photom_model,
        ):
            result, simul, contam, simul_slits = wfss_contam.contam_corr(
                output_model,
                waverange_model,
                photom_model,
                self.maximum_cores,
                orders=orders,
                magnitude_limit=self.magnitude_limit,
                oversample_factor=self.wl_oversample,
                max_pixels_per_chunk=self.max_pixels_per_chunk,
            )
        if simul is None:
            # Input model is returned as result, no intermediate models created
            output_model.meta.cal_step.wfss_contam = "SKIPPED"
            return output_model

        # Save intermediate results, if requested
        if self.save_simulated_image:
            simul_path = self.save_model(simul, suffix="simul", force=True)
            log.info(f'Full-frame simulated grism image saved to "{simul_path}"')
            simul_slits_path = self.save_model(simul_slits, suffix="simul_slits", force=True)
            log.info(f'Simulated slits saved to "{simul_slits_path}"')
        if self.save_contam_images:
            contam_path = self.save_model(contam, suffix="contam", force=True)
            log.info(f'Contamination estimates saved to "{contam_path}"')

        # Close intermediate files
        simul.close()
        simul_slits.close()
        contam.close()

        # If the step succeeded, it created a new output model, so
        # close the input if it was opened here.
        if output_model is not input_data:
            output_model.close()

        return result
