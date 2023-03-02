#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import wfss_contam


__all__ = ["WfssContamStep"]


class WfssContamStep(Step):
    """
    This Step performs contamination correction of WFSS spectra.
    """

    class_alias = "wfss_contam"

    spec = """
        save_simulated_image = boolean(default=False)  # Save full-frame simulated image
        save_contam_images = boolean(default=False)  # Save source contam estimates
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none')
        skip = boolean(default=True)
    """

    reference_file_types = ['photom', 'wavelengthrange']

    def process(self, input_model, *args, **kwargs):

        with datamodels.open(input_model) as dm:

            max_cores = self.maximum_cores

            # Get the wavelengthrange ref file
            waverange_ref = self.get_reference_file(dm, 'wavelengthrange')
            self.log.info(f'Using WAVELENGTHRANGE reference file {waverange_ref}')
            waverange_model = datamodels.WavelengthrangeModel(waverange_ref)

            # Get the photom ref file
            photom_ref = self.get_reference_file(dm, 'photom')
            self.log.info(f'Using PHOTOM reference file {photom_ref}')
            photom_model = datamodels.open(photom_ref)

            result, simul, contam = wfss_contam.contam_corr(dm,
                                                            waverange_model,
                                                            photom_model,
                                                            max_cores)

            # Save intermediate results, if requested
            if self.save_simulated_image:
                simul_path = self.save_model(simul, suffix="simul", force=True)
                self.log.info(f'Full-frame simulated grism image saved to "{simul_path}"')
            if self.save_contam_images:
                contam_path = self.save_model(contam, suffix="contam", force=True)
                self.log.info(f'Contamination estimates saved to "{contam_path}"')

        # Return the corrected data
        return result
