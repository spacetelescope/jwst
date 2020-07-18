#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import wavecorr

__all__ = ["WavecorrStep"]


class WavecorrStep(Step):
    """
    This step applies wavelength offsets to off-center NIRSpec sources.
    """

    spec = """
    """

    reference_file_types = ['wavecorr']

    def process(self, step_input):

        wavecorr_supported_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                                    'NRS_AUTOFLAT']

        # Open the input
        with datamodels.open(step_input) as input_model:

            # Check for valid exposure type
            exp_type = input_model.meta.exposure.type.upper()
            if exp_type not in wavecorr_supported_modes:
                self.log.info(f'Skipping wavecorr correction for EXP_TYPE {exp_type}')
                input_model.meta.cal_step.wavecorr = "SKIPPED"
                return input_model

            # Check for prerequisites
            if hasattr(input_model.meta.cal_step, 'assign_wcs') and input_model.meta.cal_step.assign_wcs == 'SKIPPED':
                self.log.warning("assign_wcs was skipped")
                self.log.warning("Wavecorr step will be skipped")
                input_model.meta.cal_step.wavecorr = 'SKIPPED'
                return input_model

            # Check for existence of WCS
            if isinstance(input_model, datamodels.SlitModel):
                if not (hasattr(input_model.meta, 'wcs') and input_model.meta.wcs is not None):
                    raise AttributeError("Input model does not have a WCS object; assign_wcs should "
                                         "be run before wavecorr.")
            else:
                if not (hasattr(input_model.slits[0].meta, 'wcs') and input_model.slits[0].meta.wcs is not None):
                    raise AttributeError("Input model does not have a WCS object; assign_wcs should "
                                         "be run before wavecorr.")

            # Get the reference file
            reffile = self.get_reference_file(input_model, 'wavecorr')
            self.log.info(f'Using WAVECORR reference file {reffile}')
            if reffile == 'N/A':
                self.log.warning('No WAVECORR reference file found')
                self.log.warning('Wavecorr step will be skipped')
                input_model.meta.cal_step.wavecorr = 'SKIPPED'
                return input_model

            # Apply the correction
            output_model = wavecorr.do_correction(input_model, reffile)

        return output_model
