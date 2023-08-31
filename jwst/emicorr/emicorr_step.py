#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import emicorr
import os


# These are the exposure types to correct for EMI
exp_types2correct = [
    'MIR_4QPM'
]


__all__ = ["EmiCorrStep"]


class EmiCorrStep(Step):
    """
    EmiCorrStep: Apply MIRI EMI correction to a science image.
    """

    class_alias = "emicorr"

    spec = """
        save_intermediate_results = boolean(default=False)
        user_supplied_reffile = string(default=None)  # User-supplied reference file
    """

    reference_file_types = ['emicorr']

    rfile_suffix = 'emi_ref_waves'

    def process(self, input):
        with datamodels.open(input) as input_model:

            # Catch the cases to skip
            exp_type = input_model.meta.exposure.type.upper()
            if exp_type not in exp_types2correct:
                log.info(f'EMI correction not implemented for EXP_TYPE {exp_type}.')
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # Setup outlier detection parameters
            pars = {
                'save_intermediate_results': self.save_intermediate_results,
                'user_supplied_reffile': self.user_supplied_reffile
            }

            # Get the reference file
            save_onthefly_reffile = False
            if self.user_supplied_reffile is None:
                emicorr_ref_filename = self.get_reference_file(input_model, 'emicorr')
                # Create the reference file
                if emicorr_ref_filename == 'N/A':
                    save_onthefly_reffile = True
                    emicorr_ref_filename = Step._make_output_path(self, suffix=rfile_suffix)
                    emicorr.mk_reffile_waveform(input_model, emicorr_ref_filename)
                    self.log.info('Using on-the-fly reference file: %s', emicorr_ref_filename)
                else:
                    self.log.info('Using CRDS reference file: %s', emicorr_ref_filename)

            else:
                emicorr_ref_filename = self.user_supplied_reffile
                log.info('Using user-supplied reference file: '.format(emicorr_ref_filename))

            # Load the reference file
            emicorr_model = datamodels.EmiModel(emicorr_ref_filename)

            # Do the correction
            output_model = emicorr.do_correction(input_model, emicorr_model, **pars)
            output_model.meta.cal_step.emicorr = 'COMPLETE'

            # close and remove the reference file created on-the-fly
            emicorr_model.close()
            if save_onthefly_reffile and not self.save_intermediate_results:
                if os.path.isfile(emicorr_ref_filename):
                    os.remove(emicorr_ref_filename)
                    self.log.debug('On-the-fly reference file deleted: %s', emicorr_ref_filename)

        return output_model
