#! /usr/bin/env python



### TESTING my local datamodels branch
import os, sys
datamdl_branch = os.path.dirname('/Users/pena/Documents/SCSB/stdatamodels/stdatamodels')
print(datamdl_branch)
sys.path.insert(1, datamdl_branch)
print('\n *** tesing new datamodels branch *** \n')
###


from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import emicorr
import traceback


# These are the exposure types to correct for EMI
exp_types2correct = [
    'MIR_4QPM',
    'MIR_IMAGE',
    'MIR_LRS-SLITLESS',
    'MIR_LYOT'
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

    def process(self, input):
        with datamodels.open(input) as input_model:

            # Catch the cases to skip
            exp_type = input_model.meta.exposure.type.upper()
            if exp_type not in exp_types2correct:
                self.log.info('EMI correction not implemented for EXP_TYPE: {}'.format(exp_type))
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # Setup outlier detection parameters
            pars = {
                'save_intermediate_results': self.save_intermediate_results,
                'user_supplied_reffile': self.user_supplied_reffile
            }

            # Get the reference file
            save_onthefly_reffile = None
            if self.user_supplied_reffile is None:
                try:
                    emicorr_ref_filename = self.get_reference_file(input_model, 'emicorr')
                    # Create the reference file only of told to save outputs, else correct on-the-fly
                    if emicorr_ref_filename == 'N/A':
                        emicorr_model = None
                        if self.save_intermediate_results:
                            emicorr_ref_filename = Step._make_output_path(self, suffix='emi_ref_')
                            save_onthefly_reffile = emicorr_ref_filename
                    else:
                        self.log.info('Using CRDS reference file: {}'.format(emicorr_ref_filename))
                        emicorr_model = datamodels.EmiModel(emicorr_ref_filename)
                except Exception:
                    # CRDS rules not there yet
                    #print(traceback.print_exc())
                    self.log.info('CRDS rules for emicorr reference file not implemented yet. Creating on-the-fly reference file.')
                    emicorr_model = None

            else:
                self.log.info('Using user-supplied reference file: {}'.format(self.user_supplied_reffile))
                emicorr_model = datamodels.EmiModel(self.user_supplied_reffile)

            # Do the correction
            #save_onthefly_reffile = Step._make_output_path(self, suffix='emi_ref_')
            output_model = emicorr.do_correction(input_model, emicorr_model, save_onthefly_reffile, **pars)
            if isinstance(output_model, str):
                # in this case output_model=subarray_readpatt configuration
                self.log.info('No correction match for this configuration: {}'.format(output_model))
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # close and remove the reference file created on-the-fly
            output_model.meta.cal_step.emicorr = 'COMPLETE'

        return output_model
