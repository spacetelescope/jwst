#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import emicorr


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

            # Setup outlier detection parameters
            pars = {
                'save_intermediate_results': self.save_intermediate_results,
                'user_supplied_reffile': self.user_supplied_reffile
            }

            # Get the reference file
            save_onthefly_reffile = False
            if self.user_supplied_reffile is None:
                emicorr_ref_file = self.get_reference_file(input_model, 'emicorr')
                if emicorr_ref_file == 'N/A':
                    save_onthefly_reffile = True
            else:
                emicorr_ref_file = self.user_supplied_reffile

            # Read/create the reference file
            self.emicorr_ref_file = emicorr.reffile_waveform(input_model, emicorr_ref_file)
            self.log.info('Using emicorr reference file: %s', self.emicorr_ref_file)

            # Load the reference file
            emicorr_model = datamodels.EmiModel(self.emicorr_ref_file)

            # Do the correction
            output_model = emicorr.do_correction(input_model, emicorr_model, **pars)

            # save the reference file created on-the-fly
            if save_onthefly_reffile:
                reffile_path = emicorr_model.save(suffix=rfile_suffix, force=True)
                self.log.info('On-the-fly reference file saved as: %s', reffile_path)

            emicorr_model.close()

        return output_model
