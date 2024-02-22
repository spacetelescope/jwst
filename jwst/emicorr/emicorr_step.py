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
        user_supplied_reffile = string(default=None)  # ASDF user-supplied reference file
        nints_to_phase = integer(default=None)  # Number of integrations to phase
        nbins = integer(default=None)  # Number of bins in one phased wave
        scale_reference = boolean(default=True)  # If True, the reference wavelength will be scaled to the data's phase amplitude
        skip = boolean(default=True)
    """

    reference_file_types = ['emicorr']

    def process(self, input):
        with datamodels.open(input) as input_model:

            # Catch the cases to skip, i.e. all instruments other than MIRI
            instrument = input_model.meta.instrument.name
            if instrument != 'MIRI':
                self.log.warning('EMI correction not implemented for instrument: {}'.format(instrument))
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # Setup parameters
            pars = {
                'save_intermediate_results': self.save_intermediate_results,
                'user_supplied_reffile': self.user_supplied_reffile,
                'nints_to_phase': self.nints_to_phase,
                'nbins': self.nbins,
                'scale_reference': self.scale_reference
            }

            # Get the reference file
            save_onthefly_reffile, emicorr_ref_filename, emicorr_model = None, None, None
            if self.user_supplied_reffile is None:
                try:
                    emicorr_ref_filename = self.get_reference_file(input_model, 'emicorr')
                    # Create the reference file only if told to save outputs, else correct on-the-fly
                    if emicorr_ref_filename == 'N/A':
                        emicorr_ref_filename = None
                    else:
                        self.log.info('Using CRDS reference file: {}'.format(emicorr_ref_filename))
                        emicorr_model = datamodels.EmiModel(emicorr_ref_filename)
                except Exception:
                    # No reference file in CRDS
                    self.log.info('No CRDS emicorr reference file found. Creating on-the-fly reference file.')

            else:
                self.log.info('Using user-supplied reference file: {}'.format(self.user_supplied_reffile))
                emicorr_model = datamodels.EmiModel(self.user_supplied_reffile)

            # Do the correction
            if self.save_intermediate_results:
                if emicorr_ref_filename is None and self.user_supplied_reffile is None:
                    # get the same full path as input file to save the on-the-fly reference file
                    emicorr_ref_filename = Step._make_output_path(self, suffix='emi_ref_waves')
                    save_onthefly_reffile = emicorr_ref_filename
                else:
                    save_onthefly_reffile = None
            output_model = emicorr.do_correction(input_model, emicorr_model, save_onthefly_reffile, **pars)
            if isinstance(output_model, str) or output_model is None:
                # in this case output_model=subarray_readpatt configuration
                self.log.warning('No correction match for this configuration')
                self.log.warning('Step skipped')
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # close and remove the reference file created on-the-fly
            output_model.meta.cal_step.emicorr = 'COMPLETE'

        return output_model
