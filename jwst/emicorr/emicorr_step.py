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
        skip = boolean(default=True)  # Skip the step
        onthefly_corr_freq = float_list(default=None)  # Frequencies to use for correction
        use_n_cycles = integer(default=3)  # Use N cycles to calculate the phase, to use all integrations set to None
    """

    reference_file_types = ['emicorr']

    def process(self, input):
        # Open the input data model
        with datamodels.open(input) as input_model:

            # Catch the cases to skip
            instrument = input_model.meta.instrument.name
            if instrument != 'MIRI':
                self.log.warning('EMI correction not implemented for instrument: {}'.format(instrument))
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            readpatt = input_model.meta.exposure.readpatt
            allowed_readpatts = ['FAST', 'FASTR1', 'SLOW', 'SLOWR1']
            if readpatt.upper() not in allowed_readpatts:
                self.log.warning('EMI correction not implemented for read pattern: {}'.format(readpatt))
                input_model.meta.cal_step.emicorr = 'SKIPPED'
                return input_model

            # Setup parameters
            pars = {
                'save_intermediate_results': self.save_intermediate_results,
                'user_supplied_reffile': self.user_supplied_reffile,
                'nints_to_phase': self.nints_to_phase,
                'nbins': self.nbins,
                'scale_reference': self.scale_reference,
                'onthefly_corr_freq': self.onthefly_corr_freq,
                'use_n_cycles': self.use_n_cycles
            }

            # Get the reference file
            save_onthefly_reffile, emicorr_ref_filename, emicorr_model = None, None, None
            if self.onthefly_corr_freq is not None:
                emicorr_ref_filename = None
                self.log.info('Correcting with reference file created on-the-fly.')

            elif self.user_supplied_reffile is None:
                emicorr_ref_filename = self.get_reference_file(input_model, 'emicorr')
                # Skip the spep if no reference file is found
                if emicorr_ref_filename == 'N/A':
                    self.log.warning('No reference file found.')
                    self.log.warning('EMICORR step will be skipped')
                    input_model.meta.cal_step.emicorr = 'SKIPPED'
                    return input_model
                else:
                    self.log.info('Using CRDS reference file: {}'.format(emicorr_ref_filename))
                    emicorr_model = datamodels.EmiModel(emicorr_ref_filename)

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
