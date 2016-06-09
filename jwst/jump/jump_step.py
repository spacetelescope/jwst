#! /usr/bin/env python

from ..stpipe import Step, cmdline
import ..datamodels
from .jump import detect_jumps


class JumpStep( Step ):
    """
    JumpStep: Performs CR/jump detection on each ramp integration within an
    exposure. The 2-point difference method is applied on a first pass to all
    of the data. The y-intercept method is then applied in a second pass to
    only those pixels that have signal levels in the readout noise regime.
    """

    spec = """
        rejection_threshold = float(default=4.0,min=0) # CR rejection threshold
        do_yintercept = boolean(default=False) # do y-intercept method?
        yint_threshold = float(default=1.0,min=0) # y-intercept signal threshold
    """

    reference_file_types = ['gain','readnoise']

    def process(self, input):

        with models.open(input) as input_model:

            # Check for an input model with NGROUPS <=2
            if input_model.meta.exposure.ngroups<=2:
                self.log.warn('Can not apply jump detection when NGROUPS<=2;')
                self.log.warn('Jump step will be skipped')

                result = input_model.copy()
                result.meta.cal_step.jump = 'SKIPPED'
                return result

            # Retrieve the parameter values
            rej_thresh = self.rejection_threshold
            do_yint    = self.do_yintercept
            sig_thresh = self.yint_threshold
            self.log.info('CR rejection threshold = %g sigma', rej_thresh)
            if do_yint:
                self.log.info('Y-intercept signal threshold = %g', sig_thresh)

            # Get the gain and readnoise reference files
            gain_filename = self.get_reference_file( input_model, 'gain')
            self.log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = models.GainModel( gain_filename )

            readnoise_filename = self.get_reference_file( input_model,
                                                          'readnoise')
            self.log.info('Using READNOISE reference file: %s',
                          readnoise_filename)
            readnoise_model = models.ReadnoiseModel( readnoise_filename )

            # Call the jump detection routine
            result = detect_jumps( input_model, gain_model, readnoise_model,
                                   rej_thresh, do_yint, sig_thresh )

            gain_model.close()
            readnoise_model.close()


        result.meta.cal_step.jump = 'COMPLETE'

        return result


if __name__ == '__main__':
    cmdline.step_script( JumpStep )
