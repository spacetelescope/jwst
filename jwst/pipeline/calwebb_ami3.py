#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from ..associations import Association
from .. import datamodels


# step imports
from ..ami import ami_analyze_step
from ..ami import ami_average_step
from ..ami import ami_normalize_step


__version__ = "1.0"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Ami3Pipeline(Pipeline):
    """

    Ami3Pipeline: Apply all level-3 calibration steps to an
    association of level-2b AMI exposures. Included steps are:
    ami_analyze (fringe detection)
    ami_average (average results of fringe detection)
    ami_normalize (normalize results by reference target)

    """

    spec = """
    """

    # Define alias to steps
    step_defs = {'ami_analyze': ami_analyze_step.AmiAnalyzeStep,
                 'ami_average': ami_average_step.AmiAverageStep,
                 'ami_normalize': ami_normalize_step.AmiNormalizeStep
                 }

    def process(self, input):

        log.info('Starting calwebb_ami3 ...')

        # Load the input association table
        with open(input, 'r') as input_fh:
            asn = Association.load(input_fh)

        # We assume there's one final product defined by the association
        prod = asn['products'][0]

        # Construct lists of all the PSF and science target members
        psf_files = []
        targ_files = []
        for member in prod['members']:
            if member['exptype'].upper() == 'PSF':
                psf_files.append(member['expname'])
                log.debug(' psf_file {0} = {1}'.format(len(psf_files),
                          member['expname']))
            if member['exptype'].upper() == 'SCIENCE':
                targ_files.append(member['expname'])
                log.debug(' targ_file {0} = {1}'.format(len(targ_files),
                          member['expname']))

        # Make sure we found some science target members
        if len(targ_files) == 0:
            log.error('No science target members found in association table')
            log.error('Calwebb_ami3 processing will be aborted')
            return

        # If there aren't any PSF images, we can't do the normalize step
        if len(psf_files) == 0:
            log.info('No PSF reference members found in association table')
            log.info('ami_normalize step will be skipped')

        # Loop over all the images in the assocation
        for member in prod['members']:
            input_file = member['expname']

            # Do the LG analysis for this image
            log.debug(' Do LG processing for member %s', input_file)
            result = self.ami_analyze(input_file)

            # Save the LG analysis results to a file
            output_file = mk_filename(self.output_dir, input_file, 'lg')
            self.log.debug(' Saving LG results to %s', output_file)
            result.save(output_file)

        # Average the PSF image results
        psf_avg = None
        if len(psf_files) > 0:
            self.log.debug(' Calling ami_average for PSF results ...')
            psf_avg = self.ami_average(psf_files)

            # Save the results to a file
            output_file = mk_filename(self.output_dir, psf_files[0], 'lgavg')
            self.log.debug(' Saving average PSF results to %s', output_file)
            psf_avg.save(output_file)

        # Average the science target image results
        if len(targ_files) > 0:
            self.log.debug(' Calling ami_average for target results ...')
            targ_avg = self.ami_average(targ_files)

            # Save the results to a file
            output_file = mk_filename(self.output_dir, targ_files[0], 'lgavg')
            self.log.debug(' Saving average target results to %s', output_file)
            targ_avg.save(output_file)

        # Now that all LGAVG products have been produced, do normalization of
        # the target results by reference results, if reference results exist
        if psf_avg is not None:

            result = self.ami_normalize(targ_avg, psf_avg)

            # Save the result
            output_file = mk_filename(self.output_dir, prod['name'], 'lgnorm')
            self.log.info(' Saving result to %s', output_file)
            result.save(output_file)
            result.close()

        # We're done
        log.info('... ending calwebb_ami3')

        return


def mk_filename(output_dir, filename, suffix):

    # If the user specified an output_dir, replace any existing
    # path with output_dir
    if output_dir is not None:
        dirname, filename = os.path.split(filename)
        filename = os.path.join(output_dir, filename)

    # Now replace the existing suffix with the new one
    base, ext = os.path.splitext(filename)
    return base[:base.rfind('_')] + '_' + suffix + ext
