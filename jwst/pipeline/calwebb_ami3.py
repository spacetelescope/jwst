#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from ..associations import load_asn
from .. import datamodels


# step imports
from ..ami import ami_analyze_step
from ..ami import ami_average_step
from ..ami import ami_normalize_step


__version__ = "1.2"

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
        save_averages = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'ami_analyze': ami_analyze_step.AmiAnalyzeStep,
                 'ami_average': ami_average_step.AmiAverageStep,
                 'ami_normalize': ami_normalize_step.AmiNormalizeStep
                 }

    def process(self, input):

        log.info('Starting calwebb_ami3')

        # Load the input association table
        with open(input, 'r') as input_fh:
            asn = load_asn(input_fh)

        # We assume there's one final product defined by the association
        prod = asn['products'][0]

        # Construct lists of all the PSF and science target members
        # so that we know what we've been given to work with
        psf_files = [m['expname'] for m in prod['members'] if m['exptype'].upper() =='PSF']
        targ_files = [m['expname'] for m in prod['members'] if m['exptype'].upper() =='SCIENCE']

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
        psf_files = []
        targ_files = []
        for member in prod['members']:
            input_file = member['expname']

            # Do the LG analysis for this image
            log.debug('Do LG processing for member %s', input_file)
            result = self.ami_analyze(input_file)

            # Save the LG analysis results to a file
            output_file = mk_filename(self.output_dir, input_file, 'ami')
            self.log.info('Saving LG results to %s', output_file)
            result.save(output_file)

            # Save the result file name for input to ami_average
            if member['exptype'].upper() == 'PSF':
                psf_files.append(output_file)
            if member['exptype'].upper() == 'SCIENCE':
                targ_files.append(output_file)

        # Average the reference PSF image results
        psf_avg = None
        if len(psf_files) > 0:
            self.log.debug('Calling ami_average for PSF results ...')
            psf_avg = self.ami_average(psf_files)

            # Save the results to a file, if requested
            if self.save_averages:
                output_file = mk_prodname(self.output_dir, prod['psf_name'], 'amiavg')
                self.log.info('Saving averaged PSF results to %s', output_file)
                psf_avg.save(output_file)

        # Average the science target image results
        if len(targ_files) > 0:
            self.log.debug(' Calling ami_average for target results ...')
            targ_avg = self.ami_average(targ_files)

            # Save the results to a file, if requested
            if self.save_averages:
                output_file = mk_prodname(self.output_dir, prod['name'], 'amiavg')
                self.log.info('Saving averaged target results to %s', output_file)
                targ_avg.save(output_file)

        # Now that all LGAVG products have been produced, do normalization of
        # the target results by the reference results, if reference results exist
        if psf_avg is not None:

            result = self.ami_normalize(targ_avg, psf_avg)

            # Save the result
            output_file = mk_prodname(self.output_dir, prod['name'], 'aminorm')
            self.log.info('Saving normalized result to %s', output_file)
            result.save(output_file)
            result.close()

        # We're done
        log.info('... ending calwebb_ami3')

        return


def mk_filename(output_dir, filename, suffix):

    """
    Build a file name to use when saving results.

    An existing input file name is used as a template. A user-specified
    output directory path is prepended to the root of the input file name.
    The last product type suffix contained in the input file name is
    replaced with the specified new suffix. Any existing file name
    extension (e.g. ".fits") is preserved.

    Args:
        output_dir (str): The output_dir requested by the user
        filename (str): The input file name, to be reworked
        suffix (str): The desired file type suffix for the new file name

    Returns:
        string: The new file name

    Examples:
        For output_dir='/my/path', filename='jw12345_nrca_cal.fits', and
        suffix='i2d', the returned file name will be
        '/my/path/jw12345_nrca_i2d.fits'
    """

    # If the user specified an output_dir, replace any existing
    # path with output_dir
    if output_dir is not None:
        dirname, filename = os.path.split(filename)
        filename = os.path.join(output_dir, filename)

    # Now replace the existing suffix with the new one
    base, ext = os.path.splitext(filename)
    return base[:base.rfind('_')] + '_' + suffix + ext


def mk_prodname(output_dir, filename, suffix):

    """
    Build a file name based on an ASN product name template.

    The input ASN product name is used as a template. A user-specified
    output directory path is prepended to the root of the product name.
    The input product type suffix is appended to the root of the input
    product name, preserving any existing file name extension
    (e.g. ".fits").

    Args:
        output_dir (str): The output_dir requested by the user
        filename (str): The input file name, to be reworked
        suffix (str): The desired file type suffix for the new file name

    Returns:
        string: The new file name

    Examples:
        For output_dir='/my/path', filename='jw12345_nrca_cal.fits', and
        suffix='i2d', the returned file name will be
        '/my/path/jw12345_nrca_cal_i2d.fits'
    """

    # If the user specified an output_dir, replace any existing
    # path with output_dir
    if output_dir is not None:
        dirname, filename = os.path.split(filename)
        filename = os.path.join(output_dir, filename)

    # Now append the new suffix to the root name, preserving
    # any existing extension
    base, ext = os.path.splitext(filename)
    if len(ext) == 0:
        ext = ".fits"
    return base + '_' + suffix + ext
