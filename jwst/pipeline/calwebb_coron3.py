#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from ..associations import load_asn
from .. import datamodels
from ..resample import blend

# step imports
from ..coron import stack_refs_step
from ..coron import align_refs_step
from ..coron import klip_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_step


__version__ = "0.7.1"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Coron3Pipeline(Pipeline):
    """
    Coron3Pipeline: Apply all level-3 calibration steps to a
    coronagraphic association of exposures. Included steps are:
    stack_refs (assemble reference PSF inputs)
    align_refs (align reference PSFs to target images)
    klip (PSF subtraction using the KLIP algorithm)
    outlier_detection (flag outliers)
    resample (image combination and resampling)
    """

    spec = """
    """

    # Define aliases to steps
    step_defs = {'stack_refs': stack_refs_step.StackRefsStep,
                 'align_refs': align_refs_step.AlignRefsStep,
                 'klip': klip_step.KlipStep,
                 'outlier_detection': outlier_detection_step.OutlierDetectionStep,
                 'resample': resample_step.ResampleStep
                 }

    def process(self, input):

        log.info('Starting calwebb_coron3 ...')

        # Load the input association table
        with open(input, 'r') as input_fh:
            asn = load_asn(input_fh)

        # We assume there's one final product defined by the association
        prod = asn['products'][0]

        # Construct lists of all the PSF and science target members
        psf_files = [m['expname'] for m in prod['members'] if m['exptype'].upper() == 'PSF']
        targ_files = [m['expname'] for m in prod['members'] if m['exptype'].upper() == 'SCIENCE']

        # Make sure we found some PSF and target members
        if len(psf_files) == 0:
            log.error('No reference PSF members found in association table')
            log.error('Calwebb_coron3 processing will be aborted')
            return

        if len(targ_files) == 0:
            log.error('No science target members found in association table')
            log.error('Calwebb_coron3 processing will be aborted')
            return

        # Assemble all the input psf files into a single ModelContainer
        psf_models = datamodels.ModelContainer()
        for i in range(len(psf_files)):
            psf_input = datamodels.CubeModel(psf_files[i])
            psf_models.append(psf_input)
            psf_input.close()

        # Call the stack_refs step to stack all the PSF images into
        # a single CubeModel
        psf_stack = self.stack_refs(psf_models)
        psf_models.close()

        # Save the resulting PSF stack
        output_file = mk_prodname(self.output_dir, prod['name'], 'psfstack')
        log.info('Saving psfstack file %s', output_file)
        psf_stack.save(output_file)

        # Call the sequence of steps align_refs, klip, and outlier_detection
        # once for each input target exposure
        resample_input = datamodels.ModelContainer()
        for target_file in targ_files:

            # Call align_refs
            log.debug('Calling align_refs for member %s', target_file)
            psf_aligned = self.align_refs(target_file, psf_stack)

            # Save the alignment results
            filename = mk_filename(self.output_dir, target_file, 'psfalign')
            log.info('Saving psfalign file %s', filename)
            psf_aligned.save(filename)

            # Call KLIP
            log.debug('Calling klip for member %s', target_file)
            #psf_sub, psf_fit = self.klip(target_file, psf_aligned)
            psf_sub = self.klip(target_file, psf_aligned)
            psf_aligned.close()

            # Save the psf subtraction results
            filename = mk_filename(self.output_dir, target_file, 'psfsub')
            log.info('Saving psfsub file %s', filename)
            psf_sub.save(filename)

            # Create a ModelContainer of the psf_sub results to send to
            # outlier_detection
            log.debug('Building ModelContainer of klip results')
            target_models = datamodels.ModelContainer()
            for i in range(psf_sub.data.shape[0]):
                image = datamodels.ImageModel(data=psf_sub.data[i],
                        err=psf_sub.err[i], dq=psf_sub.dq[i])
                image.update(psf_sub)
                image.meta.wcs = psf_sub.meta.wcs
                target_models.append(image)

            # Call outlier_detection
            self.outlier_detection.resample_data=False
            target_models = self.outlier_detection(target_models)

            # Create a level-2c output product
            log.debug('Creating and saving Level-2C result')
            lev2c_name = mk_filename(self.output_dir, target_file, 'calints-'+asn['asn_id'])
            lev2c_model = psf_sub.copy()
            # Update/replace L2B product DQ array with L2C results
            for i in range(len(target_models)):
                lev2c_model.dq[i] = target_models[i].dq
            lev2c_model.meta.cal_step.outlier_detection = 'COMPLETE'
            lev2c_model.save(lev2c_name)

            # Append results from this target exposure to resample input model
            for i in range(len(target_models)):
                resample_input.append(target_models[i])

        # Call the resample step to combine all the psf-subtracted target images
        result = self.resample(resample_input)
        output_file = mk_prodname(self.output_dir, prod['name'], 'i2d')

        if result == resample_input:
        # Resampling was skipped, 
        #  yet we need to return a DrizProductModel, so...
            log.warning('Creating fake resample results until step is available')
            result = datamodels.DrizProductModel(data=resample_input[0].data,
                                             con=resample_input[0].dq,
                                             wht=resample_input[0].err)
            result.update(resample_input[0])
            # The resample step blends headers already...
            log.debug('Blending metadata for {}'.format(output_file))
            blend.blendfitsdata(targ_files, result)
            result.meta.asn.pool_name = asn['asn_pool']
            result.meta.asn.table_name = input

        result.meta.cal_step.outlier_detection = 'COMPLETE'
        result.meta.cal_step.resample = 'COMPLETE'
        result.meta.model_type = 'DrizProductModel'

        # Save the final result
        log.info('Saving final result to %s', output_file)
        result.save(output_file)
        result.close()

        # We're done
        log.info('... ending calwebb_coron3')

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
