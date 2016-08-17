#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from ..associations import Association
from .. import datamodels

# step imports
from ..coron import stack_refs_step
from ..coron import align_refs_step
from ..coron import klip_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_step


__version__ = "0.7.0"

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

    # Define alias to steps
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
            asn = Association.load(input_fh)

        # We assume there's one final product defined by the association
        prod = asn['products'][0]

        # Construct lists of all the PSF and science target members
        psf_files = []
        targ_files = []
        for member in prod['members']:
            if member['exptype'].upper() == 'PSF':
                psf_files.append(member['expname'])
                log.debug(' psf_file {0} = {1}'.format(len(psf_files), member['expname']))
            if member['exptype'].upper() == 'SCIENCE':
                targ_files.append(member['expname'])
                log.debug(' targ_file {0} = {1}'.format(len(targ_files), member['expname']))

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
            input = datamodels.CubeModel(psf_files[i])
            psf_models.append(input)
            input.close()

        # Call the stack_refs step to stack all the PSF images into
        # a single CubeModel
        psf_stack = self.stack_refs(psf_models)
        psf_models.close()

        # Save the resulting PSF stack
        filename = mk_filename(self.output_dir, prod['name'], 'psfstack')
        log.info('Saving psfstack file %s', filename)
        psf_stack.save(filename)

        # Call the sequence of steps align_refs, klip, and outlier_detection
        # once for each input target exposure
        resample_input = datamodels.ModelContainer()
        for target_file in targ_files:

            target_model = datamodels.CubeModel(target_file)

            # Call align_refs
            log.debug(' Calling align_refs for member %s', target_file)
            psf_aligned = self.align_refs(target_model, psf_stack)
    
            # Save the alignment results
            filename = mk_filename(self.output_dir, target_file, 'psfalign')
            log.info('Saving psfalign file %s', filename)
            psf_aligned.save(filename)

            # Call KLIP
            log.debug(' Calling klip for member %s', target_file)
            psf_sub, psf_fit = self.klip(target_model, psf_aligned)
            target_model.close()
            psf_aligned.close()

            # Save the psf subtraction results
            filename = mk_filename(self.output_dir, target_file, 'psfsub')
            log.info('Saving psfsub file %s', filename)
            psf_sub.save(filename)
         
            # Create a ModelContainer of the psf_sub results to send to
            # outlier_detection
            log.debug(' Building ModelContainer of klip results')
            outlier_input = datamodels.ModelContainer()
            for i in range(psf_sub.data.shape[0]):
                image = datamodels.ImageModel(data=psf_sub.data[i],
                        err=psf_sub.err[i], dq=psf_sub.dq[i])
                outlier_input.append(image)

            # Call outlier_detection
            targ_cr = self.outlier_detection(outlier_input)
            outlier_input.close()

            # Append results from this target exposure to resample input model
            for i in range(len(targ_cr)):
                resample_input.append(targ_cr[i])

        # Call the resample step to combine all the psf-subtracted target images
        result = self.reample(resample_input)

        # Save the final result
        filename = prod['name']
        if self.output_dir is not None:
            dirname, filename = os.path.split(filename)
            filename = os.path.join(self.output_dir, filename)
        self.log.info(' Saving final result to %s', filename)
        result.save(filename)
        result.close()

        # We're done
        log.info('... ending calwebb_coron3')

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

