#!/usr/bin/env python
import os.path as op

from ..stpipe import Pipeline
from .. import datamodels
from ..model_blender import blendmeta

# step imports
from ..coron import stack_refs_step
from ..coron import align_refs_step
from ..coron import klip_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_step

__all__ = ['Coron3Pipeline']


class Coron3Pipeline(Pipeline):
    """Class for defining Coron3Pipeline.

    Coron3Pipeline: Apply all level-3 calibration steps to a
    coronagraphic association of exposures. Included steps are:

    #. stack_refs (assemble reference PSF inputs)
    #. align_refs (align reference PSFs to target images)
    #. klip (PSF subtraction using the KLIP algorithm)
    #. outlier_detection (flag outliers)
    #. resample (image combination and resampling)

    """

    spec = """
        suffix = string(default='i2d')
    """

    # Define aliases to steps
    step_defs = {
        'stack_refs': stack_refs_step.StackRefsStep,
        'align_refs': align_refs_step.AlignRefsStep,
        'klip': klip_step.KlipStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample': resample_step.ResampleStep
    }

    prefetch_references = False

    def process(self, user_input):
        """Primary method for performing pipeline."""
        self.log.info('Starting calwebb_coron3 ...')

        # Load the input association table
        asn = self.load_as_level3_asn(user_input)
        acid = asn.get('asn_id', '')

        # We assume there's one final product defined by the association
        prod = asn['products'][0]
        self.output_file = prod.get('name', self.output_file)

        # Setup required output products and formats
        self.outlier_detection.suffix = f'{acid}_crfints'
        self.outlier_detection.save_results = self.save_results
        self.resample.blendheaders = False

        # Save the original outlier_detection.skip setting from the
        # input, because it may get toggled off within loops for
        # processing individual inputs
        skip_outlier_detection = self.outlier_detection.skip

        # Construct lists of all the PSF and science target members
        psf_files = [m['expname'] for m in prod['members']
                     if m['exptype'].upper() == 'PSF']
        targ_files = [m['expname'] for m in prod['members']
                      if m['exptype'].upper() == 'SCIENCE']

        # Make sure we found some PSF and target members
        if len(psf_files) == 0:
            err_str1 = 'No reference PSF members found in association table.'
            self.log.error(err_str1)
            self.log.error('Calwebb_coron3 processing will be aborted')
            return

        if len(targ_files) == 0:
            err_str1 = 'No science target members found in association table'
            self.log.error(err_str1)
            self.log.error('Calwebb_coron3 processing will be aborted')
            return

        for member in psf_files + targ_files:
            self.prefetch(member)

        # Assemble all the input psf files into a single ModelContainer
        psf_models = datamodels.ModelContainer()
        for i in range(len(psf_files)):
            psf_input = datamodels.CubeModel(psf_files[i])
            psf_models.append(psf_input)
            psf_input.close()

        # Perform outlier detection on the PSFs.
        if not skip_outlier_detection:
            for model in psf_models:
                self.outlier_detection(model)
                # step may have been skipped for this model;
                # turn back on for next model
                self.outlier_detection.skip = False
        else:
            self.log.info('Outlier detection skipped for PSF\'s')

        # Stack all the PSF images into a single CubeModel
        psf_stack = self.stack_refs(psf_models)
        psf_models.close()

        # Save the resulting PSF stack
        psf_stack.meta.filetype = 'psf stack'
        self.save_model(psf_stack, suffix='psfstack')

        # Call the sequence of steps outlier_detection, align_refs, and klip
        # once for each input target exposure
        resample_input = datamodels.ModelContainer()
        for target_file in targ_files:
            with datamodels.open(target_file) as target:

                # Remove outliers from the target.
                if not skip_outlier_detection:
                    target = self.outlier_detection(target)
                    # step may have been skipped for this model;
                    # turn back on for next model
                    self.outlier_detection.skip = False

                # Call align_refs
                psf_aligned = self.align_refs(target, psf_stack)

                # Save the alignment results
                psf_aligned.meta.filetype = 'psf aligned'
                self.save_model(
                    psf_aligned, output_file=target_file,
                    suffix='psfalign', acid=acid
                )

                # Call KLIP
                psf_sub = self.klip(target, psf_aligned)
                psf_aligned.close()

                # Save the psf subtraction results
                psf_sub.meta.filetype = 'psf subtracted'
                self.save_model(
                    psf_sub, output_file=target_file,
                    suffix='psfsub', acid=acid
                )

                # Split out the integrations into separate models
                # in a ModelContainer to pass to `resample`
                for model in psf_sub.to_container():
                    resample_input.append(model)

        # Call the resample step to combine all psf-subtracted target images
        result = self.resample(resample_input)

        # Blend the science headers
        try:
            completed = result.meta.cal_step.resample
        except AttributeError:
            self.log.debug('Could not determine whether resample was completed. Presuming not.')
            completed = 'SKIPPED'
        if completed == 'COMPLETE':
            self.log.debug(f'Blending metadata for {result}')
            blendmeta.blendmodels(result, inputs=targ_files)

        try:
            result.meta.asn.pool_name = asn['asn_pool']
            result.meta.asn.table_name = op.basename(user_input)
        except AttributeError:
            self.log.debug(f'Cannot set association information on final result {result}')

        # Save the final result
        self.save_model(result, suffix=self.suffix)

        # We're done
        self.log.info('...ending calwebb_coron3')

        return
