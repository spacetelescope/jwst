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

    def process(self, input):
        """Primary method for performing pipeline."""
        self.log.info('Starting calwebb_coron3 ...')

        # Load the input association table
        asn = self.load_as_level3_asn(input)
        acid = asn.get('asn_id', '')

        # We assume there's one final product defined by the association
        prod = asn['products'][0]
        self.output_file = prod.get('name', self.output_file)

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

        # Stack all the PSF images into a single CubeModel
        psf_stack = self.stack_refs(psf_models)
        psf_models.close()

        # Save the resulting PSF stack
        self.save_model(psf_stack, suffix='psfstack')

        # Call the sequence of steps align_refs, klip, and outlier_detection
        # once for each input target exposure
        resample_input = datamodels.ModelContainer()
        for target_file in targ_files:

            # Call align_refs
            self.log.debug('Calling align_refs for member %s', target_file)
            psf_aligned = self.align_refs(target_file, psf_stack)

            # Save the alignment results
            self.save_model(
                psf_aligned, output_file=target_file,
                suffix='psfalign', acid=acid
            )

            # Call KLIP
            self.log.debug('Calling klip for member %s', target_file)
            psf_sub = self.klip(target_file, psf_aligned)
            psf_aligned.close()

            # Save the psf subtraction results
            self.save_model(
                psf_sub, output_file=target_file,
                suffix='psfsub', acid=acid
            )

            # Create a ModelContainer of the psf_sub results to send to
            # outlier_detection
            self.log.debug('Building ModelContainer of klip results')
            target_models = datamodels.ModelContainer()
            for i in range(psf_sub.data.shape[0]):
                image = datamodels.ImageModel(data=psf_sub.data[i],
                                              err=psf_sub.err[i],
                                              dq=psf_sub.dq[i])
                image.update(psf_sub)
                image.meta.wcs = psf_sub.meta.wcs
                target_models.append(image)

            # Call outlier_detection
            target_models = self.outlier_detection(target_models)

            # Create Level 2c products
            if target_models[0].meta.cal_step.outlier_detection == 'COMPLETE':
                err_str1 = "Creating Level 2c output with updated DQ arrays..."
                self.log.info(err_str1)
                lev2c_model = psf_sub.copy()
                # Replace Level 2b product DQ array with Level 2c DQ array
                for i in range(len(target_models)):
                    lev2c_model.dq[i] = target_models[i].dq
                lev2c_model.meta.cal_step.outlier_detection = 'COMPLETE'
                self.save_model(lev2c_model, output_file=target_file,
                                suffix='crfints', acid=acid)

            # Append results from this target exposure to resample input model
            for i in range(len(target_models)):
                resample_input.append(target_models[i])

        # Call the resample step to combine all psf-subtracted target images
        result = self.resample(resample_input)

        if result == resample_input:
            # Resampling was skipped
            #     yet we need to return an output resampled ImageModel, so...
            #warn1 = 'Creating fake resample results until step is available'
            #self.log.warning(warn1)
            #result = datamodels.ImageModel(data=resample_input[0].data,
            #                                     con=resample_input[0].dq,
            #                                     wht=resample_input[0].err)
            #result.update(resample_input[0])
            # The resample step blends headers already...
            self.log.debug('Blending metadata for {}'.format(
                result.meta.filename))
            blendmeta.blendmodels(result, inputs=targ_files)

        result.meta.asn.pool_name = asn['asn_pool']
        result.meta.asn.table_name = op.basename(input)

        # Save the final result
        self.save_model(result, suffix=self.suffix)

        # We're done
        self.log.info('...ending calwebb_coron3')

        return
