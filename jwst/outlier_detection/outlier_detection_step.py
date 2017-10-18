from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import outlier_detection
from . import outlier_detection_scaled
from . import outlier_detection_ifu
from . import outlier_detection_spec
from ..resample import resample, resample_spec

# Categorize all supported versions of outlier_detection
outlier_registry = {'imaging' : outlier_detection.OutlierDetection,
                    'scaled' : outlier_detection_scaled.OutlierDetectionScaled,
                    'ifu' : outlier_detection_ifu.OutlierDetectionIFU,
                    'slitspec' : outlier_detection_spec.OutlierDetectionSpec
                    }

# Categorize all supported modes
IMAGE_MODES = ['NRC_IMAGE', 'MIR_IMAGE', 'NRS_IMAGE', 'NIS_IMAGE', 'FGS_IMAGE']
SLIT_SPEC_MODES = ['NRC_GRISM', 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NIS_WFSS']
TSO_SPEC_MODES = ['NIS_SOSS','MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NRS_BRIGHTOBJ']
IFU_SPEC_MODES = ['NRS_IFU', 'MIR_MRS']
TSO_IMAGE_MODES = ['NRC_TSIMAGE']
CORON_IMAGE_MODES = ['NRC_CORON', 'MIR_LYOT', 'MIR_4QPM']


class OutlierDetectionStep(Step):
    """
    Flag outlier bad pixels and cosmic rays in the DQ array of each input image

    Input images can listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input : asn file or ModelContainer
        Single filename association table, or a datamodels.ModelContainer.
    """
    # The members of spec needs to be a super-set of all parameters needed
    # by the various versions of the outlier_detection algorithms, and each 
    # version will pick and choose what they need while ignoring the rest.
    spec = """
        wht_type = option('exptime','error',None,default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square') # drizzle kernel
        fillval = string(default='INDEF')
        nlow = integer(default=0)
        nhigh = integer(default=0)
        maskpt = float(default=0.7)
        grow = integer(default=1)
        snr = string(default='4.0 3.0')
        scale = string(default='0.5 0.4')
        backg = float(default=0.0)
        save_intermediate_results = boolean(default=False)
        resample_data = boolean(default=True)
        good_bits = integer(default=4)
        scale_detection = boolean(default=False)
    """
    reference_file_types = ['gain', 'readnoise']
    prefetch_references = False

    def process(self, input):

        with datamodels.open(input) as input_models:
            if not isinstance(input_models, datamodels.ModelContainer):
                self.log.warning("Input is not a ModelContainer.")
                self.log.warning("Outlier detection step will be skipped.")
                result = input_models.copy()
                result.meta.cal_step.outlier_detection = "SKIPPED"
                return result

            self.log.info("Performing outlier detection on {} inputs".format(len(input_models)))
            self.input_models = input_models

            reffiles = {}
            reffiles['gain'] = self._build_reffile_container('gain')
            reffiles['readnoise'] = self._build_reffile_container('readnoise')

            pars = {
                'wht_type': self.wht_type,
                'pixfrac': self.pixfrac,
                'kernel': self.kernel,
                'fillval': self.fillval,
                'nlow': self.nlow,
                'nhigh': self.nhigh,
                'maskpt': self.maskpt,
                'grow': self.grow,
                'snr': self.snr,
                'scale': self.scale,
                'backg': self.backg,
                'save_intermediate_results': self.save_intermediate_results,
                'resample_data': self.resample_data,
                'good_bits': self.good_bits
                }
            # Add logic here to select which version of OutlierDetection 
            # needs to be used depending on the input data
            exptype = input_models[0].meta.exposure.type 
            
            if exptype in IMAGE_MODES:
                # default mode: imaging with resampling
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 'i2d'
            elif exptype in TSO_SPEC_MODES: 
                # algorithm selected for TSO data (no resampling)
                pars['resample_data'] = False # force resampling off...
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 'i2d'
            elif exptype in TSO_IMAGE_MODES+CORON_IMAGE_MODES:
                # algorithm selected for TSO data (no resampling)
                pars['resample_data'] = False # force resampling off...
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 's2d'                
            elif exptype in TSO_IMAGE_MODES and self.scale_detection:
                # selected scaled algorithm for TSO data 
                detection_step = outlier_registry['scaled']
                pars['resample_suffix'] = 'i2d'
            elif exptype in SLIT_SPEC_MODES:
                detection_step = outlier_registry['slit']
                pars['resample_suffix'] = 's2d'
            elif exptype in IFU_SPEC_MODES:  
                # select algorithm for IFU data
                detection_step = outlier_registry['ifu']
                pars['resample_suffix'] = 's3d'
            else:
                self.log.error("Outlier detection failed for unknown/unsupported exposure type: {}".format(exptype))
                result = input_models.copy()
                result.meta.cal_step.outlier_detection = "SKIPPED"
                return result
                
            # Set up outlier detection, then do detection
            step = detection_step(self.input_models,
                reffiles=reffiles, **pars)
            step.do_detection()

            for model in self.input_models:
                model.meta.cal_step.outlier_detection = 'COMPLETE'

            return self.input_models

    def _build_reffile_container(self, reftype):
        """
        Return a ModelContainer of reference file models.

        Parameters
        ----------

        input_models: ModelContainer
            the science data, ImageModels in a ModelContainer

        reftype: string
            type of reference file

        Returns
        -------

        a ModelContainer with corresponding reference files for each input model
        """

        reffile_to_model = {'gain': datamodels.GainModel,
            'readnoise': datamodels.ReadnoiseModel}

        reffiles = [self.get_reference_file(im, reftype) for im in self.input_models]
        self.log.debug("Using {} reffile(s):".format(reftype.upper()))
        for r in set(reffiles):
            self.log.debug("    {}".format(r))

        # Check if all the ref files are the same.  If so build it by reading
        # the reference file just once.
        if len(set(reffiles)) <= 1:
            length = len(self.input_models)
            ref_list = [reffile_to_model.get(reftype)(reffiles[0])] * length
        else:
            ref_list = [reffile_to_model.get(reftype)(ref) for ref in reffiles]

        return datamodels.ModelContainer(ref_list)


if __name__ == '__main__':
    cmdline.step_script(OutlierDetectionStep)
