"""Public common step definition for OutlierDetection processing."""
from functools import partial

from ..stpipe import Step
from .. import datamodels

from . import outlier_detection
from . import outlier_detection_scaled
from . import outlier_detection_ifu
from . import outlier_detection_spec

# Categorize all supported versions of outlier_detection
outlier_registry = {'imaging': outlier_detection.OutlierDetection,
                    'scaled': outlier_detection_scaled.OutlierDetectionScaled,
                    'ifu': outlier_detection_ifu.OutlierDetectionIFU,
                    'slitspec': outlier_detection_spec.OutlierDetectionSpec
                    }

# Categorize all supported modes
IMAGE_MODES = ['NRC_IMAGE', 'MIR_IMAGE', 'NRS_IMAGE', 'NIS_IMAGE', 'FGS_IMAGE']
SLIT_SPEC_MODES = ['NRC_WFSS', 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT',
                   'NRS_MSASPEC', 'NIS_WFSS']
TSO_SPEC_MODES = ['NIS_SOSS', 'MIR_LRS-SLITLESS', 'NRC_TSGRISM',
                  'NRS_BRIGHTOBJ']
IFU_SPEC_MODES = ['NRS_IFU', 'MIR_MRS']
TSO_IMAGE_MODES = ['NRC_TSIMAGE']
CORON_IMAGE_MODES = ['NRC_CORON', 'MIR_LYOT', 'MIR_4QPM']

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(Step):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already opened
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
        weight_type = option('exptime','error',None,default='exptime')
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
        good_bits = integer(default=6)
        scale_detection = boolean(default=False)
        search_output_file = boolean(default=False)
    """

    def process(self, input):
        """Perform outlier detection processing on input data."""
        with datamodels.open(input) as input_models:
            self.input_models = input_models
            if not isinstance(self.input_models, datamodels.ModelContainer):
                self.input_container = False
            else:
                self.input_container = True

            # Setup output path naming if associations are involved.
            asn_id = None
            try:
                asn_id = self.input_models.meta.asn_table.asn_id
            except (AttributeError, KeyError):
                pass
            if asn_id is None:
                asn_id = self.search_attr('asn_id')
            if asn_id is not None:
                _make_output_path = self.search_attr(
                    '_make_output_path', parent_first=True
                )
                self._make_output_path = partial(
                    _make_output_path,
                    asn_id=asn_id
                )

            # Setup outlier detection parameters
            pars = {
                'weight_type': self.weight_type,
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
                'good_bits': self.good_bits,
                'make_output_path': self.make_output_path,
            }

            # Add logic here to select which version of OutlierDetection
            # needs to be used depending on the input data
            if self.input_container:
                exptype = self.input_models[0].meta.exposure.type
            else:
                exptype = self.input_models.meta.exposure.type
            self.check_input()

            if exptype in IMAGE_MODES:
                # default mode: imaging with resampling
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 'i2d'
            elif exptype in TSO_SPEC_MODES:
                # algorithm selected for TSO data (no resampling)
                pars['resample_data'] = False  # force resampling off...
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 's2d'
            elif exptype in TSO_IMAGE_MODES+CORON_IMAGE_MODES and \
                    not self.scale_detection:
                # algorithm selected for TSO data (no resampling)
                pars['resample_data'] = False  # force resampling off...
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 'i2d'
            elif exptype in TSO_IMAGE_MODES and self.scale_detection:
                # selected scaled algorithm for TSO data
                detection_step = outlier_registry['scaled']
                pars['resample_suffix'] = 'i2d'
            elif exptype in SLIT_SPEC_MODES:
                detection_step = outlier_registry['slitspec']
                pars['resample_suffix'] = 's2d'
            elif exptype in IFU_SPEC_MODES:
                # select algorithm for IFU data
                detection_step = outlier_registry['ifu']
                pars['resample_suffix'] = 's3d'
            else:
                self.log.error("Outlier detection failed for \
                                unknown/unsupported exposure type: {}".format(
                                exptype))
                self.valid_input = False

            if not self.valid_input:
                result = input_models
                for input in result:
                    input.meta.cal_step.outlier_detection = "SKIPPED"
                    self.skip = True
                return result

            self.log.debug("Using {} class for outlier_detection".format(
                           detection_step.__name__))

            reffiles = {}

            # Set up outlier detection, then do detection
            step = detection_step(self.input_models, reffiles=reffiles, **pars)
            step.do_detection()

            if self.input_container:
                for model in self.input_models:
                    model.meta.cal_step.outlier_detection = 'COMPLETE'
            else:
                self.input_models.meta.cal_step.outlier_detection = 'COMPLETE'

            return self.input_models

    def _build_reffile_container(self, reftype):
        """Return a ModelContainer of reference file models.

        This method builds a ModelContainer object which contains all the
        reference files needed for processing the inputs.

        Parameters
        ----------
        input_models: ModelContainer
            the science data, ImageModels in a ModelContainer

        reftype: string
            type of reference file

        Returns
        -------
        ModelContainer with corresponding reference files for each input model

        """
        reffile_to_model = {'gain': datamodels.GainModel,
                            'readnoise': datamodels.ReadnoiseModel}
        if self.input_container:
            reffiles = [self.get_reference_file(im, reftype)
                        for im in self.input_models]
            length = len(self.input_models)
        else:
            reffiles = [self.get_reference_file(self.input_models, reftype)]
            length = 1

        self.log.debug("Using {} reffile(s):".format(reftype.upper()))
        for r in set(reffiles):
            self.log.debug("    {}".format(r))

        # Check if all the ref files are the same.  If so build it by reading
        # the reference file just once.
        if len(set(reffiles)) <= 1:
            ref_list = [reffile_to_model.get(reftype)(reffiles[0])] * length
        else:
            ref_list = [reffile_to_model.get(reftype)(ref) for ref in reffiles]

        return datamodels.ModelContainer(ref_list)

    def check_input(self):
        """Use this method to determine whether input is valid or not."""
        if self.input_container:
            self._check_input_container()
        else:
            self._check_input_cube()

    def _check_input_container(self):
        """Check to see whether input is the expected ModelContainer object."""
        ninputs = len(self.input_models)
        if not isinstance(self.input_models, datamodels.ModelContainer):
            self.log.warning("Input is not a ModelContainer.")
            self.log.warning("Outlier detection step will be skipped.")
            self.valid_input = False
        elif ninputs < 2:
            self.log.warning("Input only contains %d exposures." % (ninputs))
            self.log.warning("Outlier detection step will be skipped.")
            self.valid_input = False
        else:
            self.valid_input = True
            self.log.info("Performing outlier detection on %d inputs" %
                          (ninputs))

    def _check_input_cube(self):
        """Check to see whether input is the expected CubeModel object."""
        ninputs = self.input_models.shape[0]
        if not isinstance(self.input_models, datamodels.CubeModel):
            self.log.warning("Input is not the expected CubeModel.")
            self.log.warning("Outlier detection step will be skipped.")
            self.valid_input = False
        elif ninputs < 2:
            self.log.warning("Input only contains %d exposures." % (ninputs))
            self.log.warning("Outlier detection step will be skipped.")
            self.valid_input = False
        else:
            self.valid_input = True
            self.log.info("Performing outlier detection with %d inputs" %
                          (ninputs))
