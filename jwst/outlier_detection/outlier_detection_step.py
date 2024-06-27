"""Public common step definition for OutlierDetection processing."""
from functools import partial

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.stpipe import Step
from jwst.lib.pipe_utils import is_tso

from jwst.outlier_detection import outlier_detection
from jwst.outlier_detection import outlier_detection_ifu
from jwst.outlier_detection import outlier_detection_spec
from jwst.outlier_detection import outlier_detection_tso

# Categorize all supported versions of outlier_detection
outlier_registry = {'imaging': outlier_detection.OutlierDetection,
                    'ifu': outlier_detection_ifu.OutlierDetectionIFU,
                    'slitspec': outlier_detection_spec.OutlierDetectionSpec,
                    'tso': outlier_detection_tso.OutlierDetectionTSO
                    }

# Categorize all supported modes
IMAGE_MODES = ['NRC_IMAGE', 'MIR_IMAGE', 'NRS_IMAGE', 'NIS_IMAGE', 'FGS_IMAGE']
SLIT_SPEC_MODES = ['NRC_WFSS', 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT',
                   'NRS_MSASPEC', 'NIS_WFSS']
TSO_SPEC_MODES = ['NIS_SOSS', 'MIR_LRS-SLITLESS', 'NRC_TSGRISM',
                  'NRS_BRIGHTOBJ']
IFU_SPEC_MODES = ['NRS_IFU', 'MIR_MRS']
TSO_IMAGE_MODES = ['NRC_TSIMAGE']  # missing MIR_IMAGE with TSOVIST=True, not really addable
CORON_IMAGE_MODES = ['NRC_CORON', 'MIR_LYOT', 'MIR_4QPM']

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(Step):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input_data : asn file or ~jwst.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

    class_alias = "outlier_detection"

    # The members of spec needs to be a super-set of all parameters needed
    # by the various versions of the outlier_detection algorithms, and each
    # version will pick and choose what they need while ignoring the rest.
    # There are effectively N modes which use the following parameters:
    #    "-" = parameter is disabled
    #    " " = parameter is ignored
    #    "+" = parameter is used
    #    "!" = used only if another parameter is set
    #    "?" = idk?
    #                           | Image | Coron | TSO | Spec(IFU) | Spec |
    # --------------------------------------------------------------------
    # weight_type               |   !   |       |  +  |           |   !  |
    # pixfrac                   |   !   |       |     |           |   !  |
    # kernel                    |   !   |       |     |           |   !  |
    # fillval                   |   !   |       |     |           |   !  |
    # maskpt                    |       |       |  +  |           |      |
    # snr                       |   +   |   +   |  +  |           |      |
    # scale                     |   +   |   +   |  +  |           |      |
    # backg                     |   +   |   +   |  +  |           |      |
    # kernel_size               |       |       |     |     +     |      |
    # threshold_percent         |       |       |     |     +     |      |
    # rolling_window_width      |       |       |  +  |           |      |
    # ifu_second_check          |       |       |     |     +     |      |
    # save_intermediate_results |   +   |   +   |  +  |     +     |   +  |
    # resample_data             |   +   |   -   |  -  |           |   +  |
    # good_bits                 |   +   |   +   |  +  |     +     |   +  |
    # search_output_file *2     |       |       |     |           |      |
    # allowed_memory            |   !   |       |     |           |   !  |
    # in_memory *1              |   +   |   +   |  +  |     +     |   +  |
    #
    # *1 in_memory is used in the generic step when opening the input. But
    #    not all modes use the parameter internally.
    #
    # *2 This is used internally in stpipe. The table will not document this
    #    feature as it's difficult to know when it will and won't be used.
    #
    # For each mode the input is:
    # Image : ModelContainer of ImageModel
    # Coron : CubeModel (either a psf or target)
    # TSO   : CubeModel
    # IFU   : ModelContainer of IFUImageModel
    # Spec  : SourceModelContainer of SlitModel (created from MultiSlitModel)
    #         or ModelContainer of ImageModel (possible for MIR_LRS-FIXEDSLIT)
    #         or ModelContainer of SlitModel (for...?)
    #
    # So the only mode that benefits from `in_memory` is Image.
    #
    # For Coron (which uses only the non-resampled version of Image) the fake
    # ModelContainer made from the input CubeModel is all in memory so
    # `_remove_file` is skipped (for the open models) for both the
    # "drizzled_models" and "blot_models".
    #
    # For TSO `in_memory` isn't used (even for the calls to the parent
    # `compute_weight_threshold` and `save_median`).
    #
    # For IFU `in_memory` isn't used (and the algorithm is very different using
    # almost nothing from the parent class).
    #
    # For Spec(non-IFU) `in_memory` is used. However the input
    # `SourceModelContainer` doesn't use the "save_open" "return_open"
    # arguments when made from a `MultiSlitModel` (as in spec3). So the
    # input is a list of "open" `SlitModel`s. The call to `ResampleSpec` does
    # use `in_memory` and writes out the resampled `SlitModels` returning
    # a `ModelContainer` (if `resample=True`). `create_median` and
    # `blot_median` are used from the parent class and use `in_memory` in
    # the same way as the Image mode.

    spec = """
        weight_type = option('ivm','exptime',default='ivm')
        pixfrac = float(default=1.0)
        kernel = string(default='square') # drizzle kernel
        fillval = string(default='INDEF')
        nlow = integer(default=0)  # DEPRECATED this setting has no effect and will be removed
        nhigh = integer(default=0)  # DEPRECATED this setting has no effect and will be removed
        maskpt = float(default=0.7)
        snr = string(default='5.0 4.0')
        scale = string(default='1.2 0.7')
        backg = float(default=0.0)
        kernel_size = string(default='7 7')
        threshold_percent = float(default=99.8)
        rolling_window_width = integer(default=25)
        ifu_second_check = boolean(default=False)
        save_intermediate_results = boolean(default=False)
        resample_data = boolean(default=True)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
        search_output_file = boolean(default=False)
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=False)
    """

    def process(self, input_data):
        """Perform outlier detection processing on input data."""

        with datamodels.open(input_data, save_open=self.in_memory) as input_models:
            if not isinstance(input_models, ModelContainer):
                self.input_container = False
            else:
                self.input_container = True
            # Setup output path naming if associations are involved.
            asn_id = None
            try:
                asn_id = input_models.meta.asn_table.asn_id
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
                'weight_type': self.weight_type,  # for calling the resample step
                'wht_type': self.weight_type,  # for calling the resample class directly
                'pixfrac': self.pixfrac,
                'kernel': self.kernel,
                'fillval': self.fillval,
                'nlow': self.nlow,
                'nhigh': self.nhigh,
                'maskpt': self.maskpt,
                'snr': self.snr,
                'scale': self.scale,
                'backg': self.backg,
                'kernel_size': self.kernel_size,
                'threshold_percent': self.threshold_percent,
                'rolling_window_width': self.rolling_window_width,
                'ifu_second_check': self.ifu_second_check,
                'allowed_memory': self.allowed_memory,
                'in_memory': self.in_memory,
                'save_intermediate_results': self.save_intermediate_results,
                'resample_data': self.resample_data,
                'good_bits': self.good_bits,
                'make_output_path': self.make_output_path,
            }

            # Select which version of OutlierDetection
            # needs to be used depending on the input data
            if self.input_container:
                single_model = input_models[0]
            else:
                single_model = input_models
            exptype = single_model.meta.exposure.type
            self.check_input(input_models)

            if is_tso(single_model):
                # force resampling off and use rolling median
                pars['resample_data'] = False
                detection_step = outlier_registry['tso']
            elif exptype in CORON_IMAGE_MODES:
                # force resampling off but use same workflow as imaging
                pars['resample_data'] = False
                detection_step = outlier_registry['imaging']
            elif exptype in IMAGE_MODES:
                # imaging with resampling
                detection_step = outlier_registry['imaging']
                pars['resample_suffix'] = 'i2d'
            elif exptype in SLIT_SPEC_MODES:
                detection_step = outlier_registry['slitspec']
                pars['resample_suffix'] = 's2d'
            elif exptype in IFU_SPEC_MODES:
                # select algorithm for IFU data
                detection_step = outlier_registry['ifu']
            else:
                self.log.error("Outlier detection failed for unknown/unsupported ",
                               f"exposure type: {exptype}")
                self.valid_input = False

            if not self.valid_input:
                if self.input_container:
                    for model in input_models:
                        model.meta.cal_step.outlier_detection = "SKIPPED"
                else:
                    input_models.meta.cal_step.outlier_detection = "SKIPPED"
                return input_models

            self.log.debug(f"Using {detection_step.__name__} class for outlier_detection")

            # Set up outlier detection, then do detection
            step = detection_step(input_models, asn_id=asn_id, **pars)
            step.do_detection(input_models)

            state = 'COMPLETE'
            if self.input_container:
                for model in input_models:
                    model.meta.cal_step.outlier_detection = state
            else:
                input_models.meta.cal_step.outlier_detection = state
            return input_models


    def check_input(self, input_models):
        """Use this method to determine whether input is valid or not."""
        if self.input_container:
            self._check_input_container(input_models)
        else:
            self._check_input_cube(input_models)

    def _check_input_container(self, input_models):
        """Check to see whether input is the expected ModelContainer object."""
        ninputs = len(input_models)
        if not isinstance(input_models, ModelContainer):
            self.log.warning("Input is not a ModelContainer")
            self.log.warning("Outlier detection step will be skipped")
            self.valid_input = False
        elif ninputs < 2:
            self.log.warning(f"Input only contains {ninputs} exposure")
            self.log.warning("Outlier detection step will be skipped")
            self.valid_input = False
        else:
            self.valid_input = True
            self.log.info(f"Performing outlier detection on {ninputs} inputs")

    def _check_input_cube(self, input_models):
        """Check to see whether input is the expected CubeModel object."""
        ninputs = input_models.shape[0]
        if type(input_models) not in [datamodels.CubeModel, datamodels.SlitModel]:
            self.log.warning("Input is not the expected CubeModel")
            self.log.warning("Outlier detection step will be skipped")
            self.valid_input = False
        elif ninputs < 2:
            self.log.warning(f"Input only contains {ninputs} integration")
            self.log.warning("Outlier detection step will be skipped")
            self.valid_input = False
        else:
            self.valid_input = True
            self.log.info(f"Performing outlier detection with {ninputs} inputs")
