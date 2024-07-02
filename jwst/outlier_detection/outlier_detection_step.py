"""Public common step definition for OutlierDetection processing."""
from functools import partial

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.stpipe import Step
from jwst.lib.pipe_utils import is_tso

from . import ifu, imaging, tso, spec

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
            if isinstance(input_models, ModelContainer):
                ninputs = len(input_models)
                if ninputs < 2:
                    self.log.warning(f"Input only contains {ninputs} exposure")
                    self.log.warning("Outlier detection step will be skipped")
                    return self._set_status(input_models, "SKIPPED")
            elif isinstance(input_models, (datamodels.CubeModel, datamodels.SlitModel)):
                ninputs = input_models.shape[0]
                if ninputs < 2:
                    self.log.warning(f"Input only contains {ninputs} integration")
                    self.log.warning("Outlier detection step will be skipped")
                    return self._set_status(input_models, "SKIPPED")
            else:
                self.log.warning("Input {input_models} is not supported")
                self.log.warning("Outlier detection step will be skipped")
                return self._set_status(input_models, "SKIPPED")
            self.log.info(f"Performing outlier detection with {ninputs} inputs")

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

            # Select which version of OutlierDetection
            # needs to be used depending on the input data
            if isinstance(input_models, ModelContainer):
                single_model = input_models[0]
            else:
                single_model = input_models
            exptype = single_model.meta.exposure.type

            snr1, snr2 = [float(v) for v in self.snr.split()]
            scale1, scale2 = [float(v) for v in self.scale.split()]

            if is_tso(single_model):
                tso.detect_outliers(
                    input_models,
                    self.save_intermediate_results,
                    self.good_bits,
                    self.maskpt,
                    self.rolling_window_width,
                    snr1,
                    snr2,
                    scale1,
                    scale2,
                    self.backg,
                    asn_id,
                    self.make_output_path,
                )
            elif exptype in CORON_IMAGE_MODES:
                imaging.detect_outliers(
                    input_models,
                    self.save_intermediate_results,
                    self.good_bits,
                    self.maskpt,
                    snr1,
                    snr2,
                    scale1,
                    scale2,  # TODO turn off the things that aren't actually used
                    self.backg,
                    False,  # force resampling off but use the same workflow as imaging
                    self.weight_type,
                    self.pixfrac,
                    self.kernel,
                    self.fillval,
                    self.allowed_memory,
                    self.in_memory,
                    asn_id,
                    self.make_output_path,
                )
            elif exptype in IMAGE_MODES:
                imaging.detect_outliers(
                    input_models,
                    self.save_intermediate_results,
                    self.good_bits,
                    self.maskpt,
                    snr1,
                    snr2,
                    scale1,
                    scale2,
                    self.backg,
                    self.resample_data,
                    self.weight_type,
                    self.pixfrac,
                    self.kernel,
                    self.fillval,
                    self.allowed_memory,
                    self.in_memory,
                    asn_id,
                    self.make_output_path,
                )
            elif exptype in SLIT_SPEC_MODES:
                spec.detect_outliers(
                    input_models,
                    self.save_intermediate_results,
                    self.good_bits,
                    self.maskpt,
                    snr1,
                    snr2,
                    scale1,
                    scale2,
                    self.backg,
                    self.resample_data,
                    self.weight_type,
                    self.pixfrac,
                    self.kernel,
                    self.fillval,
                    self.in_memory,
                    asn_id,
                    self.make_output_path,
                )
            elif exptype in IFU_SPEC_MODES:
                ifu.detect_outliers(
                    input_models,
                    self.save_intermediate_results,
                    self.kernel_size,
                    self.ifu_second_check,
                    self.threshold_percent,
                    self.make_output_path,
                )
            else:
                self.log.error("Outlier detection failed for unknown/unsupported ",
                               f"exposure type: {exptype}")
                return self._set_status(input_models, "SKIPPED")

            return self._set_status(input_models, "COMPLETE")

    def _set_status(self, input_models, status):
        if isinstance(input_models, ModelContainer):
            for model in input_models:
                model.meta.cal_step.outlier_detection = status
        else:
            input_models.meta.cal_step.outlier_detection = status
        return input_models
