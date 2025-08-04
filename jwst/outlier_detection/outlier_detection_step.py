"""Public common step definition for OutlierDetection processing."""

from functools import partial

from stdatamodels import filetype
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.lib.pipe_utils import is_tso
from jwst.outlier_detection import coron, ifu, imaging, spec, tso
from jwst.stpipe import Step
from jwst.stpipe.utilities import query_step_status, record_step_status

# Categorize all supported modes
IMAGE_MODES = ["NRC_IMAGE", "MIR_IMAGE", "NRS_IMAGE", "NIS_IMAGE", "FGS_IMAGE"]
SLIT_SPEC_MODES = ["NRC_WFSS", "MIR_LRS-FIXEDSLIT", "NRS_FIXEDSLIT", "NRS_MSASPEC", "NIS_WFSS"]
TSO_SPEC_MODES = ["NIS_SOSS", "MIR_LRS-SLITLESS", "NRC_TSGRISM", "NRS_BRIGHTOBJ"]
IFU_SPEC_MODES = ["NRS_IFU", "MIR_MRS"]
TSO_IMAGE_MODES = ["NRC_TSIMAGE"]  # missing MIR_IMAGE with TSOVIST=True, not really addable
CORON_IMAGE_MODES = ["NRC_CORON", "MIR_LYOT", "MIR_4QPM"]

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(Step):
    """
    Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or dictionary,
    or already opened with a ModelContainer or ModelLibrary.
    DQ arrays are modified in place.
    SCI, ERR, VAR_RNOISE, VAR_FLAT, and VAR_POISSON arrays are updated with
    NaN values matching the DQ flags.
    """

    class_alias = "outlier_detection"

    spec = """
        weight_type = option('ivm','exptime',default='ivm')
        pixfrac = float(min=0.0, max=1.0, default=1.0)  # Pixel shrinkage factor
        kernel = option('square','point','turbo',default='square')  # Flux distribution kernel
        fillval = string(default='NAN')
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
        in_memory = boolean(default=True) # in_memory flag ignored if run within the pipeline; set at pipeline level instead
    """  # noqa: E501

    def process(self, input_data):
        """
        Perform outlier detection processing on input data.

        Parameters
        ----------
        input_data : asn file, ~jwst.datamodels.ModelContainer, or ~jwst.datamodels.ModelLibrary
            The input association.
            For imaging modes a ModelLibrary is expected, whereas for spectroscopic modes a
            ModelContainer is expected.

        Returns
        -------
        result_models : ~jwst.datamodels.ModelContainer or ~jwst.datamodels.ModelLibrary
            The modified input data with DQ flags set for detected outliers.
        """
        # Open the input data, making a copy as needed.
        result_models = self._open_models(input_data)

        # determine the "mode" (if not set by the pipeline)
        mode = self._guess_mode(input_data)
        if mode is None:
            record_step_status(result_models, "outlier_detection", False)
            return result_models

        self.log.info(f"Outlier Detection mode: {mode}")

        # determine the asn_id (if not set by the pipeline)
        self._get_asn_id(input_data)

        snr1, snr2 = [float(v) for v in self.snr.split()]
        scale1, scale2 = [float(v) for v in self.scale.split()]

        if mode == "tso":
            result_models = tso.detect_outliers(
                result_models,
                self.save_intermediate_results,
                self.good_bits,
                self.maskpt,
                self.rolling_window_width,
                snr1,
                self.make_output_path,
            )
        elif mode == "coron":
            result_models = coron.detect_outliers(
                result_models,
                self.save_intermediate_results,
                self.good_bits,
                self.maskpt,
                snr1,
                self.make_output_path,
            )
        elif mode == "imaging":
            result_models = imaging.detect_outliers(
                result_models,
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
                self.make_output_path,
            )
        elif mode == "spec":
            result_models = spec.detect_outliers(
                result_models,
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
                self.make_output_path,
            )
        elif mode == "ifu":
            result_models = ifu.detect_outliers(
                result_models,
                self.save_intermediate_results,
                self.kernel_size,
                self.ifu_second_check,
                self.threshold_percent,
                self.make_output_path,
            )
        else:
            self.log.error(f"Outlier detection failed for unknown/unsupported mode: {mode}")
            record_step_status(result_models, "outlier_detection", False)

        if query_step_status(result_models, "outlier_detection") != "SKIPPED":
            record_step_status(result_models, "outlier_detection", True)
        return result_models

    def _guess_mode(self, input_models):
        # The pipelines should set this mode or ideally these should
        # be separate steps (but that would require new crds reference files).
        if hasattr(self, "mode"):
            return self.mode

        # guess mode from input type
        if isinstance(input_models, str):
            if filetype.check(input_models) == "asn":
                input_models = datamodels.open(input_models, asn_n_members=1)
            else:
                input_models = datamodels.open(input_models)
        elif isinstance(input_models, (dict, list)):
            input_models = datamodels.open(input_models, asn_n_members=1)

        # Select which version of OutlierDetection
        # needs to be used depending on the input data
        if isinstance(input_models, ModelContainer):
            exptype = input_models[0].meta.exposure.type
        elif isinstance(input_models, ModelLibrary):
            with input_models:
                single_meta = input_models.read_metadata(0)
            exptype = single_meta["meta.exposure.type"]
        else:
            exptype = input_models.meta.exposure.type
            # only need to check for TSO type here because container and library are not
            # expected inputs for TSO type
            if is_tso(input_models):
                return "tso"

        if exptype in CORON_IMAGE_MODES:
            return "coron"
        if exptype in IMAGE_MODES:
            return "imaging"
        if exptype in SLIT_SPEC_MODES:
            return "spec"
        if exptype in IFU_SPEC_MODES:
            return "ifu"

        self.log.error(f"Outlier detection failed for unknown/unsupported exposure type: {exptype}")
        return None

    def _get_asn_id(self, input_models):
        """Update make_output_path to include the association ID in the output path."""
        # handle if input_models isn't open
        if isinstance(input_models, (str, dict)):
            input_models = datamodels.open(input_models, asn_n_members=1)

        # Setup output path naming if associations are involved.
        try:
            if isinstance(input_models, ModelLibrary):
                asn_id = input_models.asn["asn_id"]
            elif isinstance(input_models, ModelContainer):
                asn_id = input_models.asn_table["asn_id"]
            else:
                asn_id = input_models.meta.asn_table.asn_id
        except (AttributeError, KeyError):
            asn_id = None

        if asn_id is None:
            asn_id = self.search_attr("asn_id")
        if asn_id is not None:
            _make_output_path = self.search_attr("_make_output_path", parent_first=True)

            self._make_output_path = partial(_make_output_path, asn_id=asn_id)
        self.log.info(f"Outlier Detection asn_id: {asn_id}")
        return

    def _open_models(self, input_models):
        """
        Open the input data, making a copy if necessary.

        If the input data is a filename or path, it is opened
        and the open model is returned.

        If it is a list of models, it is opened as a ModelContainer.
        In this case, or if the input is a simple datamodel or a
        ModelContainer, a deep copy of the model/container is returned,
        in order to avoid modifying the input models.

        If the input is a ModelLibrary, it is simply returned, in order
        to avoid making unnecessary copies for performance-critical
        use cases.

        Parameters
        ----------
        input_models : str, list, JwstDataModel, ModelContainer, or ModelLibrary
            Input data to open.

        Returns
        -------
        JwstDataModel, ModelContainer, or ModelLibrary
            The opened datamodel(s).
        """
        # Check whether input contains datamodels
        make_copy = False
        if isinstance(input_models, list):
            is_datamodel = [isinstance(m, datamodels.JwstDataModel) for m in input_models]
            if any(is_datamodel):
                make_copy = True
        elif isinstance(input_models, (datamodels.JwstDataModel, ModelContainer)):
            make_copy = True

        if not isinstance(input_models, (datamodels.JwstDataModel, ModelLibrary, ModelContainer)):
            # Input might be a filename or path.
            # In that case, open it.
            input_models = datamodels.open(input_models)
        if make_copy:
            # For regular models, make a copy to avoid modifying the input.
            # Leave libraries alone for memory management reasons.
            input_models = input_models.copy()
        return input_models
