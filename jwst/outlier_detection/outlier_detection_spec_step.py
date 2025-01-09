"""Perform outlier detection on spectra."""

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.utilities import record_step_status
from jwst.stpipe import Step

from ..resample import resample_spec
from .utils import (flag_crs_in_models,
                    flag_crs_in_models_with_resampling,
                    median_with_resampling,
                    median_without_resampling,
                    OutlierDetectionStepBase)

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutlierDetectionSpecStep"]

class OutlierDetectionSpecStep(Step, OutlierDetectionStepBase):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.
    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.
    Parameters
    -----------
    input_data : ~jwst.datamodels.ModelContainer or str
        ModelContainer or association file that can initialize a ModelContainer.
    """

    class_alias = "outlier_detection_spec"

    spec = """
        weight_type = option('ivm','exptime',default='ivm')
        pixfrac = float(default=1.0)
        kernel = string(default='square') # drizzle kernel
        fillval = string(default='INDEF')
        maskpt = float(default=0.7)
        snr = string(default='5.0 4.0')
        scale = string(default='1.2 0.7')
        backg = float(default=0.0)
        save_intermediate_results = boolean(default=False)
        resample_data = boolean(default=True)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
    """

    def process(self, input_models):
        """Perform outlier detection processing on input data."""

        # determine the asn_id (if not set by the pipeline)
        asn_id = self._get_asn_id(input_models)
        self.log.info(f"Outlier Detection asn_id: {asn_id}")

        snr1, snr2 = [float(v) for v in self.snr.split()]
        scale1, scale2 = [float(v) for v in self.scale.split()]

        if not isinstance(input_models, ModelContainer):
            input_models = ModelContainer(input_models)

        if len(input_models) < 2:
            log.warning(f"Input only contains {len(input_models)} exposures")
            log.warning("Outlier detection will be skipped")
            record_step_status(input_models, "outlier_detection", False)
            return input_models

        # convert to library for resample
        # for compatibility with image3 pipeline
        library = ModelLibrary(input_models, on_disk=False)

        if self.resample_data is True:
            # Start by creating resampled/mosaic images for
            #  each group of exposures
            resamp = resample_spec.ResampleSpecData(
                input_models,
                single=True,
                blendheaders=False,
                wht_type=self.weight_type,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=self.fillval,
                good_bits=self.good_bits,
            )

            median_data, median_wcs, median_err = median_with_resampling(
                library,
                resamp,
                self.maskpt,
                save_intermediate_results=self.save_intermediate_results,
                make_output_path=self.make_output_path,
                return_error=True,
            )
        else:
            median_data, median_wcs, median_err = median_without_resampling(
                library,
                self.maskpt,
                self.weight_type,
                self.good_bits,
                save_intermediate_results=self.save_intermediate_results,
                make_output_path=self.make_output_path,
                return_error=True
            )

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        if self.resample_data:
            flag_crs_in_models_with_resampling(
                input_models,
                median_data,
                median_wcs,
                snr1,
                snr2,
                scale1,
                scale2,
                self.backg,
                median_err=median_err,
                save_blot=self.save_intermediate_results,
                make_output_path=self.make_output_path
            )
        else:
            flag_crs_in_models(
                input_models,
                median_data,
                snr1,
                median_err=median_err,
            )
        return self._set_status(input_models, True)
