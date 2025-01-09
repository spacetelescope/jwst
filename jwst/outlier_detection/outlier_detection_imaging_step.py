"""
Submodule for performing outlier detection on imaging data.
"""

import logging

from jwst.datamodels import ModelLibrary
from jwst.resample import resample
from jwst.stpipe.utilities import record_step_status
from jwst.stpipe import Step

from .utils import (flag_model_crs,
                    flag_resampled_model_crs,
                    median_without_resampling,
                    median_with_resampling,
                    OutlierDetectionStepBase)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionImagingStep"]


class OutlierDetectionImagingStep(Step, OutlierDetectionStepBase):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.
    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.
    Parameters
    -----------
    input_data : asn file or ~jwst.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.
    """

    class_alias = "outlier_detection_imaging"

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
        in_memory = boolean(default=False)
    """

    def process(self, input_models):
        """Perform outlier detection processing on input data."""

        # determine the asn_id (if not set by the pipeline)
        asn_id = self._get_asn_id(input_models)
        self.log.info(f"Outlier Detection asn_id: {asn_id}")

        snr1, snr2 = [float(v) for v in self.snr.split()]
        scale1, scale2 = [float(v) for v in self.scale.split()]

        if not isinstance(input_models, ModelLibrary):
            input_models = ModelLibrary(input_models, on_disk=not self.in_memory)

        if len(input_models) < 2:
            log.warning(f"Input only contains {len(input_models)} exposures")
            log.warning("Outlier detection will be skipped")
            record_step_status(input_models, "outlier_detection", False)
            return input_models
            
        if self.resample_data:
            resamp = resample.ResampleData(
                input_models,
                single=True,
                blendheaders=False,
                wht_type=self.weight_type,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=self.fillval,
                good_bits=self.good_bits,
            )
            median_data, median_wcs = median_with_resampling(
                input_models,
                resamp,
                self.maskpt,
                save_intermediate_results=self.save_intermediate_results,
                make_output_path=self.make_output_path,
            )
        else:
            median_data, median_wcs = median_without_resampling(
                input_models,
                self.maskpt,
                self.weight_type,
                self.good_bits,
                save_intermediate_results=self.save_intermediate_results,
                make_output_path=self.make_output_path,
            )


        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        with input_models:
            for image in input_models:
                if self.resample_data:
                    flag_resampled_model_crs(
                        image,
                        median_data,
                        median_wcs,
                        snr1,
                        snr2,
                        scale1,
                        scale2,
                        self.backg,
                        save_blot=self.save_intermediate_results,
                        make_output_path=self.make_output_path,
                    )
                else:
                    flag_model_crs(image, median_data, snr1)
                input_models.shelve(image, modify=True)

        return self._set_status(input_models, True)
