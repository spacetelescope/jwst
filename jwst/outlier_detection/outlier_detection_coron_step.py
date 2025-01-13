"""
Submodule for performing outlier detection on coronagraphy data.
"""

import logging

import numpy as np

from stdatamodels.jwst import datamodels
from jwst.stpipe import Step

from jwst.resample.resample_utils import build_mask

from .utils import create_cube_median, flag_model_crs, OutlierDetectionStepBase
from ._fileio import save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionCoronStep"]


class OutlierDetectionCoronStep(Step, OutlierDetectionStepBase):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.
    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.
    Parameters
    -----------
    input_model : ~jwst.datamodels.CubeModel
        CubeModel or filename pointing to a CubeModel
    """

    class_alias = "outlier_detection_coron"

    spec = """
        maskpt = float(default=0.7)
        snr = float(default=5.0)
        save_intermediate_results = boolean(default=False)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
        suffix = string(default="crfints")
    """

    def process(self, input_model):
        """Perform outlier detection processing on input data."""

        # determine the asn_id (if not set by the pipeline)
        asn_id = self._get_asn_id(input_model)
        self.log.info(f"Outlier Detection asn_id: {asn_id}")

        if not isinstance(input_model, datamodels.JwstDataModel):
            input_model = datamodels.open(input_model)

        if not isinstance(input_model, datamodels.CubeModel):
            raise TypeError(f"Input must be a CubeModel: {input_model}")

        # FIXME weight_type could now be used here. Similar to tso data coron
        # data was previously losing var_rnoise due to the conversion from a cube
        # to a ModelContainer (which makes the default ivm weight ignore var_rnoise).
        # Now that it's handled as a cube we could use the var_rnoise.
        input_model.wht = build_mask(input_model.dq, self.good_bits).astype(np.float32)

        # Perform median combination on set of drizzled mosaics
        median_data = create_cube_median(input_model, self.maskpt)

        if self.save_intermediate_results:
            # make a median model
            median_model = datamodels.ImageModel(median_data)
            median_model.update(input_model)
            median_model.meta.wcs = input_model.meta.wcs

            save_median(median_model, self.make_output_path)
            del median_model

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        flag_model_crs(
            input_model,
            median_data,
            self.snr,
        )
        return self._set_status(input_model, True)
