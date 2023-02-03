"""Step interface for performing scaled outlier_detection."""
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import outlier_detection_scaled


class OutlierDetectionScaledStep(Step):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input : asn file or ~jwst.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

    class_alias = "outlier_detection_scaled"

    spec = """
        weight_type = option('ivm','exptime',default='ivm')
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
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
    """

    def process(self, input):
        """Step interface to running outlier_detection."""
        with datamodels.open(input) as input_models:

            if not isinstance(input_models, datamodels.CubeModel):
                self.log.warning("Input is not a CubeModel.")
                self.log.warning("Outlier detection step will be skipped.")
                result = input_models.copy()
                result.meta.cal_step.outlier_detection = "SKIPPED"
                return result

            self.input_models = input_models

            reffiles = {}

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
                'good_bits': self.good_bits
            }

            # Setup for creating file names
            pars['make_output_path'] = self.make_output_path

            # Set up outlier detection, then do detection
            step = outlier_detection_scaled.OutlierDetectionScaled(
                self.input_models,
                reffiles=reffiles,
                **pars)
            step.do_detection()

            self.input_models.meta.cal_step.outlier_detection = 'COMPLETE'

            return self.input_models
