"""Step defined for outlier detection for stacked observations."""

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from ..stpipe import Step
from . import outlier_detection


class OutlierDetectionStackStep(Step):
    """Class definition for stacked outlier detection.

    Flag outlier bad pixels and cosmic rays in the DQ array of each input image
    of a stack of exposures, which in the case of TSO data are from the same
    data cube.

    Input images can listed in an input association file or already opened
    with a ModelContainer.

    DQ arrays are modified in place.

    By default, resampling has been disabled.  The 'resample_data' attribute
    can be reset to 'True' to turn on resampling if desired for the data.

    Parameters
    -----------
    input : asn file or ~jwst.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

    class_alias = "outlier_detection_stack"

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
        resample_data = boolean(default=False)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
    """

    def process(self, input):
        """Step interface for performing outlier_detection processing."""
        with datamodels.open(input) as input_models:

            if not isinstance(input_models, ModelContainer):
                self.log.warning("Input is not a ModelContainer.")
                self.log.warning("Outlier detection stack step will \
                                  be skipped.")
                result = input_models.copy()
                result.meta.cal_step.outlier_detection = "SKIPPED"
                return result

            self.log.info("Performing outlier detection on stack of \
                          {} inputs".format(len(input_models)))
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
                'resample_data': self.resample_data,
                'good_bits': self.good_bits
            }

            # Set up outlier detection, then do detection
            step = outlier_detection.OutlierDetection(self.input_models,
                                                      reffiles=reffiles,
                                                      **pars)
            step.do_detection()

            for model in self.input_models:
                model.meta.cal_step.outlier_detection = 'COMPLETE'

            return self.input_models
