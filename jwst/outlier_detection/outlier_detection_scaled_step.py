"""Step interface for performing scaled outlier_detection."""

from ..stpipe import Step
from .. import datamodels
from . import outlier_detection_scaled


class OutlierDetectionScaledStep(Step):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input : asn file or ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

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
        good_bits = integer(default=6)
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

    def _build_reffile_container(self, reftype):
        """Return a ModelContainer of reference file models.

        Parameters
        ----------
        input_models: ModelContainer
            the science data, ImageModels in a ModelContainer

        reftype: string
            type of reference file

        Returns
        -------
        a ModelContainer with corresponding reference files for
            each input model

        """
        reffile_to_model = {'gain': datamodels.GainModel,
                            'readnoise': datamodels.ReadnoiseModel}
        reffile_model = reffile_to_model[reftype]

        reffiles = [self.input_models.meta.ref_file.instance[reftype]['name']]

        self.log.debug("Using {} reffile(s):".format(reftype.upper()))
        for r in set(reffiles):
            self.log.debug("    {}".format(r))

        # Use get_reference_file method to insure latest reference file
        # always gets used...especially since only one name will ever be needed
        ref_list = [reffile_model(self.get_reference_file(
                                  self.input_models, reftype))]

        return datamodels.ModelContainer(ref_list)
