from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import outlier_detection


class OutlierDetectionStep(Step):
    """
    OutlierDetectionStep: Flag outlier pixels (bad pixels, cosmic-ray hits,...)
    in the DQ arrays of each input image of an ASN.  This step relies on the
    same basic algorithm (and even some of the core code) from DrizzlePac.

    Parameters
    -----------
    input : str or model
        Single filename for either a single image or an association table.  This
        would then get used to create a AsnModel(?) object for this step.
    """

    spec = """
        wht_type = option('exptime','error',None,default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        nlow = integer(default=0)
        nhigh = integer(default=0)
        hthresh = float(default=-1.0)
        lthresh = float(default=-1.0)
        nsigma = string(default='4.0 3.0')
        maskpt = float(default=0.7)
        grow = integer(default=1)
        ctegrow = integer(default=0)
        snr = string(default='4.0 3.0')
        scale = string(default='0.5 0.4')
        backg = float(default=0.0)
    """
    reference_file_types = ['gain', 'readnoise']

    def process(self, input, to_file=False):

        with datamodels.open(input) as input_models:

            if not isinstance(input_models, datamodels.ModelContainer):
                self.log.warning("Input is not a ModelContainer.")
                self.log.warning("Outlier detection step will be skipped.")
                result = input_models.copy()
                result.meta.cal_step.outlier_detection = "SKIPPED"
                return result

            self.input_models = input_models

            reffiles= {}
            reffiles['gain'] = self._build_reffile_container('gain')
            reffiles['readnoise'] = self._build_reffile_container('readnoise')

            pars = {
                'wht_type': self.wht_type,
                'pixfrac': self.pixfrac,
                'kernel': self.kernel,
                'fillval': self.fillval,
                'nlow': self.nlow,
                'nhigh': self.nhigh,
                'hthresh': self.hthresh,
                'lthresh': self.lthresh,
                'nsigma': self.nsigma,
                'maskpt': self.maskpt,
                'grow': self.grow,
                'ctegrow': self.ctegrow,
                'snr': self.snr,
                'scale': self.scale,
                'backg': self.backg}

            # Set up outlier detection, then do detection
            step = outlier_detection.OutlierDetection(self.input_models,
                reffiles=reffiles, to_file=to_file, **pars)
            step.do_detection()

            return self.input_models

    def _build_reffile_container(self, reftype):
        """
        Return a ModelContainer of reference file models.

        Parameters
        ----------

        input_models: ModelContainer
            the science data, ImageModels in a ModelContainer

        reftype: string
            type of reference file

        Returns
        -------

        a ModelContainer with corresponding reference files for each input model
        """

        reffile_to_model = {'gain': datamodels.GainModel,
            'readnoise': datamodels.ReadnoiseModel}

        reffiles = [self.get_reference_file(im, reftype) for im in self.input_models]
        self.log.debug("Using {} reffile(s):".format(reftype.upper()))
        for r in set(reffiles):
            self.log.debug("    {}".format(r))

        # Check if all the ref files are the same.  If so build it by reading
        # the reference file just once.
        if len(set(reffiles)) <= 1:
            length = len(self.input_models)
            ref_list = [reffile_to_model.get(reftype)(reffiles[0])] * length
        else:
            ref_list = [reffile_to_model.get(reftype)(ref) for ref in reffiles]
        return datamodels.ModelContainer(ref_list)



if __name__ == '__main__':
    cmdline.step_script(OutlierDetectionStep)
