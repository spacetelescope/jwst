from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import outlier_detection

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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
        nsigma = string(default='4 3')
        maskpt = float(default=0.7)
        grow = integer(default=1)
        ctegrow = integer(default=0)
        snr = string(default='4.0 3.0')
        scale = string(default='0.5 0.4')
        backg = float(default=0.0)
    """
    reference_file_types = ['gain', 'readnoise'] # No ref file for Build6...

    def process(self, input, to_file=False):

        self.input_models = datamodels.open(input)

        self.ref_filename = {}
        self.ref_filename['gain'] = self.build_reffile_container('gain')
        self.ref_filename['readnoise'] = self.build_reffile_container('readnoise')

        # Call the resampling routine
        self.step = outlier_detection.OutlierDetection(self.input_models,
                                to_file=to_file,
                                ref_filename=self.ref_filename)
        self.step.do_detection()

        return self.input_models

    def build_reffile_container(self, reftype):
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
        reffiles = [self.get_reference_file(im, reftype) for im in self.input_models]

        # Check if all the ref files are the same.  If so build it by reading
        # the reference file just once.
        if len(set(reffiles)) <= 1:
            length = len(self.input_models)
            ref_list = [datamodels.open(reffiles[0])] * length
            #ref_list = [ref_model(reffiles[0])] * length
        else:
            ref_list = [datamodels.open(ref) for ref in reffiles]
            #ref_list = [ref_model(ref) for ref in reffiles]
        return datamodels.ModelContainer(ref_list) #, build_groups=False)



if __name__ == '__main__':
    cmdline.step_script(OutlierDetectionStep)
