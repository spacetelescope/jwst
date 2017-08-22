from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import outlier_detection


class OutlierDetectionStep(Step):
    """
    Flag outlier bad pixels and cosmic rays in the DQ array of each input image

    Input images can listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input : asn file or ModelContainer
        Single filename association table, or a datamodels.ModelContainer.
    """

    spec = """
        wht_type = option('exptime','error',None,default='exptime')
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
        resample_data = boolean(default=True)
        good_bits = integer(default=4)
    """
    reference_file_types = ['gain', 'readnoise']

    def process(self, input):

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
                reffiles=reffiles, **pars)
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
        reffile_model = reffile_to_model[reftype]           

        reffiles = [im.meta.ref_file.instance[reftype]['name'] for im in self.input_models]
        self.log.debug("Using {} reffile(s):".format(reftype.upper()))
        for r in set(reffiles):
            self.log.debug("    {}".format(r))

        # Check if all the ref files are the same.  If so build it by reading
        # the reference file just once.
        if len(set(reffiles)) <= 1:
            length = len(self.input_models)
            ref_list = [reffile_model(self.reference_uri_to_cache_path(reffiles[0]))]*length
        else:
            ref_list = [reffile_model(self.reference_uri_to_cache_path(ref)) for ref in reffiles]
        return datamodels.ModelContainer(ref_list)



if __name__ == '__main__':
    cmdline.step_script(OutlierDetectionStep)
