#! /usr/bin/env python

from ..stpipe import Step, cmdline
import ..datamodels
from . import resample


class ResampleStep(Step):
    """
    ResampleStep: Uses the drizzle process to resample (geometric correction)
    a single 2D image or resample and combine a set of 2D images specified as
    an association.

    Parameters
    -----------
    input : str or model
        Single filename for either a single image or an association table.  This
        would then get used to create a AsnModel(?) object for this step.
    """

    spec = """
        single = boolean(default=False)
        wht_type = option('exptime','error',None,default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        good_bits = integer(default=-1)
    """
    reference_file_types=['drizpars']

    def process(self, input):

        input_models = datamodels.open(input)
        if type(input_models) != type(datamodels.ModelContainer()): # single exposure
            s = datamodels.ModelContainer()
            s.append(input_models)
            s.assign_group_ids()
            s.meta.resample.output = s.group_names[0].replace('_resamp','_drz')
            
            self.input_models = s
        else:
            self.input_models = input_models
        
        # identify what reference file has been associated with these inputs
        try:
            self.ref_filename = self.get_reference_file(self.input_models[0], 'drizpars' )
        except:
            # This is only in place for initial testing of this code, prior
            # to ref file being included in CRDS
            self.ref_filename = self.input_models[0].meta.storage.get_fits_header('PRIMARY')['r_resamp']
        self.log.info('Reference file to use: %s ', self.ref_filename)

        # Call the resampling routine
        self.step = resample.ResampleData(self.input_models, self.ref_filename,
                                single=self.single, wht_type=self.wht_type,
                                pixfrac=self.pixfrac, kernel=self.kernel,
                                fillval=self.fillval, good_bits=self.good_bits)
        self.step.do_drizzle()

        #self.input_models.close()
        
        if len(self.step.output_models) == 1:
            output_model = self.step.output_models[0]
        else:
            output_model = self.step.output_models
        return output_model


if __name__ == '__main__':
    cmdline.step_script(ResampleStep)
