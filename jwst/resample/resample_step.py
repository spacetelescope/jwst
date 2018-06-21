from ..stpipe import Step
from .. import datamodels
from . import resample
from ..assign_wcs.util import update_s_region


class ResampleStep(Step):
    """
    ResampleStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : DataModel or Association
        Single filename for either a single image or an association table.
    """

    spec = """
        single = boolean(default=False)
        wht_type = option('exptime','error',None,default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        good_bits = integer(default=4)
        blendheaders = boolean(default=True)
    """
    reference_file_types = ['drizpars']

    def process(self, input):

        input = datamodels.open(input)

        # If single input, wrap in a ModelContainer
        if not isinstance(input, datamodels.ModelContainer):
            input_models = datamodels.ModelContainer([input])
            input_models.meta.resample.output = input.meta.filename
            self.blendheaders = False
        else:
            input_models = input

        for reftype in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], reftype)

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, ref_filename=ref_filename,
            single=self.single, wht_type=self.wht_type, pixfrac=self.pixfrac,
            kernel=self.kernel, fillval=self.fillval, good_bits=self.good_bits,
            blendheaders=self.blendheaders)
        resamp.do_drizzle()

        for model in resamp.output_models:
            model.meta.cal_step.resample = "COMPLETE"
            update_s_region(model)
            model.meta.asn.pool_name = input_models.meta.pool_name
            model.meta.asn.table_name = input_models.meta.table_name


        if len(resamp.output_models) == 1:
            result = resamp.output_models[0]
        else:
            result = resamp.output_models

        return result
