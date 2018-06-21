from ..stpipe import Step, cmdline
from .. import datamodels
from . import resample_spec
from ..exp_to_source import multislit_to_container
from ..assign_wcs.util import update_s_region


class ResampleSpecStep(Step):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : DataModel, Association
    """

    spec = """
        single = boolean(default=False)
        wht_type = option('exptime', 'error', None, default='exptime')
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

        # Multislits get converted to a ModelContainer per slit
        if all([isinstance(i, datamodels.MultiSlitModel) for i in input_models]):
            container_dict = multislit_to_container(input_models)
            output = datamodels.MultiProductModel()
            output.update(input_models[0])
            for k, v in container_dict.items():
                input_models = v

                # Set up the resampling object as part of this step
                resamp = resample_spec.ResampleSpecData(input_models,
                    ref_filename, single=self.single,
                    wht_type=self.wht_type, pixfrac=self.pixfrac,
                    kernel=self.kernel, fillval=self.fillval,
                    good_bits=self.good_bits)
                # Do the resampling
                resamp.do_drizzle()

                for model in resamp.output_models:
                    update_s_region(model)

                if len(resamp.output_models) == 1:
                    out_slit = resamp.output_models[0]
                    output.products.append(out_slit)
                    output.products[-1].bunit_data = input_models[0].meta.bunit_data
                else:
                    out_slit = resamp.output_models
            result = output
            result.meta.cal_step.resample = "COMPLETE"
            result.meta.asn.pool_name = input_models.meta.pool_name
            result.meta.asn.table_name = input_models.meta.table_name
        else:
            # Set up the resampling object as part of this step
            resamp = resample_spec.ResampleSpecData(input_models,
                ref_filename, single=self.single, wht_type=self.wht_type,
                pixfrac=self.pixfrac, kernel=self.kernel,
                fillval=self.fillval, good_bits=self.good_bits)
            # Do the resampling
            resamp.do_drizzle()

            for model in resamp.output_models:
                model.meta.cal_step.resample = "COMPLETE"
                update_s_region(model)
                model.meta.asn.pool_name = input_models.meta.pool_name
                model.meta.asn.table_name = input_models.meta.table_name

            # Return either the single resampled datamodel, or the container
            # of datamodels.
            if len(resamp.output_models) == 1:
                result = resamp.output_models[0]
            else:
                result = resamp.output_models

        return result

