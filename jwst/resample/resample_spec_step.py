from ..stpipe import Step
from .. import datamodels
from . import resample_spec, ResampleStep
from ..exp_to_source import multislit_to_container
from ..assign_wcs.util import update_s_region


__all__ = ["ResampleSpecStep"]


class ResampleSpecStep(ResampleStep):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : DataModel, Association
    """

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

        if ref_filename != 'N/A':
            self.log.info('Drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:
            # Deal with NIRSpec which currently has no default drizpars reffile
            kwargs = self._set_spec_defaults()

        # Multislits get converted to a ModelContainer per slit
        if all([isinstance(i, datamodels.MultiSlitModel) for i in input_models]):
            container_dict = multislit_to_container(input_models)
            output = datamodels.MultiProductModel()
            output.update(input_models[0])
            for k, v in container_dict.items():
                input_models = v

                # Set up the resampling object as part of this step
                resamp = resample_spec.ResampleSpecData(input_models, **kwargs)
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
                ref_filename, single=self.single, weight_type=self.weight_type,
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
