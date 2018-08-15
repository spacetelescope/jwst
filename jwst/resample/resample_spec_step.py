from .. import datamodels
from ..datamodels import MultiSlitModel, ModelContainer
from . import resample_spec, ResampleStep
from ..exp_to_source import multislit_to_container
from ..assign_wcs.util import update_s_region_spectral


__all__ = ["ResampleSpecStep"]


class ResampleSpecStep(ResampleStep):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : `~jwst.datamodels.MultSlitModel`, `~jwst.datamodels.ModelContainer`, Association
        A singe datamodel, a container of datamodels, or an association file
    """

    def process(self, input):
        input = datamodels.open(input)

        # If single DataModel input, wrap in a ModelContainer
        if not isinstance(input, ModelContainer):
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
            self.log.info("No NIRSpec DIRZPARS reffile")
            kwargs = self._set_spec_defaults()

        self.drizpars = kwargs

        if isinstance(input_models[0], MultiSlitModel):
            # result is a MultiProductModel
            result = self._process_multislit(input_models)
        else:
            # result is a DrizProductModel
            result = self._process_slit(input_models)
        return result


    def _process_multislit(self, input_models):
        """
        Resample MultiSlit data

        Parameters
        ----------
        input : `~jwst.datamodels.ModelContainer`
            A container of `~jwst.datamodels.MultiSlitModel`

        Returns
        -------
        result : `~jwst.datamodels.MultiProductModel`
            The resampled output, one per source
        """
        containers = multislit_to_container(input_models)
        result = datamodels.MultiProductModel()
        result.update(input_models[0])
        for container in containers.values():
            resamp = resample_spec.ResampleSpecData(container, **self.drizpars)
            drizzled_models = resamp.do_drizzle()

            for model in drizzled_models:
                model.meta.cal_step.resample = "COMPLETE"
                model.meta.asn.pool_name = input_models.meta.pool_name
                model.meta.asn.table_name = input_models.meta.table_name
                update_s_region_spectral(model)

            # Everything resampled to single output model
            if len(drizzled_models) == 1:
                result.products.append(drizzled_models[0])
                result.products[-1].bunit_data = container[0].meta.bunit_data
            else:
                # For singlely-resampled images
                for model in drizzled_models:
                    result.products.append(model)
                    result.products[-1].bunit_data = container[0].meta.bunit_data

        return result


    def _process_slit(self, input_models):
        """
        Resample Slit data

        Parameters
        ----------
        input : `~jwst.datamodels.ModelContainer`
            A container of `~jwst.datamodels.ImageModel`
            or `~jwst.datamodels.SlitModel`

        Returns
        -------
        result : `~jwst.datamodels.DrizProductModel`
            The resampled output, one per source
        """
        resamp = resample_spec.ResampleSpecData(input_models, **self.drizpars)
        drizzled_models = resamp.do_drizzle()

        for model in drizzled_models:
            model.meta.cal_step.resample = "COMPLETE"
            model.meta.asn.pool_name = input_models.meta.pool_name
            model.meta.asn.table_name = input_models.meta.table_name
            update_s_region_spectral(model)

        # Return either the single resampled datamodel, or the container
        # of datamodels.
        if len(drizzled_models) == 1:
            result = drizzled_models[0]
        else:
            result = drizzled_models

        return result
