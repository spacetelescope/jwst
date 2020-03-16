__all__ = ["ResampleSpecStep"]

from .. import datamodels
from ..datamodels import MultiSlitModel, ModelContainer, ImageModel
from . import resample_spec, ResampleStep
from ..exp_to_source import multislit_to_container
from ..assign_wcs.util import update_s_region_spectral


class ResampleSpecStep(ResampleStep):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : `~jwst.datamodels.MultiSlitModel`, `~jwst.datamodels.ModelContainer`, Association
        A singe datamodel, a container of datamodels, or an association file
    """

    # Spec is all the same except for the suffix
    spec = """
    """

    def process(self, input):
        input = datamodels.open(input)

        if isinstance(input, ImageModel):
            slit_model = datamodels.SlitModel()
            input.meta.bunit_err = None  # this prevents error array populated in slit model
            slit_model.update(input)
            slit_model.data = input.data
            input = slit_model
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
            result = self._process_multislit(input_models)

        elif len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError("Input {} is not a 2D image.".format(input_models[0]))
        else:
            # result is a SlitModel
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
        result : `~jwst.datamodels.MultiSlitModel`
            The resampled output, one per source
        """
        containers = multislit_to_container(input_models)
        result = datamodels.MultiSlitModel()
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
                result.slits.append(drizzled_models[0])
                result.slits[-1].bunit_data = container[0].meta.bunit_data
            else:
                # When each input is resampled to its own output
                for model in drizzled_models:
                    result.slits.append(model)
                    result.slits[-1].bunit_data = container[0].meta.bunit_data
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
        result : `~jwst.datamodels.ImageModel`
            The resampled output, one per source
        """

#        result = datamodels.SlitModel()
#        result.update(input_models[0])
        resamp = resample_spec.ResampleSpecData(input_models, **self.drizpars)
        drizzled_models = resamp.do_drizzle()
        
#        drizzled_models[0].update(input_models[0]) #creates a empty err array
#        drizzled_models[0].update(input_models[0],only='PRIMARY')
#        drizzled_models[0].update(input_models[0],only='SCI')

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
