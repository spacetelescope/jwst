__all__ = ["ResampleSpecStep"]

from .. import datamodels
from ..datamodels import MultiSlitModel, ModelContainer, ImageModel
from . import resample_spec, ResampleStep
from ..exp_to_source import multislit_to_container
from ..assign_wcs.util import update_s_region_spectral


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = '~DO_NOT_USE+NON_SCIENCE'


class ResampleSpecStep(ResampleStep):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : `~jwst.datamodels.MultiSlitModel`, `~jwst.datamodels.ModelContainer`, Association
        A singe datamodel, a container of datamodels, or an association file
    """

    def process(self, input):

        # Define input_new, because if input is ImageModel, it will
        # get recreated as a SlitModel
        input_new = datamodels.open(input)

        if isinstance(input_new, ImageModel):
            slit_model = datamodels.SlitModel()
            slit_model.update(input_new, only="PRIMARY")
            slit_model.update(input_new, only="SCI")
            slit_model.meta.wcs = input_new.meta.wcs
            slit_model.data = input_new.data
            input_new = slit_model

        # If single DataModel input, wrap in a ModelContainer
        if not isinstance(input_new, ModelContainer):
            input_models = datamodels.ModelContainer([input_new])
            input_models.meta.resample.output = input_new.meta.filename
            self.blendheaders = False
        else:
            input_models = input_new

        # Get the drizpars reference file
        for reftype in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], reftype)

        if ref_filename != 'N/A':
            self.log.info('Drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:
            # Deal with NIRSpec, which currently has no default drizpars reffile
            self.log.info("No NIRSpec DIRZPARS reffile")
            kwargs = self._set_spec_defaults()
            kwargs['blendheaders'] = self.blendheaders

        # Update user-supplied kwargs
        kwargs['allowed_memory'] = self.allowed_memory

        # Call resampling
        self.drizpars = kwargs
        if isinstance(input_models[0], MultiSlitModel):
            result = self._process_multislit(input_models)

        elif len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError("Input {} is not a 2D image.".format(input_models[0]))

        else:
            # result is a SlitModel
            result = self._process_slit(input_models)

        # Update ASNTABLE in output
        result.meta.asn.table_name = input_models[0].meta.asn.table_name
        result.meta.filetype = 'resampled'

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

        result.update(input_models[0], only="PRIMARY")
        result.update(input_models[0], only="SCI")

        for container in containers.values():
            resamp = resample_spec.ResampleSpecData(container, **self.drizpars)
            drizzled_models = resamp.do_drizzle()

            for model in drizzled_models:
                model.meta.cal_step.resample = "COMPLETE"
                model.meta.asn.pool_name = input_models.meta.pool_name
                model.meta.asn.table_name = input_models.meta.table_name

                # Delete the BUNIT keyword for the ERR extension, so that datamodels
                # doesn't create an empty ERR extension (just for that keyword)
                if hasattr(model.meta, "bunit_err") and model.meta.bunit_err is not None:
                    del model.meta.bunit_err

                update_s_region_spectral(model)

            # Everything resampled to single output model
            if len(drizzled_models) == 1:
                result.slits.append(drizzled_models[0])
                if container[0].meta.bunit_data is not None:
                    result.slits[-1].meta.bunit_data = container[0].meta.bunit_data
            else:
                # When each input is resampled to its own output
                for model in drizzled_models:
                    result.slits.append(model)
                    result.slits[-1].meta.bunit_data = container[0].meta.bunit_data
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

        resamp = resample_spec.ResampleSpecData(input_models, **self.drizpars)

        if input_models[0].meta.exposure.type == "MIR_LRS-FIXEDSLIT":
            bb = input_models[0].meta.wcs.bounding_box
            ((x1, x2), (y1, y2)) = bb
            xmin = int(min(x1, x2))
            ymin = int(min(y1, y2))
            xmax = int(max(x1, x2))
            ymax = int(max(y1, y2))
            drizzled_models = resamp.do_drizzle(xmin=xmin, xmax=xmax,
                                                ymin=ymin, ymax=ymax)
        else:
            drizzled_models = resamp.do_drizzle()

        result = drizzled_models[0]
        result.meta.cal_step.resample = "COMPLETE"
        result.meta.asn.pool_name = input_models.meta.pool_name
        result.meta.asn.table_name = input_models.meta.table_name

        # Delete BUNIT keyword for ERR extension to prevent datamodels from
        # creating an empty ERR extension (just for the keyword)
        if hasattr(result.meta, "bunit_err") and result.meta.bunit_err is not None:
            del result.meta.bunit_err

        update_s_region_spectral(result)
        result.meta.bunit_data = drizzled_models[0].meta.bunit_data
        return result
