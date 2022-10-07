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

    class_alias = "resample_spec"

    def process(self, input):
        self.wht_type = self.weight_type
        input_new = datamodels.open(input)

        # Convert ImageModel to SlitModel (needed for MIRI LRS)
        if isinstance(input_new, ImageModel):
            input_new = datamodels.SlitModel(input_new)

        if isinstance(input_new, ModelContainer):
            input_models = input_new
            try:
                output = input_models.meta.asn_table.products[0].name
            except AttributeError:
                # NIRSpec MOS data goes through this path, as the container
                # is only ModelContainer-like, and doesn't have an asn_table
                # attribute attached.  Output name handling gets done in
                # _process_multislit() via the update method
                # TODO: the container-like object should retain asn_table
                output = None
        else:
            input_models = datamodels.ModelContainer([input_new])
            output = input_new.meta.filename
            self.blendheaders = False

        # Get the drizpars reference file
        for reftype in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], reftype)

        if ref_filename != 'N/A':
            self.log.info('Drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:
            # Deal with NIRSpec, which currently has no default drizpars reffile
            self.log.info("No DRIZPARS reffile")
            kwargs = self._set_spec_defaults()
            kwargs['blendheaders'] = self.blendheaders

        kwargs['allowed_memory'] = self.allowed_memory
        kwargs['output'] = output

        # Issue a warning about the use of exptime weighting
        if self.wht_type == 'exptime':
            self.log.warning("Use of EXPTIME weighting will result in incorrect")
            self.log.warning("propagated errors in the resampled product")

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
        result.meta.asn.pool_name = input_models[0].meta.asn.pool_name

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
                self.update_slit_metadata(model)
                update_s_region_spectral(model)
                result.slits.append(model)

        result.meta.cal_step.resample = "COMPLETE"
        result.meta.asn.pool_name = input_models.asn_pool_name
        result.meta.asn.table_name = input_models.asn_table_name
        result.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
        result.meta.resample.pixfrac = self.pixfrac

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
        result : `~jwst.datamodels.SlitModel`
            The resampled output
        """

        resamp = resample_spec.ResampleSpecData(input_models, **self.drizpars)

        drizzled_models = resamp.do_drizzle()

        result = drizzled_models[0]
        result.meta.cal_step.resample = "COMPLETE"
        result.meta.asn.pool_name = input_models.asn_pool_name
        result.meta.asn.table_name = input_models.asn_table_name
        result.meta.bunit_data = drizzled_models[0].meta.bunit_data
        result.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
        result.meta.resample.pixfrac = self.pixfrac
        self.update_slit_metadata(result)
        update_s_region_spectral(result)

        return result

    def update_slit_metadata(self, model):
        """
        Update slit attributes in the resampled slit image.

        This is needed because model.slit attributes are not in model.meta, so
        the normal update() method doesn't work with them. Updates output_model
        in-place.
        """
        for attr in ['name', 'xstart', 'xsize', 'ystart', 'ysize',
                     'slitlet_id', 'source_id', 'source_name', 'source_alias',
                     'stellarity', 'source_type', 'source_xpos', 'source_ypos',
                     'dispersion_direction', 'shutter_state']:
            try:
                val = getattr(self.input_models[-1], attr)
            except AttributeError:
                pass
            else:
                if val is not None:
                    setattr(model, attr, val)
