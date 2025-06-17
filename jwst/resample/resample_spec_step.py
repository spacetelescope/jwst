__all__ = ["ResampleSpecStep"]

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import MultiSlitModel, ImageModel

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.lib.wcs_utils import get_wavelengths
from jwst.resample.resample_utils import load_custom_wcs, find_miri_lrs_sregion

from . import resample_spec, ResampleStep
from jwst.exp_to_source import multislit_to_container
from jwst.assign_wcs.util import update_s_region_spectral
from jwst.stpipe import Step


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = "~DO_NOT_USE+NON_SCIENCE"


class ResampleSpecStep(Step):
    """Resample spectral data onto a regular grid using the drizzle algorithm."""

    class_alias = "resample_spec"

    spec = """
        pixfrac = float(min=0.0, max=1.0, default=1.0)  # Pixel shrinkage factor
        kernel = option('square', 'point', default='square')  # Flux distribution kernel
        fillval = string(default='NAN')  # Output value for pixels with no weight or flux
        weight_type = option('ivm', 'exptime', None, default='ivm')  # Input image weighting type
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        pixel_scale_ratio = float(default=1.0)  # Ratio of input to output spatial pixel scale
        pixel_scale = float(default=None)  # Spatial pixel scale in arcsec
        output_wcs = string(default='')  # Custom output WCS
        single = boolean(default=False)  # Resample each input to its own output grid
        blendheaders = boolean(default=True)  # Blend metadata from inputs into output
        in_memory = boolean(default=True)  # Keep images in memory
    """  # noqa: E501

    def process(self, input_data):
        """
        Run the resample step on the input data.

        Parameters
        ----------
        input_data : MultiSlitModel, ModelContainer, str
            A single datamodel, a container of datamodels, or an association file.

        Returns
        -------
        SlitModel or MultiSlitModel
            The resampled output, one slit per source.
        """
        input_new = datamodels.open(input_data)

        # Check if input_new is a MultiSlitModel
        model_is_msm = isinstance(input_new, MultiSlitModel)

        #  If input is a 3D rateints MultiSlitModel (unsupported) skip the step
        if model_is_msm and len((input_new[0]).shape) == 3:
            self.log.warning("Resample spec step will be skipped")
            input_new.meta.cal_step.resample = "SKIPPED"

            return input_new

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
            input_models = ModelContainer([input_new])
            output = input_new.meta.filename
            self.blendheaders = False

        # Setup drizzle-related parameters
        kwargs = self.get_drizpars()
        kwargs["output"] = output
        self.drizpars = kwargs

        # Call resampling
        if isinstance(input_models[0], MultiSlitModel):
            result = self._process_multislit(input_models)

        elif len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise TypeError(f"Input {input_models[0]} is not a 2D image.")

        else:
            # result is a SlitModel
            result = self._process_slit(input_models)

        # Update ASNTABLE in output
        result.meta.cal_step.resample = "COMPLETE"
        result.meta.asn.table_name = input_models[0].meta.asn.table_name
        result.meta.asn.pool_name = input_models[0].meta.asn.pool_name

        # populate the result wavelength attribute for MultiSlitModel
        if isinstance(result, MultiSlitModel):
            for slit_idx, _slit in enumerate(result.slits):
                wl_array = get_wavelengths(result.slits[slit_idx])
                result.slits[slit_idx].wavelength = wl_array
        else:
            # populate the result wavelength attribute for SlitModel
            wl_array = get_wavelengths(result)
            result.wavelength = wl_array

        return result

    def _process_multislit(self, input_models):
        """
        Resample MultiSlit data.

        Parameters
        ----------
        input_models : `~jwst.datamodels.ModelContainer`
            A container of `~jwst.datamodels.MultiSlitModel`

        Returns
        -------
        result : `~jwst.datamodels.MultiSlitModel`
            The resampled output, one per source
        """
        containers = multislit_to_container(input_models)
        result = datamodels.MultiSlitModel()

        result.update(input_models[0])

        pscale_ratio = None
        for container in containers.values():
            # Make sure all input models have consistent NaN and DO_NOT_USE values
            for model in container:
                match_nans_and_flags(model)

            # Call the resampling routine
            if self.single:
                resamp = resample_spec.ResampleSpec(
                    container, enable_var=False, compute_err="driz_err", **self.drizpars
                )
                drizzled_library = resamp.resample_many_to_many(in_memory=self.in_memory)
            else:
                resamp = resample_spec.ResampleSpec(
                    container, enable_var=True, compute_err="from_var", **self.drizpars
                )
                drizzled_library = resamp.resample_many_to_one()
                drizzled_library = ModelLibrary(
                    [
                        drizzled_library,
                    ],
                    on_disk=False,
                )

            with drizzled_library:
                for i, model in enumerate(drizzled_library):
                    self.update_slit_metadata(model)
                    update_s_region_spectral(model)
                    result.slits.append(model)
                    drizzled_library.shelve(model, i, modify=False)
            del drizzled_library

            # Keep the first computed pixel scale ratio for storage
            if self.pixel_scale is not None and pscale_ratio is None:
                pscale_ratio = resamp.pixel_scale_ratio

        if self.pixel_scale is None or pscale_ratio is None:
            result.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
        else:
            result.meta.resample.pixel_scale_ratio = pscale_ratio
        result.meta.resample.pixfrac = self.pixfrac

        return result

    def get_drizpars(self):
        """
        Load all drizzle-related parameter values into kwargs list.

        Returns
        -------
        kwargs : dict
            Dictionary of drizzle parameters
        """
        # Define the keys pulled from step parameters
        kwargs = {
            "pixfrac": self.pixfrac,
            "kernel": self.kernel,
            "fillval": self.fillval,
            "weight_type": self.weight_type,
            "good_bits": GOOD_BITS,
            "blendheaders": self.blendheaders,
        }

        # Custom output WCS parameters
        wcs_pars = {}
        wcs_pars["output_shape"] = ResampleStep.check_list_pars(
            self.output_shape, "output_shape", min_vals=[1, 1]
        )
        kwargs["output_wcs"] = load_custom_wcs(self.output_wcs, wcs_pars["output_shape"])
        wcs_pars["pixel_scale"] = self.pixel_scale
        wcs_pars["pixel_scale_ratio"] = self.pixel_scale_ratio

        # Report values to processing log
        for k, v in kwargs.items():
            self.log.debug(f"   {k}={v}")

        kwargs["wcs_pars"] = wcs_pars

        return kwargs

    def _process_slit(self, input_models):
        """
        Resample Slit data.

        Parameters
        ----------
        input_models : `~jwst.datamodels.ModelContainer`
            A container of `~jwst.datamodels.ImageModel`
            or `~jwst.datamodels.SlitModel`

        Returns
        -------
        result : `~jwst.datamodels.SlitModel`
            The resampled output
        """
        # Make sure all input models have consistent NaN and DO_NOT_USE values
        for model in input_models:
            match_nans_and_flags(model)

        # Call the resampling routine
        if self.single:
            resamp = resample_spec.ResampleSpec(
                input_models, enable_var=False, compute_err="driz_err", **self.drizpars
            )
            drizzled_library = resamp.resample_many_to_many(in_memory=self.in_memory)
            with drizzled_library:
                result = drizzled_library.borrow(0)
                drizzled_library.shelve(result, 0, modify=False)
            del drizzled_library

        else:
            resamp = resample_spec.ResampleSpec(
                input_models, enable_var=True, compute_err="from_var", **self.drizpars
            )
            result = resamp.resample_many_to_one()

        # result.meta.bunit_data = input_models[0].meta.bunit_data
        if self.pixel_scale is None:
            result.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
        else:
            result.meta.resample.pixel_scale_ratio = resamp.pixel_scale_ratio
        result.meta.resample.pixfrac = self.pixfrac

        self.update_slit_metadata(result)
        if result.meta.exposure.type.lower() == "mir_lrs-fixedslit":
            s_region_model1 = input_models[0].meta.wcsinfo.s_region
            s_region = find_miri_lrs_sregion(s_region_model1, result.meta.wcs)
            result.meta.wcsinfo.s_region = s_region
            self.log.info(f"Updating S_REGION: {s_region}.")
        else:
            update_s_region_spectral(result)
        return result

    def update_slit_metadata(self, model):
        """
        Update slit attributes in the resampled slit image.

        This is needed because model.slit attributes are not in model.meta, so
        the normal update() method doesn't work with them. Updates output_model
        in-place.
        """
        for attr in [
            "name",
            "xstart",
            "xsize",
            "ystart",
            "ysize",
            "slitlet_id",
            "source_id",
            "source_name",
            "source_alias",
            "stellarity",
            "source_type",
            "source_xpos",
            "source_ypos",
            "dispersion_direction",
            "shutter_state",
        ]:
            try:
                val = getattr(self.input_models[-1], attr)
            except AttributeError:
                pass
            else:
                if val is not None:
                    setattr(model, attr, val)
