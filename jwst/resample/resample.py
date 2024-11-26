import logging
import os
import warnings
import json

import numpy as np
from drizzle.resample import Drizzle
from spherical_geometry.polygon import SphericalPolygon

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary
from jwst.associations.asn_from_list import asn_from_list

from jwst.model_blender.blender import ModelBlender
from jwst.resample import resample_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleData"]


class ResampleData:
    """
    This is the controlling routine for the resampling process.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model, such as pixfrac,
         weight type, exposure time (if relevant), and kernel, and merges
         them with any user-provided values.
      2. Creates output WCS based on input images and define mapping function
         between all input arrays and the output array.
      3. Updates output data model with output arrays from drizzle, including
         a record of metadata from all input models.
    """

    def __init__(self, input_models, output=None, single=False, blendheaders=True,
                 pixfrac=1.0, kernel="square", fillval="NAN", wht_type="ivm",
                 good_bits=0, pscale_ratio=1.0, pscale=None, **kwargs):
        """
        Parameters
        ----------
        input_models : library of objects
            library of data models, one for each input image

        output : str
            filename for output

        kwargs : dict
            Other parameters.

            .. note::
                ``output_shape`` is in the ``x, y`` order.

            .. note::
                ``in_memory`` controls whether or not the resampled
                array from ``resample_many_to_many()``
                should be kept in memory or written out to disk and
                deleted from memory. Default value is `True` to keep
                all products in memory.
        """
        self.output_dir = None
        self.output_filename = output
        if output is not None and '.fits' not in str(output):
            self.output_dir = output
            self.output_filename = None
        self.intermediate_suffix = 'outlier_i2d'

        self.pscale_ratio = pscale_ratio
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = wht_type
        self.good_bits = good_bits
        self.in_memory = kwargs.get('in_memory', True)
        self.input_pixscale0 = None  # computed pixel scale of the first image (deg)
        self._recalc_pscale_ratio = pscale is not None

        log.info(f"Driz parameter kernel: {self.kernel}")
        log.info(f"Driz parameter pixfrac: {self.pixfrac}")
        log.info(f"Driz parameter fillval: {self.fillval}")
        log.info(f"Driz parameter weight_type: {self.weight_type}")

        output_wcs = kwargs.get('output_wcs', None)
        output_shape = kwargs.get('output_shape', None)
        crpix = kwargs.get('crpix', None)
        crval = kwargs.get('crval', None)
        rotation = kwargs.get('rotation', None)

        self.asn_id = kwargs.get('asn_id', None)

        if pscale is not None:
            log.info(f'Output pixel scale: {pscale} arcsec.')
            pscale /= 3600.0
        else:
            log.info(f'Output pixel scale ratio: {pscale_ratio}')

        if output_wcs:
            # Use user-supplied reference WCS for the resampled image:
            self.output_wcs = output_wcs
            if output_shape is not None:
                self.output_wcs.array_shape = output_shape[::-1]

            if output_wcs.pixel_area is None:
                self.output_pix_area = compute_image_pixel_area(self.output_wcs)
                if self.output_pix_area is None:
                    raise ValueError(
                        "Unable to compute output pixel area from 'output_wcs'."
                    )
            else:
                self.output_pix_area = output_wcs.pixel_area

        else:
            # Define output WCS based on all inputs, including a reference WCS:
            self.output_wcs = resample_utils.make_output_wcs(
                input_models,
                ref_wcs=output_wcs,
                pscale_ratio=self.pscale_ratio,
                pscale=pscale,
                rotation=rotation,
                shape=None if output_shape is None else output_shape[::-1],
                crpix=crpix,
                crval=crval
            )
            # Estimate output pixel area in Sr. NOTE: in principle we could
            # use the same algorithm as for when output_wcs is provided by the
            # user.
            tr = self.output_wcs.pipeline[0].transform
            self.output_pix_area = (
                np.deg2rad(tr['cdelt1'].factor.value) *
                np.deg2rad(tr['cdelt2'].factor.value)
            )

        self.output_array_shape = tuple(self.output_wcs.array_shape)
        log.debug(f"Output mosaic size: {self.output_array_shape}")

        if pscale is None:
            pscale = np.rad2deg(np.sqrt(self.output_pix_area))
            log.info(f'Computed output pixel scale: {3600.0 * pscale} arcsec.')

        self.pscale = pscale  # in deg

        log.debug(f"Output mosaic size: {tuple(self.output_wcs.pixel_shape)}")


    def _create_output_model(self, ref_input_model=None):
        """ Create a new blank model and update it's meta with info from ``ref_input_model``. """
        output_model = datamodels.ImageModel(None)  # tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self.output_wcs
        if self.output_pix_area is not None:
            output_model.meta.photometry.pixelarea_steradians = self.output_pix_area
            output_model.meta.photometry.pixelarea_arcsecsq = (
                self.output_pix_area * np.rad2deg(3600)**2
            )
        return output_model

    def do_drizzle(self, input_models):
        """Pick the correct drizzling mode based on self.single
        """
        if self.single:
            return self.resample_many_to_many(input_models)
        else:
            return self.resample_many_to_one(input_models)

    def _get_intensity_scale(self, img):
        """
        Compute an intensity scale from the input and output pixel area.

        For imaging data, the scaling is used to account for differences
        between the nominal pixel area and the average pixel area for
        the input data.

        For spectral data, the scaling is used to account for flux
        conservation with non-unity pixel scale ratios, when the
        data units are flux density.

        Parameters
        ----------
        img : DataModel
            The input data model.

        Returns
        -------
        iscale : float
            The scale to apply to the input data before drizzling.
        """
        input_pixflux_area = img.meta.photometry.pixelarea_steradians
        if input_pixflux_area:
            if 'SPECTRAL' in img.meta.wcs.output_frame.axes_type:
                # Use the nominal area as is
                input_pixel_area = input_pixflux_area

                # If input image is in flux density units, correct the
                # flux for the user-specified change to the spatial dimension
                if resample_utils.is_flux_density(img.meta.bunit_data):
                    input_pixel_area *= self.pscale_ratio
            else:
                img.meta.wcs.array_shape = img.data.shape
                input_pixel_area = compute_image_pixel_area(img.meta.wcs)
                if input_pixel_area is None:
                    raise ValueError(
                        "Unable to compute input pixel area from WCS of input "
                        f"image {repr(img.meta.filename)}."
                    )
                if self.input_pixscale0 is None:
                    self.input_pixscale0 = np.rad2deg(
                        np.sqrt(input_pixel_area)
                    )
                    if self._recalc_pscale_ratio:
                        self.pscale_ratio = self.pscale / self.input_pixscale0
            iscale = np.sqrt(input_pixflux_area / input_pixel_area)
        else:
            iscale = 1.0
        return iscale

    def resample_group(self, input_models, indices, compute_error=False):
        """Apply resample_many_to_many for one group

        Parameters
        ----------
        input_models : ModelLibrary

        indices : list

        compute_error : bool, optional
            If True, an approximate error image will be resampled
            alongside the science image.
        """
        output_model = None

        # Initialize the output with the wcs
        driz = Drizzle(
            out_shape=self.output_array_shape,
            kernel=self.kernel,
            fillval=self.fillval,
            disable_ctx=True,
        )
        # Also make a temporary model to hold error data
        if compute_error:
            driz_error = Drizzle(
                out_shape=self.output_array_shape,
                kernel=self.kernel,
                fillval=self.fillval,
                disable_ctx=True,
            )
        log.info(f"{len(indices)} exposures to drizzle together")
        for index in indices:
            img = input_models.borrow(index)
            if output_model is None:
                output_model = self._create_output_model(
                    ref_input_model=img
                )
                # Determine output file type from input exposure filenames
                # Use this for defining the output filename
                indx = img.meta.filename.rfind('.')
                output_type = img.meta.filename[indx:]
                output_root = '_'.join(img.meta.filename.replace(
                    output_type,
                    ''
                ).split('_')[:-1])
                output_model.meta.filename = (
                    f'{output_root}_'
                    f'{self.intermediate_suffix}{output_type}'
                )
                copy_asn_info_from_library(input_models, output_model)

            if isinstance(img, datamodels.SlitModel):
                # must call this explicitly to populate area extension
                # although the existence of this extension may not be necessary
                img.area = img.area

            iscale = self._get_intensity_scale(img)
            log.debug(f'Using intensity scale iscale={iscale}')

            inwht = resample_utils.build_driz_weight(
                img,
                weight_type=self.weight_type,
                good_bits=self.good_bits
            )

            # apply sky subtraction
            blevel = img.meta.background.level
            if not img.meta.background.subtracted and blevel is not None:
                data = img.data - blevel
            else:
                data = img.data

            xmin, xmax, ymin, ymax = resample_utils._resample_range(
                data.shape,
                img.meta.wcs.bounding_box
            )
            pixmap = resample_utils.calc_gwcs_pixmap(
                img.meta.wcs,
                self.output_wcs,
                img.data.shape,
            )

            driz.add_image(
                data=data,
                exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default was 1.0
                pixmap=pixmap,
                scale=iscale,
                weight_map=inwht,
                wht_scale=1.0,  # hard-coded for JWST count-rate data
                pixfrac=self.pixfrac,
                in_units="cps",  # GWCSDrizzle.add_image param default
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
            )
            del data

            # make an approximate error image by drizzling it
            # in the same way the image is handled
            if compute_error:
                driz_error.add_image(
                    data=img.err,
                    exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default
                    pixmap=pixmap,
                    scale=iscale,
                    weight_map=inwht,
                    wht_scale=1.0,  # hard-coded for JWST count-rate data
                    pixfrac=self.pixfrac,
                    in_units="cps",  # GWCSDrizzle.add_image param default
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
            input_models.shelve(img, index, modify=False)
            del img

        output_model.data = driz.out_img
        output_model.wht = driz.out_wht
        del driz
        # copy the drizzled error into the output model
        if compute_error:
            output_model.err = driz_error.out_img
            del driz_error

        return output_model

    def resample_many_to_many(self, input_models):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection
        """
        output_models = []
        for group_id, indices in input_models.group_indices.items():

            output_model = self.resample_group(input_models, indices)

            if not self.in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.filename
                if self.output_dir is not None:
                    output_name = os.path.join(self.output_dir, output_name)
                output_model.save(output_name)
                log.info(f"Saved model in {output_name}")
                output_models.append(output_name)
            else:
                output_models.append(output_model)

        if self.in_memory:
            # build ModelLibrary as a list of in-memory models
            return ModelLibrary(output_models, on_disk=False)
        else:
            # build ModelLibrary as an association from the output files
            # this saves memory if there are multiple groups
            asn = asn_from_list(output_models, product_name='outlier_i2d')
            asn_dict = json.loads(asn.dump()[1])  # serializes the asn and converts to dict
            return ModelLibrary(asn_dict, on_disk=True)

    def resample_many_to_one(self, input_models):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        output_model = None

        if self.blendheaders:
            blender = ModelBlender(
                blend_ignore_attrs=[
                    'meta.photometry.pixelarea_steradians',
                    'meta.photometry.pixelarea_arcsecsq',
                    'meta.filename',
                ]
            )

        driz = Drizzle(
            out_shape=self.output_array_shape,
            kernel=self.kernel,
            fillval=self.fillval,
            max_ctx_id=len(input_models),
            disable_ctx=False,
        )
        self._init_variance_arrays()
        self._init_exptime_counters()

        log.info("Resampling science and variance data")

        leading_group_idx = [v[0] for v in input_models.group_indices.values()]

        with input_models:
            for idx, img in enumerate(input_models):
                if output_model is None:
                    output_model = self._create_output_model(
                        ref_input_model=img
                    )
                    # Determine output file type from input exposure filenames
                    # Use this for defining the output filename
                    output_model.meta.filename = self.output_filename
                    output_model.meta.resample.weight_type = self.weight_type
                    output_model.meta.resample.pointings = len(input_models.group_names)

                    # copy over asn information
                    copy_asn_info_from_library(input_models, output_model)

                if idx in leading_group_idx:
                    self._update_exptime(img)

                if self.blendheaders:
                    blender.accumulate(img)

                iscale = self._get_intensity_scale(img)
                log.debug(f'Using intensity scale iscale={iscale}')

                inwht = resample_utils.build_driz_weight(
                    img,
                    weight_type=self.weight_type,
                    good_bits=self.good_bits,
                )

                # apply sky subtraction
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    data = img.data - blevel
                else:
                    data = img.data.copy()

                in_image_limits = resample_utils._resample_range(
                    data.shape,
                    img.meta.wcs.bounding_box
                )
                xmin, xmax, ymin, ymax = in_image_limits

                pixmap = resample_utils.calc_gwcs_pixmap(
                    img.meta.wcs,
                    output_model.meta.wcs,
                    data.shape,
                )

                driz.add_image(
                    data=data,
                    exptime=img.meta.exposure.exposure_time,  # GWCSDrizzle.add_image param default
                    pixmap=pixmap,
                    scale=iscale,
                    weight_map=inwht,
                    wht_scale=1.0,  # hard-coded for JWST count-rate data
                    pixfrac=self.pixfrac,
                    in_units="cps",  # GWCSDrizzle.add_image param default
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
                # Resample variance arrays in input_models to output_model
                self._resample_variance_arrays(
                    model=img,
                    iscale=iscale,
                    inwht=inwht,
                    pixmap=pixmap,
                    in_image_limits=in_image_limits,
                    output_shape=self.output_array_shape,
                )

                del data, inwht

                input_models.shelve(img)

        # Since the context array is dynamic, it must be re-assigned
        # back to the product's `con` attribute.
        output_model.data = driz.out_img
        output_model.wht = driz.out_wht
        if driz.out_ctx is not None:
            output_model.con = driz.out_ctx

        del driz

        # compute final variances:
        self._compute_resample_variance_totals(output_model)

        var_components = [
            output_model.var_rnoise,
            output_model.var_poisson,
            output_model.var_flat
        ]
        output_model.err = np.sqrt(np.nansum(var_components, axis=0))

        # nansum returns zero for input that is all NaN -
        # set those values to NaN instead
        all_nan = np.all(np.isnan(var_components), axis=0)
        output_model.err[all_nan] = np.nan

        if self.blendheaders:
            blender.finalize_model(output_model)

        self._get_exptime_totals(output_model)

        return ModelLibrary([output_model,], on_disk=False)

    def _init_variance_arrays(self):
        shape = self.output_array_shape
        self._weighted_rn_var = np.full(shape, np.nan, dtype=np.float32)
        self._weighted_pn_var = np.full(shape, np.nan, dtype=np.float32)
        self._weighted_flat_var = np.full(shape, np.nan, dtype=np.float32)
        self._total_weight_rn_var = np.zeros(shape, dtype=np.float32)
        self._total_weight_pn_var = np.zeros(shape, dtype=np.float32)
        self._total_weight_flat_var = np.zeros(shape, dtype=np.float32)

    def _resample_variance_arrays(self, model, iscale, inwht, pixmap,
                                  in_image_limits, output_shape):
        xmin, xmax, ymin, ymax = in_image_limits

        # Do the read noise variance first, so it can be
        # used for weights if needed
        rn_var = self._resample_one_variance_array(
            "var_rnoise",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        # Find valid weighting values in the variance
        if rn_var is not None:
            mask = (rn_var > 0) & np.isfinite(rn_var)
        else:
            mask = np.full_like(rn_var, False)

        # Set the weight for the image from the weight type
        weight = np.ones(output_shape)
        if self.weight_type == "ivm" and rn_var is not None:
            weight[mask] = rn_var[mask] ** -1
        elif self.weight_type == "exptime":
            if resample_utils.check_for_tmeasure(model):
                weight[:] = model.meta.exposure.measurement_time
            else:
                weight[:] = model.meta.exposure.exposure_time

        # Weight and add the readnoise variance
        # Note: floating point overflow is an issue if variance weights
        # are used - it can't be squared before multiplication
        if rn_var is not None:
            mask = (rn_var >= 0) & np.isfinite(rn_var) & (weight > 0)
            self._weighted_rn_var[mask] = np.nansum(
                [
                    self._weighted_rn_var[mask],
                    rn_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_rn_var[mask] += weight[mask]

        # Now do poisson and flat variance, updating only valid new values
        # (zero is a valid value; negative, inf, or NaN are not)
        pn_var = self._resample_one_variance_array(
            "var_poisson",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )
        if pn_var is not None:
            mask = (pn_var >= 0) & np.isfinite(pn_var) & (weight > 0)
            self._weighted_pn_var[mask] = np.nansum(
                [
                    self._weighted_pn_var[mask],
                    pn_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_pn_var[mask] += weight[mask]

        flat_var = self._resample_one_variance_array(
            "var_flat",
            input_model=model,
            iscale=iscale,
            inwht=inwht,
            pixmap=pixmap,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )
        if flat_var is not None:
            mask = (flat_var >= 0) & np.isfinite(flat_var) & (weight > 0)
            self._weighted_flat_var[mask] = np.nansum(
                [
                    self._weighted_flat_var[mask],
                    flat_var[mask] * weight[mask] * weight[mask]
                ],
                axis=0
            )
            self._total_weight_flat_var[mask] += weight[mask]

    def _compute_resample_variance_totals(self, output_model):
        # Divide by the total weights, squared, and set in the output model.
        # Zero weight and missing values are NaN in the output.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "invalid value*", RuntimeWarning)
            warnings.filterwarnings("ignore", "divide by zero*", RuntimeWarning)

            output_variance = (
                self._weighted_rn_var / self._total_weight_rn_var /
                self._total_weight_rn_var
            )
            setattr(output_model, "var_rnoise", output_variance)

            output_variance = (
                self._weighted_pn_var / self._total_weight_pn_var /
                self._total_weight_pn_var
            )
            setattr(output_model, "var_poisson", output_variance)

            output_variance = (
                self._weighted_flat_var / self._total_weight_flat_var /
                self._total_weight_flat_var
            )
            setattr(output_model, "var_flat", output_variance)

        del (
            self._weighted_rn_var,
            self._weighted_pn_var,
            self._weighted_flat_var,
            self._total_weight_rn_var,
            self._total_weight_pn_var,
            self._total_weight_flat_var,
        )

    def _resample_one_variance_array(self, name, input_model, iscale,
                                     inwht, pixmap,
                                     xmin=None, xmax=None, ymin=None, ymax=None):
        """Resample one variance image from an input model.

        The error image is passed to drizzle instead of the variance, to
        better match kernel overlap and user weights to the data, in the
        pixel averaging process. The drizzled error image is squared before
        returning.
        """
        variance = getattr(input_model, name)
        if variance is None or variance.size == 0:
            log.debug(
                f"No data for '{name}' for model "
                f"{repr(input_model.meta.filename)}. Skipping ..."
            )
            return

        elif variance.shape != input_model.data.shape:
            log.warning(
                f"Data shape mismatch for '{name}' for model "
                f"{repr(input_model.meta.filename)}. Skipping ..."
            )
            return

        output_shape = self.output_array_shape

        # Resample the error array. Fill "unpopulated" pixels with NaNs.
        driz = Drizzle(
            out_shape=output_shape,
            kernel=self.kernel,
            fillval=np.nan,
            disable_ctx=True
        )

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {variance.shape}")
        log.debug(f"Output Sci shape: {output_shape}")

        # Call 'drizzle' to perform image combination
        log.info(f"Drizzling {variance.shape} --> {output_shape}")

        driz.add_image(
            data=np.sqrt(variance),
            exptime=input_model.meta.exposure.exposure_time,
            pixmap=pixmap,
            scale=iscale,
            weight_map=inwht,
            wht_scale=1.0,  # hard-coded for JWST count-rate data
            pixfrac=self.pixfrac,
            in_units="cps",
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        return driz.out_img ** 2

    def _init_exptime_counters(self):
        self._total_exposure_time = 0.
        self._exposure_times = {'start': [], 'end': []}
        self._duration = 0.0
        self._total_measurement_time = 0.0
        self._measurement_time_failures = []

    def _update_exptime(self, model):
        self._total_exposure_time += model.meta.exposure.exposure_time
        if not resample_utils.check_for_tmeasure(model):
            self._measurement_time_failures.append(1)
        else:
            self._total_measurement_time += model.meta.exposure.measurement_time
            self._measurement_time_failures.append(0)
        self._exposure_times['start'].append(model.meta.exposure.start_time)
        self._exposure_times['end'].append(model.meta.exposure.end_time)
        self._duration += model.meta.exposure.duration

    def _get_exptime_totals(self, output_model):
        # Update some basic exposure time values based on output_model
        output_model.meta.exposure.exposure_time = self._total_exposure_time
        if not any(self._measurement_time_failures):
            output_model.meta.exposure.measurement_time = self._total_measurement_time
        output_model.meta.exposure.start_time = min(self._exposure_times['start'])
        output_model.meta.exposure.end_time = max(self._exposure_times['end'])

        # Update other exposure time keywords:
        # XPOSURE (identical to the total effective exposure time, EFFEXPTM)
        xposure = self._total_exposure_time
        output_model.meta.exposure.effective_exposure_time = xposure
        # DURATION (identical to TELAPSE, elapsed time)
        output_model.meta.exposure.duration = self._duration
        output_model.meta.exposure.elapsed_exposure_time = self._duration

        del (
            self._total_exposure_time,
            self._exposure_times,
            self._duration,
            self._total_measurement_time,
            self._measurement_time_failures,
        )

    def update_exposure_times(self, output_model, input_models):
        """Modify exposure time metadata in-place"""
        self._init_exptime_counters()
        with input_models:
            for _, indices in input_models.group_indices.items():
                model = input_models.borrow(indices[0])
                self._update_exptime(model)
        self._get_exptime_totals(output_model)


def _get_boundary_points(xmin, xmax, ymin, ymax, dx=None, dy=None, shrink=0):
    """
    xmin, xmax, ymin, ymax - integer coordinates of pixel boundaries
    step - distance between points along an edge
    shrink - number of pixels by which to reduce `shape`

    Returns a list of points and the area of the rectangle
    """
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    if dx is None:
        dx = nx
    if dy is None:
        dy = ny

    if nx - 2 * shrink < 1 or ny - 2 * shrink < 1:
        raise ValueError("Image size is too small.")

    sx = max(1, int(np.ceil(nx / dx)))
    sy = max(1, int(np.ceil(ny / dy)))

    xmin += shrink
    xmax -= shrink
    ymin += shrink
    ymax -= shrink

    size = 2 * sx + 2 * sy
    x = np.empty(size)
    y = np.empty(size)

    b = np.s_[0:sx]  # bottom edge
    r = np.s_[sx:sx + sy]  # right edge
    t = np.s_[sx + sy:2 * sx + sy]  # top edge
    l = np.s_[2 * sx + sy:2 * sx + 2 * sy]  # left

    x[b] = np.linspace(xmin, xmax, sx, False)
    y[b] = ymin
    x[r] = xmax
    y[r] = np.linspace(ymin, ymax, sy, False)
    x[t] = np.linspace(xmax, xmin, sx, False)
    y[t] = ymax
    x[l] = xmin
    y[l] = np.linspace(ymax, ymin, sy, False)

    area = (xmax - xmin) * (ymax - ymin)
    center = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))

    return x, y, area, center, b, r, t, l


def compute_image_pixel_area(wcs):
    """ Computes pixel area in steradians.
    """
    if wcs.array_shape is None:
        raise ValueError("WCS must have array_shape attribute set.")

    valid_polygon = False
    spatial_idx = np.where(np.array(wcs.output_frame.axes_type) == 'SPATIAL')[0]

    ny, nx = wcs.array_shape
    ((xmin, xmax), (ymin, ymax)) = wcs.bounding_box

    xmin = max(0, int(xmin + 0.5))
    xmax = min(nx - 1, int(xmax - 0.5))
    ymin = max(0, int(ymin + 0.5))
    ymax = min(ny - 1, int(ymax - 0.5))
    if xmin > xmax:
        (xmin, xmax) = (xmax, xmin)
    if ymin > ymax:
        (ymin, ymax) = (ymax, ymin)

    k = 0
    dxy = [1, -1, -1, 1]
    ra, dec, center = np.nan, np.nan, (np.nan, np.nan)
    while xmin < xmax and ymin < ymax:
        try:
            x, y, image_area, center, b, r, t, l = _get_boundary_points(
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                dx=min((xmax - xmin) // 4, 15),
                dy=min((ymax - ymin) // 4, 15)
            )
        except ValueError:
            return None

        world = wcs(x, y)
        ra = world[spatial_idx[0]]
        dec = world[spatial_idx[1]]

        limits = [ymin, xmax, ymax, xmin]

        for j in range(4):
            sl = [b, r, t, l][k]
            if not (np.all(np.isfinite(ra[sl])) and
                    np.all(np.isfinite(dec[sl]))):
                limits[k] += dxy[k]
                k = (k + 1) % 4
                break
            k = (k + 1) % 4
        else:
            valid_polygon = True
            break

        ymin, xmax, ymax, xmin = limits

    if not valid_polygon:
        return None

    world = wcs(*center)
    wcenter = (world[spatial_idx[0]], world[spatial_idx[1]])

    sky_area = SphericalPolygon.from_radec(ra, dec, center=wcenter).area()
    if sky_area > 2 * np.pi:
        log.warning(
            "Unexpectedly large computed sky area for an image. "
            "Setting area to: 4*Pi - area"
        )
        sky_area = 4 * np.pi - sky_area
    pix_area = sky_area / image_area

    return pix_area


def copy_asn_info_from_library(library, output_model):
    """
    Transfer association information from the input library to the output model.

    Parameters
    ----------
    library : ModelLibrary
        The input library of data models.

    output_model : DataModel
        The output data model to which the association information will be copied.
    """
    if not hasattr(library, "asn"):
        # No ASN table, occurs when input comes from ModelContainer in spectroscopic modes
        # in this case do nothing; the asn info will be passed along later
        # by code inside ResampleSpecStep
        return
    if (asn_pool := library.asn.get("asn_pool", None)) is not None:
        output_model.meta.asn.pool_name = asn_pool
    if (
        asn_table_name := library.asn.get("table_name", None)
    ) is not None:
        output_model.meta.asn.table_name = asn_table_name
