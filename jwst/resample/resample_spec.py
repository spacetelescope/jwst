import logging
import math
import warnings

import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import (
    Const1D,
    Linear1D,
    Mapping,
    Pix2Sky_TAN,
    RotateNative2Celestial,
    Tabular1D,
)
from astropy.modeling.fitting import LinearLSQFitter
from astropy.stats import sigma_clip

from astropy.utils.exceptions import AstropyUserWarning
from gwcs import wcstools, WCS
from gwcs import coordinate_frames as cf

from stdatamodels.jwst import datamodels

from jwst.assign_wcs.util import compute_scale, wcs_bbox_from_shape, wrap_ra
from jwst.resample import resample_utils
from jwst.resample.resample import ResampleImage
from jwst.datamodels import ModelLibrary

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["ResampleSpec"]


class ResampleSpec(ResampleImage):
    """
    Resample spectral data.

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

    def __init__(self, input_models, good_bits=0, output_wcs=None, wcs_pars=None, **kwargs):
        """
        Initialize the ResampleSpec object.

        Parameters
        ----------
        input_models : list
            List of data models, one for each input image
        good_bits : int
            Bit values that should be considered good when creating a mask
        output_wcs : dict
            Output WCS parameters
        wcs_pars : dict
            Additional parameters for WCS
        **kwargs : dict
            Additional parameters to be passed into `ResampleImage.__init__()`.
            See the docstring of that method for more details.
        """
        shape = None
        pixel_scale = None
        pixel_area = None
        pixel_scale_ratio = 1.0

        if isinstance(output_wcs, dict):
            output_wcs_dict = {k: v for k, v in output_wcs.items() if k != "wcs"}
            output_wcs = output_wcs["wcs"]
            pixel_scale = output_wcs_dict.get("pixel_scale")
            pixel_area = output_wcs_dict.get("pixel_area")
            pixel_scale_ratio = output_wcs_dict.get("pixel_scale_ratio", 1.0)
            shape = output_wcs.array_shape
        else:
            output_wcs_dict = None
            if output_wcs is None and wcs_pars is not None:
                shape = wcs_pars.get("output_shape")
                pixel_scale = wcs_pars.get("pixel_scale")
                pixel_scale_ratio = wcs_pars.get("pixel_scale_ratio", 1.0)

        if pixel_scale is None and pixel_area is not None:
            pixel_scale = math.sqrt(pixel_area)
        elif pixel_scale is not None and pixel_area is None:
            pixel_area = pixel_scale**2

        # Get an average input pixel scale for parameter calculations
        disp_axis = input_models[0].meta.wcsinfo.dispersion_direction
        input_pixscale0 = 3600.0 * compute_spectral_pixel_scale(
            input_models[0].meta.wcs, disp_axis=disp_axis
        )
        if np.isnan(input_pixscale0):
            log.warning("Input pixel scale could not be determined.")
            if pixel_scale is not None:
                log.warning(
                    "Output pixel scale setting is not supported without an "
                    "input pixel scale. Setting pixel_scale=None."
                )
                pixel_scale = None
                pixel_area = None

        nominal_area = input_models[0].meta.photometry.pixelarea_steradians
        if nominal_area is None:
            log.warning("Nominal pixel area not set in input data.")
            log.warning(
                "Setting output pixel scale is not supported without an "
                "input pixel scale. Setting pixel_scale=None."
            )
            pixel_scale = None
            pixel_area = None

        if output_wcs:
            # Use user-supplied reference WCS for the resampled image:
            if pixel_area is None:
                if nominal_area is None:
                    log.warning("Unable to compute output pixel area from 'output_wcs'.")
                    output_pix_area = None
                else:
                    # Compare input and output spatial scale to update nominal area
                    output_pscale = 3600.0 * compute_spectral_pixel_scale(
                        output_wcs, disp_axis=disp_axis
                    )
                    if np.isnan(output_pscale) or np.isnan(input_pixscale0):
                        log.warning("Output pixel scale could not be determined.")
                        output_pix_area = None
                    else:
                        log.debug(
                            f"Setting output pixel area from the approximate "
                            f"output spatial scale: {output_pscale}"
                        )
                        output_pix_area = output_pscale * nominal_area / input_pixscale0

            else:
                log.debug(f"Using output pixel area: {pixel_area}")
                output_pix_area = pixel_area

            # Set the pixel scale ratio for scaling reasons
            if output_pix_area is None:
                pixel_scale_ratio = 1.0
            else:
                pixel_scale_ratio = nominal_area / output_pix_area

            # Set the output shape if specified
            if shape is not None:
                output_wcs.array_shape = shape
        else:
            if pixel_scale is not None and nominal_area is not None:
                log.info(f"Specified output pixel scale: {pixel_scale} arcsec.")

                # Set the pscale ratio from the input pixel scale
                # (pixel scale ratio is output / input)
                if pixel_scale_ratio != 1.0:
                    log.warning(
                        "Ignoring input pixel_scale_ratio in favor of explicit pixel_scale."
                    )
                pixel_scale_ratio = input_pixscale0 / pixel_scale
                log.info(f"Computed output pixel scale ratio: {pixel_scale_ratio:.5g}")

            # Define output WCS based on all inputs, including a reference WCS.
            # These functions internally use pixel_scale_ratio to accommodate
            # user settings.
            # Any other customizations (crpix, crval, rotation) are ignored.
            if resample_utils.is_sky_like(input_models[0].meta.wcs.output_frame):
                if input_models[0].meta.instrument.name != "NIRSPEC":
                    output_wcs = self.build_interpolated_output_wcs(
                        input_models, pixel_scale_ratio=pixel_scale_ratio
                    )
                else:
                    output_wcs = self.build_nirspec_output_wcs(
                        input_models, good_bits=good_bits, pixel_scale_ratio=pixel_scale_ratio
                    )
            else:
                output_wcs = self.build_nirspec_lamp_output_wcs(
                    input_models, pixel_scale_ratio=pixel_scale_ratio
                )

            # Use the nominal output pixel area in sr if available,
            # scaling for user-set pixel_scale ratio if needed.
            if nominal_area is not None:
                # Note that there is only one spatial dimension so the
                # pixel_scale_ratio is not squared.
                output_pix_area = nominal_area / pixel_scale_ratio
            else:
                output_pix_area = None

        self._spec_output_pix_area = output_pix_area

        if pixel_scale is None:
            log.info(f"Specified output pixel scale ratio: {pixel_scale_ratio}.")
            pixel_scale = 3600.0 * compute_spectral_pixel_scale(output_wcs, disp_axis=disp_axis)
            log.info(f"Computed output pixel scale: {pixel_scale:.5g} arcsec.")

        if output_wcs_dict is None:
            output_wcs_dict = {}

        output_wcs_dict["wcs"] = output_wcs
        output_wcs_dict["pixel_scale"] = pixel_scale
        output_wcs_dict["pixel_scale_ratio"] = pixel_scale_ratio

        library = ModelLibrary(input_models, on_disk=False)

        super().__init__(
            library, good_bits=good_bits, output_wcs=output_wcs_dict, wcs_pars=None, **kwargs
        )
        self.intermediate_suffix = "outlier_s2d"

    def create_output_jwst_model(self, ref_input_model=None):
        """
        Create a new blank model and update its meta with info from ``ref_input_model``.

        Parameters
        ----------
        ref_input_model : `~jwst.datamodels.JwstDataModel`, optional
            The reference input model from which to copy meta data.

        Returns
        -------
        SlitModel
            A new blank model with updated meta data.
        """
        output_model = datamodels.SlitModel(None)
        # update meta data and wcs
        if ref_input_model is not None:
            output_model.update(ref_input_model)
        output_model.meta.wcs = self.output_wcs
        return output_model

    def update_output_model(self, model, info_dict):
        """
        Add spectroscopy-specific meta information to the output model.

        Parameters
        ----------
        model : SlitModel
            The output model to be updated.
        info_dict : dict
            A dictionary containing information about the resampling process.
        """
        super().update_output_model(model, info_dict)
        if self._spec_output_pix_area is None:
            model.meta.photometry.pixelarea_steradians = None
            model.meta.photometry.pixelarea_arcsecsq = None
        else:
            model.meta.photometry.pixelarea_steradians = self._spec_output_pix_area
            model.meta.photometry.pixelarea_arcsecsq = (
                self._spec_output_pix_area * np.rad2deg(3600) ** 2
            )

        # TODO: this is helpful info that should be stored in products.
        #       Not storing this at this time in order to reduce the number of
        #       failures in the regression tests.
        # model.meta.resample.pixel_scale_ratio
        # model.meta.resample.pixfrac
        # model.meta.resample.weight_type
        # model.meta.resample.pointings
        # model.meta.cal_step.resample

    def build_nirspec_output_wcs(
        self, input_models, refmodel=None, good_bits=None, pixel_scale_ratio=1.0
    ):
        """
        Create a spatial/spectral WCS covering the footprint of the input.

        Creates the output frame by linearly fitting RA, Dec along the slit
        and producing a lookup table to interpolate wavelengths in the
        dispersion direction.

        For NIRSpec, the output WCS must also provide slit coordinates
        to support source location in the spectral extraction step. To do so,
        this step creates a lookup table for virtual slit coordinates,
        corresponding to the slit y-position at the center of the array
        in the input reference model.  Slit x-position is set to zero for
        all pixels.

        Frames available in the output WCS are:

            - `detector`: image x, y
            - `slit_frame`: slit x, slit y, wavelength
            - `world`: RA, Dec, wavelength

        Parameters
        ----------
        refmodel : `~jwst.datamodels.JwstDataModel`, optional
            The reference input image from which the fiducial WCS is created.
            If not specified, the first image in input_models. If the
            first model is empty (all-NaN or all-zero), the first non-empty
            model is used.

        Returns
        -------
        output_wcs : `~gwcs.WCS`
            A gwcs WCS object defining the output frame WCS.
        """
        all_wcs = [m.meta.wcs for m in input_models if m is not refmodel]
        if refmodel:
            all_wcs.insert(0, refmodel.meta.wcs)
        else:
            # Use the first model with a reasonable amount of good data
            # as the reference model
            for model in input_models:
                dq_mask = resample_utils.build_mask(model.dq, good_bits)
                good = np.isfinite(model.data) & (model.data != 0) & dq_mask
                if np.sum(good) > 100 and refmodel is None:
                    refmodel = model
                    break

            # If no good data was found, use the first model.
            if refmodel is None:
                refmodel = input_models[0]

        # Make a copy of the data array for internal manipulation
        refmodel_data = refmodel.data.copy()

        # Renormalize to the minimum value, for best results when
        # computing the weighted mean below
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
            refmodel_data -= np.nanmin(refmodel_data)

        # Save the wcs of the reference model
        refwcs = refmodel.meta.wcs

        # Set up the transforms that are needed
        s2d = refwcs.get_transform("slit_frame", "detector")
        d2s = refwcs.get_transform("detector", "slit_frame")
        if "moving_target" in refwcs.available_frames:
            s2w = refwcs.get_transform("slit_frame", "moving_target")
            w2s = refwcs.get_transform("moving_target", "slit_frame")
        else:
            s2w = refwcs.get_transform("slit_frame", "world")
            w2s = refwcs.get_transform("world", "slit_frame")

        # Estimate position of the target without relying on the meta.target:
        # compute the mean spatial and wavelength coords weighted
        # by the spectral intensity
        bbox = refwcs.bounding_box
        grid = wcstools.grid_from_bounding_box(bbox)
        _, s, lam = np.array(d2s(*grid))

        # Find invalid values
        good = np.isfinite(s) & np.isfinite(lam) & np.isfinite(refmodel_data)
        refmodel_data[~good] = np.nan

        # Reject the worst outliers in the data
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=AstropyUserWarning, message=".*automatically clipped.*"
            )
            weights = sigma_clip(refmodel_data, masked=True, sigma=100.0)
        weights = np.ma.filled(weights, fill_value=0.0)
        if not np.all(weights == 0.0):
            wmean_s = np.average(s[good], weights=weights[good])
        else:
            wmean_s = np.nanmean(s)
        wmean_l = np.nanmean(lam)

        # Transform the weighted means into target RA/Dec
        # (at the center of the slit in x)
        targ_ra, targ_dec, _ = s2w(0, wmean_s, wmean_l)
        sx, sy = s2d(0, wmean_s, wmean_l)
        log.debug(f"Fiducial RA, Dec, wavelength: {targ_ra}, {targ_dec}, {wmean_l}")
        log.debug(f"Index at fiducial center: x={sx}, y={sy}")

        # Estimate spatial sampling from the reference model
        # at the center of the array
        lam_center_idx = int(np.mean(bbox, axis=1)[0])
        log.debug(f"Center of dispersion axis: {lam_center_idx}")
        grid_center = grid[0][:, lam_center_idx], grid[1][:, lam_center_idx]
        ra_ref, dec_ref, _ = np.array(refwcs(*grid_center))

        # Convert ra and dec to tangent projection
        tan = Pix2Sky_TAN()
        native2celestial = RotateNative2Celestial(targ_ra, targ_dec, 180)
        undist2sky = tan | native2celestial
        x_tan, y_tan = undist2sky.inverse(ra_ref, dec_ref)
        is_nan = np.isnan(x_tan) | np.isnan(y_tan)
        x_tan = x_tan[~is_nan]
        y_tan = y_tan[~is_nan]

        # Estimate the spatial sampling from the tangent projection
        # offset from center
        fitter = LinearLSQFitter()
        fit_model = Linear1D()

        xstop = x_tan.shape[0] * pixel_scale_ratio
        x_idx = np.linspace(0, xstop, x_tan.shape[0], endpoint=False)
        ystop = y_tan.shape[0] * pixel_scale_ratio
        y_idx = np.linspace(0, ystop, y_tan.shape[0], endpoint=False)
        pix_to_xtan = fitter(fit_model, x_idx, x_tan)
        pix_to_ytan = fitter(fit_model, y_idx, y_tan)

        # Check whether sampling is more along RA or along Dec
        swap_xy = abs(pix_to_xtan.slope) < abs(pix_to_ytan.slope)
        log.debug(f"Swap xy: {swap_xy}")

        # Get output wavelengths from all data
        ref_lam = _find_nirspec_output_sampling_wavelengths(all_wcs)
        n_lam = len(ref_lam)
        if not n_lam:
            raise ValueError("Not enough data to construct output WCS.")

        # Find the spatial extent in x/y tangent
        min_tan_x, max_tan_x, min_tan_y, max_tan_y = self._max_spatial_extent(
            all_wcs, undist2sky.inverse
        )
        diff_y = np.abs(max_tan_y - min_tan_y)
        diff_x = np.abs(max_tan_x - min_tan_x)

        pix_to_tan_slope_y = np.abs(pix_to_ytan.slope)
        slope_sign_y = np.sign(pix_to_ytan.slope)
        pix_to_tan_slope_x = np.abs(pix_to_xtan.slope)
        slope_sign_x = np.sign(pix_to_xtan.slope)

        # Image size in spatial dimension from the maximum slope
        # and tangent offset span, plus one pixel to make sure
        # we catch all the data
        if swap_xy:
            ny = int(np.ceil(diff_y / pix_to_tan_slope_y)) + 1
        else:
            ny = int(np.ceil(diff_x / pix_to_tan_slope_x)) + 1

        # Correct the intercept for the new minimum value.
        # Also account for integer pixel size to make sure the
        # data is centered in the array.
        offset_y = ny / 2 * pix_to_tan_slope_y - diff_y / 2
        offset_x = ny / 2 * pix_to_tan_slope_x - diff_x / 2

        if slope_sign_y > 0:
            zero_value_y = min_tan_y
        else:
            zero_value_y = max_tan_y

        if slope_sign_x > 0:
            zero_value_x = min_tan_x
        else:
            zero_value_x = max_tan_x

        pix_to_ytan.intercept = zero_value_y - slope_sign_y * offset_y
        pix_to_xtan.intercept = zero_value_x - slope_sign_x * offset_x

        # Now set up the final transforms

        # For wavelengths, extrapolate 1/2 pixel at the edges and
        # make tabular model w/inverse
        pixel_coord = list(range(n_lam))
        if len(pixel_coord) > 1:
            # left:
            slope = (ref_lam[1] - ref_lam[0]) / pixel_coord[1]
            ref_lam.insert(0, -0.5 * slope + ref_lam[0])
            pixel_coord.insert(0, -0.5)

            # right:
            slope = (ref_lam[-1] - ref_lam[-2]) / (pixel_coord[-1] - pixel_coord[-2])
            ref_lam.append(slope * (pixel_coord[-1] + 0.5) + ref_lam[-2])
            pixel_coord.append(pixel_coord[-1] + 0.5)
        else:
            ref_lam = 3 * ref_lam
            pixel_coord = [-0.5, 0, 0.5]
        wavelength_transform = Tabular1D(
            points=pixel_coord, lookup_table=ref_lam, bounds_error=False, fill_value=np.nan
        )

        # For spatial coordinates, map detector pixels to tangent offset,
        # then to world coordinates (RA, Dec, wavelength in um).
        # Make sure the inverse returns the axis with the larger slope,
        # in case the smaller one is close to zero
        mapping = Mapping((1, 1, 0))
        if swap_xy:
            mapping.inverse = Mapping((2, 1))
        else:
            mapping.inverse = Mapping((2, 0))
        pix2world = mapping | (pix_to_xtan & pix_to_ytan | undist2sky) & wavelength_transform

        # For NIRSpec, slit coordinates are still needed to locate the
        # planned source position.  Since the slit is now rectified,
        # return the central slit coords for all x, converting from pixels
        # to world coordinates, then back to slit units.
        slit_center = w2s(*pix2world(np.full(ny, n_lam // 2), np.arange(ny)))[1]

        # Make a 1D lookup table for all ny.
        # Allow linear extrapolation at the edges.
        slit_transform = Tabular1D(
            points=np.arange(ny), lookup_table=slit_center, bounds_error=False, fill_value=None
        )

        # In the transform, the first slit coordinate is always set to 0
        # to represent the "horizontal" center of the slit
        # (if we imagine the slit to be vertical in the usual
        # X-Y 2D cartesian frame).
        mapping = Mapping((0, 1, 0))
        inv_mapping = Mapping((2, 1), n_inputs=3)
        inv_mapping.inverse = mapping
        mapping.inverse = inv_mapping
        zero_model = Const1D(0)
        zero_model.inverse = zero_model

        # Final detector to slit transform (x, y -> slit_x, slit_y, wavelength)
        det2slit = mapping | zero_model & slit_transform & wavelength_transform

        # The slit to world coordinates is just the inverse of the slit transform,
        # piped back into the pixel to world transform
        slit2world = det2slit.inverse | pix2world

        # Create coordinate frames: detector, slit_frame, and world
        det = cf.Frame2D(name="detector", axes_order=(0, 1))
        slit_spatial = cf.Frame2D(
            name="slit_spatial", axes_order=(0, 1), unit=("", ""), axes_names=("x_slit", "y_slit")
        )
        spec = cf.SpectralFrame(
            name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
        )
        slit_frame = cf.CompositeFrame([slit_spatial, spec], name="slit_frame")
        sky = cf.CelestialFrame(name="sky", axes_order=(0, 1), reference_frame=coord.ICRS())
        world = cf.CompositeFrame([sky, spec], name="world")

        pipeline = [(det, det2slit), (slit_frame, slit2world), (world, None)]
        output_wcs = WCS(pipeline)

        # Compute bounding box and output array shape.
        data_size = (ny, n_lam)
        output_wcs.bounding_box = wcs_bbox_from_shape(data_size)
        output_wcs.array_shape = data_size

        return output_wcs

    def _max_spatial_extent(self, wcs_list, transform):
        """
        Compute spatial coordinate limits for all nods in the tangent plane.

        Parameters
        ----------
        wcs_list : list
            List of WCS objects for all nods.
        transform : callable
            Function to convert RA, Dec to tangent plane coordinates.

        Returns
        -------
        limits_x : tuple
            Minimum and maximum x values.
        limits_y : tuple
            Minimum and maximum y values.
        """
        limits_x = [np.inf, -np.inf]
        limits_y = [np.inf, -np.inf]
        for wcs in wcs_list:
            x, y = wcstools.grid_from_bounding_box(wcs.bounding_box)
            ra, dec, _ = wcs(x, y)

            good = np.logical_and(np.isfinite(ra), np.isfinite(dec))
            ra = ra[good]
            dec = dec[good]

            xtan, ytan = transform(ra, dec)
            for tan_all, limits in zip([xtan, ytan], [limits_x, limits_y], strict=True):
                min_tan = np.min(tan_all)
                max_tan = np.max(tan_all)

                if min_tan < limits[0]:
                    limits[0] = min_tan
                if max_tan > limits[1]:
                    limits[1] = max_tan

        return *limits_x, *limits_y

    def build_interpolated_output_wcs(self, input_models, pixel_scale_ratio=1.0):
        """
        Create a spatial/spectral WCS output frame using all the input models.

        Creates output frame by linearly fitting RA, Dec along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Frames available in the output WCS are:

            - `detector`: image x, y
            - `world`: RA, Dec, wavelength

        Parameters
        ----------
        input_models : list
            List of data models, one for each input image
        pixel_scale_ratio : float
            The ratio of the input pixel scale to the output pixel scale.

        Returns
        -------
        output_wcs : `~gwcs.WCS` object
            A gwcs WCS object defining the output frame WCS
        """
        # for each input model convert slit x,y to ra,dec,lam
        # use first input model to set spatial scale
        # use center of appended ra and dec arrays to set up
        # center of final ra,dec
        # append all ra,dec, wavelength array for each slit
        # use first model to initialize wavelength array
        # append wavelengths that fall outside the endpoint of
        # of wavelength array when looping over additional data

        all_wavelength = []
        all_ra_slit = []
        all_dec_slit = []
        xstop = 0

        all_wcs = [m.meta.wcs for m in input_models]
        for im, model in enumerate(input_models):
            wcs = model.meta.wcs
            bbox = wcs.bounding_box
            grid = wcstools.grid_from_bounding_box(bbox)
            ra, dec, lam = np.array(wcs(*grid))

            # Handle vertical (MIRI).  The following 2 variables are
            # 0 or 1, i.e. zero-indexed in x,y WCS order
            spectral_axis = find_dispersion_axis(model)
            spatial_axis = spectral_axis ^ 1

            # Compute the wavelength array, trimming NaNs from the ends
            # In many cases, a whole slice is NaNs, so ignore those warnings
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")
                wavelength_array = np.nanmedian(lam, axis=spectral_axis)
            wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

            # We need to estimate the spatial sampling to use for the output WCS.
            # It is assumed the spatial sampling is the same for all the input
            # models. So we can use the first input model to set the spatial
            # sampling.

            # Steps to do this for first input model:
            # 1. Find the middle of the spectrum in wavelength
            # 2. Pull out the ra and dec at the center of the slit.
            # 3. Find the mean ra,dec and the center of the slit this will
            #    represent the tangent point
            # 4. Convert ra,dec -> tangent plane projection: x_tan,y_tan
            # 5. using x_tan, y_tan perform a linear fit to find spatial sampling
            # first input model sets initializes wavelength array and defines
            # the spatial scale of the output wcs
            if im == 0:
                all_wavelength = np.append(all_wavelength, wavelength_array)

                # find the center ra and dec for this slit at central wavelength
                lam_center_index = int((bbox[spectral_axis][1] - bbox[spectral_axis][0]) / 2)
                if spatial_axis == 0:
                    # MIRI LRS spectral = 1, the spatial axis = 0
                    ra_slice = ra[lam_center_index, :]
                    dec_slice = dec[lam_center_index, :]
                else:
                    ra_slice = ra[:, lam_center_index]
                    dec_slice = dec[:, lam_center_index]

                # wrap RA if near zero
                ra_center_pt = np.nanmean(wrap_ra(ra_slice))
                dec_center_pt = np.nanmean(dec_slice)

                # convert ra and dec to tangent projection
                tan = Pix2Sky_TAN()
                native2celestial = RotateNative2Celestial(ra_center_pt, dec_center_pt, 180)
                undist2sky1 = tan | native2celestial

                # Filter out RuntimeWarnings due to computed NaNs in the WCS
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    # at this center of slit find x,y tangent projection - x_tan, y_tan
                    x_tan, y_tan = undist2sky1.inverse(ra, dec)

                # pull out data from center
                if spectral_axis == 0:
                    x_tan_array = x_tan.T[lam_center_index]
                    y_tan_array = y_tan.T[lam_center_index]
                else:  # MIRI LRS Spectral Axis = 1, the WCS x axis is spatial
                    x_tan_array = x_tan[lam_center_index]
                    y_tan_array = y_tan[lam_center_index]

                x_tan_array = x_tan_array[~np.isnan(x_tan_array)]
                y_tan_array = y_tan_array[~np.isnan(y_tan_array)]

                # estimate the spatial sampling
                fitter = LinearLSQFitter()
                fit_model = Linear1D()
                xstop = x_tan_array.shape[0] * pixel_scale_ratio
                x_idx = np.linspace(0, xstop, x_tan_array.shape[0], endpoint=False)
                ystop = y_tan_array.shape[0] * pixel_scale_ratio
                y_idx = np.linspace(0, ystop, y_tan_array.shape[0], endpoint=False)
                pix_to_xtan = fitter(fit_model, x_idx, x_tan_array)
                pix_to_ytan = fitter(fit_model, y_idx, y_tan_array)

            # append all ra and dec values to use later to find min and max
            # ra and dec
            ra_use = ra[~np.isnan(ra)].flatten()
            dec_use = dec[~np.isnan(dec)].flatten()
            all_ra_slit = np.append(all_ra_slit, ra_use)
            all_dec_slit = np.append(all_dec_slit, dec_use)

            # now check wavelength array to see if we need to add to it
            this_minw = np.min(wavelength_array)
            this_maxw = np.max(wavelength_array)
            all_minw = np.min(all_wavelength)
            all_maxw = np.max(all_wavelength)
            if this_minw < all_minw:
                addpts = wavelength_array[wavelength_array < all_minw]
                all_wavelength = np.append(all_wavelength, addpts)
            if this_maxw > all_maxw:
                addpts = wavelength_array[wavelength_array > all_maxw]
                all_wavelength = np.append(all_wavelength, addpts)

        # done looping over set of models
        all_ra = np.hstack(all_ra_slit)
        all_dec = np.hstack(all_dec_slit)
        all_wave = np.hstack(all_wavelength)
        all_wave = all_wave[~np.isnan(all_wave)]
        all_wave = np.sort(all_wave, axis=None)
        # Tabular interpolation model, pixels -> lambda
        wavelength_array = np.unique(all_wave)
        # Check if the data is MIRI LRS FIXED Slit. If it is then
        # the wavelength array needs to be flipped so that the resampled
        # dispersion direction matches the dispersion direction on the detector.
        if input_models[0].meta.exposure.type == "MIR_LRS-FIXEDSLIT":
            wavelength_array = np.flip(wavelength_array, axis=None)

        step = 1
        stop = wavelength_array.shape[0]
        points = np.arange(0, stop, step)
        pix_to_wavelength = Tabular1D(
            points=points,
            lookup_table=wavelength_array,
            bounds_error=False,
            fill_value=None,
            name="pix2wavelength",
        )

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is descending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(0, stop, step)

        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(
            points=points,
            lookup_table=lookup_table,
            bounds_error=False,
            fill_value=None,
            name="wavelength2pix",
        )

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))

        # Sometimes the slit is perpendicular to the RA or Dec axis.
        # For example, if the slit is perpendicular to RA, that means
        # the slope of pix_to_xtan will be nearly zero, so make sure
        # mapping.inverse uses pix_to_ytan.inverse.  The auto definition
        # of mapping.inverse is to use the 2nd spatial coordinate, i.e. Dec.

        swap_xy = abs(pix_to_xtan.slope) < abs(pix_to_ytan.slope)
        if swap_xy:
            # Account for vertical or horizontal dispersion on detector
            mapping.inverse = Mapping((2, 1) if spatial_axis else (1, 2))

        # The final transform
        # redefine the ra, dec center tangent point to include all data

        # check if all_ra crosses 0 degrees - this makes it hard to
        # define the min and max ra correctly
        all_ra = wrap_ra(all_ra)
        ra_min = np.amin(all_ra)
        ra_max = np.amax(all_ra)
        ra_center_final = (ra_max + ra_min) / 2.0
        dec_min = np.amin(all_dec)
        dec_max = np.amax(all_dec)
        dec_center_final = (dec_max + dec_min) / 2.0

        tan = Pix2Sky_TAN()
        if len(input_models) == 1:  # single model use ra_center_pt to be consistent
            # with how resample was done before
            ra_center_final = ra_center_pt
            dec_center_final = dec_center_pt

        native2celestial = RotateNative2Celestial(ra_center_final, dec_center_final, 180)
        undist2sky = tan | native2celestial

        ## Use all the wcs
        min_tan_x, max_tan_x, min_tan_y, max_tan_y = self._max_spatial_extent(
            all_wcs, undist2sky.inverse
        )
        diff_y = np.abs(max_tan_y - min_tan_y)
        diff_x = np.abs(max_tan_x - min_tan_x)

        pix_to_tan_slope_y = np.abs(pix_to_ytan.slope)
        slope_sign_y = np.sign(pix_to_ytan.slope)
        pix_to_tan_slope_x = np.abs(pix_to_xtan.slope)
        slope_sign_x = np.sign(pix_to_xtan.slope)

        if swap_xy:
            ny = int(np.ceil(diff_y / pix_to_tan_slope_y))
        else:
            ny = int(np.ceil(diff_x / pix_to_tan_slope_x))

        offset_y = 0.5 * (ny - 1) * pix_to_tan_slope_y
        offset_x = 0.5 * (ny - 1) * pix_to_tan_slope_x

        pix_to_ytan.intercept = -slope_sign_y * offset_y
        pix_to_xtan.intercept = -slope_sign_x * offset_x

        # define the output wcs
        transform = mapping | (pix_to_xtan & pix_to_ytan | undist2sky) & pix_to_wavelength

        det = cf.Frame2D(name="detector", axes_order=(0, 1))
        sky = cf.CelestialFrame(name="sky", axes_order=(0, 1), reference_frame=coord.ICRS())
        spec = cf.SpectralFrame(
            name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
        )
        world = cf.CompositeFrame([sky, spec], name="world")

        pipeline = [(det, transform), (world, None)]

        output_wcs = WCS(pipeline)

        # compute the output array size in WCS axes order, i.e. (x, y)
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = int(np.ceil(len(wavelength_array)))
        output_array_size[spatial_axis] = ny

        # turn the size into a numpy shape in (y, x) order
        output_wcs.array_shape = output_array_size[::-1]
        output_wcs.pixel_shape = output_array_size
        bounding_box = wcs_bbox_from_shape(output_array_size[::-1])
        output_wcs.bounding_box = bounding_box
        return output_wcs

    def build_nirspec_lamp_output_wcs(self, input_models, pixel_scale_ratio):
        """
        Create a spatial/spectral WCS output frame for NIRSpec lamp mode.

        Creates output frame by linearly fitting x_msa, y_msa along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Frames available in the output WCS are:

            - `detector`: image x, y
            - `world`: MSA x, MSA y, wavelength

        Parameters
        ----------
        input_models : list
            List of data models, one for each input image
        pixel_scale_ratio : float
            The ratio of the input pixel scale to the output pixel scale.

        Returns
        -------
        output_wcs : `~gwcs.WCS` object
            A gwcs WCS object defining the output frame WCS.
        """
        model = input_models[0]
        wcs = model.meta.wcs
        bbox = wcs.bounding_box
        grid = wcstools.grid_from_bounding_box(bbox)
        x_msa, y_msa, lam = np.array(wcs(*grid))
        # Handle vertical (MIRI) or horizontal (NIRSpec) dispersion.  The
        # following 2 variables are 0 or 1, i.e. zero-indexed in x,y WCS order
        spectral_axis = find_dispersion_axis(model)
        spatial_axis = spectral_axis ^ 1

        # Compute the wavelength array, trimming NaNs from the ends
        # In many cases, a whole slice is NaNs, so ignore those warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            wavelength_array = np.nanmedian(lam, axis=spectral_axis)
            wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

        # Find the center ra and dec for this slit at central wavelength
        lam_center_index = int((bbox[spectral_axis][1] - bbox[spectral_axis][0]) / 2)
        x_msa_array = x_msa.T[lam_center_index]
        y_msa_array = y_msa.T[lam_center_index]
        x_msa_array = x_msa_array[~np.isnan(x_msa_array)]
        y_msa_array = y_msa_array[~np.isnan(y_msa_array)]

        # Estimate and fit the spatial sampling
        fitter = LinearLSQFitter()
        fit_model = Linear1D()
        xstop = x_msa_array.shape[0] * pixel_scale_ratio
        x_idx = np.linspace(0, xstop, x_msa_array.shape[0], endpoint=False)
        ystop = y_msa_array.shape[0] * pixel_scale_ratio
        y_idx = np.linspace(0, ystop, y_msa_array.shape[0], endpoint=False)
        pix_to_x_msa = fitter(fit_model, x_idx, x_msa_array)
        pix_to_y_msa = fitter(fit_model, y_idx, y_msa_array)

        step = 1
        stop = wavelength_array.shape[0]
        points = np.arange(0, stop, step)
        pix_to_wavelength = Tabular1D(
            points=points,
            lookup_table=wavelength_array,
            bounds_error=False,
            fill_value=None,
            name="pix2wavelength",
        )

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is descending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(0, stop, step)

        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(
            points=points,
            lookup_table=lookup_table,
            bounds_error=False,
            fill_value=None,
            name="wavelength2pix",
        )

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))
        mapping.inverse = Mapping((2, 1))

        # The final transform
        # define the output wcs
        transform = mapping | pix_to_x_msa & pix_to_y_msa & pix_to_wavelength

        det = cf.Frame2D(name="detector", axes_order=(0, 1))
        sky = cf.Frame2D(name=f"resampled_{model.meta.wcs.output_frame.name}", axes_order=(0, 1))
        spec = cf.SpectralFrame(
            name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
        )
        world = cf.CompositeFrame([sky, spec], name="world")

        pipeline = [(det, transform), (world, None)]

        output_wcs = WCS(pipeline)

        # Compute the output array size and bounding box
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = len(wavelength_array)
        x_size = len(x_msa_array)
        output_array_size[spatial_axis] = int(np.ceil(x_size * pixel_scale_ratio))

        # turn the size into a numpy shape in (y, x) order
        output_wcs.array_shape = output_array_size[::-1]
        output_wcs.pixel_shape = output_array_size
        bounding_box = wcs_bbox_from_shape(output_array_size[::-1])
        output_wcs.bounding_box = bounding_box

        return output_wcs


def find_dispersion_axis(refmodel):
    """
    Find the dispersion axis (0-indexed) of the given 2D wavelength array.

    Parameters
    ----------
    refmodel : `~jwst.datamodels.DataModel`
        The input data model.

    Returns
    -------
    dispaxis : int
        The dispersion axis (0-indexed).
    """
    dispaxis = refmodel.meta.wcsinfo.dispersion_direction
    # Change from 1 --> X and 2 --> Y to 0 --> X and 1 --> Y.
    return dispaxis - 1


def _find_nirspec_output_sampling_wavelengths(wcs_list):
    refwcs = wcs_list[0]
    bbox = refwcs.bounding_box

    grid = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lambdas = refwcs(*grid)

    ref_lam = sorted(np.nanmedian(lambdas[:, np.any(np.isfinite(lambdas), axis=0)], axis=0))
    lam1 = ref_lam[0]
    lam2 = ref_lam[-1]

    min_delta = np.fabs(np.ediff1d(ref_lam).min())

    image_lam = []
    for w in wcs_list[1:]:
        bbox = w.bounding_box
        grid = wcstools.grid_from_bounding_box(bbox)
        ra, dec, lambdas = w(*grid)
        lam = sorted(np.nanmedian(lambdas[:, np.any(np.isfinite(lambdas), axis=0)], axis=0))
        image_lam.append((lam, np.min(lam), np.max(lam)))
        min_delta = min(min_delta, np.fabs(np.ediff1d(ref_lam).min()))

    # The code below is optimized for the case when wavelength is an increasing
    # function of the pixel index along the X-axis. It will not work correctly
    # if this assumption does not hold.

    # Estimate overlaps between ranges and decide in which order to combine
    # them:

    while image_lam:
        best_overlap = -np.inf
        best_wcs = 0
        for k, (_lam, lmin, lmax) in enumerate(image_lam):
            overlap = min(lam2, lmax) - max(lam1, lmin)
            if best_overlap < overlap:
                best_overlap = overlap
                best_wcs = k

        lam, lmin, lmax = image_lam.pop(best_wcs)
        if lmax < lam1:
            ref_lam = lam + ref_lam
            lam1 = lmin
        elif lmin > lam2:
            ref_lam.extend(lam)
            lam2 = lmax
        else:
            lam_ar = np.array(lam)
            if lmin < lam1:
                idx = np.flatnonzero(lam_ar < lam1)
                ref_lam = lam_ar[idx].tolist() + ref_lam
                lam1 = ref_lam[0]
            if lmax > lam2:
                idx = np.flatnonzero(lam_ar > lam2)
                ref_lam = ref_lam + lam_ar[idx].tolist()
                lam2 = ref_lam[-1]

    # In the resampled WCS, if two wavelengths are closer to each other
    # than 1/10 of the minimum difference between two wavelengths,
    # remove one of the points.
    ediff = np.fabs(np.ediff1d(ref_lam))
    idx = np.flatnonzero(ediff < max(0.1 * min_delta, 1e2 * np.finfo(1.0).eps))
    for i in idx[::-1]:
        del ref_lam[i]

    return ref_lam


def compute_spectral_pixel_scale(wcs, fiducial=None, disp_axis=1):
    """
    Compute an approximate spatial pixel scale for spectral data.

    Parameters
    ----------
    wcs : gwcs.WCS
        Spatial/spectral WCS.
    fiducial : tuple of float, optional
        (RA, Dec, wavelength) taken as the fiducial reference. If
        not specified, the center of the array is used.
    disp_axis : int
        Dispersion axis for the data. Assumes the same convention
        as `wcsinfo.dispersion_direction` (1 for NIRSpec, 2 for MIRI).

    Returns
    -------
    pixel_scale : float
        The spatial scale in degrees.
    """
    # Get the coordinates for the center of the array
    if fiducial is None:
        center_x, center_y = np.mean(wcs.bounding_box, axis=1)
        fiducial = wcs(center_x, center_y)

    pixel_scale = compute_scale(wcs, fiducial, disp_axis=disp_axis)
    return float(pixel_scale)
