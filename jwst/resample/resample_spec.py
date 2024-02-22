import logging
import warnings

import numpy as np
from scipy.optimize import minimize_scalar
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import (
    Mapping, Tabular1D, Linear1D, Pix2Sky_TAN, RotateNative2Celestial, Const1D
)
from astropy.modeling.fitting import LinearLSQFitter
from gwcs import wcstools, WCS
from gwcs import coordinate_frames as cf
from gwcs.geometry import SphericalToCartesian

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from ..assign_wcs.util import wrap_ra
from . import resample_utils
from .resample import ResampleData


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_S2C = SphericalToCartesian()


class ResampleSpecData(ResampleData):
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

    def __init__(self, input_models, output=None, single=False, blendheaders=False,
                 pixfrac=1.0, kernel="square", fillval=0, wht_type="ivm",
                 good_bits=0, pscale_ratio=1.0, pscale=None, **kwargs):
        """
        Parameters
        ----------
        input_models : list of objects
            list of data models, one for each input image

        output : str
            filename for output

        kwargs : dict
            Other parameters
        """
        self.input_models = input_models

        self.output_filename = output
        self.pscale_ratio = pscale_ratio
        self.pscale = pscale
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = wht_type
        self.good_bits = good_bits
        self.in_memory = kwargs.get('in_memory', True)

        output_wcs = kwargs.get('output_wcs', None)
        output_shape = kwargs.get('output_shape', None)

        self.input_pixscale0 = None  # computed pixel scale of the first image (deg)
        self._recalc_pscale_ratio = pscale is not None

        # Define output WCS based on all inputs, including a reference WCS
        if output_wcs is None:
            if resample_utils.is_sky_like(
                self.input_models[0].meta.wcs.output_frame
            ):
                if self.input_models[0].meta.instrument.name != "NIRSPEC":
                    self.output_wcs = self.build_interpolated_output_wcs()
                else:
                    self.output_wcs = self.build_nirspec_output_wcs()
            else:
                self.output_wcs = self.build_nirspec_lamp_output_wcs()
        else:
            self.output_wcs = output_wcs
            if output_shape is not None:
                self.output_wcs.array_shape = output_shape[::-1]

        self.blank_output = datamodels.SlitModel(tuple(self.output_wcs.array_shape))
        self.blank_output.update(self.input_models[0])
        self.blank_output.meta.wcs = self.output_wcs
        self.output_models = ModelContainer()

        log.info(f"Driz parameter kernal: {self.kernel}")
        log.info(f"Driz parameter pixfrac: {self.pixfrac}")
        log.info(f"Driz parameter fillval: {self.fillval}")
        log.info(f"Driz parameter weight_type: {self.weight_type}")

    def build_nirspec_output_wcs(self, refmodel=None):
        """
        Create a spatial/spectral WCS covering footprint of the input
        """
        all_wcs = [m.meta.wcs for m in self.input_models if m is not refmodel]
        if refmodel:
            all_wcs.insert(0, refmodel.meta.wcs)
        else:
            refmodel = self.input_models[0]

        # make a copy of the data array for internal manipulation
        refmodel_data = refmodel.data.copy()
        # renormalize to the minimum value, for best results when
        # computing the weighted mean below
        refmodel_data -= np.nanmin(refmodel_data)

        # save the wcs of the reference model
        refwcs = refmodel.meta.wcs

        # setup the transforms that are needed
        s2d = refwcs.get_transform('slit_frame', 'detector')
        d2s = refwcs.get_transform('detector', 'slit_frame')
        s2w = refwcs.get_transform('slit_frame', 'world')

        # estimate position of the target without relying on the meta.target:
        # compute the mean spatial and wavelength coords weighted
        # by the spectral intensity
        bbox = refwcs.bounding_box
        grid = wcstools.grid_from_bounding_box(bbox)
        _, s, lam = np.array(d2s(*grid))
        sd = s * refmodel_data
        ld = lam * refmodel_data
        good_s = np.isfinite(sd)
        if np.any(good_s):
            total = np.sum(refmodel_data[good_s])
            wmean_s = np.sum(sd[good_s]) / total
            wmean_l = np.sum(ld[good_s]) / total
        else:
            wmean_s = 0.5 * (refmodel.slit_ymax - refmodel.slit_ymin)
            wmean_l = d2s(*np.mean(bbox, axis=1))[2]

        # transform the weighted means into target RA/Dec
        targ_ra, targ_dec, _ = s2w(0, wmean_s, wmean_l)

        ref_lam = _find_nirspec_output_sampling_wavelengths(
            all_wcs,
            targ_ra, targ_dec
        )
        ref_lam = np.array(ref_lam)
        n_lam = ref_lam.size
        if not n_lam:
            raise ValueError("Not enough data to construct output WCS.")

        x_slit = np.zeros(n_lam)
        lam = 1e-6 * ref_lam

        # Find the spatial pixel scale:
        y_slit_min, y_slit_max = self._max_virtual_slit_extent(all_wcs, targ_ra, targ_dec)

        nsampl = 50
        xy_min = s2d(
            nsampl * [0],
            nsampl * [y_slit_min],
            lam[(tuple((i * n_lam) // nsampl for i in range(nsampl)), )]
        )
        xy_max = s2d(
            nsampl * [0],
            nsampl * [y_slit_max],
            lam[(tuple((i * n_lam) // nsampl for i in range(nsampl)), )]
        )

        good = np.logical_and(np.isfinite(xy_min), np.isfinite(xy_max))
        if not np.any(good):
            raise ValueError("Error estimating output WCS pixel scale.")

        xy1 = s2d(x_slit, np.full(n_lam, refmodel.slit_ymin), lam)
        xy2 = s2d(x_slit, np.full(n_lam, refmodel.slit_ymax), lam)
        xylen = np.nanmax(np.linalg.norm(np.array(xy1) - np.array(xy2), axis=0)) + 1
        pscale = (refmodel.slit_ymax - refmodel.slit_ymin) / xylen

        # compute image span along Y-axis (length of the slit in the detector plane)
        # det_slit_span = np.linalg.norm(np.subtract(xy_max, xy_min))
        det_slit_span = np.nanmax(np.linalg.norm(np.subtract(xy_max, xy_min), axis=0))
        ny = int(np.ceil(det_slit_span * self.pscale_ratio + 0.5)) + 1

        border = 0.5 * (ny - det_slit_span * self.pscale_ratio) - 0.5

        if xy_min[1][1] < xy_max[1][1]:
            y_slit_model = Linear1D(
                slope=pscale / self.pscale_ratio,
                intercept=y_slit_min - border * pscale * self.pscale_ratio
            )
        else:
            y_slit_model = Linear1D(
                slope=-pscale / self.pscale_ratio,
                intercept=y_slit_max + border * pscale * self.pscale_ratio
            )

        # extrapolate 1/2 pixel at the edges and make tabular model w/inverse:
        lam = lam.tolist()
        pixel_coord = list(range(n_lam))

        if len(pixel_coord) > 1:
            # left:
            slope = (lam[1] - lam[0]) / pixel_coord[1]
            lam.insert(0, -0.5 * slope + lam[0])
            pixel_coord.insert(0, -0.5)
            # right:
            slope = (lam[-1] - lam[-2]) / (pixel_coord[-1] - pixel_coord[-2])
            lam.append(slope * (pixel_coord[-1] + 0.5) + lam[-2])
            pixel_coord.append(pixel_coord[-1] + 0.5)

        else:
            lam = 3 * lam
            pixel_coord = [-0.5, 0, 0.5]

        wavelength_transform = Tabular1D(points=pixel_coord,
                                         lookup_table=lam,
                                         bounds_error=False, fill_value=np.nan)
        wavelength_transform.inverse = Tabular1D(points=lam,
                                                 lookup_table=pixel_coord,
                                                 bounds_error=False,
                                                 fill_value=np.nan)
        self.data_size = (ny, len(ref_lam))

        # Construct the final transform.
        # First coordinate is set to 0 to represent the "horizontal" center
        # of the slit (if we imagine slit to be vertical in the usual X-Y 2D
        # cartesian frame):
        mapping = Mapping((0, 1, 0))
        inv_mapping = Mapping((2, 1))
        inv_mapping.inverse = mapping
        mapping.inverse = inv_mapping
        zero_model = Const1D(0)
        zero_model.inverse = zero_model
        det2slit = mapping | zero_model & y_slit_model & wavelength_transform

        # Create coordinate frames
        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        slit_spatial = cf.Frame2D(name='slit_spatial', axes_order=(0, 1),
                                  unit=("", ""), axes_names=('x_slit', 'y_slit'))
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
                                unit=(u.micron,), axes_names=('wavelength',))
        slit_frame = cf.CompositeFrame([slit_spatial, spec], name='slit_frame')
        sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                                reference_frame=coord.ICRS())
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, det2slit), (slit_frame, s2w), (world, None)]
        output_wcs = WCS(pipeline)

        # Compute bounding box and output array shape.  Add one to the y (slit)
        # height to account for the half pixel at top and bottom due to pixel
        # coordinates being centers of pixels
        bounding_box = resample_utils.wcs_bbox_from_shape(self.data_size)
        output_wcs.bounding_box = bounding_box
        output_wcs.array_shape = self.data_size

        return output_wcs

    def _max_virtual_slit_extent(self, wcs_list, target_ra, target_dec):
        """
        Compute min & max slit coordinates for all nods in the "virtual"
        slit frame.

        NOTE: this code, potentially, might have troubles dealing
              with large dithers such that ``target_ra`` and ``target_dec``
              may not be converted to slit frame (i.e., result in ``NaN``).

              A more sophisticated algorithm may be needed to "stitch" large
              dithers. But then distortions may come into play.
        """
        y_slit_min = np.inf
        y_slit_max = -np.inf

        t0 = 0

        for wcs in wcs_list:
            d2s = wcs.get_transform('detector', 'slit_frame')
            w2s = wcs.get_transform('world', 'slit_frame')

            x, y = wcstools.grid_from_bounding_box(wcs.bounding_box)
            ra, dec, lam = wcs(x, y)

            good = np.logical_and(np.isfinite(ra), np.isfinite(dec))
            x = x[good]
            y = y[good]
            lm = lam[good]

            _, yslit, _ = d2s(x, y)

            # position of the target in the slit relative to its position
            # for the refence image:
            ts = w2s(target_ra, target_dec, np.mean(lm))[1] - t0

            if wcs is wcs_list[0]:
                t0 = ts
                ts = 0

            y_slit_min_i = np.min(yslit) - ts
            y_slit_max_i = np.max(yslit) - ts

            if y_slit_min_i < y_slit_min:
                y_slit_min = y_slit_min_i

            if y_slit_max_i > y_slit_max:
                y_slit_max = y_slit_max_i

        return y_slit_min, y_slit_max

    def build_interpolated_output_wcs(self, refmodel=None):
        """
        Create a spatial/spectral WCS output frame using all the input models

        Creates output frame by linearly fitting RA, Dec along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Parameters
        ----------
        refmodel : `~jwst.datamodels.JwstDataModel`
            The reference input image from which the fiducial WCS is created.
            If not specified, the first image in self.input_models is used.

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

        for im, model in enumerate(self.input_models):
            wcs = model.meta.wcs
            bbox = wcs.bounding_box
            grid = wcstools.grid_from_bounding_box(bbox)
            ra, dec, lam = np.array(wcs(*grid))
            # Handle vertical (MIRI) or horizontal (NIRSpec) dispersion.  The
            # following 2 variables are 0 or 1, i.e. zero-indexed in x,y WCS order
            spectral_axis = find_dispersion_axis(model)
            spatial_axis = spectral_axis ^ 1

            # Compute the wavelength array, trimming NaNs from the ends
            # In many cases, a whole slice is NaNs, so ignore those warnings
            warnings.simplefilter("ignore")
            wavelength_array = np.nanmedian(lam, axis=spectral_axis)
            warnings.resetwarnings()
            wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

            # We need to estimate the spatial sampling to use for the output WCS.
            # Tt is assumed the spatial sampling is the same for all the input
            # models. So we can use the first input model to set the spatial
            # sampling.

            # Steps to do this for first input model:
            # 1. find the middle of the spectrum in wavelength
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
                lam_center_index = int((bbox[spectral_axis][1] -
                                        bbox[spectral_axis][0]) / 2)
                if spatial_axis == 0:  # MIRI LRS, the WCS x axis is spatial
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
                warnings.simplefilter("ignore")
                # at this center of slit find x,y tangent projection - x_tan, y_tan
                x_tan, y_tan = undist2sky1.inverse(ra, dec)
                warnings.resetwarnings()

                # pull out data from center
                if spectral_axis == 0:  # MIRI LRS, the WCS x axis is spatial
                    x_tan_array = x_tan.T[lam_center_index]
                    y_tan_array = y_tan.T[lam_center_index]
                else:
                    x_tan_array = x_tan[lam_center_index]
                    y_tan_array = y_tan[lam_center_index]

                x_tan_array = x_tan_array[~np.isnan(x_tan_array)]
                y_tan_array = y_tan_array[~np.isnan(y_tan_array)]

                # estimate the spatial sampling
                fitter = LinearLSQFitter()
                fit_model = Linear1D()
                xstop = x_tan_array.shape[0] / self.pscale_ratio
                xstep = 1 / self.pscale_ratio
                ystop = y_tan_array.shape[0] / self.pscale_ratio
                ystep = 1 / self.pscale_ratio
                pix_to_xtan = fitter(fit_model, np.arange(0, xstop, xstep), x_tan_array)
                pix_to_ytan = fitter(fit_model, np.arange(0, ystop, ystep), y_tan_array)

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
        if self.input_models[0].meta.exposure.type == 'MIR_LRS-FIXEDSLIT':
            wavelength_array = np.flip(wavelength_array, axis=None)

        step = 1 / self.pscale_ratio
        stop = wavelength_array.shape[0] / self.pscale_ratio
        points = np.arange(0, stop, step)
        pix_to_wavelength = Tabular1D(points=points,
                                      lookup_table=wavelength_array,
                                      bounds_error=False, fill_value=None,
                                      name='pix2wavelength')

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is descending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(0, stop, step)

        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(points=points,
                                              lookup_table=lookup_table,
                                              bounds_error=False, fill_value=None,
                                              name='wavelength2pix')

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))

        # Sometimes the slit is perpendicular to the RA or Dec axis.
        # For example, if the slit is perpendicular to RA, that means
        # the slope of pix_to_xtan will be nearly zero, so make sure
        # mapping.inverse uses pix_to_ytan.inverse.  The auto definition
        # of mapping.inverse is to use the 2nd spatial coordinate, i.e. Dec.

        swap_xy = np.isclose(pix_to_xtan.slope, 0, atol=1e-8)
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
        if len(self.input_models) == 1:  # single model use ra_center_pt to be consistent
            # with how resample was done before
            ra_center_final = ra_center_pt
            dec_center_final = dec_center_pt

        native2celestial = RotateNative2Celestial(ra_center_final, dec_center_final, 180)
        undist2sky = tan | native2celestial
        # find the spatial size of the output - same in x,y
        if swap_xy:
            _, x_tan_all = undist2sky.inverse(all_ra, all_dec)
            pix_to_tan_slope = pix_to_ytan.slope
        else:
            x_tan_all, _ = undist2sky.inverse(all_ra, all_dec)
            pix_to_tan_slope = pix_to_xtan.slope

        x_min = np.amin(x_tan_all)
        x_max = np.amax(x_tan_all)
        x_size = int(np.ceil((x_max - x_min) / np.absolute(pix_to_tan_slope)))
        if swap_xy:
            pix_to_ytan.intercept = -0.5 * (x_size - 1) * pix_to_ytan.slope
        else:
            pix_to_xtan.intercept = -0.5 * (x_size - 1) * pix_to_xtan.slope

        # single model use size of x_tan_array
        # to be consistent with method before
        if len(self.input_models) == 1:
            x_size = len(x_tan_array)

        # define the output wcs
        transform = mapping | (pix_to_xtan & pix_to_ytan | undist2sky) & pix_to_wavelength

        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                                reference_frame=coord.ICRS())
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
                                unit=(u.micron,), axes_names=('wavelength',))
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, transform),
                    (world, None)]

        output_wcs = WCS(pipeline)

        # compute the output array size in WCS axes order, i.e. (x, y)
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = int(np.ceil(len(wavelength_array) / self.pscale_ratio))
        output_array_size[spatial_axis] = int(np.ceil(x_size / self.pscale_ratio))

        # turn the size into a numpy shape in (y, x) order
        output_wcs.array_shape = output_array_size[::-1]
        output_wcs.pixel_shape = output_array_size
        bounding_box = resample_utils.wcs_bbox_from_shape(output_array_size[::-1])
        output_wcs.bounding_box = bounding_box

        return output_wcs

    def build_nirspec_lamp_output_wcs(self):
        """
        Create a spatial/spectral WCS output frame for NIRSpec lamp mode

        Creates output frame by linearly fitting x_msa, y_msa along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Returns
        -------
        output_wcs : `~gwcs.WCS` object
            A gwcs WCS object defining the output frame WCS
        """
        model = self.input_models[0]
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
        warnings.simplefilter("ignore")
        wavelength_array = np.nanmedian(lam, axis=spectral_axis)
        warnings.resetwarnings()
        wavelength_array = wavelength_array[~np.isnan(wavelength_array)]

        # Find the center ra and dec for this slit at central wavelength
        lam_center_index = int((bbox[spectral_axis][1] -
                                bbox[spectral_axis][0]) / 2)
        x_msa_array = x_msa.T[lam_center_index]
        y_msa_array = y_msa.T[lam_center_index]
        x_msa_array = x_msa_array[~np.isnan(x_msa_array)]
        y_msa_array = y_msa_array[~np.isnan(y_msa_array)]

        # Estimate and fit the spatial sampling
        fitter = LinearLSQFitter()
        fit_model = Linear1D()
        xstop = x_msa_array.shape[0] / self.pscale_ratio
        xstep = 1 / self.pscale_ratio
        ystop = y_msa_array.shape[0] / self.pscale_ratio
        ystep = 1 / self.pscale_ratio
        pix_to_x_msa = fitter(fit_model, np.arange(0, xstop, xstep), x_msa_array)
        pix_to_y_msa = fitter(fit_model, np.arange(0, ystop, ystep), y_msa_array)

        step = 1 / self.pscale_ratio
        stop = wavelength_array.shape[0] / self.pscale_ratio
        points = np.arange(0, stop, step)
        pix_to_wavelength = Tabular1D(points=points,
                                      lookup_table=wavelength_array,
                                      bounds_error=False, fill_value=None,
                                      name='pix2wavelength')

        # Tabular models need an inverse explicitly defined.
        # If the wavelength array is descending instead of ascending, both
        # points and lookup_table need to be reversed in the inverse transform
        # for scipy.interpolate to work properly
        points = wavelength_array
        lookup_table = np.arange(0, stop, step)

        if not np.all(np.diff(wavelength_array) > 0):
            points = points[::-1]
            lookup_table = lookup_table[::-1]
        pix_to_wavelength.inverse = Tabular1D(points=points,
                                              lookup_table=lookup_table,
                                              bounds_error=False, fill_value=None,
                                              name='wavelength2pix')

        # For the input mapping, duplicate the spatial coordinate
        mapping = Mapping((spatial_axis, spatial_axis, spectral_axis))
        mapping.inverse = Mapping((2, 1))

        # The final transform
        # define the output wcs
        transform = mapping | pix_to_x_msa & pix_to_y_msa & pix_to_wavelength

        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        sky = cf.Frame2D(name=f'resampled_{model.meta.wcs.output_frame.name}', axes_order=(0, 1))
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,),
                                unit=(u.micron,), axes_names=('wavelength',))
        world = cf.CompositeFrame([sky, spec], name='world')

        pipeline = [(det, transform),
                    (world, None)]

        output_wcs = WCS(pipeline)

        # Compute the output array size and bounding box
        output_array_size = [0, 0]
        output_array_size[spectral_axis] = int(np.ceil(len(wavelength_array) / self.pscale_ratio))
        x_size = len(x_msa_array)
        output_array_size[spatial_axis] = int(np.ceil(x_size / self.pscale_ratio))
        # turn the size into a numpy shape in (y, x) order
        output_wcs.array_shape = output_array_size[::-1]
        output_wcs.pixel_shape = output_array_size
        bounding_box = resample_utils.wcs_bbox_from_shape(output_array_size[::-1])
        output_wcs.bounding_box = bounding_box

        return output_wcs


def find_dispersion_axis(refmodel):
    """
    Find the dispersion axis (0-indexed) of the given 2D wavelength array
    """
    dispaxis = refmodel.meta.wcsinfo.dispersion_direction
    # Change from 1 --> X and 2 --> Y to 0 --> X and 1 --> Y.
    return dispaxis - 1


def _spherical_sep(j, k, wcs, xyz_ref):
    """
    Objective function that computes the angle between two points
    on the sphere for small separations.
    """
    ra, dec, _ = wcs(k, j, with_bounding_box=False)
    return 1 - np.dot(_S2C(ra, dec), xyz_ref)


def _find_nirspec_output_sampling_wavelengths(wcs_list, targ_ra, targ_dec, mode='median'):
    assert mode in ['median', 'fast', 'accurate']
    refwcs = wcs_list[0]
    bbox = refwcs.bounding_box

    grid = wcstools.grid_from_bounding_box(bbox)
    ra, dec, lambdas = refwcs(*grid)

    if mode == 'median':
        ref_lam = sorted(np.nanmedian(lambdas[:, np.any(np.isfinite(lambdas), axis=0)], axis=0))
    else:
        ref_lam, _, _ = _find_nirspec_sampling_wavelengths(
            refwcs,
            targ_ra, targ_dec,
            ra, dec,
            fast=mode == 'fast'
        )

    lam1 = ref_lam[0]
    lam2 = ref_lam[-1]

    min_delta = np.fabs(np.ediff1d(ref_lam).min())

    image_lam = []
    for w in wcs_list[1:]:
        bbox = w.bounding_box
        grid = wcstools.grid_from_bounding_box(bbox)
        ra, dec, lambdas = w(*grid)
        if mode == 'median':
            lam = sorted(np.nanmedian(lambdas[:, np.any(np.isfinite(lambdas), axis=0)], axis=0))
        else:
            lam, _, _ = _find_nirspec_sampling_wavelengths(
                w,
                targ_ra, targ_dec,
                ra, dec,
                fast=mode == 'fast'
            )
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
        for k, (lam, lmin, lmax) in enumerate(image_lam):
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


def _find_nirspec_sampling_wavelengths(wcs, ra0, dec0, ra, dec, fast=True):
    xdet = []
    lms = []
    ys = []
    skipped = []

    eps = 10 * np.finfo(float).eps

    xyz_ref = _S2C(ra0, dec0)
    ymax, xmax = ra.shape
    good = np.logical_and(np.isfinite(ra), np.isfinite(dec))

    j0 = 0

    for k in range(xmax):
        if not any(good[:, k]):
            if xdet:
                skipped.append(k)
            continue

        idx = np.flatnonzero(good[:, k]).tolist()
        if j0 in idx:
            i = idx.index(j0)
        else:
            i = 0
            j0 = idx[0]

        dmin = _spherical_sep(j0, k, wcs, xyz_ref)

        for j in idx[i + 1:]:
            d = _spherical_sep(j, k, wcs, xyz_ref)
            if d < dmin:
                dmin = d
                j0 = j
            elif d > dmin:
                break

        for j in idx[max(i - 1, 0):None if i else 0:-1]:
            d = _spherical_sep(j, k, wcs, xyz_ref)
            if d < dmin:
                dmin = d
                j0 = j
            elif d > dmin:
                break

        if j0 == 0 or not good[j0 - 1, k]:
            j1 = j0 - 0.49999
        else:
            j1 = j0 - 0.99999

        if j0 == ymax - 1 or not good[j0 + 1, k]:
            j2 = j0 + 0.49999
        else:
            j2 = j0 + 0.99999

        if fast:
            # parabolic minimization:
            f0 = dmin
            f1 = _spherical_sep(j1, k, wcs, xyz_ref)
            if not np.isfinite(f1):
                # give another try with 1/2 step:
                j1 = 0.5 * (j1 + j0)
                f1 = _spherical_sep(j1, k, wcs, xyz_ref)

            f2 = _spherical_sep(j2, k, wcs, xyz_ref)
            if not np.isfinite(f2):
                # give another try with 1/2 step:
                j2 = 0.5 * (j2 + j0)
                f2 = _spherical_sep(j2, k, wcs, xyz_ref)

            if np.isfinite(f1) and np.isfinite(f2):
                dn = (j0 - j1) * (f0 - f2) - (j0 - j2) * (f0 - f1)
                if np.abs(dn) < eps:
                    jmin = j0
                else:
                    jmin = j0 - 0.5 * ((j0 - j1)**2 * (f0 - f2) -
                                       (j0 - j2)**2 * (f0 - f1)) / dn
                    jmin = max(min(jmin, j2), j1)
            else:
                jmin = j0

        else:
            r = minimize_scalar(
                _spherical_sep,
                method='golden',
                bracket=(j1, j2),
                args=(k, wcs, xyz_ref),
                tol=None,
                options={'maxiter': 10, 'xtol': 1e-2 / (j0 + 1)}
            )
            jmin = r['x']
            if np.isfinite(jmin):
                jmin = max(min(jmin, j2), j1)
            else:
                jmin = j0

        targ_lam = wcs(k, jmin)[-1]
        if not np.isfinite(targ_lam):
            targ_lam = wcs(k, j0)[-1]

        if not np.isfinite(targ_lam):
            if xdet:
                skipped.append(k)
            continue

        lms.append(targ_lam)
        ys.append(jmin)
        xdet.append(k)

    skipped = [s for s in skipped if s <= xdet[-1]]
    if skipped and skipped[0] < xdet[-1]:
        # there are columns with all pixels having invalid world,
        # coordinates. Fill the gaps using linear interpolation.
        raise NotImplementedError(
            "Support for discontinuous sampling was not implemented."
        )

    return lms, xdet, ys
