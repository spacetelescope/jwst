import logging
import warnings

import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import (Mapping, Tabular1D, Linear1D,
                                     Pix2Sky_TAN, RotateNative2Celestial)
from astropy.modeling.fitting import LinearLSQFitter
from gwcs import wcstools, WCS
from gwcs import coordinate_frames as cf

from ..assign_wcs.util import wrap_ra
from .. import datamodels
from . import resample_utils
from .resample import ResampleData


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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
                 pixfrac=1.0, kernel="square", fillval=0, weight_type="ivm",
                 good_bits=0, pscale_ratio=1.0, **kwargs):
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
        self.single = single
        self.blendheaders = blendheaders
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = weight_type
        self.good_bits = good_bits

        # Define output WCS based on all inputs, including a reference WCS
        self.output_wcs = self.build_interpolated_output_wcs()
        self.blank_output = datamodels.SlitModel(self.data_size)
        self.blank_output.update(self.input_models[0])
        self.blank_output.meta.wcs = self.output_wcs
        self.output_models = datamodels.ModelContainer()

    def build_interpolated_output_wcs(self, refmodel=None):
        """
        Create a spatial/spectral WCS output frame using all the input models

        Creates output frame by linearly fitting RA, Dec along the slit and
        producing a lookup table to interpolate wavelengths in the dispersion
        direction.

        Parameters
        ----------
        refmodel : `~jwst.datamodels.DataModel`
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
        # use first model to initialize wavelenth array
        # append wavelengths that fall outside the endpoint of
        # of wavelength array when looping over additional data

        all_wavelength = []
        all_ra_slit = []
        all_dec_slit = []

        for im, model in enumerate(self.input_models):
            wcs = model.meta.wcs
            bb = wcs.bounding_box
            grid = wcstools.grid_from_bounding_box(bb)
            ra, dec, lam = np.array(wcs(*grid))
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
            # first input model sets intializes wavelength array and defines
            # the spatial scale of the output wcs
            if im == 0:
                all_wavelength = np.append(all_wavelength, wavelength_array)

                # find the center ra and dec for this slit at central wavelength
                lam_center_index = int((bb[spectral_axis][1] -
                                        bb[spectral_axis][0]) / 2)
                if spatial_axis == 0:  # MIRI LRS, the WCS x axis is spatial
                    ra_slice = ra[lam_center_index, :]
                    dec_slice = dec[lam_center_index, :]
                else:
                    ra_slice = ra[:, lam_center_index]
                    dec_slice = dec[:, lam_center_index]
                # wrap RA if near zero
                ra_center_pt = np.nanmean(wrap_ra(ra_slice))
                dec_center_pt = np.nanmean(dec_slice)

                if resample_utils.is_sky_like(model.meta.wcs.output_frame):
                    # convert ra and dec to tangent projection
                    tan = Pix2Sky_TAN()
                    native2celestial = RotateNative2Celestial(ra_center_pt, dec_center_pt, 180)
                    undist2sky1 = tan | native2celestial
                    # Filter out RuntimeWarnings due to computed NaNs in the WCS
                    warnings.simplefilter("ignore")
                    # at this center of slit find x,y tangent projection - x_tan, y_tan
                    x_tan, y_tan = undist2sky1.inverse(ra, dec)
                    warnings.resetwarnings()
                else:
                    # for non sky-like output frames, no need to do tangent plane projections
                    # but we still use the same variables
                    x_tan, y_tan = ra, dec

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
        # dispersion direction matches the disperion direction on the detector.
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
        # If the wavelength array is decending instead of ascending, both
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

        if np.isclose(pix_to_ytan.slope, 0, atol=1e-8):
            mapping_tuple = (0, 1)
            # Account for vertical or horizontal dispersion on detector
            if spatial_axis:
                mapping.inverse = Mapping(mapping_tuple[::-1])
            else:
                mapping.inverse = Mapping(mapping_tuple)

        # The final transform
        # redefine the ra, dec center tangent point to include all data

        # check if all_ra crosses 0 degress - this makes it hard to
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

        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            native2celestial = RotateNative2Celestial(ra_center_final, dec_center_final, 180)
            undist2sky = tan | native2celestial
            # find the spatial size of the output - same in x,y
            x_tan_all, _ = undist2sky.inverse(all_ra, all_dec)
        else:
            x_tan_all, _ = all_ra, all_dec
        x_min = np.amin(x_tan_all)
        x_max = np.amax(x_tan_all)
        x_size = int(np.ceil((x_max - x_min) / np.absolute(pix_to_xtan.slope)))

        # single model use size of x_tan_array
        # to be consistent with method before
        if len(self.input_models) == 1:
            x_size = len(x_tan_array)

        # define the output wcs
        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            transform = mapping | (pix_to_xtan & pix_to_ytan | undist2sky) & pix_to_wavelength
        else:
            transform = mapping | pix_to_xtan & pix_to_ytan & pix_to_wavelength

        det = cf.Frame2D(name='detector', axes_order=(0, 1))
        if resample_utils.is_sky_like(model.meta.wcs.output_frame):
            sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                                    reference_frame=coord.ICRS())
        else:
            sky = cf.Frame2D(name=f'resampled_{model.meta.wcs.output_frame.name}', axes_order=(0, 1))
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
        self.data_size = tuple(output_array_size[::-1])
        bounding_box = resample_utils.wcs_bbox_from_shape(self.data_size)
        output_wcs.bounding_box = bounding_box

        return output_wcs


def find_dispersion_axis(refmodel):
    """
    Find the dispersion axis (0-indexed) of the given 2D wavelength array
    """
    dispaxis = refmodel.meta.wcsinfo.dispersion_direction
    # Change from 1 --> X and 2 --> Y to 0 --> X and 1 --> Y.
    return dispaxis - 1
