import logging

import numpy as np
from drizzle import util
from drizzle import cdrizzle
from spherical_geometry.polygon import SphericalPolygon

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from . import gwcs_drizzle
from . import resample_utils
from ..lib.basic_utils import bytes2human
from ..model_blender import blendmeta

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutputTooLargeError", "ResampleData"]


class OutputTooLargeError(RuntimeError):
    """Raised when the output is too large for in-memory instantiation"""


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
                 pixfrac=1.0, kernel="square", fillval="INDEF", wht_type="ivm",
                 good_bits=0, pscale_ratio=1.0, pscale=None, **kwargs):
        """
        Parameters
        ----------
        input_models : list of objects
            list of data models, one for each input image

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
        self.input_models = input_models
        self.output_filename = output
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
                output_pix_area = _compute_image_pixel_area(self.output_wcs)
                if output_pix_area is None:
                    raise ValueError(
                        "Unable to compute output pixel area from 'output_wcs'."
                    )
            else:
                output_pix_area = output_wcs.pixel_area

        else:
            # Define output WCS based on all inputs, including a reference WCS:
            self.output_wcs = resample_utils.make_output_wcs(
                self.input_models,
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
            output_pix_area = (
                np.deg2rad(tr['cdelt1'].factor.value) *
                np.deg2rad(tr['cdelt2'].factor.value)
            )

        if pscale is None:
            pscale = np.rad2deg(np.sqrt(output_pix_area))
            log.info(f'Computed output pixel scale: {3600.0 * pscale} arcsec.')

        self.pscale = pscale  # in deg

        log.debug('Output mosaic size: {}'.format(self.output_wcs.array_shape))
        can_allocate, required_memory = datamodels.util.check_memory_allocation(
            self.output_wcs.array_shape, kwargs['allowed_memory'], datamodels.ImageModel
        )
        if not can_allocate:
            raise OutputTooLargeError(
                f'Combined ImageModel size {self.output_wcs.array_shape} '
                f'requires {bytes2human(required_memory)}. '
                f'Model cannot be instantiated.'
            )
        self.blank_output = datamodels.ImageModel(tuple(self.output_wcs.array_shape))

        # update meta data and wcs
        self.blank_output.update(input_models[0])
        self.blank_output.meta.wcs = self.output_wcs
        self.blank_output.meta.photometry.pixelarea_steradians = output_pix_area
        self.blank_output.meta.photometry.pixelarea_arcsecsq = (
            output_pix_area * np.rad2deg(3600)**2
        )

        self.output_models = ModelContainer(open_models=False)

    def do_drizzle(self):
        """Pick the correct drizzling mode based on self.single
        """
        if self.single:
            return self.resample_many_to_many()
        else:
            return self.resample_many_to_one()

    def blend_output_metadata(self, output_model):
        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = output_model.meta.filename

        ignore_list = [
            'meta.photometry.pixelarea_steradians',
            'meta.photometry.pixelarea_arcsecsq',
        ]

        log.info('Blending metadata for {}'.format(output_file))
        blendmeta.blendmodels(
            output_model,
            inputs=self.input_models,
            output=output_file,
            ignore=ignore_list
        )

    def resample_many_to_many(self):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure, i.e. map NRCA5 and
        NRCB5 onto the same output image, as they image different areas of the
        sky.

        Used for outlier detection
        """
        for exposure in self.input_models.models_grouped:
            output_model = self.blank_output
            # Determine output file type from input exposure filenames
            # Use this for defining the output filename
            indx = exposure[0].meta.filename.rfind('.')
            output_type = exposure[0].meta.filename[indx:]
            output_root = '_'.join(exposure[0].meta.filename.replace(
                output_type, '').split('_')[:-1])
            output_model.meta.filename = f'{output_root}_outlier_i2d{output_type}'

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                            kernel=self.kernel, fillval=self.fillval)

            log.info(f"{len(exposure)} exposures to drizzle together")
            for img in exposure:
                img = datamodels.open(img)

                input_pixflux_area = img.meta.photometry.pixelarea_steradians
                if (input_pixflux_area and
                        'SPECTRAL' not in img.meta.wcs.output_frame.axes_type):
                    img.meta.wcs.array_shape = img.data.shape
                    input_pixel_area = _compute_image_pixel_area(img.meta.wcs)
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

                # TODO: should weight_type=None here?
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

                driz.add_image(
                    data,
                    img.meta.wcs,
                    iscale=iscale,
                    inwht=inwht,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax
                )
                del data
                img.close()

            if not self.in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.filename
                output_model.save(output_name)
                log.info(f"Saved model in {output_name}")
                self.output_models.append(output_name)
            else:
                self.output_models.append(output_model.copy())
            output_model.data *= 0.
            output_model.wht *= 0.

        return self.output_models

    def resample_many_to_one(self):
        """Resample and coadd many inputs to a single output.

        Used for stage 3 resampling
        """
        output_model = self.blank_output.copy()
        output_model.meta.filename = self.output_filename
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pointings = len(self.input_models.group_names)

        if self.blendheaders:
            self.blend_output_metadata(output_model)

        # Initialize the output with the wcs
        driz = gwcs_drizzle.GWCSDrizzle(output_model, pixfrac=self.pixfrac,
                                        kernel=self.kernel, fillval=self.fillval)

        log.info("Resampling science data")
        for img in self.input_models:
            input_pixflux_area = img.meta.photometry.pixelarea_steradians
            if (input_pixflux_area and
                    'SPECTRAL' not in img.meta.wcs.output_frame.axes_type):
                img.meta.wcs.array_shape = img.data.shape
                input_pixel_area = _compute_image_pixel_area(img.meta.wcs)
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

            inwht = resample_utils.build_driz_weight(img,
                                                     weight_type=self.weight_type,
                                                     good_bits=self.good_bits)
            # apply sky subtraction
            blevel = img.meta.background.level
            if not img.meta.background.subtracted and blevel is not None:
                data = img.data - blevel
            else:
                data = img.data.copy()

            xmin, xmax, ymin, ymax = resample_utils._resample_range(
                data.shape,
                img.meta.wcs.bounding_box
            )

            driz.add_image(
                data,
                img.meta.wcs,
                iscale=iscale,
                inwht=inwht,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax
            )
            del data, inwht

        # Resample variances array in self.input_models to output_model
        self.resample_variance_array("var_rnoise", output_model, iscale=iscale)
        self.resample_variance_array("var_poisson", output_model, iscale=iscale)
        self.resample_variance_array("var_flat", output_model, iscale=iscale)
        output_model.err = np.sqrt(
            np.nansum(
                [
                    output_model.var_rnoise,
                    output_model.var_poisson,
                    output_model.var_flat
                ],
                axis=0
            )
        )

        self.update_exposure_times(output_model)
        self.output_models.append(output_model)

        return self.output_models

    def resample_variance_array(self, name, output_model, iscale):
        """Resample variance arrays from self.input_models to the output_model

        Resample the ``name`` variance array to the same name in output_model,
        using a cumulative sum.

        This modifies output_model in-place.
        """
        output_wcs = output_model.meta.wcs
        inverse_variance_sum = np.full_like(output_model.data, np.nan)

        log.info(f"Resampling {name}")
        for model in self.input_models:
            variance = getattr(model, name)
            if variance is None or variance.size == 0:
                log.debug(
                    f"No data for '{name}' for model "
                    f"{repr(model.meta.filename)}. Skipping ..."
                )
                continue

            elif variance.shape != model.data.shape:
                log.warning(
                    f"Data shape mismatch for '{name}' for model "
                    f"{repr(model.meta.filename)}. Skipping ..."
                )
                continue

            # Make input weight map of unity where there is science data
            inwht = resample_utils.build_driz_weight(
                model,
                weight_type=None,
                good_bits=self.good_bits
            )

            resampled_variance = np.zeros_like(output_model.data)
            outwht = np.zeros_like(output_model.data)
            outcon = np.zeros_like(output_model.con)

            xmin, xmax, ymin, ymax = resample_utils._resample_range(
                variance.shape,
                model.meta.wcs.bounding_box
            )

            # Resample the variance array. Fill "unpopulated" pixels with NaNs.
            self.drizzle_arrays(
                variance,
                inwht,
                model.meta.wcs,
                output_wcs,
                resampled_variance,
                outwht,
                outcon,
                iscale=iscale**2,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=np.nan,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax
            )

            # Add the inverse of the resampled variance to a running sum.
            # Update only pixels (in the running sum) with valid new values:
            mask = resampled_variance > 0

            inverse_variance_sum[mask] = np.nansum(
                [inverse_variance_sum[mask], np.reciprocal(resampled_variance[mask])],
                axis=0
            )

        # We now have a sum of the inverse resampled variances.  We need the
        # inverse of that to get back to units of variance.
        output_variance = np.reciprocal(inverse_variance_sum)

        setattr(output_model, name, output_variance)

    def update_exposure_times(self, output_model):
        """Modify exposure time metadata in-place"""
        total_exposure_time = 0.
        exposure_times = {'start': [], 'end': []}
        duration = 0.0
        for exposure in self.input_models.models_grouped:
            total_exposure_time += exposure[0].meta.exposure.exposure_time
            exposure_times['start'].append(exposure[0].meta.exposure.start_time)
            exposure_times['end'].append(exposure[0].meta.exposure.end_time)
            duration += exposure[0].meta.exposure.duration

        # Update some basic exposure time values based on output_model
        output_model.meta.exposure.exposure_time = total_exposure_time
        output_model.meta.exposure.start_time = min(exposure_times['start'])
        output_model.meta.exposure.end_time = max(exposure_times['end'])
        output_model.meta.resample.product_exposure_time = total_exposure_time

        # Update other exposure time keywords:
        # XPOSURE (identical to the total effective exposure time, EFFEXPTM)
        xposure = total_exposure_time
        output_model.meta.exposure.effective_exposure_time = xposure
        # DURATION (identical to TELAPSE, elapsed time)
        output_model.meta.exposure.duration = duration
        output_model.meta.exposure.elapsed_exposure_time = duration

    @staticmethod
    def drizzle_arrays(insci, inwht, input_wcs, output_wcs, outsci, outwht,
                       outcon, uniqid=1, xmin=0, xmax=0, ymin=0, ymax=0,
                       iscale=1.0, pixfrac=1.0, kernel='square',
                       fillval="INDEF", wtscale=1.0):
        """
        Low level routine for performing 'drizzle' operation on one image.

        The interface is compatible with STScI code. All images are Python
        ndarrays, instead of filenames. File handling (input and output) is
        performed by the calling routine.

        Parameters
        ----------

        insci : 2d array
            A 2d numpy array containing the input image to be drizzled.

        inwht : 2d array
            A 2d numpy array containing the pixel by pixel weighting.
            Must have the same dimensions as insci. If none is supplied,
            the weighting is set to one.

        input_wcs : gwcs.WCS object
            The world coordinate system of the input image.

        output_wcs : gwcs.WCS object
            The world coordinate system of the output image.

        outsci : 2d array
            A 2d numpy array containing the output image produced by
            drizzling. On the first call it should be set to zero.
            Subsequent calls it will hold the intermediate results.  This
            is modified in-place.

        outwht : 2d array
            A 2d numpy array containing the output counts. On the first
            call it should be set to zero. On subsequent calls it will
            hold the intermediate results.  This is modified in-place.

        outcon : 2d or 3d array, optional
            A 2d or 3d numpy array holding a bitmap of which image was an input
            for each output pixel. Should be integer zero on first call.
            Subsequent calls hold intermediate results.  This is modified
            in-place.

        uniqid : int, optional
            The id number of the input image. Should be one the first time
            this function is called and incremented by one on each subsequent
            call.

        xmin : int, optional
            This and the following three parameters set a bounding rectangle
            on the input image. Only pixels on the input image inside this
            rectangle will have their flux added to the output image. Xmin
            sets the minimum value of the x dimension. The x dimension is the
            dimension that varies quickest on the image. All four parameters
            are zero based, counting starts at zero.

        xmax : int, optional
            Sets the maximum value of the x dimension on the bounding box
            of the input image. If ``xmax = 0``, no maximum will
            be set in the x dimension (all pixels in a row of the input image
            will be resampled).

        ymin : int, optional
            Sets the minimum value in the y dimension on the bounding box. The
            y dimension varies less rapidly than the x and represents the line
            index on the input image.

        ymax : int, optional
            Sets the maximum value in the y dimension. If ``ymax = 0``,
            no maximum will be set in the y dimension (all pixels in a column
            of the input image will be resampled).

        iscale : float, optional
            A scale factor to be applied to pixel intensities of the
            input image before resampling.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel: str, optional
            The name of the kernel used to combine the input. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: "square", "gaussian", "point", "tophat", "turbo", "lanczos2",
            and "lanczos3". The square kernel is the default.

        fillval: str, optional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of INDEF does not set a value.

        Returns
        -------
        A tuple with three values: a version string, the number of pixels
        on the input image that do not overlap the output image, and the
        number of complete lines on the input image that do not overlap the
        output input image.

        """

        # Insure that the fillval parameter gets properly interpreted for use with tdriz
        if util.is_blank(str(fillval)):
            fillval = 'INDEF'
        else:
            fillval = str(fillval)

        if insci.dtype > np.float32:
            insci = insci.astype(np.float32)

        # Add input weight image if it was not passed in
        if inwht is None:
            inwht = np.ones_like(insci)

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int((uniqid - 1) / 32)

        # Check if the context image has this many planes
        if outcon.ndim == 3:
            nplanes = outcon.shape[0]
        elif outcon.ndim == 2:
            nplanes = 1
        else:
            nplanes = 0

        if nplanes <= planeid:
            raise IndexError("Not enough planes in drizzle context image")

        # Alias context image to the requested plane if 3d
        if outcon.ndim == 3:
            outcon = outcon[planeid]

        # Compute the mapping between the input and output pixel coordinates
        # for use in drizzle.cdrizzle.tdriz
        pixmap = resample_utils.calc_gwcs_pixmap(input_wcs, output_wcs, insci.shape)

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {insci.shape}")
        log.debug(f"Output Sci shape: {outsci.shape}")

        log.info(f"Drizzling {insci.shape} --> {outsci.shape}")

        _vers, _nmiss, _nskip = cdrizzle.tdriz(
            insci, inwht, pixmap,
            outsci, outwht, outcon,
            uniqid=uniqid,
            xmin=xmin, xmax=xmax,
            ymin=ymin, ymax=ymax,
            scale=iscale,
            pixfrac=pixfrac,
            kernel=kernel,
            in_units="cps",
            expscale=1.0,
            wtscale=wtscale,
            fillstr=fillval
        )


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


def _compute_image_pixel_area(wcs):
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
                ymin, xmax, ymax, xmin = limits
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
