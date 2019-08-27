from distutils.version import LooseVersion
import logging
import math

import numpy as np
from astropy import units as u
import photutils
from photutils import CircularAperture, CircularAnnulus, \
                      RectangularAperture, aperture_photometry

from .. import datamodels
from ..datamodels import dqflags
from . import spec_wcs
from . import util

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# These values are used to indicate whether the input reference file
# (if any) is JSON or IMAGE.
FILE_TYPE_JSON = "JSON"
FILE_TYPE_IMAGE = "IMAGE"

# This is to prevent calling offset_from_offset multiple times for
# multi-integration data.
OFFSET_NOT_ASSIGNED_YET = "not assigned yet"

# This is intended to be larger than any possible distance (in pixels)
# between the target and any point in the image; used by locn_from_wcs().
HUGE_DIST = 1.e10

def ifu_extract1d(input_model, ref_dict, source_type, subtract_background,
                  apply_nod_offset):
    """Extract a 1-D spectrum from an IFU cube.

    Parameters
    ----------
    input_model : JWST data model for an IFU cube (IFUCubeModel)
        The input model.

    ref_dict : dict
        The contents of the reference file.

    source_type : string
        "point" or "extended"

    subtract_background : bool or None
        User supplied flag indicating whether the background should be subtracted.
        If None, the value in the extract_1d reference file will be used.
        If not None, this parameter overrides the value in the
        extract_1d reference file.

    apply_nod_offset : bool or None
        If True, the target and background positions specified in the
        reference file (or the default position, if there is no reference
        file) will be shifted to account for nod and/or dither offset.

    Returns
    -------
    output_model : MultiSpecModel
        This will contain the extracted spectrum.
    """

    if not isinstance(input_model, datamodels.IFUCubeModel):
        log.error("Expected an IFU cube.")
        raise RuntimeError("Expected an IFU cube.")

    source_type = source_type.lower()
    if source_type != "point" and source_type != "extended":
        log.warning("source_type was '%s', setting to 'point'.", source_type)
        source_type = "point"
    else:
        log.info("source_type = %s", source_type)

    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)

    slitname = input_model.meta.exposure.type
    if slitname is None:
        slitname = "ANY"

    extract_params = get_extract_parameters(ref_dict, slitname)
    if subtract_background is not None:
        if subtract_background and source_type == "extended":
            subtract_background = False
            log.info("Turning off background subtraction because "
                     "the source is extended.")
        extract_params['subtract_background'] = subtract_background

    if extract_params:
        if extract_params['ref_file_type'] == FILE_TYPE_JSON:
            (ra, dec, wavelength, temp_flux, background, npixels, dq) = \
                    extract_ifu(input_model, source_type, extract_params)
        else:                                   # FILE_TYPE_IMAGE
            (ra, dec, wavelength, temp_flux, background, npixels, dq) = \
                    image_extract_ifu(input_model, source_type, extract_params)
    else:
        log.critical('Missing extraction parameters.')
        raise ValueError('Missing extraction parameters.')

    # Check the units.  We expect either "MJy / sr" or "mJy / arcsec^2".
    try:
        bunit = input_model.meta.bunit_data
    except AttributeError:
        bunit = None
    megajanskys = True                          # initial value
    if bunit is not None:
        if bunit.find("M") < 0 or bunit.find("arcsec") >= 0:
            megajanskys = False

    # Convert the sum to an average, for surface brightness.
    npixels_temp = np.where(npixels > 0., npixels, 1.)
    surf_bright = temp_flux / npixels_temp
    background /= npixels_temp
    del npixels_temp
    if not megajanskys:
        # Not MJy / sr?  Then assume the units are mJy / arcsec^2.
        # Target spectrum in units of surface brightness.
        sb = surf_bright * u.mJy / (u.arcsec**2)
        surf_bright = sb.to(u.MJy / u.steradian).value
        # Background spectrum.
        bkg = background * u.mJy / (u.arcsec**2)
        background = bkg.to(u.MJy / u.steradian).value

    # Compute the solid angle of a pixel in steradians.
    pixel_solid_angle = util.pixel_area(input_model.meta.wcs,
                                        input_model.data.shape,
                                        False)
    if pixel_solid_angle is None:
        log.warning("Pixel solid angle could not be determined")
        pixel_solid_angle = 1.
    if megajanskys:
        # Convert flux from MJy / steradian to Jy.
        flux = temp_flux * pixel_solid_angle * 1.e6
    else:
        # Convert flux from mJy / arcsec**2 to Jy.
        psa = pixel_solid_angle * u.steradian
        pixel_solid_angle = psa.to(u.arcsec**2).value
        flux = temp_flux * pixel_solid_angle * 1.e-3
    del temp_flux

    error = np.zeros_like(flux)
    sb_error = np.zeros_like(flux)
    berror = np.zeros_like(flux)
    spec_dtype = datamodels.SpecModel().spec_table.dtype
    otab = np.array(list(zip(wavelength,
                             flux, error, surf_bright, sb_error,
                             dq, background, berror, npixels)),
                    dtype=spec_dtype)
    spec = datamodels.SpecModel(spec_table=otab)
    spec.meta.wcs = spec_wcs.create_spectral_wcs(ra, dec, wavelength)
    spec.spec_table.columns['wavelength'].unit = 'um'
    spec.spec_table.columns['flux'].unit = "Jy"
    spec.spec_table.columns['error'].unit = "Jy"
    spec.spec_table.columns['surf_bright'].unit = "MJy/sr"
    spec.spec_table.columns['sb_error'].unit = "MJy/sr"
    spec.spec_table.columns['background'].unit = "MJy/sr"
    spec.spec_table.columns['berror'].unit = "MJy/sr"
    spec.slit_ra = ra
    spec.slit_dec = dec
    if slitname is not None and slitname != "ANY":
        spec.name = slitname
    output_model.spec.append(spec)

    # See output_model.spec[0].meta.wcs instead.
    output_model.meta.wcs = None

    return output_model


def get_extract_parameters(ref_dict, slitname):
    """Read extraction parameters for an IFU.

    Parameters
    ----------
    ref_dict : dict
        The contents of the reference file.

    slitname : str
        The name of the slit, or "ANY".

    Returns
    -------
    dict
        The extraction parameters.
    """

    extract_params = {}

    if ('ref_file_type' not in ref_dict or
        ref_dict['ref_file_type'] == FILE_TYPE_JSON):
            extract_params['ref_file_type'] = FILE_TYPE_JSON
            for aper in ref_dict['apertures']:
                if 'id' in aper and aper['id'] != "dummy" and \
                   (aper['id'] == slitname or aper['id'] == "ANY" or
                    slitname == "ANY"):
                    region_type = aper.get("region_type", "target")
                    if region_type == "target":
                        extract_params['x_center'] = aper.get('x_center')
                        extract_params['y_center'] = aper.get('y_center')
                        extract_params['method'] = aper.get('method', 'exact')
                        extract_params['subpixels'] = aper.get('subpixels', 5)
                        extract_params['radius'] = aper.get('radius')
                        extract_params['subtract_background'] = \
                              aper.get('subtract_background', False)
                        extract_params['inner_bkg'] = aper.get('inner_bkg')
                        extract_params['outer_bkg'] = aper.get('outer_bkg')
                        extract_params['width'] = aper.get('width')
                        extract_params['height'] = aper.get('height')
                        # theta is in degrees (converted to radians later)
                        extract_params['theta'] = aper.get('theta', 0.)
                    break

    elif ref_dict['ref_file_type'] == FILE_TYPE_IMAGE:
        extract_params['ref_file_type'] = FILE_TYPE_IMAGE
        foundit = False
        for im in ref_dict['ref_model'].images:
            if (im.name is None or im.name == "ANY" or slitname == "ANY" or
                im.name == slitname):
                    extract_params['ref_image'] = im
                    foundit = True
                    break

        if not foundit:
            log.error("No match for slit name %s in reference image", slitname)
            raise RuntimeError("Specify slit name or use 'any' in ref image.")

    else:
        log.error("Reference file type %s not recognized",
                  ref_dict['ref_file_type'])
        raise RuntimeError("Reference file must be JSON or a FITS image.")

    return extract_params


def extract_ifu(input_model, source_type, extract_params):
    """This function does the extraction.

    Parameters
    ----------
    input_model : IFUCubeModel
        The input model.

    source_type : string
        "point" or "extended"

    extract_params : dict
        The extraction parameters for aperture photometry.

    Returns
    -------
    ra, dec : float
        ra and dec are the right ascension and declination respectively
        at the nominal center of the image.

    wavelength : ndarray, 1-D
        The wavelength in micrometers at each plane of the IFU cube.

    temp_flux : ndarray, 1-D
        The sum of the data values in the extraction aperture minus the
        sum of the data values in the background region (scaled by the
        ratio of areas), for each plane.
        The data values are in units of surface brightness, so this value
        isn't really the flux, it's an intermediate value.  Dividing by
        `npixels` (to compute the average) will give the value for the
        `surf_bright` (surface brightness) column, and multiplying by
        the solid angle of a pixel will give the flux for a point source.

    background : ndarray, 1-D
        The background count rate that was subtracted from the total
        source data values to get `temp_flux`.

    npixels : ndarray, 1-D, float64
        For each slice, this is the number of pixels that were added
        together to get `temp_flux`.

    dq : ndarray, 1-D, uint32
        The data quality array.
    """

    data = input_model.data
    shape = data.shape
    if len(shape) != 3:
        log.error("Expected a 3-D IFU cube; dimension is %d.", len(shape))
        raise RuntimeError("The IFU cube should be 3-D.")

    # We need to allocate temp_flux, background, npixels, and dq arrays
    # no matter what.  We may need to divide by npixels, so the default
    # is 1 rather than 0.
    temp_flux = np.zeros(shape[0], dtype=np.float64)
    background = np.zeros(shape[0], dtype=np.float64)
    npixels = np.ones(shape[0], dtype=np.float64)

    dq = np.zeros(shape[0], dtype=np.uint32)

    # For an extended target, the entire aperture will be extracted, so
    # it makes no sense to shift the extraction location.
    if source_type != "extended":
        ra_targ = input_model.meta.target.ra
        dec_targ = input_model.meta.target.dec
        locn = locn_from_wcs(input_model, ra_targ, dec_targ)
        if locn is None or np.isnan(locn[0]):
            log.warning("Couldn't determine pixel location from WCS, so "
                        "nod/dither correction will not be applied.")
            x_center = extract_params['x_center']
            y_center = extract_params['y_center']
            if x_center is None:
                x_center = float(shape[-1]) / 2.
            else:
                x_center = float(x_center)
            if y_center is None:
                y_center = float(shape[-2]) / 2.
            else:
                y_center = float(y_center)
        else:
            (x_center, y_center) = locn
            log.info("Using x_center = %g, y_center = %g, based on "
                     "TARG_RA and TARG_DEC.", x_center, y_center)

    method = extract_params['method']
    # subpixels is only needed if method = 'subpixel'.
    subpixels = extract_params['subpixels']

    subtract_background = extract_params['subtract_background']
    smaller_axis = float(min(shape[-2], shape[-1]))     # for defaults
    radius = None
    inner_bkg = None
    outer_bkg = None

    if source_type == 'extended':
        # Ignore any input parameters, and extract the whole image.
        width = float(shape[-1])
        height = float(shape[-2])
        x_center = width / 2. - 0.5
        y_center = height / 2. - 0.5
        theta = 0.
        subtract_background = False
    else:
        radius = extract_params['radius']
        if radius is None:
            radius = smaller_axis / 4.
        if subtract_background:
            inner_bkg = extract_params['inner_bkg']
            if inner_bkg is None:
                inner_bkg = radius
            outer_bkg = extract_params['outer_bkg']
            if outer_bkg is None:
                outer_bkg = min(inner_bkg * math.sqrt(2.),
                                smaller_axis / 2. - 1.)
            if inner_bkg <= 0. or outer_bkg <= 0. or inner_bkg >= outer_bkg:
                log.debug("Turning background subtraction off, due to "
                          "the values of inner_bkg and outer_bkg.")
                subtract_background = False
        width = None
        height = None
        theta = None

    log.debug("IFU 1-D extraction parameters:")
    log.debug("  x_center = %s", str(x_center))
    log.debug("  y_center = %s", str(y_center))
    if source_type == 'point':
        log.debug("  radius = %s", str(radius))
        log.debug("  subtract_background = %s", str(subtract_background))
        log.debug("  inner_bkg = %s", str(inner_bkg))
        log.debug("  outer_bkg = %s", str(outer_bkg))
        log.debug("  method = %s", method)
        if method == "subpixel":
            log.debug("  subpixels = %s", str(subpixels))
    else:
        log.debug("  width = %s", str(width))
        log.debug("  height = %s", str(height))
        log.debug("  theta = %s degrees", str(theta))
        log.debug("  subtract_background = %s", str(subtract_background))
        log.debug("  method = %s", method)
        if method == "subpixel":
            log.debug("  subpixels = %s", str(subpixels))

    x0 = float(shape[2]) / 2.
    y0 = float(shape[1]) / 2.
    (ra, dec, wavelength) = get_coordinates(input_model, x0, y0)

    position = (x_center, y_center)
    if source_type == 'point':
        aperture = CircularAperture(position, r=radius)
    else:
        aperture = RectangularAperture(position, width, height, theta)

    if subtract_background and inner_bkg is not None and outer_bkg is not None:
        annulus = CircularAnnulus(position, r_in=inner_bkg, r_out=outer_bkg)
    else:
        annulus = None

    # Compute the area of the aperture and possibly also of the annulus.
    normalization = 1.
    temp = np.ones(shape[-2:], dtype=np.float64)
    phot_table = aperture_photometry(temp, aperture,
                                     method=method, subpixels=subpixels)
    aperture_area = float(phot_table['aperture_sum'][0])
    if LooseVersion(photutils.__version__) >= '0.7':
        log.debug("aperture.area = %g; aperture_area = %g",
                  aperture.area, aperture_area)
    else:
        log.debug("aperture.area() = %g; aperture_area = %g",
                  aperture.area(), aperture_area)

    if subtract_background and annulus is not None:
        # Compute the area of the annulus.
        phot_table = aperture_photometry(temp, annulus,
                                         method=method, subpixels=subpixels)
        annulus_area = float(phot_table['aperture_sum'][0])
        if LooseVersion(photutils.__version__) >= '0.7':
            log.debug("annulus.area = %g; annulus_area = %g",
                      annulus.area, annulus_area)
        else:
            log.debug("annulus.area() = %g; annulus_area = %g",
                      annulus.area(), annulus_area)
        if annulus_area > 0.:
            normalization = aperture_area / annulus_area
        else:
            log.warning("Background annulus has no area, so background "
                        "subtraction will be turned off.")
            subtract_background = False
    del temp

    npixels[:] = aperture_area
    for k in range(shape[0]):
        phot_table = aperture_photometry(data[k, :, :], aperture,
                                         method=method, subpixels=subpixels)
        temp_flux[k] = float(phot_table['aperture_sum'][0])
        if subtract_background:
            bkg_table = aperture_photometry(data[k, :, :], annulus,
                                            method=method, subpixels=subpixels)
            background[k] = float(bkg_table['aperture_sum'][0])
            temp_flux[k] = temp_flux[k] - background[k] * normalization

    # Check for NaNs in the wavelength array, flag them in the dq array,
    # and truncate the arrays if NaNs are found at endpoints (unless the
    # entire array is NaN).
    (wavelength, temp_flux, background, npixels, dq) = \
        nans_in_wavelength(wavelength, temp_flux, background, npixels, dq)

    return (ra, dec, wavelength, temp_flux, background, npixels, dq)


def locn_from_wcs(input_model, ra_targ, dec_targ):
    """Get the location of the spectrum, based on the WCS.

    Parameters
    ----------
    input_model : data model
        The input science model.

    ra_targ, dec_targ : float or None
        The right ascension and declination of the target (degrees).

    Returns
    -------
    locn : tuple of two int, or None
        X and Y coordinates (in that order) of the pixel that has right
        ascension and declination coordinates closest to the target
        location.  The spectral extraction region should be centered here.
    """

    if ra_targ is None or dec_targ is None:
        log.warning("TARG_RA and/or TARG_DEC not found; can't compute "
                    "pixel location of target.")
        locn = None
    else:
        shape = input_model.data.shape
        grid = np.indices(shape[-2:])
        z = np.zeros(shape[-2:], dtype=np.float64) + shape[0] // 2
        # The arguments are the X, Y, and Z pixel coordinates, and the
        # output arrays will be 2-D.
        (ra_i, dec_i, wl) = input_model.meta.wcs(grid[1], grid[0], z)
        cart = celestial_to_cartesian(ra_i, dec_i)
        v = celestial_to_cartesian(ra_targ, dec_targ)   # a single vector
        diff = cart - v
        # We want the pixel with the minimum distance from v, but the pixel
        # with the minimum value of distance squared will be the same.
        dist2 = (diff**2).sum(axis=-1)
        nan_mask = np.isnan(wl)
        dist2[..., :] = np.where(nan_mask, HUGE_DIST, dist2[..., :])
        del nan_mask
        k = np.argmin(dist2.ravel())
        (j, i) = divmod(k, dist2.shape[1])      # y, x coordinates

        if i <= 0 or j <= 0 or i >= shape[-1] - 1 or j >= shape[-2] - 1:
            log.warning("WCS implies the target is at or beyond the edge "
                        "of the image; this location will not be used.")
            locn = None
        else:
            locn = (i, j)                       # x, y coordinates

    return locn


def celestial_to_cartesian(ra, dec):
    """Convert celestial coordinates to Cartesian.

    Parameters
    ----------
    ra, dec : ndarray or float
        The right ascension and declination (degrees).  Both `ra` and `dec`
        should be arrays of the same shape, or they should both be float.

    Returns
    -------
    cart : ndarray
        If `ra` and `dec` are float, `cart` with be a 3-element array.
        If `ra` and `dec` are arrays, `cart` will be an array with shape
        ra.shape + (3,).
        For each element of `ra` (or `dec`), the last axis of `cart` will
        give the Cartesian coordinates of a unit vector in the direction
        `ra`, `dec`.  The elements of the vector in Cartesian coordinates
        are in the order x, y, z, where x is the direction toward the
        vernal equinox, y is the direction toward right ascension = 90
        degrees (6 hours) and declination = 0, and z is toward the north
        celestial pole.
    """

    if hasattr(ra, 'shape'):
        shape = ra.shape + (3,)
    else:
        shape = (3,)

    cart = np.zeros(shape, dtype=np.float64)
    cart[..., 2] = np.sin(dec * np.pi / 180.)
    r_xy = np.cos(dec * np.pi / 180.)
    cart[..., 1] = r_xy * np.sin(ra * np.pi / 180.)
    cart[..., 0] = r_xy * np.cos(ra * np.pi / 180.)

    return cart


def image_extract_ifu(input_model, source_type, extract_params):
    """Extraction using a reference image.

    Extended summary
    ----------------
    One of the requirements for this step is that for an extended target,
    the entire aperture is supposed to be extracted (with no background
    subtraction).  It doesn't make any sense to use an image reference file
    to extract the entire aperture; a trivially simple JSON reference file
    would do.  Therefore, we assume that if the user specified a reference
    file in image format, the user actually wanted that reference file
    to be used, so we will ignore the requirement in this case.

    Parameters
    ----------
    input_model : IFUCubeModel
        The input model.

    source_type : string
        "point" or "extended"

    extract_params : dict
        The extraction parameters.  One of these is an open file handle
        for an image that specifies which pixels should be included as
        source and which (if any) should be used as background.

    Returns
    -------
    ra, dec : float
        ra and dec are the right ascension and declination respectively
        at the centroid of the target region in the reference image.

    wavelength : ndarray, 1-D
        The wavelength in micrometers at each plane of the IFU cube.

    temp_flux : ndarray, 1-D
        The sum of the data values in the extraction aperture minus the
        sum of the data values in the background region (scaled by the
        ratio of areas), for each plane.
        The data values are in units of surface brightness, so this value
        isn't really the flux, it's an intermediate value.  Dividing by
        `npixels` (to compute the average) will give the value for the
        `surf_bright` (surface brightness) column, and multiplying by
        the solid angle of a pixel will give the flux for a point source.

    background : ndarray, 1-D
        The background count rate that was subtracted from the total
        source data values to get `temp_flux`.

    npixels : ndarray, 1-D, float64
        For each slice, this is the number of pixels that were added
        together to get `temp_flux`.

    dq : ndarray, 1-D, uint32
        The data quality array.
    """

    data = input_model.data
    # The axes are (z, y, x) in the sense that (x, y) are the ordinary
    # axes of the image for one plane, i.e. at one wavelength.  The
    # wavelengths vary along the z-axis.
    shape = data.shape

    # The dispersion direction is the first axis.
    temp_flux = np.zeros(shape[0], dtype=np.float64)
    background = np.zeros(shape[0], dtype=np.float64)
    npixels = np.ones(shape[0], dtype=np.float64)
    n_bkg = np.ones(shape[0], dtype=np.float64)

    dq = np.zeros(shape[0], dtype=np.uint32)

    ref_image = extract_params['ref_image']
    ref = ref_image.data
    if 'subtract_background' in extract_params:
        subtract_background = extract_params['subtract_background']
    else:
        subtract_background = None

    (mask_target, mask_bkg) = separate_target_and_background(ref)

    if subtract_background is not None:
        if subtract_background:
            if mask_bkg is None:
                log.info("Skipping background subtraction because "
                         "background regions are not defined.")
        else:
            if mask_bkg is not None:
                log.info("Background subtraction was turned off "
                         "- skipping it.")
            mask_bkg = None

    ra_targ = input_model.meta.target.ra
    dec_targ = input_model.meta.target.dec
    locn = locn_from_wcs(input_model, ra_targ, dec_targ)
    if locn is not None:
        log.info("Target location is x_center = %g, y_center = %g, "
                 "based on TARG_RA and TARG_DEC.", locn[0], locn[1])

    # Use the centroid of mask_target as the point where the target
    # would be without any nod/dither correction.
    (y0, x0) = im_centroid(data, mask_target)
    log.debug("Target location based on reference image is X = %g, Y = %g",
              x0, y0)

    if locn is None or np.isnan(locn[0]):
        log.warning("Couldn't determine pixel location from WCS, so "
                    "nod/dither correction will not be applied.")
    else:
        (x_center, y_center) = locn
        # Shift the reference image so it will be centered at locn.
        # Only shift by a whole number of pixels.
        delta_x = int(round(x_center - x0))         # must be integer
        delta_y = int(round(y_center - y0))
        log.debug("Shifting reference image by %g in X and %g in Y",
                  delta_x, delta_y)
        temp = shift_ref_image(mask_target, delta_y, delta_x)
        if temp is not None:
            mask_target = temp
            del temp
            if mask_bkg is not None:
                mask_bkg = shift_ref_image(mask_bkg, delta_y, delta_x)
            # Since we have shifted mask_target and mask_bkg to
            # x_center and y_center, update x0 and y0.
            x0 = x_center
            y0 = y_center

    # Extract the data.
    # First add up the values along the x direction, then add up the
    # values along the y direction.
    gross = (data * mask_target).sum(axis=2, dtype=np.float64).sum(axis=1)

    # Compute the number of pixels that were added together to get gross.
    normalization = 1.
    temp = np.ones_like(data)
    npixels[:] = (temp * mask_target).sum(axis=2, dtype=np.float64).sum(axis=1)
    if mask_bkg is not None:
        n_bkg[:] = (temp * mask_bkg).sum(axis=2, dtype=np.float64).sum(axis=1)
        n_bkg = np.where(n_bkg <= 0., 1., n_bkg)
        normalization = npixels / n_bkg
    del temp

    # Extract the background.
    if mask_bkg is not None:
        background = (data * mask_bkg).sum(axis=2, dtype=np.float64).sum(axis=1)
        background *= normalization
        temp_flux = gross - background
    else:
        background = np.zeros_like(gross)
        temp_flux = gross.copy()

    # Compute the ra, dec, and wavelength at the pixels of a column through
    # the IFU cube at the target location.  ra and dec should be constant
    # (so they're scalars), but wavelength will vary from plane to plane.

    log.debug("IFU 1-D extraction parameters (using reference image):")
    log.debug("  x_center = %s", str(x0))
    log.debug("  y_center = %s", str(y0))
    log.debug("  subtract_background parameter = %s", str(subtract_background))
    if mask_bkg is not None:
        log.debug("    background will be subtracted")
    else:
        log.debug("    background will not be subtracted")

    (ra, dec, wavelength) = get_coordinates(input_model, x0, y0)

    # Check for NaNs in the wavelength array, flag them in the dq array,
    # and truncate the arrays if NaNs are found at endpoints (unless the
    # entire array is NaN).
    (wavelength, temp_flux, background, npixels, dq) = \
        nans_in_wavelength(wavelength, temp_flux, background, npixels, dq)

    return (ra, dec, wavelength, temp_flux, background, npixels, dq)


def get_coordinates(input_model, x0, y0):
    """Get celestial coordinates and wavelengths.

    Parameters
    ----------
    input_model : IFUCubeModel
        The input model.

    x0, y0 : float
        The pixel number at which the coordinates should be determined.
        If the reference file is JSON format, this point will be the
        nominal center of the image.  For an image reference file this
        will be the centroid of the pixels that were flagged as source.

    Returns
    -------
    ra, dec : float
        ra and dec are the right ascension and declination respectively
        at pixel (0, y0, x0).

    wavelength : ndarray, 1-D
        The wavelength in micrometers at each pixel.
    """

    if hasattr(input_model.meta, 'wcs'):
        wcs = input_model.meta.wcs
    else:
        log.warning("WCS function not found in input.")
        wcs = None

    nelem = input_model.data.shape[0]

    if wcs is not None:
        x_array = np.empty(nelem, dtype=np.float64)
        x_array.fill(x0)
        y_array = np.empty(nelem, dtype=np.float64)
        y_array.fill(y0)
        z_array = np.arange(nelem, dtype=np.float64) # for wavelengths
        ra, dec, wavelength = wcs(x_array, y_array, z_array)
        # ra and dec should be constant, so nelem // 2 is an arbitrary element.
        ra = ra[nelem // 2]
        dec = dec[nelem // 2]
    else:
        (ra, dec) = (0., 0.)
        wavelength = np.arange(1, nelem + 1, dtype=np.float64)

    return (ra, dec, wavelength)


def nans_in_wavelength(wavelength, flux, background, npixels, dq):
    """Check for NaNs in the wavelength array.

    If NaNs are found in the wavelength array, flag them in the dq array,
    and truncate the arrays at either or both ends if NaNs are found at
    endpoints (unless the entire array is NaN).

    Parameters
    ----------
    wavelength : ndarray, 1-D, float64
        The wavelength in micrometers at each pixel.

    flux : ndarray, 1-D, float64
        The flux minus the background at each pixel.

    background : ndarray, 1-D, float64
        The background count rate that was subtracted from the total
        source count rate to get `flux`.

    npixels : ndarray, 1-D, float64
        For each slice, this is the number of pixels that were added
        together to get `flux`.

    dq : ndarray, 1-D, uint32
        The data quality array.

    Returns
    -------
    wavelength : ndarray, 1-D, float64

    flux : ndarray, 1-D, float64

    background : ndarray, 1-D, float64

    npixels : ndarray, 1-D, float64

    dq : ndarray, 1-D, uint32
    """

    nelem = np.size(wavelength)
    if nelem == 0:
        log.warning("Output arrays are empty!")
        return (wavelength, flux, background, npixels, dq)

    nan_mask = np.isnan(wavelength)
    n_nan = nan_mask.sum(dtype=np.intp)
    if n_nan == nelem:
        log.warning("Wavelength array is all NaN!")
        dq = np.bitwise_or(dq[:], dqflags.pixel['DO_NOT_USE'])
        return (wavelength, flux, background, npixels, dq)

    if n_nan > 0:
        log.warning("%d NaNs in wavelength array.", n_nan)
        dq[nan_mask] = np.bitwise_or(dq[nan_mask], dqflags.pixel['DO_NOT_USE'])
        not_nan = np.logical_not(nan_mask)
        flag = np.where(not_nan)
        if len(flag[0]) > 0:
            n_trimmed = flag[0][0] + nelem - (flag[0][-1] + 1)
            if n_trimmed > 0:
                slc = slice(flag[0][0], flag[0][-1] + 1)
                wavelength = wavelength[slc]
                flux = flux[slc]
                background = background[slc]
                npixels = npixels[slc]
                dq = dq[slc]
                log.info("Output arrays have been trimmed by %d elements",
                         n_trimmed)

    return (wavelength, flux, background, npixels, dq)


def separate_target_and_background(ref):
    """Create masks for target and background.

    Parameters
    ----------
    ref : ndarray, 2-D or 3-D
        This is the reference image data array.  This should be the same
        shape as one plane of the science data (or the same shape as the
        entire 3-D science data array).  A value of 1 in a pixel indicates
        that the pixel should be included when computing the target
        spectrum.  A value of 0 means the pixel is not part of either the
        target or background.  A value of -1 means the pixel should be
        included as part of the background region.

    Returns
    -------
    mask_target : ndarray, 2-D or 3-D
        This is an array of the same type and shape as the reference
        image, but with values of only 0 or 1.  A value of 1 indicates
        that the corresponding pixel of the science data array should be
        included when adding up values to make the 1-D spectrum, and a
        value of 0 means that it should not be included.

    mask_bkg : ndarray, 2-D or 3-D, or None.
        This is like `mask_target` (i.e. values are 0 or 1) but for
        background regions.  A value of 1 in `mask_bkg` indicates a pixel
        that should be included as part of the background.  If there is no
        pixel in the reference image with a value of -1, `mask_bkg` will
        be set to None.
    """

    mask_target = np.where(ref == 1., 1., 0.)

    if np.any(ref == -1.):
        mask_bkg = np.where(ref == -1., 1., 0.)
    else:
        mask_bkg = None

    return (mask_target, mask_bkg)


def im_centroid(data, mask_target):
    """Compute the mean location of the target.

    Parameters
    ----------
    data : ndarray, 3-D
        This is the science image data array.

    mask_target : ndarray, 2-D or 3-D
        This is an array of the same type and shape as one plane of the
        science image (or the same type and shape of the entire 3-D science
        image), but with values of 0 or 1, where 1 indicates a pixel within
        the source.

    Returns
    -------
    y0, x0 : tuple of two float
        The centroid of pixels flagged as source.
    """

    # Collapse the science data along the dispersion direction to get a
    # 2-D image of the IFU field of view.  Multiplying by mask_target
    # zeros out all pixels that are not regarded as part of the target
    # (or targets).
    if len(mask_target.shape) == 2:
        data_2d = data.sum(axis=0, dtype=np.float64) * mask_target
    else:
        data_2d = (data * mask_target).sum(axis=0, dtype=np.float64)
    if data_2d.sum() == 0.:
        log.warning("Couldn't compute image centroid.")
        shape = data_2d.shape
        y0 = shape[0] / 2.
        x0 = shape[0] / 2.
        return (y0, x0)

    x_profile = data_2d.sum(axis=0, dtype=np.float64)
    x = np.arange(data_2d.shape[1], dtype=np.float64)
    s_x = (x_profile * x).sum()
    x0 = s_x / x_profile.sum()

    y_profile = data_2d.sum(axis=1, dtype=np.float64)
    y = np.arange(data_2d.shape[0], dtype=np.float64)
    s_y = (y_profile * y).sum()
    y0 = s_y / y_profile.sum()

    return (y0, x0)


def shift_ref_image(mask, delta_y, delta_x, fill=0):
    """Apply nod/dither offset to target or background for ref image.

    Parameters
    ----------
    mask : ndarray, 2-D or 3-D
        This is either the target mask or the background mask, which was
        created from an image reference file.
        It is assumed that all pixels in the science data are good.  If
        that is not correct, `shift_ref_image` may be called with a data
        quality array as the first argument, in order to obtain a shifted
        data quality array.  In this case, `fill` should be set to a
        positive value, e.g. 1.

    delta_y, delta_x : int
        These are the shifts to be applied to the vertical and horizontal
        axes respectively.

    fill : int or float
        The output array will be initialized to this value.  This should
        be 0 (the default) if `mask` is a target or background mask, but
        it should be set to 1 (or some other positive value) if `mask` is
        a data quality array.

    Returns
    -------
    temp : ndarray, same type and shape as `mask`
        A copy of `mask`, but shifted by `delta_y` and `delta_x`.
    """

    if delta_x == 0 and delta_y == 0:
        return mask.copy()

    shape = mask.shape
    if abs(delta_y) >= shape[-2] or abs(delta_x) >= shape[-1]:
        log.warning("Nod offset %d or %d is too large, skipping ...",
                    delta_y, delta_x)
        return None

    if delta_y > 0:
        islice_y = slice(0, -delta_y)
        oslice_y = slice(delta_y, shape[-2])
    elif delta_y < 0:
        islice_y = slice(-delta_y, shape[-2])
        oslice_y = slice(0, delta_y)
    else:
        islice_y = slice(0, shape[-2])
        oslice_y = slice(0, shape[-2])

    if delta_x > 0:
        islice_x = slice(0, -delta_x)
        oslice_x = slice(delta_x, shape[-1])
    elif delta_x < 0:
        islice_x = slice(-delta_x, shape[-1])
        oslice_x = slice(0, delta_x)
    else:
        islice_x = slice(0, shape[-1])
        oslice_x = slice(0, shape[-1])

    temp = np.zeros_like(mask) + fill
    temp[..., oslice_y, oslice_x] = mask[..., islice_y, islice_x]

    return temp
