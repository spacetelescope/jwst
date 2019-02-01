import logging
from collections import namedtuple
import copy
import json
import math

import numpy as np
from astropy.modeling import polynomial
from .. import datamodels
from ..datamodels import dqflags
from .. assign_wcs import niriss        # for specifying spectral order number
from .. transforms import models as trmodels
from .. lib import pipe_utils
from . import extract1d
from . import ifu
from . import spec_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

WFSS_EXPTYPES = ['NIS_WFSS', 'NRC_WFSS', 'NRC_GRISM', 'NRC_TSGRISM']
"""Exposure types to be regarded as wide-field slitless spectroscopy."""

# These values are used to indicate whether the input reference file
# (if any) is JSON or IMAGE.
FILE_TYPE_JSON = "JSON"
FILE_TYPE_IMAGE = "IMAGE"
FILE_TYPE_OTHER = "N/A"

# This is to prevent calling offset_from_offset multiple times for
# multi-integration data.
OFFSET_NOT_ASSIGNED_YET = "not assigned yet"

ANY = "ANY"
"""Wildcard for slit name.

Extended summary
----------------
For full-frame input data, keyword SLTNAME may not be populated, so the
slit name will be set to this string to indicate that the first slit in
the reference file should be used.
A slit name in the reference file can also be ANY, in which case the
reference information for that slit will be regarded as matching any slit
name from the input data.
"""

ANY_ORDER = 1000
"""Wildcard for spectral order number in a reference image.

Extended summary
----------------
If the reference file contains images, keyword SPORDER gives the order
number of the spectrum that would be extracted using a given image in
the reference file.
"""

HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""

# These values are assigned in get_extract_parameters, using key "match".
# If there was an aperture in the reference file for which the "id" key
# matched, that's (at least) a partial match.  If "spectral_order" also
# matched, that's an exact match.

NO_MATCH = "no match"
PARTIAL = "partial match"
EXACT = "exact match"

DUMMY = "dummy"
"""Distinguish between a slit object and not a slit object.

Extended summary
----------------
The input model can contain a sequence of slits (SlitModel objects,
2-D cutouts, with attributes), or it might be a full detector image or
subarray.  In the former case, the step will try to extract a spectrum
from each slit.  In the latter case, the input can be a CubeModel or an
ImageModel, and not all of the attributes are the same.  DUMMY is
passed to functions instead of a SlitModel object if the input model
does not contain slits.
"""

Aperture = namedtuple('Aperture', ['xstart', 'ystart', 'xstop', 'ystop'])


class Extract1dError(Exception):
    pass

class InvalidSpectralOrderNumberError(Extract1dError):
    """The spectral order number was invalid or off the detector."""
    def __init__(self, message=None):
        super().__init__()
        self.message = message


def load_ref_file(refname):
    """Open the reference file.

    Parameters
    ----------
    refname : str
        The name of the reference file.  This file is expected to be
        either a JSON file giving extraction information, or a file
        containing one or more images that are to be used as masks that
        define the extraction region and optionally background regions.

    Returns
    -------
    ref_dict : dict
        If the reference file is in JSON format, ref_dict will be the
        dictionary returned by json.load(), except that the file type
        ('JSON') will also be included with key 'ref_file_type'.
        If the reference file is an image, ref_dict will be a
        dictionary with two keys:  ref_dict['ref_file_type'] = 'IMAGE'
        and ref_dict['ref_model'].  The latter will be the open file
        handle for the jwst.datamodels object for the reference file.
    """

    if refname == "N/A":
        ref_dict = None
    else:
        # Try reading the file as JSON.
        fd = open(refname)
        try:
            ref_dict = json.load(fd)
            ref_dict['ref_file_type'] = FILE_TYPE_JSON
            fd.close()
        except UnicodeDecodeError:
            fd.close()
            # Try opening the file as a reference image.
            try:
                fd = datamodels.MultiExtract1dImageModel(refname)
                ref_dict = {'ref_file_type': FILE_TYPE_IMAGE}
                ref_dict['ref_model'] = fd      # only used for images
            except OSError:
                log.info("The reference file should be JSON or FITS.")
                log.error("Don't know how to read %s.", refname)
                raise

    return ref_dict


def get_extract_parameters(ref_dict,
                           input_model, slitname, sp_order,
                           meta, smoothing_length, bkg_order):
    """Get reference file values.

    Parameters
    ----------
    ref_dict : dict or None
        For a reference file in JSON format, `ref_dict` will be the entire
        contents of the file.  For a reference image, `ref_dict` will have
        just two entries, 'ref_file_type' (a string) and 'ref_model', a
        JWST data model for a collection of images.  If there is no
        reference file, `ref_dict` will be None.

    input_model : data model
        This can be either the input science file or one SlitModel out of
        a list of slits.

    slitname : str
        The name of the slit, or "ANY", or for NIRISS SOSS data this will
        be the subarray name.

    sp_order : int
        The spectral order number.

    meta : metadata for current slit.
        Not currently used.

    smoothing_length : int or None
        Width of a boxcar function for smoothing the background regions.
        If None, the smoothing length will be gotten from `ref_dict`, or
        it will be set to 0 (no background smoothing) if this key is
        not found in `ref_dict`.
        If `smoothing_length` is not None, that means that the user
        explicitly specified the value, so that value will be used.
        This argument is only used if background regions have been
        specified.

    bkg_order : int or None
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background.  If None, the polynomial
        order will be gotten from `ref_dict`, or it will be set to 0 if
        not found in `ref_dict`.
        A value of 0 means that a simple average of the background
        regions, column by column (or row by row), will be used.
        If `bkg_order` is not None, that means that the user explicitly
        specified the value, so that value will be used.
        This argument must be positive or zero, and it is only used if
        background regions have been specified.

    Returns
    -------
    extract_params : dict
        Information copied out of `ref_dict`.  The items will be selected
        based on `slitname` and `sp_order`.  Default values will be
        assigned if `ref_dict` is None.  For a reference image, the key
        'ref_image' gives the image model.
    """

    extract_params = {'match': NO_MATCH}        # initial value

    if ref_dict is None:
        # There is no reference file; use "reasonable" default values.
        extract_params['ref_file_type'] = FILE_TYPE_OTHER
        extract_params['match'] = EXACT
        shape = input_model.data.shape
        extract_params['spectral_order'] = sp_order
        extract_params['xstart'] = 0                    # first pixel in X
        extract_params['xstop'] = shape[-1] - 1         # last pixel in X
        extract_params['ystart'] = 0                    # first pixel in Y
        extract_params['ystop'] = shape[-2] - 1         # last pixel in Y
        extract_params['extract_width'] = None
        extract_params['src_coeff'] = None
        extract_params['bkg_coeff'] = None
        extract_params['nod_correction'] = 0
        extract_params['independent_var'] = 'pixel'
        extract_params['smoothing_length'] = 0  # because no background sub.
        extract_params['bkg_order'] = 0         # because no background sub.
        # Note that extract_params['dispaxis'] is not assigned.  This will
        # be done later by calling find_dispaxis().

    elif ref_dict['ref_file_type'] == FILE_TYPE_JSON:
        extract_params['ref_file_type'] = ref_dict['ref_file_type']
        for aper in ref_dict['apertures']:
            if 'id' in aper and aper['id'] != "dummy" and \
               (aper['id'] == slitname or aper['id'] == "ANY" or
                slitname == "ANY"):
                extract_params['match'] = PARTIAL
                # region_type is retained for backward compatibility; it is
                # not required to be present.
                # spectral_order is a secondary selection criterion.  The
                # default is the expected value, so if the key is not present
                # in the JSON file, the current aperture will be selected.
                # If the current aperture in the JSON file has
                # "spectral_order": "ANY", that aperture will be selected.
                region_type = aper.get("region_type", "target")
                if region_type != "target":
                    continue
                spectral_order = aper.get("spectral_order", sp_order)
                if spectral_order == sp_order or spectral_order == ANY:
                    extract_params['match'] = EXACT
                    extract_params['spectral_order'] = sp_order
                    disp = aper.get('dispaxis')
                    if disp is None:
                        # Will be found later by calling find_dispaxis().
                        log.warning("dispaxis not specified in reference file")
                    elif disp != HORIZONTAL and disp != VERTICAL:
                        log.error("dispaxis = %d is not valid.", disp)
                        raise ValueError('dispaxis must be 1 or 2.')
                    else:
                        extract_params['dispaxis'] = disp
                    extract_params['src_coeff'] = aper.get('src_coeff')
                    extract_params['bkg_coeff'] = aper.get('bkg_coeff')
                    extract_params['independent_var'] = \
                          aper.get('independent_var', 'pixel').lower()
                    if smoothing_length is None:
                        extract_params['smoothing_length'] = \
                              aper.get('smoothing_length', 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params['smoothing_length'] = smoothing_length
                    if bkg_order is None:
                        extract_params['bkg_order'] = aper.get('bkg_order', 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params['bkg_order'] = bkg_order
                    extract_params['xstart'] = aper.get('xstart')
                    extract_params['xstop'] = aper.get('xstop')
                    extract_params['ystart'] = aper.get('ystart')
                    extract_params['ystop'] = aper.get('ystop')
                    extract_params['extract_width'] = aper.get('extract_width')
                    extract_params['nod_correction'] = 0        # default value
                    break

    elif ref_dict['ref_file_type'] == FILE_TYPE_IMAGE:
        extract_params['ref_file_type'] = ref_dict['ref_file_type']
        foundit = False
        for im in ref_dict['ref_model'].images:
            if im.name == slitname or im.name == ANY or slitname == ANY:
                extract_params['match'] = PARTIAL
                if (im.spectral_order == sp_order or
                    im.spectral_order >= ANY_ORDER):
                    extract_params['match'] = EXACT
                    extract_params['spectral_order'] = sp_order
                    foundit = True
                    break

        if foundit:
            extract_params['ref_image'] = im
            if hasattr(im, "dispersion_axis"):
                if (im.dispersion_axis == HORIZONTAL or
                    im.dispersion_axis == VERTICAL):
                    extract_params['dispaxis'] = im.dispersion_axis
            # else dispaxis will be set by find_dispaxis()
            if smoothing_length is None:
                extract_params['smoothing_length'] = im.smoothing_length
            else:
                # The user-supplied value takes precedence.
                extract_params['smoothing_length'] = smoothing_length
            extract_params['nod_correction'] = 0

    else:
        log.error("Reference file type %s not recognized",
                  ref_dict['ref_file_type'])

    return extract_params


def find_dispaxis(input_model, slit, spectral_order, extract_params):
    """Get the location of the spectrum, based on the WCS.

    Parameters:
    -----------
    input_model : data model
        The input science file.

    slit : data model or str
        This is a slit from a `MultiSlitModel` (or similar), or "dummy"
        if the input is not an array of 2-D cutouts.
        We use the `meta.wcs` and `wavelength` attributes.

    spectral_order : int
        The spectral order number.

    extract_params : dict, may be modified in-place
        Parameters read from the 1-D extraction reference file.
        The dispersion direction will be determined by comparing the
        increment in wavelength from one pixel in `slit` to the next,
        in the horizontal and vertical directions.
        If key 'dispaxis' is present in `extract_params` and is not None,
        that value will not be modified; however, if 'dispaxis' was not in
        `extract_params` or if its value was None, then
        `extract_params['dispaxis']` will be updated with the value
        determined by this function.
    """

    if 'dispaxis' in extract_params:
        initial_value = extract_params['dispaxis']
    else:
        initial_value = None
        # There needs to be a default value for dispaxis.  If this can't be
        # updated with a valid value, we will not extract the spectrum.
        extract_params['dispaxis'] = None

    if slit == DUMMY:
        shape = input_model.data.shape[-2:]
    else:
        shape = slit.data.shape[-2:]

    wcs = None                                  # initial value
    if slit == DUMMY:
        if input_model.meta.exposure.type == "NIS_SOSS":
            if hasattr(input_model.meta, 'wcs'):
                try:
                    transform = niriss.niriss_soss_set_input(
                                        input_model, spectral_order)
                except ValueError:
                    if spectral_order == 1:
                        log.warning("Spectral order 1 not found")
                        log.warning("Can't determine dispaxis from the WCS.")
                        return
                    else:
                        log.warning("Spectral order %d not found, using 1",
                                    spectral_order)
                        transform = niriss.niriss_soss_set_input(
                                        input_model, 1)
                wcs = transform                 # not None
        elif hasattr(input_model.meta, 'wcs'):
            wcs = input_model.meta.wcs
            transform = wcs.forward_transform
    elif hasattr(slit, 'meta') and hasattr(slit.meta, 'wcs'):
        wcs = slit.meta.wcs
        transform = wcs.forward_transform

    if wcs is None:
        log.warning("find_dispaxis:  WCS not found")
        return

    n_inputs = None
    if hasattr(transform, 'n_inputs'):
        n_inputs = transform.n_inputs
    elif (hasattr(transform, 'forward_transform') and
          hasattr(transform.forward_transform, 'n_inputs')):
            n_inputs = transform.forward_transform.n_inputs
    elif (hasattr(wcs, 'forward_transform') and
          hasattr(wcs.forward_transform, 'n_inputs')):
            n_inputs = wcs.forward_transform.n_inputs
    else:
        log.warning("Can't find n_inputs, assuming n_inputs = 2")
        n_inputs = 2
    if n_inputs is None or n_inputs < 2 or n_inputs > 3:
        log.warning("n_inputs for wcs is %s; should be 2 or 3", str(n_inputs))
        log.warning("Can't determine dispaxis from the WCS.")
        return

    if slit == DUMMY:
        got_wavelength = False
    else:
        if hasattr(slit, "wavelength"):
            wl_array = slit.wavelength.copy()
        else:
            wl_array = None
        if (wl_array is None or len(wl_array) == 0 or
            wl_array.min() == 0. and wl_array.max() == 0.):
                got_wavelength = False
        else:
            got_wavelength = True

    bb = wcs.bounding_box
    if bb is None:
        x_cent = shape[-1] // 2
        y_cent = shape[-2] // 2
    else:
        x_cent = float(bb[0][0] + bb[0][1]) / 2.
        y_cent = float(bb[1][0] + bb[1][1]) / 2.
        x_cent = math.floor(x_cent)
        y_cent = math.floor(y_cent)

    # Find which is the dispersion axis, by comparing the increment in
    # wavelength from one pixel to the next, in the horizontal and
    # vertical directions.
    dispaxis = None                             # pessimistic value
    if got_wavelength:
        mask = (wl_array == 0)
        if np.any(mask):
            wl_array[mask] = np.nan
        del mask
        dwlx = wl_array[:, 1:] - wl_array[:, 0:-1]
        dwly = wl_array[1:, :] - wl_array[0:-1, :]
        dwlx = np.nanmean(dwlx)
        dwly = np.nanmean(dwly)
        log.debug("find_dispaxis, wavelength attribute:  dwlx = %s dwly = %s",
                  str(dwlx), str(dwly))

        dwlx = np.abs(dwlx)
        dwly = np.abs(dwly)
        if dwlx > dwly:
            dispaxis = HORIZONTAL
        elif dwlx < dwly:
            dispaxis = VERTICAL
    else:
        if n_inputs > 2:
            stuff = transform(x_cent, y_cent, spectral_order)
            wl_00 = stuff[2]
            stuff = transform(x_cent, y_cent + 1, spectral_order)
            wl_01 = stuff[2]
            stuff = transform(x_cent + 1, y_cent, spectral_order)
            wl_10 = stuff[2]
        else:
            stuff = transform(x_cent, y_cent)
            wl_00 = stuff[2]
            stuff = transform(x_cent, y_cent + 1)
            wl_01 = stuff[2]
            stuff = transform(x_cent + 1, y_cent)
            wl_10 = stuff[2]
        dwlx = wl_10 - wl_00
        dwly = wl_01 - wl_00
        if (np.isnan(dwlx) or np.isnan(dwly) or dwlx == dwly or
            wl_00 == 0 or wl_01 == 0 or wl_10 == 0):
                if dwlx == dwly:
                    log.warning("One-pixel offset gives dwlx = {}, dwly = {}"
                                .format(dwlx, dwly))
                else:
                    log.warning("wavelength from WCS is NaN or 0 "
                                "within the bounding box")
                log.warning("    computing differences over "
                            "the whole bounding box ...")
                if input_model.meta.exposure.type in WFSS_EXPTYPES:
                    wl_wcs = np.zeros(shape, dtype=np.float64)
                    log.debug("Starting to compute wavelengths ...")
                    if n_inputs > 2:
                        for j in range(shape[0]):
                            for i in range(shape[1]):
                                stuff = transform(i, j, spectral_order)
                                if stuff[2] == 0.:
                                    wl_wcs[j, i] = np.nan
                                else:
                                    wl_wcs[j, i] = stuff[2]
                    else:
                        for j in range(shape[0]):
                            for i in range(shape[1]):
                                stuff = transform(i, j)
                                if stuff[2] == 0.:
                                    wl_wcs[j, i] = np.nan
                                else:
                                    wl_wcs[j, i] = stuff[2]
                    log.debug("... finished computing wavelengths")
                else:
                    grid = np.indices(shape, dtype=np.float64)
                    if n_inputs > 2:
                        stuff = transform(grid[1], grid[0], spectral_order)
                    else:
                        stuff = transform(grid[1], grid[0])
                    wl_wcs = stuff[2].copy()
                    del grid, stuff
                    # Flag wavelength = 0 as invalid.
                    mask = (wl_wcs == 0)
                    if np.any(mask):
                        wl_wcs[mask] = np.nan
                    del mask
                dwlx = wl_wcs[:, 1:] - wl_wcs[:, 0:-1]
                dwly = wl_wcs[1:, :] - wl_wcs[0:-1, :]
                dwlx = np.nanmean(dwlx)
                dwly = np.nanmean(dwly)
        log.debug("find_dispaxis, using wcs:  dwlx = %s dwly = %s",
                  str(dwlx), str(dwly))
        if np.isnan(dwlx) or np.isnan(dwly):
            log.warning("dwlx and/or dwly is STILL NaN!  "
                        "Can't determine dispaxis from WCS")
        else:
            dwlx = np.abs(dwlx)
            dwly = np.abs(dwly)
            if dwlx > dwly:
                dispaxis = HORIZONTAL
            elif dwlx < dwly:
                dispaxis = VERTICAL

    if dispaxis is None:
        log.warning("Can't determine dispaxis from the WCS.")
    elif initial_value is None:
        # Only assign this value if dispaxis wasn't present in the
        # reference file, or if there was no reference file.
        extract_params['dispaxis'] = dispaxis

    if initial_value is None or initial_value == dispaxis:
        log_fcn = log.debug
    else:
        log_fcn = log.warning
    log_fcn("find_dispaxis:  dispaxis from ref file = %s, "
            "from wavelengths = %s", str(initial_value), str(dispaxis))


def log_initial_parameters(extract_params):
    """Log some of the initial extraction parameters.

    Parameters
    ----------
    extract_params : dict
        Information read from the reference file.
    """

    if not "xstart" in extract_params:
        return

    log.debug("Initial parameters:")
    log.debug("dispaxis = %s", str(extract_params["dispaxis"]))
    log.debug("spectral order = %d", extract_params["spectral_order"])
    log.debug("independent_var = %s", extract_params["independent_var"])
    log.debug("smoothing_length = %d", extract_params["smoothing_length"])
    log.debug("initial xstart = %s", str(extract_params["xstart"]))
    log.debug("initial xstop = %s", str(extract_params["xstop"]))
    log.debug("initial ystart = %s", str(extract_params["ystart"]))
    log.debug("initial ystop = %s", str(extract_params["ystop"]))
    log.debug("extract_width = %s", str(extract_params["extract_width"]))
    log.debug("initial src_coeff = %s", str(extract_params["src_coeff"]))
    log.debug("initial bkg_coeff = %s", str(extract_params["bkg_coeff"]))
    log.debug("bkg_order = %d", extract_params["bkg_order"])


def get_aperture(im_shape, wcs, verbose, extract_params):
    """Get the extraction limits xstart, xstop, ystart, ystop.

    Parameters
    ----------
    im_shape : tuple
        The shape (2-D) of the input data.  This will be for the current
        integration, if the input contains more than one integration.

    wcs : a WCS object, or None
        The wcs (if any) for the input data or slit.

    verbose : bool
        If True, log messages.

    extract_params : dict
        Parameters read from the reference file.

    Returns
    -------
    ap_ref : namedtuple or an empty dict
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    if extract_params['ref_file_type'] == FILE_TYPE_IMAGE:
        return {}

    ap_ref = aperture_from_ref(extract_params, im_shape)

    (ap_ref, truncated) = update_from_shape(ap_ref, im_shape)
    if truncated and verbose:
        log.warning("Extraction limits extended outside the input image "
                    "borders; limits have been truncated.")

    if wcs is not None:
        ap_wcs = aperture_from_wcs(wcs, verbose)
    else:
        ap_wcs = None

    # If the xstart, etc., values were not specified for the dispersion
    # direction, the extraction region should be centered within the
    # WCS bounding box (domain).
    ap_ref = update_from_wcs(ap_ref, ap_wcs, extract_params["extract_width"],
                             extract_params["dispaxis"], verbose)
    ap_ref = update_from_width(ap_ref, extract_params["extract_width"],
                               extract_params["dispaxis"])

    return ap_ref


def aperture_from_ref(extract_params, im_shape):
    """Get extraction region from reference file or image shape.

    Parameters
    ----------
    extract_params : dict
        Parameters read from the reference file.

    im_shape : tuple of int
        The last two elements are the height and width of the input image
        (slit).

    Returns
    -------
    ap_ref : namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    nx = im_shape[-1]
    ny = im_shape[-2]

    xstart = extract_params.get('xstart', None)
    xstop = extract_params.get('xstop', None)
    ystart = extract_params.get('ystart', None)
    ystop = extract_params.get('ystop', None)

    if xstart is None:
        xstart = 0
    if xstop is None:
        xstop = nx - 1                  # limits are inclusive
    if ystart is None:
        ystart = 0
    if ystop is None:
        ystop = ny - 1

    ap_ref = Aperture(xstart=xstart, xstop=xstop, ystart=ystart, ystop=ystop)

    return ap_ref


def update_from_width(ap_ref, extract_width, direction):
    """Update XD extraction limits based on extract_width.

    If extract_width was specified, that value should override
    ystop - ystart (or xstop - xstart, depending on dispersion direction).

    Parameters
    ----------
    ap_ref : namedtuple
        Contains xstart, xstop, ystart, ystop.  These are the initial
        values as read from the reference file, except that they may
        have been truncated at the image borders.

    extract_width : int or None
        The number of pixels in the cross-dispersion direction to add
        together to make a 1-D spectrum from a 2-D image.

    direction : int
        HORIZONTAL (1) if the dispersion direction is predominantly
        horizontal.  VERTICAL (2) if the dispersion direction is
        predominantly vertical.

    Returns
    -------
    ap_width : namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    if extract_width is None:
        return ap_ref

    if direction == HORIZONTAL:
        temp_width = ap_ref.ystop - ap_ref.ystart + 1
    else:
        temp_width = ap_ref.xstop - ap_ref.xstart + 1

    if extract_width == temp_width:
        return ap_ref                                   # OK as is

    # An integral value corresponds to the center of a pixel.  If the
    # extraction limits were not specified via polynomial coefficients,
    # assign_polynomial_limits will create polynomial functions using
    # values from an Aperture, and these lower and upper limits will be
    # expanded by 0.5 to give polynomials (constant functions) for the
    # lower and upper edges of the bounding pixels.
    width = float(extract_width)
    if direction == HORIZONTAL:
        lower = float(ap_ref.ystart)
        upper = float(ap_ref.ystop)
        lower = (lower + upper) / 2. - (width - 1.) / 2.
        upper = lower + (width - 1.)
        ap_width = Aperture(xstart=ap_ref.xstart, xstop=ap_ref.xstop,
                            ystart=lower, ystop=upper)
    else:
        lower = float(ap_ref.xstart)
        upper = float(ap_ref.xstop)
        lower = (lower + upper) / 2. - (width - 1.) / 2.
        upper = lower + (width - 1.)
        ap_width = Aperture(xstart=lower, xstop=upper,
                            ystart=ap_ref.ystart, ystop=ap_ref.ystop)

    return ap_width


def update_from_shape(ap, im_shape):
    """Truncate extraction region based on input image shape.

    Parameters
    ----------
    ap : namedtuple
        Extraction region.

    im_shape : tuple of int
        The last two elements are the height and width of the input image.

    Returns
    -------
    ap_shape : namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.

    truncated : bool
         True if any value was truncated at an image edge.
    """

    nx = im_shape[-1]
    ny = im_shape[-2]

    xstart = ap.xstart
    xstop = ap.xstop
    ystart = ap.ystart
    ystop = ap.ystop

    truncated = False
    if ap.xstart < 0:
        xstart = 0
        truncated = True
    if ap.xstop >= nx:
        xstop = nx - 1                          # limits are inclusive
        truncated = True
    if ap.ystart < 0:
        ystart = 0
        truncated = True
    if ap.ystop >= ny:
        ystop = ny - 1
        truncated = True

    ap_shape = Aperture(xstart=xstart, xstop=xstop,
                        ystart=ystart, ystop=ystop)

    return (ap_shape, truncated)


def aperture_from_wcs(wcs, verbose):
    """Get the limits over which the WCS is defined.

    Parameters
    ----------
    wcs : data model
        The world coordinate system interface.

    verbose : bool
        If True, log messages.

    Returns
    -------
    ap_wcs : namedtuple or None
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.  These are the
        limits copied directly from wcs.bounding_box.
    """

    got_bounding_box = False
    try:
        bounding_box = wcs.bounding_box
        got_bounding_box = True
    except AttributeError:
        if verbose:
            log.info("wcs.bounding_box not found; using wcs.domain instead.")
        bounding_box = ((wcs.domain[0]['lower'], wcs.domain[0]['upper']),
                        (wcs.domain[1]['lower'], wcs.domain[1]['upper']))

    if got_bounding_box and bounding_box is None:
        if verbose:
            log.warning("wcs.bounding_box is None")
        return None

    # bounding_box should be a tuple of tuples, each of the latter
    # consisting of (lower, upper) limits.
    if len(bounding_box) < 2:
        if verbose:
            log.warning("wcs.bounding_box has the wrong shape")
        return None

    # These limits are float, and they are inclusive.
    xstart = bounding_box[0][0]
    xstop = bounding_box[0][1]
    ystart = bounding_box[1][0]
    ystop = bounding_box[1][1]

    ap_wcs = Aperture(xstart=xstart, xstop=xstop, ystart=ystart, ystop=ystop)

    return ap_wcs


def update_from_wcs(ap_ref, ap_wcs, extract_width, direction, verbose):
    """Limit the extraction region to the WCS bounding box.

    Parameters
    ----------
    ap_ref : namedtuple
        Contains xstart, xstop, ystart, ystop.  These are the values of
        the extraction region as specified by the reference file or the
        image size.

    ap_wcs : namedtuple or None
        These are the bounding box limits.

    extract_width : int
        The number of pixels in the cross-dispersion direction to add
        together to make a 1-D spectrum from a 2-D image.

    direction : int
        HORIZONTAL (1) if the dispersion direction is predominantly
        horizontal.  VERTICAL (2) if the dispersion direction is
        predominantly vertical.

    verbose : bool
        If True, log messages.

    Returns
    -------
    ap : namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    if ap_wcs is None:
        return ap_ref

    # If the wcs limits don't pass the sanity test, ignore the bounding box.
    if not sanity_check_limits(ap_ref, ap_wcs, verbose):
        return ap_ref

    # ap_wcs has the limits over which the WCS transformation is defined;
    # take those as the outer limits over which we will extract.
    xstart = compare_start(ap_ref.xstart, ap_wcs.xstart)
    ystart = compare_start(ap_ref.ystart, ap_wcs.ystart)
    xstop = compare_stop(ap_ref.xstop, ap_wcs.xstop)
    ystop = compare_stop(ap_ref.ystop, ap_wcs.ystop)

    if extract_width is not None:
        if direction == HORIZONTAL:
            width = ystop - ystart + 1
        else:
            width = xstop - xstart + 1
        if width < extract_width:
            if verbose:
                log.warning("extract_width was truncated from %g to %g",
                            extract_width, width)

    ap = Aperture(xstart=xstart, xstop=xstop, ystart=ystart, ystop=ystop)

    return ap


def sanity_check_limits(ap_ref, ap_wcs, verbose):
    """Sanity check.

    Parameters
    ----------
    ap_ref : namedtuple
        Contains xstart, xstop, ystart, ystop.  These are the values of
        the extraction region as specified by the reference file or the
        image size.

    ap_wcs : namedtuple
        These are the bounding box limits.

    verbose : bool
        If True, log messages.

    Returns
    -------
    flag : boolean
        True if ap_ref and ap_wcs do overlap, i.e. if the sanity test passes.
    """

    if (ap_wcs.xstart >= ap_ref.xstop or ap_wcs.xstop <= ap_ref.xstart or
        ap_wcs.ystart >= ap_ref.ystop or ap_wcs.ystop <= ap_ref.ystart):
        if verbose:
            log.warning("The WCS bounding box is outside the aperture:")
            log.warning("  aperture:  %s, %s, %s, %s",
                        str(ap_ref.xstart), str(ap_ref.xstop),
                        str(ap_ref.ystart), str(ap_ref.ystop))
            log.warning("  wcs:       %s, %s, %s, %s",
                        str(ap_wcs.xstart), str(ap_wcs.xstop),
                        str(ap_wcs.ystart), str(ap_wcs.ystop))
            log.warning(" so the wcs bounding box will be ignored.")
        flag = False
    else:
        flag = True

    return flag


def compare_start(start_ref, start_wcs):
    """Compare the start limit from the aperture with the WCS lower limit.

    Extended summary
    ----------------
    The more restrictive (i.e. larger) limit is the one upon which the
    output value will be based.  If the WCS limit is larger, the value will
    be increased to an integer, on the assumption that WCS lower limits
    correspond to the lower edge of the bounding pixel.  If this value
    will actually be used for an extraction limit (i.e. if the limits were
    not already specified by polynomial coefficients), then
    `assign_polynomial_limits` will create a polynomial function using this
    value, except that it will be decreased by 0.5 to correspond to the
    lower edge of the bounding pixels.

    Parameters
    ----------
    start_ref : int or float
        xstart or ystart, as specified by the reference file or the image
        size.

    start_wcs : int or float
        The lower limit from the WCS bounding box.

    Returns
    -------
    value : int or float
        The start limit, possibly constrained by the WCS start limit.
    """

    if start_ref >= start_wcs:          # ref is inside WCS limit
        value = start_ref
    else:                               # outside (below) WCS limit
        value = math.ceil(start_wcs)

    return value


def compare_stop(stop_ref, stop_wcs):
    """Compare the stop limit from the aperture with the WCS upper limit.

    The more restrictive (i.e. smaller) limit is the one upon which the
    output value will be based.  If the WCS limit is smaller, the value
    will be truncated to an integer, on the assumption that WCS upper
    limits correspond to the upper edge of the bounding pixel.  If this
    value will actually be used for an extraction limit (i.e. if the
    limits were not already specified by polynomial coefficients), then
    `assign_polynomial_limits` will create a polynomial function using this
    value, except that it will be increased by 0.5 to correspond to the
    upper edge of the bounding pixels.

    Parameters
    ----------
    stop_ref : int or float
        xstop or ystop, as specified by the reference file or the image
        size.

    stop_wcs : int or float
        The upper limit from the WCS bounding box.

    Returns
    -------
    value : int or float
        The stop limit, possibly constrained by the WCS stop limit.
    """

    if stop_ref <= stop_wcs:            # ref is inside WCS limit
        value = stop_ref
    else:                               # outside (above) WCS limit
        value = math.floor(stop_wcs)

    return value


def create_poly(coeff):
    """Create a polynomial model from coefficients.

    Parameters
    ----------
    coeff : list of float
        The coefficients of the polynomial, constant term first, highest
        order term last.

    Returns
    -------
    `astropy.modeling.polynomial.Polynomial1D` object, or None if `coeff`
        is empty.
    """

    n = len(coeff)
    if n < 1:
        return None

    coeff_dict = {}
    for i in range(n):
        key = "c{}".format(i)
        coeff_dict[key] = coeff[i]

    return polynomial.Polynomial1D(degree=n - 1, **coeff_dict)


class ExtractBase:
    """Base class for 1-D extraction info and methods.

    Attributes
    ----------
    exp_type : str
        Exposure type.

    dispaxis : int
        Dispersion direction:  1 is horizontal, 2 is vertical.

    spectral_order : int
        Spectral order number.

    xstart : int or None
        First pixel (zero indexed) in extraction region.

    xstop : int or None
        Last pixel (zero indexed) in extraction region.

    ystart : int or None
        First pixel (zero indexed) in extraction region.

    ystop : int or None
        Last pixel (zero indexed) in extraction region.

    extract_width : int or None
        Height (in the cross-dispersion direction) of the extraction
        region.

    independent_var : str
        The polynomial functions for computing the upper and lower
        boundaries of extraction and background regions can be functions
        of pixel number or wavelength (in microns).  These options are
        distinguished by `independent_var`.

    src_coeff : list of lists of float, or None
        These are coefficients of polynomial functions that define the
        cross-dispersion limits of one or more source extraction regions
        (yes, there can be more than one extraction region).  Note that
        these are float values, and the limits can include fractions of
        pixels.
        If specified, this takes priority over `ystart`, `ystop`, and
        `extract_width`, though `xstart` and `xstop` will still be used
        for the limits in the dispersion direction.  (Interchange "x" and
        "y" if the dispersion is in the vertical direction.)

        For example, the value could be:

-           [[1, 2], [3, 4, 5], [6], [7, 8]]

        which means:

-           [1, 2] coefficients for the lower limit of the first region
-           [3, 4, 5] coefficients for the upper limit of the first region
-           [6] coefficient(s) for the lower limit of the second region
-           [7, 8] coefficients for the upper limit of the second region

        The coefficients are listed with the constant term first, highest
        order term last.  For example, [3, 4, 5] means 3 + 4 * x + 5 * x^2.

    bkg_coeff : list of lists of float, or None
        This has the same format as `src_coeff`, but the polynomials
        define one or more background regions.

    p_src : list of astropy.modeling.polynomial.Polynomial1D
        These are Astropy polynomial functions defining the limits in the
        cross-dispersion direction of the source extraction region(s).
        If `src_coeff` was specified, `p_src` will be created directly
        from `src_coeff`; otherwise, a constant function based on
        `ystart`, `ystop`, `extract_width` will be assigned to `p_src`.

    p_bkg : list of astropy.modeling.polynomial.Polynomial1D, or None
        These are Astropy polynomial functions defining the limits in the
        cross-dispersion direction of the background extraction regions.
        This list will be populated from `bkg_coeff`, if that was specified.

    smoothing_length : int
        Width of a boxcar function for smoothing the background regions.
        This argument must be an odd positive number or zero, and it is
        only used if background regions have been specified.

    bkg_order : int
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background.
        This argument must be positive or zero, and it is only used if
        background regions have been specified.

    nod_correction : float
        If not zero, this will be added to the extraction region limits
        for the cross-dispersion direction, both target and background.

    wcs : WCS object
        For computing the right ascension, declination, and wavelength at
        one or more pixels.
    """

    def __init__(self):
        self.exp_type = ""
        """
        Issue #1781.
        self.instrument_name = ""
        """
        self.dispaxis = None
        self.spectral_order = None
        self.xstart = None
        self.xstop = None
        self.ystart = None
        self.ystop = None
        self.extract_width = None
        self.independent_var = "pixel"
        self.src_coeff = None
        self.bkg_coeff = None
        self.p_src = None
        self.p_bkg = None
        self.smoothing_length = 0
        self.bkg_order = 0
        self.nod_correction = 0.
        self.wcs = None


    def update_extraction_limits(self, ap):
        pass


    def assign_polynomial_limits(self, verbose):
        pass


    def offset_from_offset(self, input_model, slit, verbose):
        """Get nod/dither pixel offset from [xy]_offset.

        Parameters
        ----------
        input_model : data model
            The input science data.

        slit : SlitModel, or "dummy"
            One slit from a MultiSlitModel (or similar), or "dummy" if
            there are no slits.

        verbose : boolean
            If True, write log messages.

        Returns
        -------
        offset : float
            The offset of the exposure from the nominal position, due to
            nod or dither.  This is the component of the offset
            perpendicular to the dispersion direction.  A positive value
            means that the spectrum is at a larger pixel number than the
            nominal location.
        """

        instrument = input_model.meta.instrument.name

        total_points = input_model.meta.dither.total_points
        if total_points is None or total_points < 2:
            if verbose:
                log.info("Total number of dither points = %s; assuming no "
                         "nod/dither offset", str(total_points))
            return 0.

        missing = False
        if 'detector' not in self.wcs.available_frames:
            if verbose:
                log.warning("detector frame is not available,")
            missing = True
        if 'v2v3' not in self.wcs.available_frames:
            if verbose:
                if missing:
                    log.warning("and v2v3 frame is not available,")
                else:
                    log.warning("v2v3 frame is not available,")
            missing = True
        if missing:
            if verbose:
                log.warning("so can't compute nod/dither offset")
            return 0.
        v2v3_detector = self.wcs.get_transform('v2v3', 'detector')

        xoffset = input_model.meta.dither.x_offset      # in arcsec
        yoffset = input_model.meta.dither.y_offset      # in arcsec
        if verbose:
            log.debug("xoffset = %s, yoffset = %s", str(xoffset), str(yoffset))
        if xoffset is None or yoffset is None:
            if verbose:
                log.warning("XOFFSET and/or YOFFSET not found; "
                            "assuming no nod/dither offset")
            return 0.

        v2ref = None
        v3ref = None
        v3idlyangle = None
        vparity = None
        wl_start = None
        wl_end = None
        if slit != DUMMY:
            if hasattr(slit.meta.wcsinfo, "v2_ref"):
                v2ref = slit.meta.wcsinfo.v2_ref                # in arcsec
            if hasattr(slit.meta.wcsinfo, "v3_ref"):
                v3ref = slit.meta.wcsinfo.v3_ref                # in arcsec
            if hasattr(slit.meta.wcsinfo, "v3yangle"):
                v3idlyangle = slit.meta.wcsinfo.v3yangle        # in degrees
            if hasattr(slit.meta.wcsinfo, "vparity"):
                vparity = slit.meta.wcsinfo.vparity
            # These wavelengths are in meters; convert to microns later.
            if hasattr(slit.meta.wcsinfo, "waverange_start"):
                wl_start = slit.meta.wcsinfo.waverange_start
            if hasattr(slit.meta.wcsinfo, "waverange_end"):
                wl_end = slit.meta.wcsinfo.waverange_end
        else:
            if hasattr(input_model.meta.wcsinfo, "v2_ref"):
                v2ref = input_model.meta.wcsinfo.v2_ref         # in arcsec
            if hasattr(input_model.meta.wcsinfo, "v3_ref"):
                v3ref = input_model.meta.wcsinfo.v3_ref         # in arcsec
            if hasattr(input_model.meta.wcsinfo, "v3yangle"):
                v3idlyangle = input_model.meta.wcsinfo.v3yangle # in degrees
            if hasattr(input_model.meta.wcsinfo, "vparity"):
                vparity = input_model.meta.wcsinfo.vparity
            if hasattr(input_model.meta.wcsinfo, "waverange_start"):
                wl_start = input_model.meta.wcsinfo.waverange_start
            if hasattr(input_model.meta.wcsinfo, "waverange_end"):
                wl_end = input_model.meta.wcsinfo.waverange_end

        if v3idlyangle is None and instrument == "NIRSPEC":
            v3idlyangle = 138.892975
            if verbose:
                log.warning("v3yangle not specified; "
                            "setting it to %.6f degrees", v3idlyangle)
        if vparity is None:
            vparity = -1.
            if verbose:
                log.warning("vparity not specified; setting it to %g",
                            vparity)

        if verbose:
            log.debug("v2ref = %s, v3ref = %s, v3yangle = %s, "
                      "vparity = %s, wl_start = %s, wl_end = %s",
                      str(v2ref), str(v3ref), str(v3idlyangle), str(vparity),
                        str(wl_start), str(wl_end))
        if (v2ref is None or v3ref is None or
            v3idlyangle is None or vparity is None or
            wl_start is None or wl_end is None):
                if verbose:
                    log.warning("Missing wcsinfo values; "
                                "can't compute nod/dither offset")
                return 0.

        idl_v23 = trmodels.IdealToV2V3(v3idlyangle, v2ref, v3ref, vparity)

        # Wavelength at the middle of the range, in microns.
        wavelength = (0.5 * (wl_end - wl_start) + wl_start) * 1.e6

        # Compute the location in V2,V3 [in arcsec]
        xv0, yv0 = idl_v23(0., 0.)
        xv, yv = idl_v23(xoffset, yoffset)
        (x0, y0) = v2v3_detector(xv0, yv0, wavelength)
        (x, y) = v2v3_detector(xv, yv, wavelength)
        if verbose:
            log.debug("x0 = %s, x = %s, y0 = %s, y = %s",
                      str(x0), str(x), str(y0), str(y))
        if x0 is None or y0 is None or x is None or y is None:
            if verbose:
                log.warning("One or more of x, y, x0, y0 is None; "
                            "can't compute nod/dither offset")
            return 0.
        # offsets from xoffset = 0, yoffset = 0
        dx = x - x0
        dy = y - y0
        if self.dispaxis == HORIZONTAL:
            offset = dy
        else:
            offset = dx
        if np.isnan(offset):
            if verbose:
                log.warning("Nod/dither offset is NaN; setting it to 0.")
            offset = 0.

        return offset


class ExtractModel(ExtractBase):
    """The extraction region was specified in a JSON file."""

    def __init__(self, input_model, slit,
                 ref_file_type=None,
                 match="unknown",
                 dispaxis=HORIZONTAL, spectral_order=1,
                 xstart=None, xstop=None, ystart=None, ystop=None,
                 extract_width=None, src_coeff=None, bkg_coeff=None,
                 independent_var="pixel",
                 smoothing_length=0, bkg_order=0, nod_correction=0.,
                 x_center=None, y_center=None,
                 inner_bkg=None, outer_bkg=None, method='subpixel'):
        """Create a polynomial model from coefficients.

        Extended summary
        ----------------
        If InvalidSpectralOrderNumberError is raised, processing of the
        current slit or spectral order should be skipped.

        Parameters
        ----------
        input_model : data model
            The input science data.

        slit : an input slit, or a dummy value if not used
            For MultiSlit or MultiProduct data, `slit` is one slit from
            a list of slits in the input.  For other types of data, `slit`
            will not be used.

        ref_file_type : str
            This indicates whether the reference file (if any) was a JSON
            file or an image.

        match : str
            An entry in the reference file (if there is one) should match
            the slit name and the spectral order number of the current
            slit (or it can be "ANY").  `match` will be "exact match" if
            both the name of the current slit and the spectral order number
            match the selected entry in the reference file, and it will be
            "partial match" if only the slit name matches.  If neither
            match, `match` will be "no match".

        dispaxis : int
            Dispersion direction:  1 is horizontal, 2 is vertical.

        spectral_order : int
            Spectral order number.

        xstart : int
            First pixel (zero indexed) in extraction region.

        xstop : int
            Last pixel (zero indexed) in extraction region.

        ystart : int
            First pixel (zero indexed) in extraction region.

        ystop : int
            Last pixel (zero indexed) in extraction region.

        extract_width : int
            Height (in the cross-dispersion direction) of the extraction
            region.

        src_coeff : list of lists of float, or None
            These are coefficients of polynomial functions that define
            the cross-dispersion limits of one or more source extraction
            regions.

        bkg_coeff : list of lists of float, or None
            This has the same format as `src_coeff`, but the polynomials
            define one or more background regions.

        independent_var : str
            This can be either "pixel" or "wavelength" to specify the
            independent variable for polynomial functions.

        smoothing_length : int
            Width of a boxcar function for smoothing the background
            regions.

        bkg_order : int
            Polynomial order for fitting to each column (or row, if the
            dispersion is vertical) of background.

        nod_correction : float
            If not zero, this will be added to the extraction region limits
            for the cross-dispersion direction, both target and background.

        x_center : float
            This is not relevant; it's only used for IFU data.

        y_center : float
            This is not relevant; it's only used for IFU data.

        inner_bkg : float
            This is not relevant; it's only used for IFU data.

        outer_bkg : float
            This is not relevant; it's only used for IFU data.

        method : str
            This is not relevant; it's only used for IFU data.
        """

        super().__init__()

        self.exp_type = input_model.meta.exposure.type

        self.dispaxis = dispaxis
        self.spectral_order = spectral_order

        # xstart, xstop, ystart, or ystop may be overridden with src_coeff,
        # they may be limited by the input image size or by the WCS bounding
        # box, or they may be modified if extract_width was specified
        # (because extract_width takes precedence).
        # If these values are specified, the limits in the cross-dispersion
        # direction should be integers, but they may later be replaced with
        # fractional values, depending on extract_width, in order to center
        # the extraction window in the originally specified xstart to xstop
        # (or ystart to ystop).
        if xstart is None:
            self.xstart = None
        elif self.dispaxis == VERTICAL:
            r = round(xstart)
            if xstart == r:
                self.xstart = xstart
            else:
                log.warning("xstart %s should have been an integer; "
                            "rounding to %s", str(xstart), str(r))
                self.xstart = r
        else:                           # dispaxis is HORIZONTAL
            self.xstart = xstart

        if xstop is None:
            self.xstop = None
        elif self.dispaxis == VERTICAL:
            r = round(xstop)
            if xstop == r:
                self.xstop = xstop
            else:
                log.warning("xstop %s should have been an integer; "
                            "rounding to %s", str(xstop), str(r))
                self.xstop = r
        else:
            self.xstop = xstop

        if ystart is None:
            self.ystart = None
        elif self.dispaxis == HORIZONTAL:
            r = round(ystart)
            if ystart == r:
                self.ystart = ystart
            else:
                log.warning("ystart %s should have been an integer; "
                            "rounding to %s", str(ystart), str(r))
                self.ystart = r
        else:
            self.ystart = ystart

        if ystop is None:
            self.ystop = None
        elif self.dispaxis == HORIZONTAL:
            r = round(ystop)
            if ystop == r:
                self.ystop = ystop
            else:
                log.warning("ystop %s should have been an integer; "
                            "rounding to %s", str(ystop), str(r))
                self.ystop = r
        else:
            self.ystop = ystop

        if extract_width is None:
            self.extract_width = None
        else:
            self.extract_width = int(round(extract_width))
        # 'wavelength' or 'pixel', the independent variable for functions
        # for lower and upper limits of source and background regions.
        self.independent_var = independent_var.lower()
        if (self.independent_var != "wavelength" and
            self.independent_var != "pixel" and
            self.independent_var != "pixels"):
            log.error("independent_var = '%s'; "
                      "specify 'wavelength' or 'pixel'", self.independent_var)
            raise RuntimeError("Invalid value for independent_var")

        # Coefficients for source (i.e. target) and background limits and
        # corresponding polynomial functions.
        self.src_coeff = copy.deepcopy(src_coeff)
        self.bkg_coeff = copy.deepcopy(bkg_coeff)

        # These functions will be assigned by assign_polynomial_limits.
        # The "p" in the attribute name indicates a polynomial function.
        self.p_src = None
        self.p_bkg = None

        if smoothing_length is None:
            smoothing_length = 0
        if smoothing_length > 0 and \
           smoothing_length // 2 * 2 == smoothing_length:
            log.warning("smoothing_length was even (%d), so incremented by 1",
                        smoothing_length)
            smoothing_length += 1               # must be odd
        self.smoothing_length = smoothing_length
        self.bkg_order = bkg_order
        self.nod_correction = nod_correction

        self.wcs = None                         # initial value
        if input_model.meta.exposure.type == "NIS_SOSS":
            if hasattr(input_model.meta, 'wcs'):
                try:
                    self.wcs = niriss.niriss_soss_set_input(
                                input_model, self.spectral_order)
                except ValueError:
                    raise InvalidSpectralOrderNumberError(
                                "Spectral order {} is not valid"
                                .format(self.spectral_order))
        elif slit == DUMMY:
            if hasattr(input_model.meta, 'wcs'):
                self.wcs = input_model.meta.wcs
        elif hasattr(slit, 'meta') and hasattr(slit.meta, 'wcs'):
            self.wcs = slit.meta.wcs
        if self.wcs is None:
            log.warning("WCS function not found in input.")


    def add_nod_correction(self, verbose):
        """Add the nod offset to the extraction location (in-place).

        Extended summary
        ----------------
        If source extraction coefficients src_coeff were specified, this
        method will add the nod offset correction to the first coefficient
        of every coefficient list; otherwise, the nod offset will be added
        to xstart & xstop or to ystart & ystop.
        If background extraction coefficients bkg_coeff were specified,
        this method will add the nod offset to the first coefficients.
        Note that background coefficients are handled independently of
        src_coeff.

        Parameters
        ----------
        verbose : bool
            If True, messages can be logged.
        """

        if self.nod_correction == 0.:
            return

        if self.src_coeff is None:
            if self.dispaxis == HORIZONTAL:
                dir = "y"
                self.ystart += self.nod_correction
                self.ystop += self.nod_correction
            else:
                dir = "x"
                self.xstart += self.nod_correction
                self.xstop += self.nod_correction
            if verbose:
                log.info("Applying nod/dither offset of %s "
                         "to %sstart and %sstop",
                         str(self.nod_correction), dir, dir)

        if self.src_coeff is not None or self.bkg_coeff is not None:
            if verbose:
                log.info("Applying nod/dither offset of %s "
                         "to polynomial coefficients",
                         str(self.nod_correction))

        if self.src_coeff is not None:
            n_lists = len(self.src_coeff)
            for i in range(n_lists):
                coeff_list = self.src_coeff[i]
                coeff_list[0] += self.nod_correction
                self.src_coeff[i] = copy.copy(coeff_list)

        if self.bkg_coeff is not None:
            n_lists = len(self.bkg_coeff)
            for i in range(n_lists):
                coeff_list = self.bkg_coeff[i]
                coeff_list[0] += self.nod_correction
                self.bkg_coeff[i] = copy.copy(coeff_list)


    def update_extraction_limits(self, ap):
        """Update start and stop limits.

        Extended summary
        ----------------
        Copy the values of xstart, etc., to the attributes.  Note, however,
        that if src_coeff was specified, that will override the values
        given by xstart, etc.
        The limits in the dispersion direction will be rounded to integer.

        Parameters
        ----------
        ap : namedtuple
        """

        self.xstart = ap.xstart
        self.xstop = ap.xstop
        self.ystart = ap.ystart
        self.ystop = ap.ystop

        if self.dispaxis == HORIZONTAL:
            self.xstart = int(round(self.xstart))
            self.xstop = int(round(self.xstop))
        else:                           # vertical
            self.ystart = int(round(self.ystart))
            self.ystop = int(round(self.ystop))


    def log_extraction_parameters(self):
        """Log the updated extraction parameters."""

        log.debug("Updated parameters:")
        log.debug("nod_correction = %s", str(self.nod_correction))
        note_x = ""
        note_y = ""
        if self.src_coeff is not None:
            # Since src_coeff was specified, that will be used instead of
            # xstart & xstop (or ystart & ystop).
            if self.dispaxis == HORIZONTAL:
                note_y = " (not actually used)"
            else:
                note_x = " (not actually used)"
        log.debug("xstart = %s%s", str(self.xstart), note_x)
        log.debug("xstop = %s%s", str(self.xstop), note_x)
        log.debug("ystart = %s%s", str(self.ystart), note_y)
        log.debug("ystop = %s%s", str(self.ystop), note_y)
        if self.src_coeff is not None:
            log.debug("src_coeff = %s", str(self.src_coeff))
        if self.bkg_coeff is not None:
            log.debug("bkg_coeff = %s", str(self.bkg_coeff))


    def assign_polynomial_limits(self, verbose):
        """Create polynomial functions for extraction limits.

        Extended summary
        ----------------
        self.src_coeff and self.bkg_coeff contain lists of polynomial
        coefficients.  These will be used to create corresponding lists of
        polynomial functions, self.p_src and self.p_bkg.  Note, however,
        that the structures of those two lists are not the same.
        The coefficients lists have this form:

-           [[1, 2], [3, 4, 5], [6], [7, 8]]

        which means:

-           [1, 2] coefficients for the lower limit of the first region
-           [3, 4, 5] coefficients for the upper limit of the first region
-           [6] coefficient(s) for the lower limit of the second region
-           [7, 8] coefficients for the upper limit of the second region

        The lists of coefficients must always be in pairs, for the lower
        and upper limits respectively, but they're not explicitly in an
        additional layer of two-element lists,

-           i.e., not like this:  [[[1, 2], [3, 4, 5]], [[6], [7, 8]]]

        That seemed unnecessarily messy and harder for the user to specify.
        For the lists of polynomial functions, on the other hand, that
        additional layer of list is used:

-           [[fcn_lower1, fcn_upper1], [fcn_lower2, fcn_upper2]]

        where:

-           fcn_lower1 is 1 + 2 * x
-           fcn_upper1 is 3 + 4 * x + 5 * x**2
-           fcn_lower2 is 6
-           fcn_upper2 is 7 + 8 * x

        Parameters
        ----------
        verbose : bool
            If True, messages can be logged.
        """

        if self.src_coeff is None:
            # Create constant functions.

            if self.dispaxis == HORIZONTAL:
                lower = float(self.ystart) - 0.5
                upper = float(self.ystop) + 0.5
            else:
                lower = float(self.xstart) - 0.5
                upper = float(self.xstop) + 0.5
            if verbose:
                log.debug("Converting extraction limits to [[%g], [%g]]",
                          lower, upper)
            self.p_src = [[create_poly([lower]), create_poly([upper])]]
        else:
            # The source extraction can include more than one region.
            n_lists = len(self.src_coeff)
            if n_lists // 2 * 2 != n_lists:
                raise RuntimeError("src_coeff must contain alternating lists "
                                   "of lower and upper limits.")
            self.p_src = []
            expect_lower = True                         # toggled in loop
            for coeff_list in self.src_coeff:
                if expect_lower:
                    lower = create_poly(coeff_list)
                else:
                    upper = create_poly(coeff_list)
                    self.p_src.append([lower, upper])
                expect_lower = not expect_lower

        if self.bkg_coeff is not None:
            n_lists = len(self.bkg_coeff)
            if n_lists // 2 * 2 != n_lists:
                raise RuntimeError("bkg_coeff must contain alternating lists "
                                   "of lower and upper limits.")
            self.p_bkg = []
            expect_lower = True                         # toggled in loop
            for coeff_list in self.bkg_coeff:
                if expect_lower:
                    lower = create_poly(coeff_list)
                else:
                    upper = create_poly(coeff_list)
                    self.p_bkg.append([lower, upper])
                expect_lower = not expect_lower


    def extract(self, data, wl_array, verbose):
        """Do the extraction.

        Extended summary
        ----------------
        This version is for the case that the reference file is a JSON
        file, or that there is no reference file.

        Parameters
        ----------
        data : ndarray, 2-D
            Data array from which the spectrum will be extracted.

        wl_array : ndarray, 2-D, or None
            Wavelengths corresponding to `data`, or None if no WAVELENGTH
            extension was found in the input file.

        verbose : bool
            If True, log messages.

        Returns
        -------
        ra, dec : float
            ra and dec are the right ascension and declination respectively
            at the nominal center of the slit.

        wavelength : ndarray, 1-D
            The wavelength in micrometers at each pixel.

        net : ndarray, 1-D
            The count rate (counts / s) minus the background at each pixel.

        background : ndarray, 1-D
            The background count rate that was subtracted from the total
            source count rate to get `net`.

        dq : ndarray, 1-D, int32
            The data quality array.
        """

        # If the wavelength attribute exists and is populated, use it
        # in preference to the wavelengths returned by the wcs function.
        if wl_array is None or len(wl_array) == 0:
            got_wavelength = False
        else:
            got_wavelength = True               # may be reset later

        # The default value is 0, so all 0 values means that the
        # wavelength attribute was not populated.
        if not got_wavelength or wl_array.min() == 0. and wl_array.max() == 0.:
            got_wavelength = False
        if got_wavelength:
            if verbose:
                log.debug("Wavelengths are from wavelength attribute.")
            # We need a 1-D array of wavelengths, one element for each
            # output table row.
            # These are slice limits.
            sx0 = int(round(self.xstart))
            sx1 = int(round(self.xstop)) + 1
            sy0 = int(round(self.ystart))
            sy1 = int(round(self.ystop)) + 1
            # Convert non-positive values to NaN, to easily ignore them.
            wl = wl_array.copy()                # don't modify wl_array
            nan_flag = np.isnan(wl)
            # To avoid a warning about invalid value encountered in less_equal.
            wl[nan_flag] = -1000.
            wl = np.where(wl <= 0., np.nan, wl)
            if self.dispaxis == HORIZONTAL:
                wavelength = np.nanmean(wl[sy0:sy1, sx0:sx1], axis=0)
            else:
                wavelength = np.nanmean(wl[sy0:sy1, sx0:sx1], axis=1)

        # Now call the wcs function to compute the celestial coordinates.
        # Also use the returned wavelengths if we weren't able to get them
        # from the wavelength attribute.

        # Used for computing the celestial coordinates.
        if self.dispaxis == HORIZONTAL:
            slice0 = int(round(self.xstart))
            slice1 = int(round(self.xstop)) + 1
            x_array = np.arange(slice0, slice1, dtype=np.float64)
            y_array = np.empty(x_array.shape, dtype=np.float64)
            y_array.fill((self.ystart + self.ystop) / 2.)
        else:
            slice0 = int(round(self.ystart))
            slice1 = int(round(self.ystop)) + 1
            y_array = np.arange(slice0, slice1, dtype=np.float64)
            x_array = np.empty(y_array.shape, dtype=np.float64)
            x_array.fill((self.xstart + self.xstop) / 2.)

        if self.wcs is not None:
            if verbose and not got_wavelength:
                log.debug("Wavelengths are from the wcs function.")
            nelem = slice1 - slice0
            if self.exp_type in WFSS_EXPTYPES:
                # We expect two (x and y) or three (x, y, spectral order).
                n_inputs = self.wcs.forward_transform.n_inputs
                ra = np.zeros(nelem, dtype=np.float64)
                dec = np.zeros(nelem, dtype=np.float64)
                # Temporary variable so as not to clobber `wavelength`.
                wcs_wl = np.zeros(nelem, dtype=np.float64)
                transform = self.wcs.forward_transform
                if n_inputs == 2:
                    for i in range(nelem):
                        ra[i], dec[i], wcs_wl[i], _ = transform(
                                        x_array[i], y_array[i])
                elif n_inputs == 3:
                    for i in range(nelem):
                        ra[i], dec[i], wcs_wl[i], _ = transform(
                                x_array[i], y_array[i], self.spectral_order)
                else:
                    if verbose:
                        log.warning("n_inputs for wcs function is %d",
                                    n_inputs)
                        log.warning("WCS function was expected to take "
                                    "either 2 or 3 arguments.")
                    ra[:] = -999.
                    dec[:] = -999.
                    wcs_wl[:] = -999.
            else:
                ra, dec, wcs_wl = self.wcs(x_array, y_array)
            # We need one right ascension and one declination, representing
            # the direction of pointing.
            mask = np.isnan(ra)
            not_nan = np.logical_not(mask)
            if np.any(not_nan):
                ra2 = ra[not_nan]
                min_ra = ra2.min()
                max_ra = ra2.max()
                ra = (min_ra + max_ra) / 2.
            else:
                if verbose:
                    log.warning("All right ascension values are NaN; "
                                "assigning dummy value -999.")
                ra = -999.
            mask = np.isnan(dec)
            not_nan = np.logical_not(mask)
            if np.any(not_nan):
                dec2 = dec[not_nan]
                min_dec = dec2.min()
                max_dec = dec2.max()
                dec = (min_dec + max_dec) / 2.
            else:
                if verbose:
                    log.warning("All declination values are NaN; "
                                "assigning dummy value -999.")
                dec = -999.

        else:
            (ra, dec, wcs_wl) = (None, None, None)

        if not got_wavelength:
            wavelength = wcs_wl                 # from wcs, or None

        # Range (slice) of pixel numbers in the dispersion direction.
        disp_range = [slice0, slice1]
        if self.dispaxis == HORIZONTAL:
            image = data
        else:
            image = np.transpose(data, (1, 0))
        if wavelength is None:
            if verbose:
                log.warning("Wavelengths could not be determined.")
            if slice0 <= 0:
                wavelength = np.arange(1, slice1 - slice0 + 1,
                                       dtype=np.float64)
            else:
                wavelength = np.arange(slice0, slice1, dtype=np.float64)

        temp_wl = wavelength.copy()
        nan_mask = np.isnan(wavelength)
        n_nan = nan_mask.sum(dtype=np.intp)
        if n_nan > 0:
            if verbose:
                log.warning("%d NaNs in wavelength array", n_nan)
            temp_wl[nan_mask] = 0.01            # because NaNs cause problems

        # src total flux, area, total weight
        (net, background) = \
        extract1d.extract1d(image, temp_wl, disp_range,
                            self.p_src, self.p_bkg, self.independent_var,
                            self.smoothing_length, self.bkg_order,
                            weights=None)
        del temp_wl

        dq = np.zeros(net.shape, dtype=np.int32)
        if n_nan > 0:
            (wavelength, net, background, dq) = \
                nans_at_endpoints(wavelength, net, background, dq, verbose)

        return (ra, dec, wavelength, net, background, dq)


class ImageExtractModel(ExtractBase):
    """This uses an image that specifies the extraction region."""

    def __init__(self, input_model, slit,
                 ref_file_type=None,
                 match="unknown",
                 spectral_order=1,
                 ref_image=None,
                 dispaxis=HORIZONTAL,
                 smoothing_length=0,
                 nod_correction=0):
        """Extract using a reference image to define the extraction and
           background regions.

        Parameters
        ----------
        input_model : data model
            The input science data.

        slit : an input slit, or a dummy value if not used
            For MultiSlit or MultiProduct data, `slit` is one slit from
            a list of slits in the input.  For other types of data, `slit`
            will not be used.

        ref_file_type : str
            This indicates whether the reference file (if any) was a JSON
            file or an image.

        match : str
            An entry in the reference file (if there is one) should match
            the slit name and the spectral order number of the current
            slit (or it can be "ANY").  `match` will be "exact match" if
            both the name of the current slit and the spectral order number
            match the selected entry in the reference file, and it will be
            "partial match" if only the slit name matches.  If neither
            match, `match` will be "no match".

        ref_image : data model
            The reference image.

        spectral_order : int
            Spectral order number.

        dispaxis : int
            Dispersion direction:  1 is horizontal, 2 is vertical.

        smoothing_length : int
            Width of a boxcar function for smoothing the background
            regions.

        nod_correction : float
            If not zero, the reference image will be shifted in the
            cross-dispersion direction by this number of pixels, rounded
            to an integer.  That is, a feature with pixel location [y0, x0]
            will be moved to [y0 + nod, x0], where `nod` is
            int(round(nod_correction)), if the dispersion direction is
            horizontal.
        """

        super().__init__()

        self.exp_type = input_model.meta.exposure.type
        """
        issue #1781
        self.instrument_name = input_model.meta.instrument.name
        """
        # ref_model contains one or more images; ref_image is the one that
        # matches the current configuration (slit name and spectral order).
        self.ref_image = ref_image
        self.spectral_order = spectral_order
        self.dispaxis = dispaxis
        self.nod_correction = nod_correction

        if smoothing_length is None:
            smoothing_length = 0
        if (smoothing_length > 0 and
            smoothing_length // 2 * 2 == smoothing_length):
            log.warning("smoothing_length was even (%d), so incremented by 1",
                        smoothing_length)
            smoothing_length += 1               # must be odd
        self.smoothing_length = smoothing_length

        if self.exp_type == "NIS_SOSS":
            if hasattr(input_model.meta, 'wcs'):
                try:
                    self.wcs = niriss.niriss_soss_set_input(
                                input_model, self.spectral_order)
                except ValueError:
                    raise InvalidSpectralOrderNumberError(
                                "Spectral order {} is not valid"
                                .format(self.spectral_order))
        elif slit == DUMMY:
            if hasattr(input_model.meta, 'wcs'):
                self.wcs = input_model.meta.wcs
        elif hasattr(slit, 'meta') and hasattr(slit.meta, 'wcs'):
            self.wcs = slit.meta.wcs
        if self.wcs is None:
            log.warning("WCS function not found in input.")


    def add_nod_correction(self, verbose):
        """Shift the reference image (in-place).

        Parameters
        ----------
        verbose : bool
            If True, messages can be logged.
        """

        if self.nod_correction == 0:
            return

        if verbose:
            log.info("Applying nod/dither offset of %s",
                     str(self.nod_correction))

        # Shift the image in the cross-dispersion direction.
        ref = self.ref_image.data.copy()
        shift = self.nod_correction
        ishift = int(round(shift))
        if ishift != shift:
            if verbose:
                log.info("Rounding nod/dither offset of %g to %d",
                         shift, ishift)
        if self.dispaxis == HORIZONTAL:
            if abs(ishift) >= ref.shape[0]:
                if verbose:
                    log.warning("Nod offset %d is too large, skipping ...",
                                ishift)
                return
            self.ref_image.data[:, :] = 0.
            if ishift > 0:
                self.ref_image.data[ishift:, :] = ref[:-ishift, :]
            else:
                ishift = -ishift
                self.ref_image.data[:-ishift, :] = ref[ishift:, :]
        else:
            if abs(ishift) >= ref.shape[1]:
                if verbose:
                    log.warning("Nod offset %d is too large, skipping ...",
                                ishift)
                return
            self.ref_image.data[:, :] = 0.
            if ishift > 0:
                self.ref_image.data[:, ishift:] = ref[:, :-ishift]
            else:
                ishift = -ishift
                self.ref_image.data[:, :-ishift] = ref[:, ishift:]


    def log_extraction_parameters(self):
        """Log the updated extraction parameters."""

        log.debug("Using a reference image that defines extraction regions.")
        log.debug("dispaxis = %d", self.dispaxis)
        log.debug("spectral order = %s", str(self.spectral_order))
        log.debug("smoothing_length = %d", self.smoothing_length)
        log.debug("nod_correction = %s", str(self.nod_correction))


    def extract(self, data, wl_array, verbose):
        """
        Do the actual extraction, for the case that the reference file
        is an image.

        Parameters
        ----------
        data : ndarray, 2-D
            Science data array.

        wl_array : ndarray, 2-D, or None
            Wavelengths corresponding to `data`, or None if no WAVELENGTH
            extension was found in the input file.

        verbose : bool
            If True, log messages.

        Returns
        -------
        ra, dec : float
            ra and dec are the right ascension and declination respectively
            at the nominal center of the slit.

        wavelength : ndarray, 1-D
            The wavelength in micrometers at each pixel.

        net : ndarray, 1-D
            The count rate (counts / s) minus the background at each pixel.

        background : ndarray, 1-D
            The background count rate that was subtracted from the total
            source count rate to get `net`.

        dq : ndarray, 1-D, int32
        """

        shape = data.shape
        # Truncate or expand reference image to match the science data.
        ref = self.match_shape(shape)

        # This is the axis along which to add up the data.
        if self.dispaxis == HORIZONTAL:
            axis = 0
        else:
            axis = 1

        # The values of these arrays will be just 0 or 1.  If ref did not
        # define any background pixels, however, mask_bkg will be None.
        (mask_target, mask_bkg) = self.separate_target_and_background(ref)

        # This is the number of pixels in the cross-dispersion direction,
        # in the target extraction region.
        n_target = mask_target.sum(axis=axis, dtype=np.float)

        # Extract the data.
        gross = (data * mask_target).sum(axis=axis, dtype=np.float)

        # Extract the background.
        if mask_bkg is not None:
            n_bkg = mask_bkg.sum(axis=axis, dtype=np.float)
            # -1 is used as a flag, and also to avoid dividing by zero.
            n_bkg = np.where(n_bkg == 0., -1., n_bkg)
            background = (data * mask_bkg).sum(axis=axis, dtype=np.float)
            # Boxcar smoothing.
            if self.smoothing_length > 1:
                background = extract1d.bxcar(background, self.smoothing_length)
            scalefactor = n_target / n_bkg
            scalefactor = np.where(n_bkg > 0., scalefactor, 0.)
            background *= scalefactor
            net = gross - background
        else:
            background = np.zeros_like(gross)
            net = gross.copy()
        del gross

        if wl_array is None or len(wl_array) == 0:
            got_wavelength = False
        else:
            got_wavelength = True               # may be reset below
        # If wl_array has all 0 values, interpret that to mean that the
        # wavelength attribute was not populated.
        if not got_wavelength or wl_array.min() == 0. and wl_array.max() == 0.:
            got_wavelength = False

        # Used for computing the celestial coordinates and the 1-D array
        # of wavelengths.
        flag = (mask_target > 0.)
        grid = np.indices(shape)
        masked_grid = flag.astype(np.float) * grid[axis]
        g_sum = masked_grid.sum(axis=axis)
        f_sum = flag.sum(axis=axis, dtype=np.float)
        f_sum_zero = np.where(f_sum <= 0.)
        f_sum[f_sum_zero] = 1.                  # to avoid dividing by zero

        spectral_trace = g_sum / f_sum
        del f_sum, g_sum, masked_grid, grid, flag

        # We want x_array and y_array to be 1-D arrays, with the X values
        # initially running from 0 at the left edge of the input cutout to
        # the right edge, and the Y values being near the middle of
        # the spectral extraction region.  So the locations
        # (x_array[i], y_array[i]) should be the spectral trace.  Near the
        # left and right edges, there might not be any non-zero values in
        # mask_target, so a slice will be extracted from both x_array and
        # y_array in order to exclude pixels that are not within the
        # extraction region.
        if self.dispaxis == HORIZONTAL:
            x_array = np.arange(shape[1], dtype=np.float)
            y_array = spectral_trace
        else:
            x_array = spectral_trace
            y_array = np.arange(shape[0], dtype=np.float)

        # Trim off the ends, if there's no data there.  Save trim_slc.
        mask = np.where(n_target > 0.)
        if len(mask[0]) > 0:
            trim_slc = slice(mask[0][0], mask[0][-1] + 1)
            net = net[trim_slc]
            background = background[trim_slc]
            n_target = n_target[trim_slc]
            x_array = x_array[trim_slc]
            y_array = y_array[trim_slc]

        if got_wavelength:
            if verbose:
                log.debug("Wavelengths are from wavelength attribute.")
            indx = np.around(x_array).astype(np.int)
            indy = np.around(y_array).astype(np.int)
            indx = np.where(indx < 0, 0, indx)
            indx = np.where(indx >= shape[1], shape[1] - 1, indx)
            indy = np.where(indy < 0, 0, indy)
            indy = np.where(indy >= shape[0], shape[0] - 1, indy)
            wavelength = wl_array[indy, indx]

        nelem = len(x_array)

        if self.wcs is not None:
            if verbose and not got_wavelength:
                log.debug("Wavelengths are from the wcs function.")
            if self.exp_type in WFSS_EXPTYPES:
                # We expect two (x and y) or three (x, y, spectral order).
                n_inputs = self.wcs.forward_transform.n_inputs
                ra = np.zeros(nelem, dtype=np.float)
                dec = np.zeros(nelem, dtype=np.float)
                # Temporary variable so as not to clobber `wavelength`.
                wcs_wl = np.zeros(nelem, dtype=np.float)
                transform = self.wcs.forward_transform
                if n_inputs == 2:
                    for i in range(nelem):
                        ra[i], dec[i], wcs_wl[i], _ = transform(
                                        x_array[i], y_array[i])
                elif n_inputs == 3:
                    for i in range(nelem):
                        ra[i], dec[i], wcs_wl[i], _ = transform(
                                x_array[i], y_array[i], self.spectral_order)
                else:
                    if verbose:
                        log.warning("n_inputs for wcs function is %d",
                                    n_inputs)
                        log.warning("WCS function was expected to take "
                                    "either 2 or 3 arguments.")
                    ra[:] = -999.
                    dec[:] = -999.
                    wcs_wl[:] = -999.
            else:
                """
                See issue #1781
                if self.instrument_name == "NIRSPEC":
                    # xxx temporary:  NIRSpec wcs is one-based.
                    ra, dec, wcs_wl = self.wcs(x_array + 1., y_array + 1.)
                else:
                    ra, dec, wcs_wl = self.wcs(x_array, y_array)
                """
                ra, dec, wcs_wl = self.wcs(x_array, y_array)
            # We need one right ascension and one declination, representing
            # the direction of pointing.
            middle = ra.shape[0] // 2           # ra and dec have same shape
            mask = np.isnan(ra)
            not_nan = np.logical_not(mask)
            if not_nan[middle]:
                if verbose:
                    log.debug("Using midpoint of spectral trace "
                              "for RA and Dec.")
                ra = ra[middle]
            else:
                if np.any(not_nan):
                    if verbose:
                        log.warning("Midpoint of coordinate array is NaN; "
                                    "using the average of non-NaN min and "
                                    "max values.")
                    ra = (np.nanmin(ra) + np.nanmax(ra)) / 2.
                else:
                    if verbose:
                        log.warning("All right ascension values are NaN; "
                                    "assigning dummy value -999.")
                    ra = -999.
            mask = np.isnan(dec)
            not_nan = np.logical_not(mask)
            if not_nan[middle]:
                dec = dec[middle]
            else:
                if np.any(not_nan):
                    dec = (np.nanmin(dec) + np.nanmax(dec)) / 2.
                else:
                    if verbose:
                        log.warning("All declination values are NaN; "
                                    "assigning dummy value -999.")
                    dec = -999.

        else:
            (ra, dec, wcs_wl) = (None, None, None)

        if not got_wavelength:
            wavelength = wcs_wl                 # from wcs, or None
        if wavelength is None:
            if self.dispaxis == HORIZONTAL:
                wavelength = np.arange(shape[1], dtype=np.float)
            else:
                wavelength = np.arange(shape[0], dtype=np.float)
            wavelength = wavelength[trim_slc]

        dq = np.zeros(net.shape, dtype=np.int32)
        nan_mask = np.isnan(wavelength)
        n_nan = nan_mask.sum(dtype=np.intp)
        if n_nan > 0:
            if verbose:
                log.warning("%d NaNs in wavelength array", n_nan)
            (wavelength, net, background, dq) = \
                nans_at_endpoints(wavelength, net, background, dq, verbose)

        return (ra, dec, wavelength, net, background, dq)


    def match_shape(self, shape):
        """Truncate or expand reference image to match the science data.

        Extended summary
        ----------------
        The science data may be 2-D or 3-D, but the reference image only
        needs to be 2-D.

        Parameters
        ----------
        shape : tuple
            The shape of the science data.

        Returns
        -------
        ndarray, 2-D
            This is either the reference image (the data array, not the
            complete data model), or an array of the same type, either
            larger or smaller than the actual reference image, but matching
            the science data array both in size and location on the
            detector.
        """

        ref = self.ref_image.data

        ref_shape = ref.shape
        if shape == ref_shape:
            return ref

        # This is the shape of the last two axes of the science data.
        buf = np.zeros((shape[-2], shape[-1]), dtype=ref.dtype)
        y_max = min(shape[-2], ref_shape[0])
        x_max = min(shape[-1], ref_shape[1])
        slc0 = slice(0, y_max)
        slc1 = slice(0, x_max)
        buf[slc0, slc1] = ref[slc0, slc1].copy()

        return buf


    def separate_target_and_background(self, ref):
        """Create masks for source and background.

        Parameters
        ----------
        ref : ndarray, 2-D
            This is the reference image as returned by `match_shape`,
            i.e. it might be a subset of the original reference image.

        Returns
        -------
        mask_target : ndarray, 2-D
            This is an array of the same type and shape as the science
            image, but with values of only 0 or 1.  A value of 1 indicates
            that the corresponding pixel in the science data array should
            be included when adding up values to make the 1-D spectrum,
            and a value of 0 means that it should not be included.

        mask_bkg : ndarray, 2-D, or None.
            This is like `mask_target` but for background regions.
            A negative value in the reference image flags a pixel that
            should be included in the background region(s).  If there is
            no pixel in the reference image with a negative value,
            `mask_bkg` will be set to None.
        """

        mask_target = np.where(ref > 0., 1., 0.)

        if np.any(ref < 0.):
            mask_bkg = np.where(ref < 0., 1., 0.)
        else:
            mask_bkg = None

        return (mask_target, mask_bkg)


def interpolate_response(wavelength, relsens, verbose):
    """Interpolate within the relative response table.

    Parameters
    ----------
    wavelength : ndarray, 1-D
        Wavelengths in the science data

    relsens : record array
        Contains two columns, 'wavelength' and 'response'.

    verbose : bool
        If True, write log messages.

    Returns
    -------
    r_factor : ndarray, 1-D
        The response, interpolated at `wavelength`, with extrapolated
        elements and zero or negative response values set to 1.  Divide
        the net count rate by r_factor to obtain the flux.
    """

    # "_relsens" indicates that the values were read from the RELSENS table.
    wl_relsens = relsens['wavelength']
    resp_relsens = relsens['response']
    MICRONS_100 = 1.e-4                 # 100 microns, in meters
    if wl_relsens.max() > 0. and wl_relsens.max() < MICRONS_100:
        if verbose:
            log.warning("Converting RELSENS wavelengths to microns.")
        wl_relsens *= 1.e6

    bad = False
    if np.any(np.isnan(wl_relsens)):
        log.error("In RELSENS, the 'wavelength' column contains NaNs.")
        bad = True
    if np.any(np.isnan(resp_relsens)):
        log.error("In RELSENS, the 'response' column contains NaNs.")
        bad = True
    if bad:
        raise ValueError("Found NaNs in RELSENS table.")

    # np.interp requires that wl_relsens be increasing.
    if wl_relsens[-1] < wl_relsens[0]:
        if verbose:
            log.warning("The wavelength column in RELSENS was decreasing.")
        wl_relsens = wl_relsens[::-1].copy()
        resp_relsens = resp_relsens[::-1].copy()

    # `r_factor` is the response, interpolated at the wavelengths in the
    # science data.  -2048 is a flag value, to check for extrapolation.
    r_factor = np.interp(wavelength, wl_relsens, resp_relsens, -2048., -2048.)
    mask = np.where(r_factor == -2048.)
    if len(mask[0]) > 0:
        if verbose:
            log.warning("Using RELSENS, %d elements were extrapolated; "
                        "these values will be set to 1.", len(mask[0]))
        r_factor[mask] = 1.
    mask = np.where(r_factor <= 0.)
    if len(mask[0]) > 0:
        if verbose:
            log.warning("Using RELSENS, %d interpolated response values "
                        "were <= 0; these values will be set to 1.",
                        len(mask[0]))
        r_factor[mask] = 1.

    return r_factor


def do_extract1d(input_model, refname, smoothing_length, bkg_order,
                 log_increment):
    """Extract 1-D spectra.

    Extended summary
    ----------------

    Parameters
    ----------
    input_model : data model
        The input science model.

    refname : str
        The name of the reference file, or "N/A".

    smoothing_length : int
        Width of a boxcar function for smoothing the background regions.

    bkg_order : int
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background.

    log_increment : int
        if `log_increment` is greater than 0 and the input data are
        multi-integration, a message will be written to the log every
        `log_increment` integrations.

    Returns
    -------
    output_model : data model
        A new MultiSpecModel containing the extracted spectra.
    """

    output_model = datamodels.MultiSpecModel()
    if hasattr(input_model, "int_times"):
        output_model.int_times = input_model.int_times.copy()
    output_model.update(input_model)

    # This will be relevant if we're asked to extract a spectrum and the
    # spectral order is zero.  That's only OK if the disperser is a prism.
    prism_mode = is_prism(input_model)

    # Read and interpret the reference file.
    ref_dict = load_ref_file(refname)

    if isinstance(input_model, datamodels.MultiSlitModel) or \
       isinstance(input_model, datamodels.MultiProductModel):

        if isinstance(input_model, datamodels.MultiSlitModel):
            slits = input_model.slits
        else:                           # MultiProductModel
            slits = input_model.products

        # Loop over the slits in the input model
        for slit in slits:
            log.info('Working on slit %s', slit.name)
            prev_offset = OFFSET_NOT_ASSIGNED_YET
            if np.size(slit.data) <= 0:
                log.info('No data for slit %s, skipping ...', slit.name)
                continue
            sp_order = get_spectral_order(slit)
            if sp_order == 0 and not prism_mode:
                log.info("Spectral order 0 is a direct image, skipping ...")
                continue
            extract_params = get_extract_parameters(
                                ref_dict,
                                slit, slit.name, sp_order,
                                input_model.meta, smoothing_length, bkg_order)
            if extract_params['match'] == NO_MATCH:
                log.critical('Missing extraction parameters.')
                raise ValueError('Missing extraction parameters.')
            elif extract_params['match'] == PARTIAL:
                log.info('Spectral order %d not found, skipping ...', sp_order)
                continue
            find_dispaxis(input_model, slit, sp_order, extract_params)
            if extract_params['dispaxis'] is None:
                log.warning("The dispersion direction couldn't be determined, "
                            "so skipping ...")
                continue

            try:
                (ra, dec, wavelength, net, background, dq,
                 prev_offset) = extract_one_slit(
                                        input_model, slit, -1,
                                        prev_offset, True, extract_params)
            except InvalidSpectralOrderNumberError as e:
                log.info(str(e) + ", skipping ...")
                continue
            got_relsens = True
            try:
                relsens = slit.relsens
            except AttributeError:
                got_relsens = False
            if got_relsens and len(relsens) == 0:
                got_relsens = False
            if got_relsens:
                r_factor = interpolate_response(wavelength, relsens, True)
                flux = net / r_factor
            else:
                log.warning("No relsens for current slit, "
                            "so can't compute flux.")
                flux = np.zeros_like(net)
            fl_error = np.ones_like(net)
            nerror = np.ones_like(net)
            berror = np.ones_like(net)
            spec = datamodels.SpecModel()
            otab = np.array(list(zip(wavelength, flux, fl_error, dq,
                                 net, nerror, background, berror)),
                            dtype=spec.spec_table.dtype)
            spec = datamodels.SpecModel(spec_table=otab)
            spec.meta.wcs = spec_wcs.create_spectral_wcs(ra, dec, wavelength)
            spec.spec_table.columns['wavelength'].unit = 'um'
            spec.spec_table.columns['flux'].unit = 'mJy'
            spec.spec_table.columns['error'].unit = 'mJy'
            spec.spec_table.columns['net'].unit = 'DN/s'
            spec.spec_table.columns['nerror'].unit = 'DN/s'
            spec.spec_table.columns['background'].unit = 'DN/s'
            spec.spec_table.columns['berror'].unit = 'DN/s'
            spec.slit_ra = ra
            spec.slit_dec = dec
            spec.spectral_order = sp_order
            copy_keyword_info(slit, slit.name, spec)
            output_model.spec.append(spec)
    else:
        slitname = input_model.meta.exposure.type
        if slitname is None:
            slitname = ANY
        if slitname == 'NIS_SOSS':
            slitname = input_model.meta.subarray.name
        log.debug('slitname=%s', slitname)

        # Loop over these spectral order numbers.
        if input_model.meta.exposure.type == "NIS_SOSS":
            # This list of spectral order numbers may need to be assigned
            # differently for other exposure types.
            spectral_order_list = [1, 2, 3]
        else:
            # For this case, we'll call get_spectral_order to get the order.
            spectral_order_list = ["not set yet"]

        if isinstance(input_model, (datamodels.ImageModel,
                                    datamodels.DrizProductModel)):
            prev_offset = OFFSET_NOT_ASSIGNED_YET
            for sp_order in spectral_order_list:
                if sp_order == "not set yet":
                    sp_order = get_spectral_order(input_model)
                if sp_order == 0 and not prism_mode:
                    log.info("Spectral order 0 is a direct image, "
                             "skipping ...")
                    continue

                extract_params = get_extract_parameters(
                                    ref_dict,
                                    input_model, slitname, sp_order,
                                    input_model.meta, smoothing_length,
                                    bkg_order)
                if extract_params['match'] == EXACT:
                    slit = DUMMY
                    find_dispaxis(input_model, slit, sp_order, extract_params)
                    if extract_params['dispaxis'] is None:
                        log.warning("The dispersion direction couldn't be "
                                    "determined, so skipping ...")
                        continue
                    try:
                        (ra, dec, wavelength, net, background, dq,
                         prev_offset) = extract_one_slit(
                                        input_model, slit, -1,
                                        prev_offset, True, extract_params)
                    except InvalidSpectralOrderNumberError as e:
                        log.info(str(e) + ", skipping ...")
                        continue
                elif extract_params['match'] == PARTIAL:
                    log.info('Spectral order %d not found, skipping ...',
                             sp_order)
                    continue
                else:
                    log.critical('Missing extraction parameters.')
                    raise ValueError('Missing extraction parameters.')
                got_relsens = True
                try:
                    relsens = input_model.relsens
                except AttributeError:
                    got_relsens = False
                if got_relsens and len(relsens) == 0:
                    got_relsens = False
                if got_relsens:
                    r_factor = interpolate_response(wavelength, relsens, True)
                    flux = net / r_factor
                else:
                    log.warning("No relsens for input file, "
                                "so can't compute flux.")
                    flux = np.zeros_like(net)
                fl_error = np.ones_like(net)
                nerror = np.ones_like(net)
                berror = np.ones_like(net)
                spec = datamodels.SpecModel()
                otab = np.array(list(zip(wavelength, flux, fl_error, dq,
                                     net, nerror, background, berror)),
                                dtype=spec.spec_table.dtype)
                spec = datamodels.SpecModel(spec_table=otab)
                spec.meta.wcs = spec_wcs.create_spectral_wcs(
                                        ra, dec, wavelength)
                spec.spec_table.columns['wavelength'].unit = 'um'
                spec.spec_table.columns['flux'].unit = 'mJy'
                spec.spec_table.columns['error'].unit = 'mJy'
                spec.spec_table.columns['net'].unit = 'DN/s'
                spec.spec_table.columns['nerror'].unit = 'DN/s'
                spec.spec_table.columns['background'].unit = 'DN/s'
                spec.spec_table.columns['berror'].unit = 'DN/s'
                spec.slit_ra = ra
                spec.slit_dec = dec
                spec.spectral_order = sp_order
                if slitname is not None and slitname != "ANY":
                    spec.name = slitname
                output_model.spec.append(spec)

        elif isinstance(input_model, (datamodels.CubeModel,
                                      datamodels.SlitModel)):

            slit = DUMMY

            # NRS_BRIGHTOBJ exposures are instances of SlitModel.
            prev_offset = OFFSET_NOT_ASSIGNED_YET
            for sp_order in spectral_order_list:
                if sp_order == "not set yet":
                    sp_order = get_spectral_order(input_model)
                    if sp_order == 0 and not prism_mode:
                        log.info("Spectral order 0 is a direct image, "
                                 "skipping ...")
                        continue

                extract_params = get_extract_parameters(
                                    ref_dict,
                                    input_model, slitname, sp_order,
                                    input_model.meta, smoothing_length,
                                    bkg_order)
                if extract_params['match'] == NO_MATCH:
                    log.critical('Missing extraction parameters.')
                    raise ValueError('Missing extraction parameters.')
                elif extract_params['match'] == PARTIAL:
                    log.warning('Spectral order %d not found, skipping ...',
                                sp_order)
                    continue
                find_dispaxis(input_model, slit, sp_order, extract_params)
                if extract_params['dispaxis'] is None:
                    log.warning("The dispersion direction couldn't be "
                                "determined, so skipping ...")
                    continue

                got_relsens = True
                try:
                    relsens = input_model.relsens
                except AttributeError:
                    got_relsens = False
                if got_relsens and len(relsens) == 0:
                    got_relsens = False
                if not got_relsens:
                    log.warning("No relsens for input file, "
                                "so can't compute flux.")

                # Loop over each integration in the input model
                verbose = True          # for just the first integration
                if input_model.data.shape[0] == 1:
                    log.info("Beginning loop, just 1 integration ...")
                else:
                    log.info("Beginning loop over %d integrations ...",
                             input_model.data.shape[0])
                for integ in range(input_model.data.shape[0]):
                    # Extract spectrum
                    try:
                        (ra, dec, wavelength, net, background, dq,
                         prev_offset) = extract_one_slit(
                                        input_model, slit, integ,
                                        prev_offset, verbose, extract_params)
                    except InvalidSpectralOrderNumberError as e:
                        log.info(str(e) + ", skipping ...")
                        break
                    if got_relsens:
                        r_factor = interpolate_response(
                                        wavelength, input_model.relsens,
                                        verbose)
                        flux = net / r_factor
                    else:
                        flux = np.zeros_like(net)
                    fl_error = np.ones_like(net)
                    nerror = np.ones_like(net)
                    berror = np.ones_like(net)
                    spec = datamodels.SpecModel()
                    otab = np.array(list(zip(wavelength, flux, fl_error, dq,
                                         net, nerror, background, berror)),
                                    dtype=spec.spec_table.dtype)
                    spec = datamodels.SpecModel(spec_table=otab)
                    spec.meta.wcs = spec_wcs.create_spectral_wcs(
                                        ra, dec, wavelength)
                    spec.spec_table.columns['wavelength'].unit = 'um'
                    spec.spec_table.columns['flux'].unit = 'mJy'
                    spec.spec_table.columns['error'].unit = 'mJy'
                    spec.spec_table.columns['net'].unit = 'DN/s'
                    spec.spec_table.columns['nerror'].unit = 'DN/s'
                    spec.spec_table.columns['background'].unit = 'DN/s'
                    spec.spec_table.columns['berror'].unit = 'DN/s'
                    spec.slit_ra = ra
                    spec.slit_dec = dec
                    spec.spectral_order = sp_order
                    output_model.spec.append(spec)

                    if (log_increment > 0 and
                        (integ + 1) % log_increment == 0):
                            if integ == 0:
                                if input_model.data.shape[0] == 1:
                                    log.info("1 integration done")
                                else:
                                    log.info("... 1 integration done")
                            elif integ == input_model.data.shape[0] - 1:
                                log.info("All %d integrations done",
                                         input_model.data.shape[0])
                            else:
                                log.info("... %d integrations done", integ + 1)
                            progress_msg_printed = True
                    else:
                            progress_msg_printed = False
                    verbose = False

                if not progress_msg_printed:
                    if input_model.data.shape[0] == 1:
                        log.info("1 integration done")
                    else:
                        log.info("All %d integrations done",
                                 input_model.data.shape[0])

        elif isinstance(input_model, datamodels.IFUCubeModel):

            try:
                source_type = input_model.meta.target.source_type.lower()
            except AttributeError:
                source_type = "unknown"
            output_model = ifu.ifu_extract1d(input_model, refname, source_type)

        else:
            log.error("The input file is not supported for this step.")
            raise RuntimeError("Can't extract a spectrum from this file.")

    # Copy the integration time information from the INT_TIMES table
    # to keywords in the output file.
    if pipe_utils.is_tso(input_model):
        populate_time_keywords(input_model, output_model)
    else:
        log.debug("Not copying from the INT_TIMES table because "
                  "this is not a TSO exposure.")

    # See output_model.spec[i].meta.wcs instead.
    output_model.meta.wcs = None

    # If the reference file is an image, explicitly close it.
    if ref_dict is not None and 'ref_model' in ref_dict:
        ref_dict['ref_model'].close()

    return output_model


def populate_time_keywords(input_model, output_model):
    """Copy the integration times keywords to header keywords.

    Parameters
    ----------
    input_model : data model
        The input science model.

    output_model : data model
        The output science model.  This may be modified in-place.
    """

    nints = input_model.meta.exposure.nints
    int_start = input_model.meta.exposure.integration_start
    if int_start is None:
        log.warning("INTSTART not found; assuming a value of 1.")
        int_start = 1
    int_start -= 1                              # zero indexed
    int_end = input_model.meta.exposure.integration_end
    if int_end is None:
        log.warning("INTEND not found; assuming a value of %d.", nints)
        int_end = nints
    int_end -= 1                                # zero indexed
    if nints > 1:
        num_integrations = int_end - int_start + 1
    else:
        num_integrations = 1

    if hasattr(input_model, 'int_times') and input_model.int_times is not None:
        nrows = len(input_model.int_times)
    else:
        nrows = 0
    if nrows < 1:
        log.warning("There is no INT_TIMES table in the input file.")
        return

    # If we have a single plane (e.g. ImageModel or MultiSlitModel),
    # we will only populate the keywords if the corresponding uncal file
    # had one integration.  If the data were or might have been segmented,
    # we use the first and last integration numbers to determine whether
    # the data were in fact averaged over integrations, and if so, we
    # should not populate the int_times-related header keywords.

    skip = False                        # initial value

    if isinstance(input_model, (datamodels.MultiSlitModel,
                                datamodels.MultiProductModel,
                                datamodels.ImageModel,
                                datamodels.DrizProductModel)):
        if num_integrations > 1:
            log.warning("Not using INT_TIMES table because the data "
                        "have been averaged over integrations.")
            skip = True
    elif isinstance(input_model, (datamodels.CubeModel,
                                  datamodels.SlitModel)):
        shape = input_model.data.shape
        if len(shape) == 2 and num_integrations > 1:
            log.warning("Not using INT_TIMES table because the data "
                        "have been averaged over integrations.")
            skip = True
        elif len(shape) != 3 or shape[0] > nrows:
            # Later, we'll check that the integration_number column actually
            # has a row corresponding to every integration in the input.
            log.warning("Not using INT_TIMES table because the data shape "
                        "is not consistent with the number of table rows.")
            skip = True
    elif isinstance(input_model, datamodels.IFUCubeModel):
        log.warning("The INT_TIMES table will be ignored for IFU data.")
        skip = True

    if skip:
        return

    int_num = input_model.int_times['integration_number']
    start_utc = input_model.int_times['int_start_MJD_UTC']
    mid_utc = input_model.int_times['int_mid_MJD_UTC']
    end_utc = input_model.int_times['int_end_MJD_UTC']
    start_tdb = input_model.int_times['int_start_BJD_TDB']
    mid_tdb = input_model.int_times['int_mid_BJD_TDB']
    end_tdb = input_model.int_times['int_end_BJD_TDB']

    # Inclusive range of integration numbers in the input data,
    # zero indexed.
    data_range = (int_start, int_end)
    # Inclusive range of integration numbers in the INT_TIMES table,
    # zero indexed.
    table_range = (int_num[0] - 1, int_num[-1] - 1)
    offset = data_range[0] - table_range[0]
    if data_range[0] < table_range[0] or data_range[1] > table_range[1]:
        log.warning("Not using the INT_TIMES table because it does not "
                    "include rows for all integrations in the data.")
        return

    log.debug("TSO data, so copying times from the INT_TIMES table.")

    if hasattr(input_model, 'data'):
        shape = input_model.data.shape
        if len(shape) == 2:
            num_integ = 1
        else:                                   # len(shape) == 3
            num_integ = shape[0]
    else:                                       # e.g. MultiSlit data
        num_integ = 1

    # This assumes that the spec attribute of output_model has already
    # been created, and spectra have been appended.
    n_output_spec = len(output_model.spec)

    # num_j is the number of spectra per integration, e.g. the number of
    # fixed-slit spectra, MSA spectra, or different spectral orders;
    # num_integ is the number of integrations.
    # The total number of output spectra is n_output_spec = num_integ * num_j
    num_j = n_output_spec // num_integ
    if n_output_spec != num_j * num_integ:      # sanity check
        log.warning("populate_time_keywords:  Don't understand "
                    "n_output_spec = %d, num_j = %d, num_integ = %d",
                    n_output_spec, num_j, num_integ)
    else:
        log.debug("Number of output spectra = %d; "
                  "number of spectra for each integration = %d; "
                  "number of integrations = %d",
                  n_output_spec, num_j, num_integ)

    # n is a counter for spectra in output_model.
    n = 0
    for j in range(num_j):                      # for each spectrum or order
        for k in range(num_integ):                  # for each integration
            row = k + offset
            spec = output_model.spec[n]             # n is incremented below
            spec.int_num = int_num[row]
            spec.time_scale = "UTC"
            spec.start_utc = start_utc[row]
            spec.mid_utc = mid_utc[row]
            spec.end_utc = end_utc[row]
            spec.start_tdb = start_tdb[row]
            spec.mid_tdb = mid_tdb[row]
            spec.end_tdb = end_tdb[row]
            n += 1


def get_spectral_order(slit):
    """Get the spectral order number.

    Parameters
    ----------
    slit : SlitModel object
        One slit from an input MultiSlitModel or similar.

    Returns
    -------
    int
        Spectral order number for `slit`.  If no information about spectral
        order is available in `wcsinfo`, a default value of 1 will be
        returned.
    """

    if hasattr(slit.meta, 'wcsinfo'):
        sp_order = slit.meta.wcsinfo.spectral_order
        if sp_order is None:
            log.warning("spectral_order is None; using 1")
            sp_order = 1
    else:
        log.warning("slit.meta doesn't have attribute wcsinfo; "
                    "setting spectral order to 1")
        sp_order = 1

    return sp_order


def is_prism(input_model):
    """Determine whether the current observing mode used a prism.

    Extended summary
    ----------------
    The reason for this test is so we can skip spectral extraction if the
    spectral order is zero and the exposure was not made using a prism.  
    In this context, therefore, a grism is not considered to be a prism.

    Parameters
    ----------
    input_model : data model
        The input science model.

    Returns
    -------
    bool
        True if the exposure used a prism; False otherwise.
    """

    detector = input_model.meta.instrument.detector
    if detector is None:
        return False

    filter = input_model.meta.instrument.filter
    if filter is None:
        filter = "NONE"
    else:
        filter = filter.upper()
    grating = input_model.meta.instrument.grating
    if grating is None:
        grating = "NONE"
    else:
        grating = grating.upper()

    prism_mode = False
    if (detector.startswith("MIR") and filter.find("P750L") >= 0 or
        detector.startswith("NRS") and grating.find("PRISM") >= 0):
            prism_mode = True

    return prism_mode


def copy_keyword_info(slit, slitname, spec):
    """Copy metadata from the input to the output spectrum.

    Parameters
    ----------
    slit : A SlitModel object
        Metadata will be copied from the input `slit` to output `spec`.

    slitname : str or None
        The name of the slit.

    spec : One element of MultiSpecModel.spec
        Metadata attributes will be updated in-place.
    """

    if slitname is not None and slitname != "ANY":
        spec.name = slitname

    if hasattr(slit, "slitlet_id"):
        spec.slitlet_id = slit.slitlet_id

    if hasattr(slit, "source_id"):
        spec.source_id = slit.source_id

    if hasattr(slit, "source_name") and slit.source_name is not None:
        spec.source_name = slit.source_name

    if hasattr(slit, "source_alias") and slit.source_alias is not None:
        spec.source_alias = slit.source_alias

    if hasattr(slit, "source_type") and slit.source_type is not None:
        spec.source_type = slit.source_type

    if hasattr(slit, "stellarity") and slit.stellarity is not None:
        spec.stellarity = slit.stellarity

    if hasattr(slit, "source_xpos"):
        spec.source_xpos = slit.source_xpos

    if hasattr(slit, "source_ypos"):
        spec.source_ypos = slit.source_ypos

    if hasattr(slit, "shutter_state"):
        spec.shutter_state = slit.shutter_state


def extract_one_slit(input_model, slit, integ,
                     prev_offset, verbose, extract_params):
    """Extract data for one slit, or spectral order, or plane.

    Parameters
    ----------
    input_model : data model
        The input science model.

    slit : one slit from a MultiSlitModel (or similar), or "dummy"
        If slit is "dummy", the data array is input_model.data; otherwise,
        the data array is slit.data.
        In the former case, if `integ` is zero or larger, the spectrum
        will be extracted from the 2-D slice input_model.data[integ].

    integ : int
        For the case that input_model is a SlitModel or a CubeModel,
        `integ` is the integration number.  If the integration number is
        not relevant (i.e. the data array is 2-D), `integ` should be -1.

    prev_offset : float or str
        When extracting from multi-integration data, the nod/dither offset
        only needs to be determined once.  `prev_offset` is either the
        previously computed offset or a value (a string) indicating that
        the offset hasn't been computed yet.  In the latter case, method
        `offset_from_offset` will be called to determine the offset.

    verbose : boolean
        If True, log more info (extraction parameters, for example).

    extract_params : dict
        Parameters read from the reference file.

    Returns
    -------
    ra, dec : float
        ra and dec are the right ascension and declination respectively
        at the nominal center of the slit.

    wavelength : ndarray, 1-D
        The wavelength in micrometers at each pixel.

    net : ndarray, 1-D
        The count rate (counts / s) minus the background at each pixel.

    background : ndarray, 1-D
        The background count rate that was subtracted from the total
        source count rate to get `net`.

    dq : ndarray, 1-D, int32
        The data quality array.

    offset : float
       The nod/dither offset in the cross-dispersion direction, either
        computed by calling `offset_from_offset` in this function, or
        copied from the input `prev_offset`.
    """

    if verbose:
        log_initial_parameters(extract_params)

    input_dq = None                             # possibly replaced below
    if integ > -1:
        data = input_model.data[integ]
        if hasattr(input_model, 'dq'):
            input_dq = input_model.dq[integ]
        try:
            wl_array = input_model.wavelength
        except AttributeError:
            wl_array = None
    elif slit == DUMMY:
        data = input_model.data
        if hasattr(input_model, 'dq'):
            input_dq = input_model.dq
        try:
            wl_array = input_model.wavelength
        except AttributeError:
            wl_array = None
    else:
        data = slit.data
        if hasattr(slit, 'dq'):
            input_dq = slit.dq
        try:
            wl_array = slit.wavelength
        except AttributeError:
            wl_array = None

    data = replace_bad_values(data, input_dq, fill=0.)

    if extract_params['ref_file_type'] == FILE_TYPE_IMAGE:
        # The reference file is an image.
        extract_model = ImageExtractModel(input_model, slit, **extract_params)
        ap = None
    else:
        # If there is a reference file (there doesn't have to be), it's in
        # JSON format.
        extract_model = ExtractModel(input_model, slit, **extract_params)
        ap = get_aperture(data.shape, extract_model.wcs,
                          verbose, extract_params)
        extract_model.update_extraction_limits(ap)

    # Only call this method for the first integration.
    if prev_offset == OFFSET_NOT_ASSIGNED_YET:
        offset = extract_model.offset_from_offset(input_model, slit, verbose)
        if offset != 0:                         # xxx should be temporary
            if verbose:
                log.debug("Computed nod/dither offset = %s, but don't "
                          "trust this yet, so assuming 0", str(offset))
            offset = 0.                         # xxx should be temporary
    else:
        offset = prev_offset
    extract_model.nod_correction = offset

    # Add the nod/dither offset to the polynomial coefficients, or shift
    # the reference image (depending on the type of reference file).
    extract_model.add_nod_correction(verbose)

    if verbose:
        extract_model.log_extraction_parameters()

    extract_model.assign_polynomial_limits(verbose)
    (ra, dec, wavelength, net, background, dq) = \
                extract_model.extract(data, wl_array, verbose)

    return (ra, dec, wavelength, net, background, dq, offset)


def replace_bad_values(data, input_dq, fill=0.):
    """Replace NaNs and values flagged with DO_NOT_USE.

    Parameters
    ----------
    data : ndarray
        The input data array.

    input_dq : ndarray or None
        If not None, this will be checked for flag value DO_NOT_USE.

    fill : float
        Pixels that are NaN in `data` or are flagged in the `input_dq`
        array (if the latter is not None) will be assigned this value.

    Returns
    -------
    ndarray
        A possibly modified copy of `data`.
    """

    mask = np.isnan(data)
    if input_dq is not None:
        bad_mask = np.bitwise_and(input_dq, dqflags.pixel['DO_NOT_USE']) > 0
        mask = np.logical_or(mask, bad_mask)

    if np.any(mask):
        mod_data = data.copy()
        mod_data[mask] = fill
        return mod_data
    else:
        return data

def nans_at_endpoints(wavelength, net, background, dq, verbose):
    """Flag NaNs in the wavelength array.

    Extended summary
    ----------------
    All five input arrays should be 1-D and have the same shape.
    If NaNs are present at endpoints of `wavelength`, the arrays will be
    trimmed to remove the NaNs.  NaNs at interior elements of `wavelength`
    will be left in place, but they will be flagged with DO_NOT_USE in the
    `dq` array.

    Parameters
    ----------
    wavelength : ndarray
        Array of wavelengths, possibly containing NaNs.

    net : ndarray
        Array of net count rates.

    background : ndarray
        Array of background values that were subtracted to get `net`.

    dq : ndarray
        Data quality array.

    verbose : bool
        If True and the arrays were trimmed, log a message.

    Returns
    -------
    wavelength, net, background, dq : ndarray
        The returned `dq` array may have NaNs flagged with DO_NOT_USE,
        and all four arrays may have been trimmed at either or both ends.
    """

    # The input arrays will not be modified in-place.
    new_wl = wavelength.copy()
    new_net = net.copy()
    new_bkg = background.copy()
    new_dq = dq.copy()
    nelem = wavelength.shape[0]

    nan_mask = np.isnan(wavelength)

    new_dq[nan_mask] = np.bitwise_or(new_dq[nan_mask],
                                     dqflags.pixel['DO_NOT_USE'])
    not_nan = np.logical_not(nan_mask)
    flag = np.where(not_nan)
    if len(flag[0]) > 0:
        n_trimmed = flag[0][0] + nelem - (flag[0][-1] + 1)
        if n_trimmed > 0:
            if verbose:
                log.info("Output arrays have been trimmed by %d elements",
                         n_trimmed)
            slc = slice(flag[0][0], flag[0][-1] + 1)
            new_wl = new_wl[slc]
            new_net = new_net[slc]
            new_bkg = new_bkg[slc]
            new_dq = new_dq[slc]
    else:
        new_dq |= dqflags.pixel['DO_NOT_USE']

    return (new_wl, new_net, new_bkg, new_dq)
