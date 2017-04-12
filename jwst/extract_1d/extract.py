from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging
from collections import namedtuple
import copy
import math

import numpy as np
import json
from astropy.modeling import polynomial
from .. import datamodels
from . import extract1d
from . import ifu

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# For full-frame input data, keyword SLTNAME may not be populated, so use
# the following string to indicate that the first slit in the reference
# file should be used.
# A slit name in the reference file can also be ANY, in which case that
# slit can be used with any slit name from the input data.
ANY = "ANY"

# Dispersion direction, predominantly horizontal or vertical.  These values
# are to be compared with keyword DISPAXIS from the input header.
HORIZONTAL = 1
VERTICAL = 2


Aperture = namedtuple('Aperture', ['xstart', 'ystart', 'xstop', 'ystop'])


def initialize_wave_model(model_name, degree):
    return getattr(polynomial, model_name)(degree)


def get_extract_parameters(refname, slitname, meta,
                           smoothing_length, bkg_order):
    extract_params = {}
    with open(refname) as f:
        ref = json.load(f)
    for aper in ref['apertures']:
        if 'id' in aper and aper['id'] != "dummy" and \
           (aper['id'] == slitname or aper['id'] == "ANY" or
            slitname == "ANY"):
            region_type = aper.get("region_type", "target")
            if region_type == "target":
                disp = aper.get('dispaxis')
                if disp is None:
                    log.warning("dispaxis not specified in %s;"
                                " assuming horizontal dispersion", refname)
                    disp = HORIZONTAL
                if disp != HORIZONTAL and disp != VERTICAL:
                    log.error("dispaxis = %d is not valid.", disp)
                    raise ValueError('dispaxis must be 1 or 2.')
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
                extract_params['nod_correction'] = get_nod_offset(aper, meta)
            break

    return extract_params


def get_nod_offset(aper, meta):
    points = meta.dither.total_points
    position = meta.dither.position_number
    nod_correction = 0.
    if points is not None:
        # Find the nod offset to the extraction slit location.
        nod_offset = "nod{0}_offset".format(points)
        if position == 2:
            position = -1
        if nod_offset in aper:
            nod_correction = position * aper[nod_offset]
    return nod_correction


def log_initial_parameters(extract_params):
    """Log some of the initial extraction parameters."""

    log.debug("dispaxis = %d", extract_params["dispaxis"])
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
    log.debug("nod_correction = %s", str(extract_params["nod_correction"]))


def get_aperture(slit, extract_params):
    """Get the extraction limits xstart, xstop, ystart, ystop.

    Parameters
    ----------
    slit: data model
        The input slit, or image.

    extract_params: dictionary
        Parameters read from the reference file.

    Returns
    -------
    ap_ref: namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    ap_ref = aperture_from_ref(extract_params, slit.data.shape)

    (ap_ref, truncated) = update_from_shape(ap_ref, slit.data.shape)
    if truncated:
        log.warning("Extraction limits extended outside the input image "
                    "borders; limits have been truncated.")

    if hasattr(slit.meta, 'wcs'):
        ap_wcs = aperture_from_wcs(slit.meta.wcs)
    else:
        ap_wcs = None

    # If the xstart, etc., values were not specified for the dispersion
    # direction, the extraction region should be centered within the
    # WCS bounding box (domain).
    ap_ref = update_from_wcs(ap_ref, ap_wcs)
    ap_ref = update_from_width(ap_ref, extract_params["extract_width"],
                               extract_params["dispaxis"])

    ap_ref = apply_nod_offset(ap_ref, extract_params["nod_correction"],
                              extract_params["dispaxis"])

    # Do this again, in case the nod offset correction was so large that
    # the extraction region would extend outside the WCS bounding box.
    ap_ref = update_from_wcs(ap_ref, ap_wcs)

    return ap_ref


def aperture_from_ref(extract_params, im_shape):
    """Get extraction region from reference file or image shape.

    Parameters
    ----------
    extract_params: dictionary
        Parameters read from the reference file.

    im_shape: tuple of int
        The last two elements are the height and width of the input image
        (slit).

    Returns
    -------
    ap_ref: namedtuple
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
    ap_ref: namedtuple
        Contains xstart, xstop, ystart, ystop.  These are the initial
        values as read from the reference file, except that they may
        have been truncated at the image borders.

    extract_width: int
        The number of pixels in the cross-dispersion direction to add
        together to make a 1-D spectrum from a 2-D image.

    direction: int
        HORIZONTAL (1) if the dispersion direction is predominantly
        horizontal.  VERTICAL (2) if the dispersion direction is
        predominantly vertical.

    Returns
    -------
    ap_width: namedtuple
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


def apply_nod_offset(ap, nod_correction, direction):
    """Add the nod offset to the aperture location.

    This function only applies the nod offset correction to the limits
    in the cross-dispersion direction that were specified with
    ystart and ystop (or xstart and xstop, depending on the dispersion
    direction).

    If the source and/or background regions were instead specified with
    src_coeff and bkg_coeff, those will need to be corrected separately.
    """

    if direction == HORIZONTAL:
        ap_corr = Aperture(xstart=ap.xstart,
                           xstop=ap.xstop,
                           ystart=ap.ystart + nod_correction,
                           ystop=ap.ystop + nod_correction)
    else:
        ap_corr = Aperture(xstart=ap.xstart + nod_correction,
                           xstop=ap.xstop + nod_correction,
                           ystart=ap.ystart,
                           ystop=ap.ystop)

    return ap_corr


def update_from_shape(ap, im_shape):
    """Truncate extraction region based on input image shape.

    Parameters
    ----------
    ap: namedtuple
        Extraction region.

    im_shape: tuple of int
        The last two elements are the height and width of the input image.

    Returns
    -------
    tuple: (ap_shape, truncated)
        ap_shape is a namedtuple with keys 'xstart', 'xstop', 'ystart',
        and 'ystop'.
        `truncated` is a boolean, True if any value was truncated at an
        image edge.
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


def aperture_from_wcs(wcs):
    """Get the limits over which the WCS is defined.

    Parameters
    ----------
    wcs: data model
        The world coordinate system interface.

    Returns
    -------
    tuple ap_wcs
        namedtuple or None
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.  These are the
        limits copied directly from wcs.bounding_box.
    """

    got_bounding_box = False
    try:
        bounding_box = wcs.bounding_box
        got_bounding_box = True
        log.debug("Using wcs.bounding_box.")
    except AttributeError:
        log.info("wcs.bounding_box not found; using wcs.domain instead.")
        bounding_box = ((wcs.domain[0]['lower'], wcs.domain[0]['upper']),
                        (wcs.domain[1]['lower'], wcs.domain[1]['upper']))

    if got_bounding_box and bounding_box is None:
        log.warning("wcs.bounding_box is None")
        return None

    # bounding_box should be a tuple of tuples, each of the latter
    # consisting of (lower, upper) limits.
    if len(bounding_box) < 2:
        log.warning("wcs.bounding_box has the wrong shape")
        return None

    # These limits are float, and they are inclusive.
    xstart = bounding_box[0][0]
    xstop = bounding_box[0][1]
    ystart = bounding_box[1][0]
    ystop = bounding_box[1][1]

    ap_wcs = Aperture(xstart=xstart, xstop=xstop, ystart=ystart, ystop=ystop)
    return ap_wcs


def update_from_wcs(ap_ref, ap_wcs):
    """Limit the extraction region to the WCS bounding box.

    Parameters
    ----------
    ap_ref: namedtuple
        Contains xstart, xstop, ystart, ystop.  These are the values of
        the extraction region as specified by the reference file or the
        image size.  The cross-dispersion limits have been shifted by the
        nod offset correction.

    ap_wcs: namedtuple
        These are the bounding box limits.

    Returns
    -------
    ap: namedtuple
        Keys are 'xstart', 'xstop', 'ystart', and 'ystop'.
    """

    if ap_wcs is None:
        return ap_ref

    # ap_wcs has the limits over which the WCS transformation is defined;
    # take those as the outer limits over which we will extract.
    xstart = compare_start(ap_ref.xstart, ap_wcs.xstart)
    ystart = compare_start(ap_ref.ystart, ap_wcs.ystart)
    xstop = compare_stop(ap_ref.xstop, ap_wcs.xstop)
    ystop = compare_stop(ap_ref.ystop, ap_wcs.ystop)

    ap = Aperture(xstart=xstart, xstop=xstop, ystart=ystart, ystop=ystop)

    return ap


def compare_start(start_ref, start_wcs):
    """Compare the start limit from the aperture with the WCS lower limit.

    The more restrictive (i.e. larger) limit is the one upon which the
    output value will be based.  If the WCS limit is larger, the value will
    be increased to an integer, on the assumption that WCS lower limits
    correspond to the lower edge of the bounding pixel.  If this value
    will actually be used for an extraction limit (i.e. if the limits were
    not already specified by polynomial coefficients), then
    assign_polynomial_limits will create a polynomial function using this
    value, except that it will be decreased by 0.5 to correspond to the
    lower edge of the bounding pixels.

    Parameters
    ----------
    start_ref: int or float
        xstart or ystart, as specified by the reference file or the image
        size, possibly shifted by the nod offset correction.

    start_wcs: int or float
        The lower limit from the WCS bounding box.

    Returns
    -------
    value: int or float
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
    assign_polynomial_limits will create a polynomial function using this
    value, except that it will be increased by 0.5 to correspond to the
    upper edge of the bounding pixels.

    Parameters
    ----------
    stop_ref: int or float
        xstop or ystop, as specified by the reference file or the image
        size, possibly shifted by the nod offset correction.

    stop_wcs: int or float
        The upper limit from the WCS bounding box.

    Returns
    -------
    value: int or float
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
    coeff: list of float
        The coefficients of the polynomial, constant term first, highest
        order term last.

    Returns
    -------
    astropy.modeling.polynomial.Polynomial1D object, or None if `coeff`
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


class ExtractModel(object):

    def __init__(self, input_model,
                 dispaxis=HORIZONTAL,
                 xstart=None, xstop=None, ystart=None, ystop=None,
                 extract_width=None, src_coeff=None, bkg_coeff=None,
                 independent_var="pixel",
                 smoothing_length=0, bkg_order=0, nod_correction=0.,
                 x_center=None, y_center=None,
                 inner_bkg=None, outer_bkg=None, method='subpixel'):

        self.dispaxis = dispaxis

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
        if hasattr(input_model, 'meta') and hasattr(input_model.meta, 'wcs'):
            self.wcs = input_model.meta.wcs
        else:
            self.wcs = None
        self._wave_model = None

        # If source extraction coefficients src_coeff were specified, this
        # method will add the nod offset correction to the first coefficient
        # of every coefficient list; otherwise, the nod offset will be added
        # to xstart & xstop or ystart & ystop in get_aperture.
        # If background extraction coefficients bkg_coeff were specified,
        # this method will add the nod offset to the first coefficients.
        # Note that background coefficients are handled independently of
        # src_coeff.
        self.add_nod_correction()


    def add_nod_correction(self):
        """Add the nod offset to src_coeff and bkg_coeff (in-place)."""

        if self.nod_correction == 0.:
            return

        if self.src_coeff is not None or self.bkg_coeff is not None:
            log.info("Applying nod offset of {0}".format(self.nod_correction))

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

        Copy the values of xstart, etc., to the attributes.  Note, however,
        that if src_coeff was specified, that will override the values
        given by xstart, etc.
        The limits in the dispersion direction will be rounded to integer.
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


    def assign_polynomial_limits(self):
        """Create polynomial functions for extraction limits.

        self.src_coeff and self.bkg_coeff contain lists of polynomial
        coefficients.  These will be used to create corresponding lists of
        polynomial functions, self.p_src and self.p_bkg.  Note, however,
        that the structures of those two lists are not the same.
        The coefficients lists have this form:
            [[1, 2], [3, 4, 5], [6], [7, 8]]
        which means:
            [1, 2] coefficients for the lower limit of the first region
            [3, 4, 5] coefficients for the upper limit of the first region
            [6] coefficient(s) for the lower limit of the second region
            [7, 8] coefficients for the upper limit of the second region
        The lists of coefficients must always be in pairs, for the lower
        and upper limits respectively, but they're not explicitly in an
        additional layer of two-element lists,
            i.e., not like this:  [[[1, 2], [3, 4, 5]], [[6], [7, 8]]]
        That seemed unnecessarily messy and harder for the user to specify.
        For the lists of polynomial functions, on the other hand, that
        additional layer of list is used:
            [[fcn_lower1, fcn_upper1], [fcn_lower2, fcn_upper2]]
        where:
            fcn_lower1 is 1 + 2 * x
            fcn_upper1 is 3 + 4 * x + 5 * x**2
            fcn_lower2 is 6
            fcn_upper2 is 7 + 8 * x
        """

        if self.src_coeff is None:
            # Create constant functions.

            if self.dispaxis == HORIZONTAL:
                lower = float(self.ystart) - 0.5
                upper = float(self.ystop) + 0.5
            else:
                lower = float(self.xstart) - 0.5
                upper = float(self.xstop) + 0.5
            log.debug("Converting extraction limits to [[%g], [%g]]",
                      lower, upper)
            self.p_src = [[create_poly([lower]), create_poly([upper])]]
        else:
            # The source extraction can include more than one region.
            n_lists = len(self.src_coeff)
            if n_lists // 2 * 2 != n_lists:
                raise RuntimeError("src_coeff must contain alternating lists"
                                   " of lower and upper limits.")
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
                raise RuntimeError("bkg_coeff must contain alternating lists"
                                   " of lower and upper limits.")
            self.p_bkg = []
            expect_lower = True                         # toggled in loop
            for coeff_list in self.bkg_coeff:
                if expect_lower:
                    lower = create_poly(coeff_list)
                else:
                    upper = create_poly(coeff_list)
                    self.p_bkg.append([lower, upper])
                expect_lower = not expect_lower


    def extract(self, data):
        """
        Do the actual extraction.

        Parameters
        ----------
        data: array_like (2-D)
            Data array for one slit of a MultiSlitModel object.

        Returns
        -------
        (wavelength, net, background)
            These are all 1-D arrays.  `wavelength` is the wavelength in
            micrometers at each pixel.  `net` is the count rate
            (counts / s) minus the background at each pixel.  `background`
            is the background count rate that was subtracted from the total
            source count rate to get `net`.
        """

        # We need integer values that are the limits of a slice in the
        # dispersion direction.
        if self.dispaxis == HORIZONTAL:
            slice0 = int(round(self.xstart))
            slice1 = int(round(self.xstop)) + 1
        else:
            slice0 = int(round(self.ystart))
            slice1 = int(round(self.ystop)) + 1

        # x_array and y_array are just used for computing the wavelengths.
        if self.dispaxis == HORIZONTAL:
            x_array = np.arange(slice0, slice1, dtype=np.float64)
            y_array = np.empty(x_array.shape, dtype=np.float64)
            y_array.fill((self.ystart + self.ystop) / 2.)
        else:
            y_array = np.arange(slice0, slice1, dtype=np.float64)
            x_array = np.empty(y_array.shape, dtype=np.float64)
            x_array.fill((self.xstart + self.xstop) / 2.)

        if self.wcs is not None:
            _, _, wavelength = self.wcs(x_array, y_array)

        elif self._wave_model is not None:
            if self.dispaxis == HORIZONTAL:
                wavelength = self._wave_model(x_array)
            else:
                wavelength = self._wave_model(y_array)

        else:
            wavelength = None
        del x_array, y_array

        # Range (slice) of pixel numbers in the dispersion direction.
        disp_range = [slice0, slice1]
        if self.dispaxis == HORIZONTAL:
            image = data
            if wavelength is None:
                wavelength = np.arange(slice0, slice1, dtype=np.float64)
        else:
            image = np.transpose(data, (1, 0))
            if wavelength is None:
                wavelength = np.arange(slice0, slice1, dtype=np.float64)

        mask = np.isnan(wavelength)
        n_nan = mask.sum(dtype=np.float64)
        if n_nan > 0:
            log.warning("%d NaNs in wavelength array; set to 0.01", n_nan)
            wavelength[mask] = 0.01         # workaround
        del mask

        # src total flux, area, total weight
        (net, background) = \
        extract1d.extract1d(image, wavelength, disp_range,
                            self.p_src, self.p_bkg, self.independent_var,
                            self.smoothing_length, self.bkg_order,
                            weights=None)

        return (wavelength, net, background)


    def __del__(self):
        self.dispaxis = None
        self.xstart = None
        self.xstop = None
        self.ystart = None
        self.ystop = None
        self.extract_width = None
        self.independent_var = None
        self.src_coeff = None
        self.bkg_coeff = None
        self.p_src = None
        self.p_bkg = None
        self.smoothing_length = None
        self.bkg_order = None
        self.nod_correction = None
        self.wcs = None
        self._wave_model = None


def interpolate_response(wavelength, relsens):
    """Interpolate within the relative response table.

    Parameters
    ----------
    wavelength: array_like, 1-D
        Wavelengths in the science data

    relsens: record array
        Contains two columns, 'wavelength' and 'response'.

    Returns
    -------
    r_factor: array_like
        The response, interpolated at `wavelength`, with extrapolated
        elements and zero or negative response values set to 1.  Divide
        the net count rate by r_factor to obtain the flux.
    """

    # "_relsens" indicates that the values were read from the RELSENS table.
    wl_relsens = relsens['wavelength']
    resp_relsens = relsens['response']
    MICRONS_100 = 1.e-4                 # 100 microns, in meters
    if wl_relsens.max() > 0. and wl_relsens.max() < MICRONS_100:
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

    # `r_factor` is the response, interpolated at the wavelengths in the
    # science data.  -2048 is a flag value, to check for extrapolation.
    r_factor = np.interp(wavelength, wl_relsens, resp_relsens, -2048., -2048.)
    mask = np.where(r_factor == -2048.)
    if len(mask[0]) > 0:
        log.warning("Using RELSENS, %d elements were extrapolated; "
                    "these values will be set to 1.", len(mask[0]))
        r_factor[mask] = 1.
    mask = np.where(r_factor <= 0.)
    if len(mask[0]) > 0:
        log.warning("Using RELSENS, %d interpolated response values "
                    "were <= 0; these values will be set to 1.", len(mask[0]))
        r_factor[mask] = 1.

    return r_factor


def do_extract1d(input_model, refname, smoothing_length, bkg_order):

    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)

    if isinstance(input_model, datamodels.MultiSlitModel) or \
       isinstance(input_model, datamodels.MultiProductModel):

        if isinstance(input_model, datamodels.MultiSlitModel):
            slits = input_model.slits
        else:                           # MultiProductModel
            slits = input_model.products

        # Loop over the slits in the input model
        for slit in slits:
            log.debug('slit name = %s' % slit.name)
            extract_params = get_extract_parameters(refname, slit.name,
                                input_model.meta, smoothing_length, bkg_order)
            wavelength, net, background = \
                extract_one_slit(slit, -1,
                                 input_model.meta,
                                 slit.name, **extract_params)
            got_relsens = True
            try:
                relsens = slit.relsens
            except AttributeError:
                got_relsens = False
            if got_relsens and len(relsens) == 0:
                got_relsens = False
            if got_relsens:
                r_factor = interpolate_response(wavelength, relsens)
                flux = net / r_factor
            else:
                log.warning("No relsens for current slit, "
                            "so can't compute flux.")
                flux = np.zeros_like(net)
            dq = np.zeros(net.shape, dtype=np.int32)
            fl_error = np.ones_like(net)
            nerror = np.ones_like(net)
            berror = np.ones_like(net)
            spec = datamodels.SpecModel()
            otab = np.array(list(zip(wavelength, flux, fl_error, dq,
                                 net, nerror, background, berror)),
                            dtype=spec.spec_table.dtype)
            spec = datamodels.SpecModel(spec_table=otab)
            output_model.spec.append(spec)
    else:
        slitname = input_model.meta.exposure.type
        if slitname is None:
            slitname = ANY
        if slitname == 'NIS_SOSS':
            slitname = input_model.meta.subarray.name
        log.debug('slitname=%s' % slitname)

        if isinstance(input_model, datamodels.ImageModel) or \
           isinstance(input_model, datamodels.DrizProductModel):
            extract_params = get_extract_parameters(refname, slitname,
                                input_model.meta, smoothing_length, bkg_order)
            if extract_params:
                wavelength, net, background = \
                        extract_one_slit(input_model, -1,
                                         input_model.meta,
                                         slitname, **extract_params)
            else:
                log.critical('Missing extraction parameters.')
                raise ValueError('Missing extraction parameters.')
            dq = np.zeros(net.shape, dtype=np.int32)
            got_relsens = True
            try:
                relsens = input_model.relsens
            except AttributeError:
                got_relsens = False
            if got_relsens and len(relsens) == 0:
                got_relsens = False
            if got_relsens:
                r_factor = interpolate_response(wavelength, relsens)
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
            output_model.spec.append(spec)

        elif isinstance(input_model, datamodels.CubeModel):

            extract_params = get_extract_parameters(refname, slitname,
                                input_model.meta, smoothing_length, bkg_order)
            if not extract_params:
                log.critical('Missing extraction parameters.')
                raise ValueError('Missing extraction parameters.')

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
            for integ in range(input_model.data.shape[0]):
                # Extract spectrum
                wavelength, net, background = \
                        extract_one_slit(input_model, integ,
                                         input_model.meta,
                                         slitname, **extract_params)
                dq = np.zeros(net.shape, dtype=np.int32)
                if got_relsens:
                    r_factor = interpolate_response(wavelength,
                                                    input_model.relsens)
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
                output_model.spec.append(spec)

        elif isinstance(input_model, datamodels.IFUCubeModel):

            try:
                source_type = input_model.meta.target.source_type.lower()
            except AttributeError:
                source_type = "unknown"
            output_model = ifu.ifu_extract1d(input_model, refname, source_type)

        else:
            log.error("The input file is not supported for this step.")
            raise RuntimeError("Can't extract a spectrum from this file.")

    return output_model


def extract_one_slit(slit, integ, meta, slitname=None, **extract_params):

    log_initial_parameters(extract_params)
    ap = get_aperture(slit, extract_params)

    extract_model = ExtractModel(slit, **extract_params)
    extract_model.update_extraction_limits(ap)
    extract_model.log_extraction_parameters()

    extract_model.assign_polynomial_limits()
    data = slit.data
    if integ > -1:
        data = slit.data[integ]
    wavelength, net, background = extract_model.extract(data)
    del extract_model

    return wavelength, net, background
