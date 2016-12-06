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
                      aper.get('independent_var', 'wavelength').lower()
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


def get_aperture(slit, meta, extract_params):

    direction = extract_params['dispaxis']

    if hasattr(slit.meta, 'wcs'):
        ap_wcs = aperture_from_wcs(slit.meta.wcs, direction)
    else:
        ap_wcs = None

    ap_shape = Aperture(0, 0, slit.data.shape[-1], slit.data.shape[-2])

    ap = reconcile_ap_limits(ap_wcs, ap_shape)

    return ap


def aperture_from_wcs(wcs, direction):

    try:
        domain = wcs.domain
    except AttributeError:
        return None

    if domain is None:
        return None

    if len(domain) < 2 or 'lower' not in domain[0] or \
                          'upper' not in domain[0] or \
                          'lower' not in domain[1] or \
                          'upper' not in domain[1]:
        return None

    # These are supposed to be integers.
    xstart = domain[0]['lower']
    xstop = domain[0]['upper']
    ystart = domain[1]['lower']
    ystop = domain[1]['upper']
    # The limits should be the limits of a slice.  Possibly modify the limits,
    # depending on domain[i]['includes_lower'] (but just in the dispersion
    # direction).
    if direction == HORIZONTAL:
        x_test = int(round(xstart))     # just in case it's not an int
        if 'includes_lower' in domain[0] and not domain[0]['includes_lower']:
            # domain[0]['lower'] is not included, so increment it.
            xstart = x_test + 1
        else:
            xstart = x_test
        x_test = int(round(xstop))
        if 'includes_upper' in domain[0] and domain[0]['includes_upper']:
            # It is included, so add one to make it an upper limit of a slice.
            xstop = x_test + 1
        else:
            xstop = x_test
    else:                               # dispersion is vertical
        y_test = int(round(ystart))
        if 'includes_lower' in domain[1] and not domain[1]['includes_lower']:
            ystart = y_test + 1
        else:
            ystart = y_test
        y_test = int(round(ystop))
        if 'includes_upper' in domain[1] and domain[1]['includes_upper']:
            ystop = y_test + 1
        else:
            ystop = y_test

    ap_wcs = Aperture(xstart, ystart, xstop, ystop)
    return ap_wcs


def reconcile_ap_limits(ap_wcs, ap_shape):

    wcs_ok = (ap_wcs is not None)       # May be reset below

    # These values are always defined, based on the input image shape.
    xstart = ap_shape.xstart
    xstop = ap_shape.xstop
    ystart = ap_shape.ystart
    ystop = ap_shape.ystop

    if wcs_ok:
        # Copy these so we can assign to them.
        wcs_xstart = ap_wcs.xstart
        wcs_xstop = ap_wcs.xstop
        wcs_ystart = ap_wcs.ystart
        wcs_ystop = ap_wcs.ystop
        # The limits in ap_wcs are with respect to the original, full-size
        # image.  If the current input image is a cutout and is outside the
        # WCS domain, flag ap_wcs as not valid.
        if ap_wcs.xstart > xstop or ap_wcs.xstop < xstart:
            wcs_xstart = None
            wcs_xstop = None
        if ap_wcs.ystart > ystop or ap_wcs.ystop < ystart:
            wcs_ystart = None
            wcs_ystop = None
        if wcs_xstart is None and wcs_xstop is None and \
           wcs_ystart is None and wcs_ystop is None:
            wcs_ok = False
            log.info("Current image is outside the WCS domain")

    if wcs_ok:
        # ap_wcs has the limits over which the WCS transformation is
        # defined; take those as the outer limits over which we will extract.
        if wcs_xstart is not None and ap_wcs.xstart > xstart:
            xstart = ap_wcs.xstart
        if wcs_xstop is not None and ap_wcs.xstop < xstop:
            xstop = ap_wcs.xstop
        if wcs_ystart is not None and ap_wcs.ystart > ystart:
            ystart = ap_wcs.ystart
        if wcs_ystop is not None and ap_wcs.ystop < ystop:
            ystop = ap_wcs.ystop

    return Aperture(xstart, ystart, xstop, ystop)


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
                 independent_var="wavelength",
                 smoothing_length=0, bkg_order=0, nod_correction=0.,
                 x_center=None, y_center=None,
                 inner_bkg=None, outer_bkg=None, method='subpixel'):

        self.dispaxis = dispaxis
        # xstart, xstop, ystart, or ystop may be overridden with src_coeff.
        if xstart is None:
            self.xstart = None
        else:
            r = round(xstart)
            if xstart == r:
                self.xstart = xstart
            else:
                log.warning("xstart %s should have been an integer; "
                            "rounding to %s", str(xstart), str(r))
                self.xstart = r
        if xstop is None:
            self.xstop = None
        else:
            r = round(xstop)
            if xstop == r:
                self.xstop = xstop
            else:
                log.warning("xstop %s should have been an integer; "
                            "rounding to %s", str(xstop), str(r))
                self.xstop = r
        if ystart is None:
            self.ystart = None
        else:
            r = round(ystart)
            if ystart == r:
                self.ystart = ystart
            else:
                log.warning("ystart %s should have been an integer; "
                            "rounding to %s", str(ystart), str(r))
                self.ystart = r
        if ystop is None:
            self.ystop = None
        else:
            r = round(ystop)
            if ystop == r:
                self.ystop = ystop
            else:
                log.warning("ystop %s should have been an integer; "
                            "rounding to %s", str(ystop), str(r))
                self.ystop = r

        if extract_width is None:
            self.extract_width = None
        else:
            self.extract_width = int(round(extract_width))
        # 'wavelength' or 'pixel', the independent variable for functions
        # for lower and upper limits of source and background regions.
        self.independent_var = independent_var

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
        # to xstart & xstop or ystart & ystop in assign_polynomial_limits.
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

    def log_extraction_parameters(self):
        """Log some of the extraction parameters."""

        log.debug("dispaxis = %d", self.dispaxis)
        log.debug("independent_var = %s", self.independent_var)
        log.debug("smoothing_length = %d", self.smoothing_length)
        log.debug("provisional xstart = %s", str(self.xstart))
        log.debug("provisional xstop = %s", str(self.xstop))
        log.debug("provisional ystart = %s", str(self.ystart))
        log.debug("provisional ystop = %s", str(self.ystop))
        log.debug("extract_width = %s", str(self.extract_width))
        log.debug("bkg_order = %d", self.bkg_order)
        log.debug("nod_correction = %s", str(self.nod_correction))
        log.debug("src_coeff = %s", str(self.src_coeff))
        log.debug("bkg_coeff = %s", str(self.bkg_coeff))

        # The following values are printed in assign_polynomial_limits,
        # because they can be assigned or modified there:
        # lower & upper, and either xstart & xstop or ystart & ystop.

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

        If src_coeff was not specified, the nod offset correction
        will be added to xstart & xstop or ystart & ystop.
        """

        if self.src_coeff is None:
            # Create constant functions.
            if self.extract_width is None:
                # These limits are inclusive, and we expect integer values.
                if self.dispaxis == HORIZONTAL:
                    width = float(round(self.ystop - self.ystart + 1))
                else:
                    width = float(round(self.xstop - self.xstart + 1))
            else:
                width = float(self.extract_width)
            # If extract_width was specified, that value should override
            # ystop - ystart (or xstop - xstart), in case of disagreement.

            if self.nod_correction != 0.:
                log.info("Applying nod offset of {0}"
                         .format(self.nod_correction))
            if self.dispaxis == HORIZONTAL:
                ystart = float(self.ystart)
                ystop = float(self.ystop)
                if self.nod_correction != 0.:
                    ystart += self.nod_correction
                    ystop += self.nod_correction
                lower = (ystart + ystop) / 2. - width / 2.
                upper = lower + width
                log.debug("xstart = %s", str(self.xstart))
                log.debug("xstop = %s", str(self.xstop))
                log.debug("lower = %s", str(lower))
                log.debug("upper = %s", str(upper))
            else:
                xstart = float(self.xstart)
                xstop = float(self.xstop)
                if self.nod_correction != 0.:
                    xstart += self.nod_correction
                    xstop += self.nod_correction
                lower = (xstart + xstop) / 2. - width / 2.
                upper = lower + width
                log.debug("lower = %s", str(lower))
                log.debug("upper = %s", str(upper))
                log.debug("ystart = %s", str(self.ystart))
                log.debug("ystop = %s", str(self.ystop))
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

    def update_extraction_limits(self, ap):
        """Possibly update start and stop limits.

        If attributes self.xstart, self.xstop, self.ystart, or self.ystop
        have not been assigned values yet, update them from `ap` (a
        named tuple).  The actual lower and upper extraction limits may be
        set to different values than these by using extract_width, or by
        assign_polynomial_limits() if src_coeff was specified.

        The limits in the dispersion direction will be rounded to integer.
        """

        if self.xstart is None:
            self.xstart = ap.xstart
        if self.xstop is None:
            self.xstop = ap.xstop
        if self.ystart is None:
            self.ystart = ap.ystart
        if self.ystop is None:
            self.ystop = ap.ystop

        if self.dispaxis == HORIZONTAL:
            self.xstart = max(self.xstart, ap.xstart)
            self.xstop = min(self.xstop, ap.xstop)
            self.xstart = int(round(self.xstart))
            self.xstop = int(round(self.xstop))
        else:                           # vertical
            self.ystart = max(self.ystart, ap.ystart)
            self.ystop = min(self.ystop, ap.ystop)
            self.ystart = int(round(self.ystart))
            self.ystop = int(round(self.ystop))

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

        # x_array and y_array are just used for computing the wavelengths.
        if self.dispaxis == HORIZONTAL:
            x_array = np.arange(self.xstart, self.xstop, dtype=np.float64)
            y_array = np.empty(x_array.shape, dtype=np.float64)
            y_array.fill((self.ystart + self.ystop) / 2.)
        else:
            y_array = np.arange(self.ystart, self.ystop, dtype=np.float64)
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

        if self.dispaxis == HORIZONTAL:
            image = data
            # Range (slice) of pixel numbers in the dispersion direction.
            disp_range = [self.xstart, self.xstop]
            if wavelength is None:
                wavelength = np.arange(self.xstart, self.xstop,
                                       dtype=np.float64)
        else:
            image = np.transpose(data, (1, 0))
            disp_range = [self.ystart, self.ystop]
            if wavelength is None:
                wavelength = np.arange(self.ystart, self.ystop,
                                       dtype=np.float64)

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
                source_type = input_model.meta.target.source_type
            except AttributeError:
                source_type = "unknown"
            output_model = ifu.ifu_extract1d(input_model, refname, source_type)

        else:
            log.error("The input file is not supported for this step.")
            raise RuntimeError("Can't extract a spectrum from this file.")

    return output_model


def extract_one_slit(slit, integ, meta, slitname=None,
                     **extract_params):

    extract_model = ExtractModel(slit, **extract_params)
    ap = get_aperture(slit, meta, extract_params)
    extract_model.update_extraction_limits(ap)
    extract_model.log_extraction_parameters()
    extract_model.assign_polynomial_limits()
    data = slit.data
    if integ > -1:
        data = slit.data[integ]
    wavelength, net, background = extract_model.extract(data)
    del extract_model

    return wavelength, net, background
