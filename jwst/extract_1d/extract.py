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
from gwcs import selector
from . import extract1d

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# For full-frame input data, keyword SLTNAME may not be populated, so use
# the following string to indicate that the first slit in the reference
# file should be used.
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
        if aper.has_key('id') and (aper['id'] == slitname or
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
                # These can be used later (for MultiSlitModel data), and
                # if they are, they will be one-indexed values.
                extract_params['slit_start1'] = None
                extract_params['slit_start2'] = None
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


def apply_nod_offset(aperture, nod_correction, dispaxis):
    """Add the nod offset (if non-zero) to the aperture location.

    Note that if the source and background regions were specified with
    the src_coeff and bkg_coeff keys, the nod correction also needs to be
    added to those locations.  This will be done by calling method
    add_nod_correction in the __init__ for ExtractModel.
    """

    if aperture is None:
        return aperture

    if nod_correction != 0.:
        log.debug("Applying nod offset of {0}".format(nod_correction))
        if dispaxis == HORIZONTAL:
            aperture = Aperture(aperture.xstart,
                                aperture.ystart + nod_correction,
                                aperture.xstop,
                                aperture.ystop + nod_correction)
        else:
            aperture = Aperture(aperture.xstart + nod_correction,
                                aperture.ystart,
                                aperture.xstop + nod_correction,
                                aperture.ystop)
    return aperture


def get_aperture(slit, meta, extract_params):

    direction = extract_params['dispaxis']

    # Copy out the values of xstart, ystart, etc., that were previously
    # read from the reference file.
    ap_ref = aperture_from_ref(extract_params)

    if hasattr(slit.meta, 'wcs'):
        ap_wcs = aperture_from_wcs(slit.meta.wcs, direction)
    else:
        ap_wcs = None

    ap_shape = Aperture(0, 0, slit.data.shape[-1], slit.data.shape[-2])

    ap = reconcile_ap_limits(ap_ref, ap_wcs, ap_shape)

    ap = apply_nod_offset(ap, extract_params['nod_correction'], direction)

    return ap


def aperture_from_ref(extract_params):

    xstart = extract_params.get('xstart', None)
    xstop = extract_params.get('xstop', None)
    ystart = extract_params.get('ystart', None)
    ystop = extract_params.get('ystop', None)

    if xstart is None and ystart is None and xstop is None and ystop is None:
        return None
    else:
        # Individual values can be None, just not all of them.
        ap_ref = Aperture(xstart, ystart, xstop, ystop)
        return ap_ref


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

    xstart = domain[0]['lower']
    xstop = domain[0]['upper']
    ystart = domain[1]['lower']
    ystop = domain[1]['upper']
    # Convert limits in the dispersion direction to integer values.
    if direction == HORIZONTAL:
        x_test = round(xstart)
        if 'includes_lower' in domain[0] and domain[0]['includes_lower']:
            if x_test < xstart:
                x_test += 1.
        else:
            if x_test <= xstart:
                x_test += 1.
        xstart = int(x_test)
        x_test = round(xstop)
        if 'includes_upper' in domain[0] and domain[0]['includes_upper']:
            if x_test > xstop:
                x_test -= 1.
        else:
            if x_test >= xstop:
                x_test -= 1.
        xstop = int(x_test)
    else:                           # dispersion is vertical
        y_test = round(ystart)
        if 'includes_lower' in domain[1] and domain[1]['includes_lower']:
            if y_test < ystart:
                y_test += 1.
        else:
            if y_test <= ystart:
                y_test += 1.
        ystart = int(y_test)
        y_test = round(ystop)
        if 'includes_upper' in domain[1] and domain[1]['includes_upper']:
            if y_test > ystop:
                y_test -= 1.
        else:
            if y_test >= ystop:
                y_test -= 1.
        ystop = int(y_test)

    ap_wcs = Aperture(xstart, ystart, xstop, ystop)
    return ap_wcs

def reconcile_ap_limits(ap_ref, ap_wcs, ap_shape):

    ref_ok = (ap_ref is not None)       # But individual elements can be None
    wcs_ok = (ap_wcs is not None)       # May be reset below

    # These values are always defined, based on the input image shape.
    xstart = ap_shape.xstart
    xstop = ap_shape.xstop
    ystart = ap_shape.ystart
    ystop = ap_shape.ystop

    # If ap_ref is populated, use it in preference to the image size, but
    # truncate at image borders, and later compare with ap_wcs as well.
    if ref_ok:
        truncated = False                       # just for info
        if ap_ref.xstart is not None:
            if ap_ref.xstart < xstart:
                truncated = True
            xstart = max(xstart, ap_ref.xstart)
        if ap_ref.xstop is not None:
            if ap_ref.xstop > xstop:
                truncated = True
            xstop = min(xstop, ap_ref.xstop)
        if ap_ref.ystart is not None:
            if ap_ref.ystart < ystart:
                truncated = True
            ystart = max(ystart, ap_ref.ystart)
        if ap_ref.ystop is not None:
            if ap_ref.ystop > ystop:
                truncated = True
            ystop = min(ystop, ap_ref.ystop)
        if truncated:
            log.info("Aperture limit(s) in reference file extended beyond"
                     " image size.")

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
        truncated = False                       # just for info
        # ap_wcs has the limits over which the WCS transformation is
        # defined; take those as the outer limits over which we will extract.
        if wcs_xstart is not None and ap_wcs.xstart > xstart:
            xstart = ap_wcs.xstart
            truncated = True
        if wcs_xstop is not None and ap_wcs.xstop < xstop:
            xstop = ap_wcs.xstop
            truncated = True
        if wcs_ystart is not None and ap_wcs.ystart > ystart:
            ystart = ap_wcs.ystart
            truncated = True
        if wcs_ystop is not None and ap_wcs.ystop < ystop:
            ystop = ap_wcs.ystop
            truncated = True
        if truncated:
            log.info("Aperture limit(s) truncated due to WCS domain")

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
                 slit_start1=None, slit_start2=None):

        self.dispaxis = dispaxis
        # possibly override with src_coeff
        if xstart is None:
            self.xstart = None
        else:
            self.xstart = int(round(xstart))
        if xstop is None:
            self.xstop = None
        else:
            self.xstop = int(round(xstop))
        if ystart is None:
            self.ystart = None
        else:
            self.ystart = int(round(ystart))
        if ystop is None:
            self.ystop = None
        else:
            self.ystop = int(round(ystop))

        self.slit_start1 = slit_start1
        self.slit_start2 = slit_start2

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

        # If src_coeff and/or bkg_coeff have been specified, add the
        # nod offset correction to the first coefficient of every
        # coefficient list.
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
            if self.extract_width is None:
                if self.dispaxis == HORIZONTAL:
                    width = round(self.ystop - self.ystart)
                else:
                    width = round(self.xstop - self.xstart)
            else:
                width = float(self.extract_width)
            # If extract_width was specified, that value should override
            # ystop - ystart (or xstop - xstart), in case of disagreement.
            if self.dispaxis == HORIZONTAL:
                ystart = float(self.ystart)
                ystop = float(self.ystop - 1)           # inclusive limit
                lower = (ystart + ystop) / 2. - width / 2.
                upper = lower + width
            else:
                xstart = float(self.xstart)
                xstop = float(self.xstop - 1)
                lower = (xstart + xstop) / 2. - width / 2.
                upper = lower + width
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
        (column, wavelength, background, countrate)
            These are all 1-D arrays.  `column` is the column (or row)
            number in the input data for each output pixel.  `wavelength`
            is the wavelength in angstroms at each pixel.  `background`
            is the background count rate that was subtracted from the
            total source count rate to get `countrate`.  `countrate` is
            the count rate (counts / s) at each pixel.
        """

        log.debug('xstart=%g, xstop=%g, ystart=%g, ystop=%g' %
                  (self.xstart, self.xstop, self.ystart, self.ystop))
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
            if self.slit_start1 is not None:
                x_array += (self.slit_start1 - 1.)
            if self.slit_start2 is not None:
                y_array += (self.slit_start2 - 1.)
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
            column = np.arange(self.xstart, self.xstop, dtype=np.float64)
        else:
            image = np.transpose(data, (1, 0))
            disp_range = [self.ystart, self.ystop]
            column = np.arange(self.ystart, self.ystop, dtype=np.float64)
        if wavelength is None:
            wavelength = column

        mask = np.isnan(wavelength)
        n_nan = mask.sum(dtype=np.float64)
        if n_nan > 0:
            log.warning("%d NaNs in wavelength array; set to 0.01", n_nan)
            wavelength[mask] = 0.01         # workaround
        del mask

        # src total flux, area, total weight
        (countrate, background) = \
        extract1d.extract1d(image, wavelength, disp_range,
                            self.p_src, self.p_bkg, self.independent_var,
                            self.smoothing_length, self.bkg_order,
                            weights=None)

        return (column, wavelength, background, countrate)

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


def do_extract1d(input_model, refname, smoothing_length, bkg_order):

    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)

    if isinstance(input_model, datamodels.MultiSlitModel):

        # Loop over the slits in the input model
        for slit in input_model.slits:
            extract_params = get_extract_parameters(refname, slit.name,
                                input_model.meta, smoothing_length, bkg_order)
            extract_params['slit_start1'] = slit.xstart         # one indexed
            extract_params['slit_start2'] = slit.ystart         # one indexed
            column, wavelength, background, countrate = \
                extract_one_slit(slit, -1,
                                 input_model.meta, refname,
                                 slit.name, **extract_params)
            spec = datamodels.SpecModel()
            otab = np.array(zip(column, wavelength, background, countrate),
                            dtype=spec.spec_table.dtype)
            spec = datamodels.SpecModel(spec_table=otab)
            output_model.spec.append(spec)
    else:
        slitname = input_model.meta.exposure.type
        if slitname == 'NIS_SOSS':
            slitname = input_model.meta.subarray.name
        log.debug('slitname=%s' % slitname)

        if isinstance(input_model, datamodels.ImageModel):
            extract_params = get_extract_parameters(refname, slitname,
                                input_model.meta, smoothing_length, bkg_order)
            if extract_params:
                column, wavelength, background, countrate = \
                        extract_one_slit(input_model, -1,
                                         input_model.meta, refname,
                                         slitname, **extract_params)
            else:
                log.critical('Missing extraction parameters.')
                raise ValueError('Missing extraction parameters.')
            spec = datamodels.SpecModel()
            otab = np.array(zip(column, wavelength, background, countrate),
                            dtype=spec.spec_table.dtype)
            spec = datamodels.SpecModel(spec_table=otab)
            output_model.spec.append(spec)

        elif isinstance(input_model, datamodels.CubeModel):

            extract_params = get_extract_parameters(refname, slitname,
                                input_model.meta, smoothing_length, bkg_order)
            if not extract_params:
                log.critical('Missing extraction parameters.')
                raise ValueError('Missing extraction parameters.')

            # Loop over each integration in the input model
            for integ in range(input_model.data.shape[0]):
                # Extract spectrum
                column, wavelength, background, countrate = \
                        extract_one_slit(input_model, integ,
                                         input_model.meta, refname,
                                         slitname, **extract_params)
                spec = datamodels.SpecModel()
                otab = np.array(zip(column, wavelength, background, countrate),
                                dtype=spec.spec_table.dtype)
                spec = datamodels.SpecModel(spec_table=otab)
                output_model.spec.append(spec)

    return output_model

def extract_one_slit(slit, integ, meta, refname, slitname=None,
                     **extract_params):

    extract_model = ExtractModel(slit, **extract_params)
    ap = get_aperture(slit, meta, extract_params)
    extract_model.update_extraction_limits(ap)
    extract_model.assign_polynomial_limits()
    data = slit.data
    if integ > -1:
        data = slit.data[integ]
    column, wavelength, background, countrate = extract_model.extract(data)
    del extract_model

    return column, wavelength, background, countrate
