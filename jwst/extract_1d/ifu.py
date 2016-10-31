import json
import logging
import math

import numpy as np
from photutils import CircularAperture, CircularAnnulus, aperture_photometry

from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def ifu_extract1d(input_model, refname, source_type, smoothing_length):
    """Extract a 1-D spectrum from an IFU cube.

    Parameters
    ----------
    input_model: IFUCubeModel
        The input model.

    refname: string
        The name of the JSON reference file.

    source_type: string
        "point" or "extended"

    smoothing_length: integer
        Number of pixels in the wavelength direction for boxcar smoothing
        of the background.  This is currently not implemented.

    Returns
    -------
    output_model: MultiSpecModel
        This will contain the extracted spectrum.
    """

    if not isinstance(input_model, datamodels.IFUCubeModel):
        log.error("Expected an IFU cube.")
        raise RuntimeError("Expected an IFU cube.")

    source_type = "extended"            # xxx test debug temporary
    if source_type != "point" and source_type != "extended":
        log.warning("source_type was '%s', setting to 'point'.", source_type)
        source_type = "point"
    log.info("source_type = %s", source_type)

    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)

    slitname = input_model.meta.exposure.type
    if slitname is None:
        slitname = "ANY"
    log.debug('slitname=%s' % slitname)

    extract_params = ifu_extract_parameters(refname, slitname,
                                            source_type, smoothing_length)

    if extract_params:
        wavelength, net, background = extract_ifu(input_model,
                                slitname, source_type, extract_params)
    else:
        log.critical('Missing extraction parameters.')
        raise ValueError('Missing extraction parameters.')

    dq = np.zeros(net.shape, dtype=np.int32)

    do_fluxcorr = True                          # default, and initial value
    try:
        relsens = input_model.relsens
    except AttributeError:
        do_fluxcorr = False
        log.warning("No relsens for input file.")

    # If there is no relsens table, copy net to flux, under the assumption
    # that the IFU cube is already in flux units.
    if do_fluxcorr:
        r_factor = ifu_interpolate_response(wavelength, relsens)
        flux = net / r_factor
    else:
        flux = net.copy()

    fl_error = np.ones_like(net)
    nerror = np.ones_like(net)
    berror = np.ones_like(net)
    spec = datamodels.SpecModel()
    otab = np.array(zip(wavelength, flux, fl_error, dq,
                        net, nerror, background, berror),
                    dtype=spec.spec_table.dtype)
    spec = datamodels.SpecModel(spec_table=otab)
    output_model.spec.append(spec)

    return output_model


def ifu_extract_parameters(refname, slitname,
                           source_type, smoothing_length=None):

    extract_params = {}
    with open(refname) as f:
        ref = json.load(f)
    for aper in ref['apertures']:
        if aper.has_key('id') and aper['id'] != "dummy" and \
           (aper['id'] == slitname or slitname == "ANY"):
            region_type = aper.get("region_type", "target")
            if region_type == "target":
                if smoothing_length is None:
                    extract_params['smoothing_length'] = \
                          aper.get('smoothing_length', 0)
                else:
                    # If the user supplied a value, use that value.
                    extract_params['smoothing_length'] = smoothing_length

                extract_params['x_center'] = aper.get('x_center')
                extract_params['y_center'] = aper.get('y_center')
                extract_params['extract_width'] = aper.get('extract_width')
                if source_type == "point":
                    extract_params['inner_bkg'] = aper.get('inner_bkg')
                    extract_params['outer_bkg'] = aper.get('outer_bkg')
                    extract_params['method'] = aper.get('method', 'subpixel')
            break

    missing = False
    if extract_params["extract_width"] is None:
        log.error("%s is missing from the reference file.", "extract_width")
        missing = True
    if extract_params["x_center"] is None:
        log.error("%s is missing from the reference file.", "x_center")
        missing = True
    if extract_params["y_center"] is None:
        log.error("%s is missing from the reference file.", "y_center")
        missing = True
    if missing:
        raise ValueError('Missing extraction parameters.')

    return extract_params


def ifu_interpolate_response(wavelength, relsens):
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

    if np.any(np.isnan(wl_relsens)) or np.any(np.isnan(resp_relsens)):
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


def extract_ifu(input_model, slitname, source_type, extract_params):
    """This function does the extraction.

    Parameters
    ----------
    input_model:

    slitname: string

    source_type: string
        "point" or "extended"

    **extract_params:

    Returns
    -------
        (wavelength, net, background)
    """

    data = input_model.data
    shape = data.shape
    if len(shape) != 3:
        log.error("Expected a 3-D IFU cube; dimension is %d.", len(shape))
        raise RuntimeError("The IFU cube should be 3-D.")

    wavelength = np.zeros(shape[0], dtype=np.float64)
    net = np.zeros(shape[0], dtype=np.float64)
    background = np.zeros(shape[0], dtype=np.float64)

    x_center = extract_params['x_center']
    y_center = extract_params['y_center']
    if x_center < 0 or x_center >= shape[2] or \
       y_center < 0 or y_center >= shape[1]:
        log.warning("Target coordinates are outside the image.")
        return (wavelength, net, background)    # all zeros

    # xxx not used yet
    smoothing_length = extract_params['smoothing_length']

    x_center = float(x_center)
    y_center = float(y_center)
    extract_width = float(extract_params['extract_width'])
    if source_type == "point":
        inner_bkg = extract_params['inner_bkg']
        outer_bkg = extract_params['outer_bkg']
        method = extract_params['method']
    else:
        inner_bkg = None
        outer_bkg = None
        method = "irrelevant"

    if inner_bkg is None or outer_bkg is None:
        subtract_background = False
    elif inner_bkg <= 0. or outer_bkg <= 0.:
        subtract_background = False
    elif inner_bkg == outer_bkg:
        subtract_background = False
    else:
        subtract_background = True
        if inner_bkg > outer_bkg:
            temp = outer_bkg
            outer_bkg = inner_bkg
            inner_bkg = temp

    try:
        wcs = input_model.meta.wcs
        got_real_wcs = True
    except AttributeError:
        wcs = input_model.get_fits_wcs()        # FITS keywords
        got_real_wcs = False

    if got_real_wcs:
        log.debug("WCS is input_model.meta.wcs")
        x_array = np.empty(shape[0], dtype=np.float64)
        x_array.fill(float(shape[2]) / 2.)
        y_array = np.empty(shape[0], dtype=np.float64)
        y_array.fill(float(shape[1]) / 2.)
        z_array = np.arange(shape[0], dtype=np.float64) # for wavelengths
        _, _, wavelength = wcs(x_array, y_array, z_array)
    else:
        log.warning("WCS is input_model.get_fits_wcs()")
        x = shape[2] / 2.
        y = shape[1] / 2.
        for k in range(shape[0]):
            wavelength[k] = float(wcs.wcs_pix2world(x, y, float(k), 0)[2])

    if source_type == "point":
        position = (x_center, y_center)
        aperture = CircularAperture(position, r=extract_width / 2.)
        if subtract_background:
            annulus = CircularAnnulus(position,
                                      r_in=inner_bkg, r_out=outer_bkg)
            normalization = aperture.area() / annulus.area()

        for k in range(shape[0]):
            phot_table = aperture_photometry(data[k, :, :], aperture,
                                             method=method)
            net[k] = float(phot_table['aperture_sum'][0])
            if subtract_background:
                bkg_table = aperture_photometry(data[k, :, :], annulus,
                                                method=method)
                background[k] = float(bkg_table['aperture_sum'][0])
                net[k] = net[k] - background[k] * normalization
    else:
        # Floating-point coordinates with zero point at lower left corner
        # of a pixel, for convenience in handling fractions of a pixel.
        x0 = x_center + 0.5
        y0 = y_center + 0.5

        # Actual lower and upper edges of extraction box.
        x_low = x0 - extract_width / 2.
        x_high = x0 + extract_width / 2.
        y_low = y0 - extract_width / 2.
        y_high = y0 + extract_width / 2.

        # Outer limits of extraction box, integer pixels.
        xo_low = math.floor(x_low)
        xo_high = math.ceil(x_high)
        yo_low = math.floor(y_low)
        yo_high = math.ceil(y_high)
        ixl = int(xo_low)
        ixh = int(xo_high)
        iyl = int(yo_low)
        iyh = int(yo_high)
        truncated = False
        if ixl < 0:
            truncated = True
            ixl = 0
            xo_low = max(0., xo_low)
            x_low = max(0., x_low)
        if iyl < 0:
            truncated = True
            iyl = 0
            yo_low = max(0., yo_low)
            y_low = max(0., y_low)
        if ixh > shape[2]:
            truncated = True
            ixh = shape[2]
            xo_high = min(float(shape[2]), xo_high)
            x_high = min(float(shape[2]), x_high)
        if iyh > shape[1]:
            truncated = True
            iyh = shape[1]
            yo_high = min(float(shape[1]), yo_high)
            y_high = min(float(shape[1]), y_high)
        if truncated:
            log.warning("Extraction box was truncated to slice "
                        "[%d:%d, %d:%d]", iy, iyh, ixl, ixh)

        # Partial weight for column or row of edge pixels.  x_low will be
        # greater than or equal to xo_low; if they're equal, the weight
        # should be 1.  Note:  These values will only be used for the case
        # that the lower and upper limits are in separate pixels.
        xw_l = 1. - (x_low - xo_low)            # left edge
        yw_l = 1. - (y_low - yo_low)            # lower edge
        xw_h = 1. - (xo_high - x_high)          # right edge
        yw_h = 1. - (yo_high - y_high)          # top edge

        nx = int(xo_high) - int(xo_low)
        ny = int(yo_high) - int(yo_low)
        nx = max(nx, 1)
        ny = max(ny, 1)
        wgt = np.ones((1, ny, nx), dtype=np.float64)
        if nx == 1:
            wgt[0, :, 0] *= (x_high - x_low)
        else: 
            wgt[0, :, 0] *= xw_l
            wgt[0, :, -1] *= xw_h
        if ny == 1:
            wgt[0, 0, :] *= (y_high - y_low)
        else:
            wgt[0, 0, :] *= yw_l
            wgt[0, -1, :] *= yw_h

        subset = data[:, iyl:iyh, ixl:ixh] * wgt
        for k in range(shape[0]):
            net[k] = subset[k, :, :].sum()

    return (wavelength, net, background)
