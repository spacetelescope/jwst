import logging
import numpy as np
from gwcs.wcstools import grid_from_bounding_box
from scipy.interpolate import interp1d
from stdatamodels.jwst.transforms.models import IdealToV2V3

from jwst.assign_wcs.util import wcs_bbox_from_shape


__all__ = ["middle_from_wcs", "location_from_wcs", "trace_from_wcs", "nod_pair_location"]

HORIZONTAL = 1
"""Horizontal dispersion axis."""
VERTICAL = 2
"""Vertical dispersion axis."""

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def middle_from_wcs(wcs, bounding_box, dispaxis):
    """
    Calculate the effective middle of the spectral region.

    Parameters
    ----------
    wcs : `~gwcs.WCS`
        WCS for the input data model, containing detector to wavelength
        transforms.
    bounding_box : tuple
        A pair of tuples, each consisting of two numbers.
        Represents the range of useful pixel values in both dimensions,
        ((xmin, xmax), (ymin, ymax)).
    dispaxis : int
        Dispersion axis.

    Returns
    -------
    middle_disp : float
        Middle pixel in the dispersion axis.
    middle_xdisp : float
        Middle pixel in the cross-dispersion axis.
    middle_wavelength : float
        Wavelength at the middle pixel.
    """
    if dispaxis == HORIZONTAL:
        # Width (height) in the cross-dispersion direction, from the start of
        # the 2-D cutout (or of the full image) to the upper limit of the bounding box.
        xd_width = int(round(bounding_box[1][1]))  # must be an int

        # Middle of the bounding_box in the dispersion direction.
        middle_disp = (bounding_box[0][0] + bounding_box[0][1]) / 2.0
        x = np.full(xd_width, middle_disp)

        # 1-D vector of cross-dispersion (y) pixel indices
        y = np.arange(xd_width, dtype=np.float64)

    else:
        # Cross-dispersion total width of bounding box; must be an int
        xd_width = int(round(bounding_box[0][1]))

        # Middle of the bounding_box in the dispersion direction.
        middle_disp = (bounding_box[1][0] + bounding_box[1][1]) / 2.0
        y = np.full(xd_width, middle_disp)

        # 1-D vector of cross-dispersion (x) pixel indices
        x = np.arange(xd_width, dtype=np.float64)

    # Get all the wavelengths at the middle dispersion element
    _, _, center_wavelengths = wcs(x, y)
    sort_idx = np.argsort(center_wavelengths)
    valid = np.isfinite(center_wavelengths[sort_idx])

    # Average to get the middle wavelength
    middle_wavelength = np.nanmean(center_wavelengths)

    # Find the effective index in cross-dispersion coordinates for the
    # averaged wavelength to get the cross-dispersion center
    if dispaxis == HORIZONTAL:
        if np.allclose(center_wavelengths, middle_wavelength):
            middle_xdisp = np.mean(y)
        else:
            middle_xdisp = np.interp(
                middle_wavelength, center_wavelengths[sort_idx][valid], y[sort_idx[valid]]
            )
    else:
        if np.allclose(center_wavelengths, middle_wavelength):
            middle_xdisp = np.mean(x)
        else:
            middle_xdisp = np.interp(
                middle_wavelength, center_wavelengths[sort_idx][valid], x[sort_idx[valid]]
            )
    return middle_disp, middle_xdisp, middle_wavelength


def location_from_wcs(input_model, slit, make_trace=True):
    """
    Get the cross-dispersion location of the spectrum, based on the WCS.

    None values will be returned if there was insufficient information
    available, e.g. if the wavelength attribute or wcs function is not
    defined.

    Parameters
    ----------
    input_model : DataModel
        The input science model containing metadata information.
    slit : DataModel or None
        One slit from a MultiSlitModel (or similar), or None.
        The WCS and target coordinates will be retrieved from `slit`
        unless `slit` is None. In that case, they will be retrieved
        from `input_model`.
    make_trace : bool, optional
        If True, the source position will be calculated for each
        dispersion element and returned in `trace`.  If False,
        None is returned.

    Returns
    -------
    middle : int or None
        Pixel coordinate in the dispersion direction within the 2-D
        cutout (or the entire input image) at the middle of the WCS
        bounding box.  This is the point at which to determine the
        nominal extraction location, in case it varies along the
        spectrum.  The offset will then be the difference between
        `location` (below) and the nominal location.
    middle_wl : float or None
        The wavelength at pixel `middle`.
    location : float or None
        Pixel coordinate in the cross-dispersion direction within the
        spectral image that is at the planned target location.
        The spectral extraction region should be centered here.
    trace : ndarray or None
        An array of source positions, one per dispersion element, corresponding
        to the location at each point in the wavelength array. If the
        input data is resampled, the trace corresponds directly to the
        location. If the trace could not be generated, or `make_trace` is
        False, None is returned.
    """
    if slit is not None:
        shape = slit.data.shape[-2:]
        wcs = slit.meta.wcs
        dispaxis = slit.meta.wcsinfo.dispersion_direction
    else:
        shape = input_model.data.shape[-2:]
        wcs = input_model.meta.wcs
        dispaxis = input_model.meta.wcsinfo.dispersion_direction

    bb = wcs.bounding_box  # ((x0, x1), (y0, y1))
    if bb is None:
        bb = wcs_bbox_from_shape(shape)
    if dispaxis == HORIZONTAL:
        lower = bb[1][0]
        upper = bb[1][1]
    else:
        lower = bb[0][0]
        upper = bb[0][1]

    # Get the wavelengths for the valid data in the sky transform,
    # average to get the middle wavelength
    middle, _, middle_wl = middle_from_wcs(wcs, bb, dispaxis)
    middle = int(np.round(middle))

    exp_type = input_model.meta.exposure.type
    trace = None
    if exp_type in ["NRS_FIXEDSLIT", "NRS_MSASPEC", "NRS_BRIGHTOBJ"]:
        log.info("Using source_xpos and source_ypos to center extraction.")
        if slit is None:
            xpos = input_model.source_xpos
            ypos = input_model.source_ypos
        else:
            xpos = slit.source_xpos
            ypos = slit.source_ypos

        slit2det = wcs.get_transform("slit_frame", "detector")
        if "gwa" in wcs.available_frames:
            # Input is not resampled, wavelengths need to be meters
            _, location = slit2det(xpos, ypos, middle_wl * 1e-6)
        else:
            _, location = slit2det(xpos, ypos, middle_wl)

        if ~np.isnan(location) and make_trace:
            trace = _nirspec_trace_from_wcs(shape, bb, wcs, xpos, ypos)

    elif exp_type == "MIR_LRS-FIXEDSLIT":
        log.info("Using dithered_ra and dithered_dec to center extraction.")
        try:
            if slit is None:
                dithra = input_model.meta.dither.dithered_ra
                dithdec = input_model.meta.dither.dithered_dec
            else:
                dithra = slit.meta.dither.dithered_ra
                dithdec = slit.meta.dither.dithered_dec
            location, _ = wcs.backward_transform(dithra, dithdec, middle_wl)

        except (AttributeError, TypeError):
            log.warning("Dithered pointing location not found in wcsinfo.")
            return None, None, None, None

        if ~np.isnan(location) and make_trace:
            trace = _miri_trace_from_wcs(shape, bb, wcs, dithra, dithdec)
    else:
        log.warning(f"Source position cannot be found for EXP_TYPE {exp_type}")
        return None, None, None, None

    if np.isnan(location):
        log.warning("Source position could not be determined from WCS.")
        return None, None, None, None

    # If the target is at the edge of the image or at the edge of the
    # non-NaN area, we can't use the WCS to find the
    # location of the target spectrum.
    if location < lower or location > upper:
        log.warning(
            f"WCS implies the target is at {location:.2f}, which is outside the bounding box,"
        )
        log.warning("so we can't get spectrum location using the WCS")
        return None, None, None, None

    return middle, middle_wl, location, trace


def _nirspec_trace_from_wcs(shape, bounding_box, wcs_ref, source_xpos, source_ypos):
    """
    Calculate NIRSpec source trace from WCS.

    The source trace is calculated by projecting the recorded source
    positions source_xpos/ypos from the NIRSpec "slit_frame" onto
    detector pixels.

    Parameters
    ----------
    shape : tuple of int
        2D shape for the full input data array, (ny, nx).
    bounding_box : tuple
        A pair of tuples, each consisting of two numbers.
        Represents the range of useful pixel values in both dimensions,
        ((xmin, xmax), (ymin, ymax)).
    wcs_ref : `~gwcs.WCS`
        WCS for the input data model, containing slit and detector
        transforms.
    source_xpos : float
        Slit position, in the x direction, for the target.
    source_ypos : float
        Slit position, in the y direction, for the target.

    Returns
    -------
    trace : ndarray of float
        Fractional pixel positions in the y (cross-dispersion direction)
        of the trace for each x (dispersion direction) pixel.
    """
    x, y = grid_from_bounding_box(bounding_box)
    nx = int(bounding_box[0][1] - bounding_box[0][0])

    # Calculate the wavelengths in the slit frame because they are in
    # meters for cal files and um for s2d files
    d2s = wcs_ref.get_transform("detector", "slit_frame")
    _, _, slit_wavelength = d2s(x, y)

    # Make an initial array of wavelengths that will cover the wavelength range of the data
    wave_vals = np.linspace(np.nanmin(slit_wavelength), np.nanmax(slit_wavelength), nx)
    # Get arrays of the source position in the slit
    pos_x = np.full(nx, source_xpos)
    pos_y = np.full(nx, source_ypos)

    # Grab the wcs transform between the slit frame where we know the
    # source position and the detector frame
    s2d = wcs_ref.get_transform("slit_frame", "detector")

    # Calculate the expected center of the source trace
    trace_x, trace_y = s2d(pos_x, pos_y, wave_vals)

    # Interpolate the trace to a regular pixel grid in the dispersion
    # direction
    interp_trace = interp1d(trace_x, trace_y, fill_value="extrapolate")

    # Get the trace position for each dispersion element
    trace = interp_trace(np.arange(nx))

    # Place the trace in the full array
    full_trace = np.full(shape[1], np.nan)
    x0 = int(np.ceil(bounding_box[0][0]))
    full_trace[x0 : x0 + nx] = trace

    return full_trace


def _miri_trace_from_wcs(shape, bounding_box, wcs_ref, source_ra, source_dec):
    """
    Calculate MIRI LRS fixed slit source trace from WCS.

    The source trace is calculated by projecting the recorded source
    positions dithered_ra/dec from the world frame onto detector pixels.

    Parameters
    ----------
    shape : tuple of int
        2D shape for the full input data array, (ny, nx).
    bounding_box : tuple
        A pair of tuples, each consisting of two numbers.
        Represents the range of useful pixel values in both dimensions,
        ((xmin, xmax), (ymin, ymax)).
    wcs_ref : `~gwcs.WCS`
        WCS for the input data model, containing sky and detector
        transforms, forward and backward.
    source_ra : float
        RA coordinate for the target.
    source_dec : float
        Dec coordinate for the target.

    Returns
    -------
    trace : ndarray of float
        Fractional pixel positions in the x (cross-dispersion direction)
        of the trace for each y (dispersion direction) pixel.
    """
    x, y = grid_from_bounding_box(bounding_box)
    ny = int(bounding_box[1][1] - bounding_box[1][0])

    # Calculate the wavelengths for the full array
    _, _, slit_wavelength = wcs_ref(x, y)

    # Make an initial array of wavelengths that will cover the wavelength range of the data
    wave_vals = np.linspace(np.nanmin(slit_wavelength), np.nanmax(slit_wavelength), ny)

    # Get arrays of the source position
    pos_ra = np.full(ny, source_ra)
    pos_dec = np.full(ny, source_dec)

    # Calculate the expected center of the source trace
    trace_x, trace_y = wcs_ref.backward_transform(pos_ra, pos_dec, wave_vals)

    # Interpolate the trace to a regular pixel grid in the dispersion
    # direction
    interp_trace = interp1d(trace_y, trace_x, fill_value="extrapolate")

    # Get the trace position for each dispersion element within the bounding box
    trace = interp_trace(np.arange(ny))

    # Place the trace in the full array
    full_trace = np.full(shape[0], np.nan)
    y0 = int(np.ceil(bounding_box[1][0]))
    full_trace[y0 : y0 + ny] = trace

    return full_trace


def trace_from_wcs(exp_type, shape, bounding_box, wcs_ref, source_x, source_y, dispaxis):
    """
    Calculate a source trace from WCS.

    The source trace is calculated by projecting a fixed source
    positions onto detector pixels, to get a source location at each
    dispersion element.  For MIRI LRS fixed slit and NIRSpec modes, this
    will be a curved trace, using the sky or slit frame as appropriate.
    For all other modes, a flat trace is returned, containing the
    cross-dispersion position at all dispersion elements.

    Parameters
    ----------
    exp_type : str
        Exposure type for the input data.
    shape : tuple of int
        2D shape for the full input data array, (ny, nx).
    bounding_box : tuple
        A pair of tuples, each consisting of two numbers.
        Represents the range of useful pixel values in both dimensions,
        ((xmin, xmax), (ymin, ymax)).
    wcs_ref : `~gwcs.WCS`
        WCS for the input data model, containing sky and detector
        transforms, forward and backward.
    source_x : float
        X pixel coordinate for the target.
    source_y : float
        Y pixel coordinate for the target.
    dispaxis : int
        Dispersion axis.

    Returns
    -------
    trace : ndarray of float
        Pixel positions in the cross-dispersion direction
        of the trace for each dispersion pixel.
    """
    if exp_type == "MIR_LRS-FIXEDSLIT":
        source_ra, source_dec, _ = wcs_ref(source_x, source_y)
        trace = _miri_trace_from_wcs(shape, bounding_box, wcs_ref, source_ra, source_dec)
    elif exp_type.startswith("NRS"):
        d2s = wcs_ref.get_transform("detector", "slit_frame")
        source_xpos, source_ypos, _ = d2s(source_x, source_y)
        trace = _nirspec_trace_from_wcs(shape, bounding_box, wcs_ref, source_xpos, source_ypos)
    else:
        # Flat trace containing the cross-dispersion position at every element
        if dispaxis == HORIZONTAL:
            trace = np.full(shape[1], np.nan)
            x0 = int(np.ceil(bounding_box[0][0]))
            nx = int(bounding_box[0][1] - bounding_box[0][0])
            trace[x0 : x0 + nx] = source_y
        else:
            trace = np.full(shape[0], np.nan)
            y0 = int(np.ceil(bounding_box[1][0]))
            ny = int(bounding_box[1][1] - bounding_box[1][0])
            trace[y0 : y0 + ny] = source_x

    return trace


def _nod_pair_from_dither(input_model, middle_wl, dispaxis):
    """
    Estimate a nod pair location from the dither offsets.

    Expected location is at the opposite spatial offset from
    the input model.  Requires 'v2v3' transform in the WCS, so
    is only available for unresampled data.

    Parameters
    ----------
    input_model : DataModel
        Model containing WCS and dither data.
    middle_wl : float
        Wavelength at the middle of the array.
    dispaxis : int
        Dispersion axis.

    Returns
    -------
    nod_location : float
        The expected location of the negative trace, in the
        cross-dispersion direction, at the middle wavelength.
    """
    if "v2v3" not in input_model.meta.wcs.available_frames:
        return np.nan

    idltov23 = IdealToV2V3(
        input_model.meta.wcsinfo.v3yangle,
        input_model.meta.wcsinfo.v2_ref,
        input_model.meta.wcsinfo.v3_ref,
        input_model.meta.wcsinfo.vparity,
    )

    if dispaxis == HORIZONTAL:
        x_offset = input_model.meta.dither.x_offset
        y_offset = -input_model.meta.dither.y_offset
    else:
        x_offset = -input_model.meta.dither.x_offset
        y_offset = input_model.meta.dither.y_offset

    dithered_v2, dithered_v3 = idltov23(x_offset, y_offset)

    # v23toworld requires a wavelength along with v2, v3, but value does not affect return
    v23toworld = input_model.meta.wcs.get_transform("v2v3", "world")
    dithered_ra, dithered_dec, _ = v23toworld(dithered_v2, dithered_v3, 0.0)

    x, y = input_model.meta.wcs.backward_transform(dithered_ra, dithered_dec, middle_wl)

    if dispaxis == HORIZONTAL:
        return y
    else:
        return x


def _nod_pair_from_slitpos(input_model, middle_wl):
    """
    Estimate a nod pair location from the source slit position.

    Expected location is at the opposite spatial position from
    the input model.  Requires 'slit_frame' transform in the WCS.
    Implemented only for NIRSpec, assuming horizontal dispersion axis.

    Parameters
    ----------
    input_model : DataModel
        Model containing WCS and dither data.
    middle_wl : float
        Wavelength at the middle of the array.

    Returns
    -------
    nod_location : float
        The expected location of the negative trace, in the
        cross-dispersion direction, at the middle wavelength.
    """
    xpos = input_model.source_xpos
    ypos = -input_model.source_ypos
    wcs = input_model.meta.wcs
    slit2det = wcs.get_transform("slit_frame", "detector")
    if "gwa" in wcs.available_frames:
        # Input is not resampled, wavelengths need to be meters
        _, location = slit2det(xpos, ypos, middle_wl * 1e-6)
    else:
        _, location = slit2det(xpos, ypos, middle_wl)
    return location


def nod_pair_location(input_model, middle_wl):
    """
    Estimate a nod pair location from the WCS.

    For MIRI, it will guess the location from the dither offsets.
    For NIRSpec, it will guess from the slit position.
    For anything else, or if the estimate fails, it will return NaN
    for the location.

    Parameters
    ----------
    input_model : DataModel
        Model containing WCS and dither data.
    middle_wl : float
        Wavelength at the middle of the array.

    Returns
    -------
    nod_location : float
        The expected location of the negative trace, in the
        cross-dispersion direction, at the middle wavelength.
    """
    exp_type = input_model.meta.exposure.type
    nod_center = np.nan
    if exp_type == "MIR_LRS-FIXEDSLIT":
        dispaxis = input_model.meta.wcsinfo.dispersion_direction
        nod_center = _nod_pair_from_dither(input_model, middle_wl, dispaxis)
    elif exp_type.startswith("NRS"):
        nod_center = _nod_pair_from_slitpos(input_model, middle_wl)

    return nod_center
