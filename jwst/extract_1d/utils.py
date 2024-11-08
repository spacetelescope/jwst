import logging
import numpy as np
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.util import wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""


def get_target_coordinates(input_model, slit):
    """Get the right ascension and declination of the target.

    For MultiSlitModel (or similar) data, each slit has the source
    right ascension and declination as attributes, and this can vary
    from one slit to another (e.g. for NIRSpec MOS, or for WFSS).  In
    this case, we want the celestial coordinates from the slit object.
    For other models, however, the celestial coordinates of the source
    are in input_model.meta.target.

    Parameters
    ----------
    input_model : data model
        The input science data model.

    slit : SlitModel or None
        One slit from a MultiSlitModel (or similar), or None if
        there are no slits.

    Returns
    -------
    targ_ra : float or None
        The right ascension of the target, or None

    targ_dec : float or None
        The declination of the target, or None
    """
    targ_ra = None
    targ_dec = None

    if slit is not None:
        # If we've been passed a slit object, get the RA/Dec
        # from the slit source attributes
        targ_ra = getattr(slit, 'source_ra', None)
        targ_dec = getattr(slit, 'source_dec', None)
    elif isinstance(input_model, datamodels.SlitModel):
        # If the input model is a single SlitModel, again
        # get the coords from the slit source attributes
        targ_ra = getattr(input_model, 'source_ra', None)
        targ_dec = getattr(input_model, 'source_dec', None)

    if targ_ra is None or targ_dec is None:
        # Otherwise get it from the generic target coords
        targ_ra = input_model.meta.target.ra
        targ_dec = input_model.meta.target.dec

    # Issue a warning if none of the methods succeeded
    if targ_ra is None or targ_dec is None:
        log.warning("Target RA and Dec could not be determined")
        targ_ra = targ_dec = None

    return targ_ra, targ_dec


def locn_from_wcs(input_model, slit, targ_ra, targ_dec):
    """Get the location of the spectrum, based on the WCS.

    Parameters
    ----------
    input_model : data model
        The input science model.

    slit : one slit from a MultiSlitModel (or similar), or None
        The WCS and target coordinates will be gotten from `slit`
        unless `slit` is None, and in that case they will be gotten
        from `input_model`.

    targ_ra : float or None
        The right ascension of the target, or None

    targ_dec : float or None
        The declination of the target, or None

    Returns
    -------
    middle : int
        Pixel coordinate in the dispersion direction within the 2-D
        cutout (or the entire input image) at the middle of the WCS
        bounding box.  This is the point at which to determine the
        nominal extraction location, in case it varies along the
        spectrum.  The offset will then be the difference between
        `locn` (below) and the nominal location.

    middle_wl : float
        The wavelength at pixel `middle`.

    locn : float
        Pixel coordinate in the cross-dispersion direction within the
        2-D cutout (or the entire input image) that has right ascension
        and declination coordinates corresponding to the target location.
        The spectral extraction region should be centered here.

    None will be returned if there was not sufficient information
    available, e.g. if the wavelength attribute or wcs function is not
    defined.
    """
    if slit is not None:
        wcs_source = slit
    else:
        wcs_source = input_model
    wcs = wcs_source.meta.wcs
    dispaxis = wcs_source.meta.wcsinfo.dispersion_direction

    bb = wcs.bounding_box  # ((x0, x1), (y0, y1))
    if bb is None:
        if slit is None:
            shape = input_model.data.shape
        else:
            shape = slit.data.shape

        bb = wcs_bbox_from_shape(shape)

    if dispaxis == HORIZONTAL:
        # Width (height) in the cross-dispersion direction, from the start of the 2-D cutout (or of the full image)
        # to the upper limit of the bounding box.
        # This may be smaller than the full width of the image, but it's all we need to consider.
        xd_width = int(round(bb[1][1]))  # must be an int
        middle = int((bb[0][0] + bb[0][1]) / 2.)  # Middle of the bounding_box in the dispersion direction.
        x = np.empty(xd_width, dtype=np.float64)
        x[:] = float(middle)
        y = np.arange(xd_width, dtype=np.float64)
        lower = bb[1][0]
        upper = bb[1][1]
    else:  # dispaxis = VERTICAL
        xd_width = int(round(bb[0][1]))  # Cross-dispersion total width of bounding box; must be an int
        middle = int((bb[1][0] + bb[1][1]) / 2.)  # Mid-point of width along dispersion direction
        x = np.arange(xd_width, dtype=np.float64)  # 1-D vector of cross-dispersion (x) pixel indices
        y = np.empty(xd_width, dtype=np.float64)  # 1-D vector all set to middle y index
        y[:] = float(middle)

        # lower and upper range in cross-dispersion direction
        lower = bb[0][0]
        upper = bb[0][1]

    # We need stuff[2], a 1-D array of wavelengths crossing the spectrum near its middle.
    fwd_transform = wcs(x, y)
    middle_wl = np.nanmean(fwd_transform[2])

    # todo - check branches and fallbacks here
    exp_type = input_model.meta.exposure.type
    if exp_type in ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ']:
        if slit is None:
            xpos = input_model.source_xpos
            ypos = input_model.source_ypos
        else:
            xpos = slit.source_xpos
            ypos = slit.source_ypos

        slit2det = wcs.get_transform('slit_frame', 'detector')
        if exp_type == 'NRS_BRIGHTOBJ':
            # Input is not resampled, wavelengths need to be meters
            x_y = slit2det(xpos, ypos, middle_wl * 1e-6)
        else:
            x_y = slit2det(xpos, ypos, middle_wl)
        log.info("Using source_xpos and source_ypos to center extraction.")

    elif exp_type == 'MIR_LRS-FIXEDSLIT':
        try:
            if slit is None:
                dithra = input_model.meta.dither.dithered_ra
                dithdec = input_model.meta.dither.dithered_dec
            else:
                dithra = slit.meta.dither.dithered_ra
                dithdec = slit.meta.dither.dithered_dec
            x_y = wcs.backward_transform(dithra, dithdec, middle_wl)
        except AttributeError:
            log.warning("Dithered pointing location not found in wcsinfo.")
            return
    else:
        log.warning(f"Source position cannot be found for EXP_TYPE {exp_type}")
        return

    # locn is the XD location of the spectrum:
    if dispaxis == HORIZONTAL:
        locn = x_y[1]
    else:
        locn = x_y[0]

    if np.isnan(locn):
        log.warning('Source position could not be determined from WCS.')
        return

    # todo - review this
    if locn < lower or locn > upper and targ_ra > 340.:
        # Try this as a temporary workaround.
        x_y = wcs.backward_transform(targ_ra - 360., targ_dec, middle_wl)

        if dispaxis == HORIZONTAL:
            temp_locn = x_y[1]
        else:
            temp_locn = x_y[0]

        if lower <= temp_locn <= upper:
            # Subtracting 360 from the right ascension worked!
            locn = temp_locn

            log.debug(f"targ_ra changed from {targ_ra} to {targ_ra - 360.}")

    # If the target is at the edge of the image or at the edge of the
    # non-NaN area, we can't use the WCS to find the
    # location of the target spectrum.
    if locn < lower or locn > upper:
        log.warning(f"WCS implies the target is at {locn:.2f}, which is outside the bounding box,")
        log.warning("so we can't get spectrum location using the WCS")
        locn = None

    return middle, middle_wl, locn
