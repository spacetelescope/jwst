import logging
import numpy as np

from scipy import ndimage, optimize
from stdatamodels.jwst.datamodels import SpecPsfModel

from jwst.extract_1d.extract1d import extract1d
from jwst.extract_1d.source_location import middle_from_wcs, nod_pair_location, trace_from_wcs

__all__ = ["psf_profile"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""

NOD_PAIR_PATTERN = ["ALONG-SLIT-NOD", "2-POINT-NOD"]


def open_psf(psf_refname, exp_type):
    """
    Open the PSF reference file.

    Parameters
    ----------
    psf_refname : str
        The name of the psf reference file.
    exp_type : str
        The exposure type of the data.

    Returns
    -------
    psf_model : SpecPsfModel
        Returns the EPSF model.
    """
    if exp_type == "MIR_LRS-FIXEDSLIT":
        # The information we read in from PSF file is:
        # center_col: psf_model.meta.psf.center_col
        # super sample factor: psf_model.meta.psf.subpix)
        # psf : psf_model.data (2d)
        # wavelength of PSF planes: psf_model.wave
        psf_model = SpecPsfModel(psf_refname)

    else:
        # So far, only MIRI LRS has a PSF datamodel defined. For any other
        # exposure type, try to use the model MIRI LRS uses to open the input model
        try:
            psf_model = SpecPsfModel(psf_refname)
        except (ValueError, AttributeError):
            raise NotImplementedError(
                f"PSF file for EXP_TYPE {exp_type} could not be read as SpecPsfModel."
            ) from None
    return psf_model


def _normalize_profile(profile, dispaxis):
    """Normalize a spatial profile along the cross-dispersion axis."""
    if dispaxis == HORIZONTAL:
        psum = np.nansum(profile, axis=0)
        nz = psum != 0
        profile[:, nz] = profile[:, nz] / psum[nz]
        profile[:, ~nz] = 0.0
    else:
        psum = np.nansum(profile, axis=1)
        nz = psum != 0
        profile[nz, :] = profile[nz, :] / psum[nz, None]
        profile[~nz, :] = 0.0
    profile[~np.isfinite(profile)] = 0.0


def _make_cutout_profile(
    xidx, yidx, psf_subpix, psf_data, dispaxis, extra_shift=0.0, nod_offset=None
):
    """
    Make a spatial profile corresponding to the data cutout.

    Input index values should already contain the shift to the trace location
    in the cross-dispersion direction and any wavelength shifts necessary
    in the dispersion direction.

    Parameters
    ----------
    xidx : ndarray of float
        Index array for x values.
    yidx : ndarray of float
        Index array for y values.
    psf_subpix : float
        Scaling factor for pixel size in the PSF data.
    psf_data : ndarray of float
        2D PSF model.
    dispaxis : int
        Dispersion axis.
    extra_shift : float, optional
        An extra shift for the primary trace location, to be added to the
        cross-dispersion indices.
    nod_offset : float, optional
        If not None, a negative trace is added to the spatial profile,
        with a cross-dispersion shift of `nod_offset`.

    Returns
    -------
    profiles : list of ndarray of float
        2D spatial profiles containing the primary trace and, optionally,
        a negative trace for a nod pair.  The profiles are normalized along
        the cross-dispersion axis.
    """
    # Add an extra spatial shift to the primary trace
    if dispaxis == HORIZONTAL:
        xmap = xidx
        ymap = yidx + extra_shift * psf_subpix
    else:
        xmap = xidx + extra_shift * psf_subpix
        ymap = yidx
    sprofile = ndimage.map_coordinates(psf_data, [ymap, xmap], order=1)
    _normalize_profile(sprofile, dispaxis)

    if nod_offset is None:
        return [sprofile]

    # Make an additional profile for the negative nod if desired
    if dispaxis == HORIZONTAL:
        ymap += psf_subpix * nod_offset
    else:
        xmap += psf_subpix * nod_offset

    nod_profile = ndimage.map_coordinates(psf_data, [ymap, xmap], order=1)
    _normalize_profile(nod_profile, dispaxis)

    return [sprofile, nod_profile * -1]


def _profile_residual(
    shifts_to_optimize, cutout, cutout_var, xidx, yidx, psf_subpix, psf_data, dispaxis, fit_bkg=True
):
    """
    Residual function to minimize for optimizing trace locations.

    Call `_make_cutout_profile` to generate a profile from input parameters.
    Call `extract1d` to generate a scene model from the data and the new profile.
    Compute a residual value from the sum of (model - cutout) ** 2 / cutout_var.

    Parameters
    ----------
    shifts_to_optimize : list of float
        The first value is used as the `extra_shift` parameter to
        `_make_cutout_profile`.  If two are provided, the second value is
        used as the `nod_offset` parameter to `_make_cutout_profile`.
        If only one value is provided, no nod offset is applied.
    cutout : ndarray
        Input data array, trimmed to the bounding box.
    cutout_var : ndarray
        Input variance, matching the data array, used for weighting the
        output residuals.
    xidx : ndarray of float
        Index array for x values.
    yidx : ndarray of float
        Index array for y values.
    psf_subpix : float
        Scaling factor for pixel size in the PSF data.
    psf_data : ndarray of float
        2D PSF model.
    dispaxis : int
        Dispersion axis.
    fit_bkg : bool, optional
        If True, background subtraction is performed during extraction.

    Returns
    -------
    float
        The residual value to minimize.
    """
    if len(shifts_to_optimize) > 1:
        nod_offset = shifts_to_optimize[1]
    else:
        nod_offset = None
    sprofiles = _make_cutout_profile(
        xidx,
        yidx,
        psf_subpix,
        psf_data,
        dispaxis,
        extra_shift=shifts_to_optimize[0],
        nod_offset=nod_offset,
    )
    extract_kwargs = {
        "extraction_type": "optimal",
        "fit_bkg": fit_bkg,
        "bkg_fit_type": "poly",
        "bkg_order": 0,
    }
    if dispaxis == HORIZONTAL:
        empty_var = np.zeros_like(cutout)
        result = extract1d(cutout, sprofiles, cutout_var, empty_var, empty_var, **extract_kwargs)
        model = result[-1]
    else:
        sprofiles = [profile.T for profile in sprofiles]
        empty_var = np.zeros_like(cutout.T)
        result = extract1d(
            cutout.T, sprofiles, cutout_var.T, empty_var, empty_var, **extract_kwargs
        )
        model = result[-1].T
    return np.nansum((model - cutout) ** 2 / cutout_var)


def psf_profile(
    input_model, trace, wl_array, psf_ref_name, optimize_shifts=True, model_nod_pair=True
):
    """
    Create a spatial profile from a PSF reference.

    Provides PSF-based profiles for point sources in slit-like data containing
    one positive trace and, optionally, one negative trace resulting from nod
    subtraction.  The location of the positive trace should be provided in the
    `trace` input parameter; the negative trace location will be guessed from
    the input metadata. If a negative trace is modeled, it is recommended that
    `optimize_shifts` also be set to True, to improve the initial guess for the
    trace location.

    Parameters
    ----------
    input_model : data model
        This can be either the input science file or one SlitModel out of
        a list of slits.
    trace : ndarray or None
        Array of source cross-dispersion position values, one for each
        dispersion element in the input model data.  If None, the source
        is assumed to be at the center of the slit.
    wl_array : ndarray
        Array of wavelength values, matching the input model data shape, for
        each pixel in the array.
    psf_ref_name : str
        PSF reference filename.
    optimize_shifts : bool, optional
        If True, the spatial location of the trace will be optimized by
        minimizing the residuals in a scene model compared to the data in
        the first integration of `input_model`.
    model_nod_pair : bool, optional
        If True, and if background subtraction has taken place, a negative
        PSF will be modeled at the mirrored spatial location of the positive
        trace.

    Returns
    -------
    profile : ndarray
        Spatial profile matching the input data.
    lower_limit : int
        Lower limit of the aperture in the cross-dispersion direction.
        For PSF profiles, this is always set to the lower edge of the bounding box,
        since the full array may have non-zero weight.
    upper_limit : int
        Upper limit of the aperture in the cross-dispersion direction.
        For PSF profiles, this is always set to the upper edge of the bounding box,
        since the full array may have non-zero weight.
    """
    # Read in reference files
    exp_type = input_model.meta.exposure.type
    psf_model = open_psf(psf_ref_name, exp_type)

    # Get the data cutout
    data_shape = input_model.data.shape[-2:]
    dispaxis = input_model.meta.wcsinfo.dispersion_direction
    wcs = input_model.meta.wcs
    bbox = wcs.bounding_box

    y0 = int(np.ceil(bbox[1][0]))
    y1 = int(np.ceil(bbox[1][1]))
    x0 = int(np.ceil(bbox[0][0]))
    x1 = int(np.ceil(bbox[0][1]))
    if input_model.data.ndim == 3:
        # use the first integration only
        cutout = input_model.data[0, y0:y1, x0:x1]
        cutout_var = input_model.var_rnoise[0, y0:y1, x0:x1]
    else:
        cutout = input_model.data[y0:y1, x0:x1]
        cutout_var = input_model.var_rnoise[y0:y1, x0:x1]
    cutout_wl = wl_array[y0:y1, x0:x1]

    # Get the nominal center of the cutout
    middle_disp, middle_xdisp, middle_wl = middle_from_wcs(wcs, bbox, dispaxis)

    # Get the effective index into the 1D PSF wavelengths from the data wavelengths
    psf_wave = psf_model.wave
    sort_idx = np.argsort(psf_wave)
    valid_wave = np.isfinite(psf_wave[sort_idx])
    wave_idx = np.interp(
        cutout_wl, psf_wave[sort_idx][valid_wave], sort_idx[valid_wave], left=np.nan, right=np.nan
    )

    if trace is None:
        # Don't try to model a negative pair if we don't have a trace to start
        if model_nod_pair:
            log.warning("Cannot model a negative nod without position information")
            model_nod_pair = False

        # Set the location to the middle of the cross-dispersion
        # all the way across the array
        location = middle_xdisp
        if dispaxis == HORIZONTAL:
            trace = trace_from_wcs(
                exp_type, data_shape, bbox, wcs, middle_disp, middle_xdisp, dispaxis
            )
        else:
            trace = trace_from_wcs(
                exp_type, data_shape, bbox, wcs, middle_xdisp, middle_disp, dispaxis
            )

    else:
        # Nominal location from the middle dispersion point
        location = trace[int(np.round(middle_disp))]

    # Trim the trace to the data cutout
    if dispaxis == HORIZONTAL:
        trace = trace[x0:x1]
    else:
        trace = trace[y0:y1]

    # Check if we need to add a negative nod pair trace
    nod_offset = None
    if model_nod_pair:
        nod_subtracted = str(input_model.meta.cal_step.bkg_subtract) == "COMPLETE"
        pattype_ok = str(input_model.meta.dither.primary_type) in NOD_PAIR_PATTERN
        if not nod_subtracted:
            log.info("Input data was not nod-subtracted. A negative trace will not be modeled.")
        elif not pattype_ok:
            log.info("Input data was not a two-point nod. A negative trace will not be modeled.")
        else:
            nod_center = nod_pair_location(input_model, middle_wl)
            if np.isnan(nod_center) or (np.abs(location - nod_center) < 2):
                log.warning("Nod center could not be estimated from the WCS.")
                log.warning("The negative nod will not be modeled.")
            else:
                if not optimize_shifts:
                    log.warning("Negative nod locations are currently approximations only.")
                    log.warning(
                        "PSF location optimization is recommended when negative nods are modeled."
                    )
                nod_offset = location - nod_center

    # Get an index grid for the data cutout
    cutout_shape = cutout.shape
    _y, _x = np.mgrid[: cutout_shape[0], : cutout_shape[1]]

    # Scale the trace location to the subsampled psf and
    # add the wavelength and spatial shifts to the coordinates to map to
    psf_subpix = psf_model.meta.psf.subpix
    psf_location = trace - bbox[0][0]
    if dispaxis == HORIZONTAL:
        psf_shift = psf_model.meta.psf.center_row - (psf_location * psf_subpix)
        xidx = wave_idx
        yidx = _y * psf_subpix + psf_shift
    else:
        psf_shift = psf_model.meta.psf.center_col - (psf_location * psf_subpix)
        xidx = _x * psf_subpix + psf_shift[:, None]
        yidx = wave_idx

    # If desired, add additional spatial shifts to the starting locations of
    # the primary trace (and negative nod pair trace if necessary)
    if optimize_shifts:
        log.info("Optimizing trace locations")
        if nod_offset is None:
            (extra_shift,) = optimize.minimize(
                _profile_residual,
                [0.0],
                (cutout, cutout_var, xidx, yidx, psf_subpix, psf_model.data, dispaxis),
                method="Nelder-Mead",
            ).x
        else:
            extra_shift, nod_offset = optimize.minimize(
                _profile_residual,
                [0.0, nod_offset],
                (cutout, cutout_var, xidx, yidx, psf_subpix, psf_model.data, dispaxis),
                method="Nelder-Mead",
            ).x
        location -= extra_shift
    else:
        extra_shift = 0.0

    log.info(f"Centering profile on spectrum at {location:.2f}, wavelength {middle_wl:.2f}")
    if nod_offset is not None:
        log.info(
            f"Also modeling a negative trace at {location - nod_offset:.2f} "
            f"(offset: {nod_offset:.2f})"
        )

    # Make a spatial profile from the shifted PSF data
    sprofiles = _make_cutout_profile(
        xidx,
        yidx,
        psf_subpix,
        psf_model.data,
        dispaxis,
        extra_shift=extra_shift,
        nod_offset=nod_offset,
    )

    # Make the output profile, matching the input data
    output_y = _y + y0
    output_x = _x + x0
    valid = (output_y >= 0) & (output_y < y1) & (output_x >= 0) & (output_x < x1)
    profiles = []
    for sprofile in sprofiles:
        profile = np.full(data_shape, 0.0)
        profile[output_y[valid], output_x[valid]] = sprofile[valid]
        profiles.append(profile)

    if dispaxis == HORIZONTAL:
        limits = (y0, y1)
    else:
        limits = (x0, x1)
    return profiles, *limits
