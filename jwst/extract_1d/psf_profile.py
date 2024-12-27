import logging
import numpy as np

from scipy.interpolate import CubicSpline
from scipy import interpolate, ndimage, optimize
from astropy.io import fits

from stdatamodels.jwst.datamodels import MiriLrsPsfModel
from stdatamodels.jwst.transforms.models import IdealToV2V3

from jwst.extract_1d.extract1d import extract1d

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""


def open_specwcs(specwcs_ref_name, exp_type):
    """Open the specwcs reference file.

    Currently only works on MIRI LRS-FIXEDSLIT exposures.

    Parameters
    ----------
    specwcs_ref_name : str
        The name of the specwcs reference file. This file contains
        information of the trace location. For MIRI LRS-FIXEDSlIT it
        is a FITS file containing the x,y center of the trace.
    exp_type : str
        The exposure type of the data.

    Returns
    -------
    trace, wave_trace, wavetab
        Center of the trace in x and y for a given wavelength.

    """
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        # use fits to read file (the datamodel does not have all that is needed)
        with fits.open(specwcs_ref_name) as ref:
            lrsdata = np.array([d for d in ref[1].data])
            # Get the zero point from the reference data.
            # The zero_point is X, Y  (which should be COLUMN, ROW)
            # These are 1-indexed in CDP-7 (i.e., SIAF convention) so must be converted to 0-indexed
            # for lrs_fixedslit
            zero_point = ref[0].header['imx'] - 1, ref[0].header['imy'] - 1

        # In the lrsdata reference table:
        # X_center, Y_center, wavelength, relative to zero_point
        xcen = lrsdata[:, 0]
        ycen = lrsdata[:, 1]
        wavetab = lrsdata[:, 2]
        trace = xcen + zero_point[0]
        wave_trace = ycen + zero_point[1]

    else:
        raise NotImplementedError(f'Specwcs files for EXP_TYPE {exp_type} '
                                  f'are not supported.')

    return trace, wave_trace, wavetab
    

def open_psf(psf_refname, exp_type):
    """Open the PSF reference file.

    Parameters
    ----------
    psf_refname : str
        The name of the psf reference file. 
    exp_type : str
        The exposure type of the data.

    Returns
    -------
    psf_model : MiriLrsPsfModel
        Currently only works on MIRI LRS-FIXEDSLIT exposures.
        Returns the EPSF model.

    """
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        # The information we read in from PSF file is:
        # center_col: psf_model.meta.psf.center_col
        # super sample factor: psf_model.meta.psf.subpix)
        # psf : psf_model.data (2d)
        # wavelength of PSF planes: psf_model.wave
        psf_model = MiriLrsPsfModel(psf_refname)

    else:
        # So far, only MIRI LRS has a PSF datamodel defined. Try to use it
        # for the input model.
        try:
            psf_model = MiriLrsPsfModel(psf_refname)
        except (ValueError, AttributeError):
            raise NotImplementedError(f'PSF file for EXP_TYPE {exp_type} '
                                      f'could not be read as MiriLrsPsfModel.') from None
    return psf_model 


def _normalize_profile(profile, dispaxis):
    """Normalize a spatial profile along the cross-dispersion axis."""
    if dispaxis == HORIZONTAL:
        psum = np.sum(profile, axis=0)
        profile[:, psum != 0] = profile[:, psum != 0] / psum[psum != 0]
        profile[:, psum == 0] = 0.0
    else:
        psum = np.sum(profile, axis=1)
        profile[psum != 0, :] = profile[psum != 0, :] / psum[psum != 0, None]
        profile[psum == 0, :] = 0.0
    profile[~np.isfinite(profile)] = 0.0


def _make_cutout_profile(xidx, yidx, psf_subpix, psf_data, dispaxis,
                         extra_shift=0.0, nod_offset=None):
    """Make a spatial profile corresponding to the data cutout.

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


def _profile_residual(param, cutout, cutout_var, xidx, yidx, psf_subpix, psf_data, dispaxis):
    """Residual function to minimize for optimizing trace locations."""
    sprofiles = _make_cutout_profile(xidx, yidx, psf_subpix, psf_data, dispaxis,
                                     extra_shift=param[0], nod_offset=param[1])
    extract_kwargs = {'extraction_type': 'optimal',
                      'fit_bkg': True,
                      'bkg_fit_type': 'poly',
                      'bkg_order': 0}
    if dispaxis == HORIZONTAL:
        empty_var = np.zeros_like(cutout)
        result = extract1d(cutout, sprofiles, cutout_var, empty_var, empty_var,
                           **extract_kwargs)
        model = result[-1]
    else:
        sprofiles = [profile.T for profile in sprofiles]
        empty_var = np.zeros_like(cutout.T)
        result = extract1d(cutout.T, sprofiles, cutout_var.T, empty_var, empty_var,
                           **extract_kwargs)
        model = result[-1].T
    return np.nansum((model - cutout) ** 2 / cutout_var)


def nod_pair_location(input_model, middle_wl, dispaxis):
    """Estimate a nod pair location from the WCS.

    Expected location is at the opposite spatial offset from
    the input model.

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
    if 'v2v3' not in input_model.meta.wcs.available_frames:
        return np.nan

    idltov23 = IdealToV2V3(
        input_model.meta.wcsinfo.v3yangle,
        input_model.meta.wcsinfo.v2_ref, input_model.meta.wcsinfo.v3_ref,
        input_model.meta.wcsinfo.vparity
    )

    if dispaxis == HORIZONTAL:
        x_offset = input_model.meta.dither.x_offset
        y_offset = -input_model.meta.dither.y_offset
    else:
        x_offset = -input_model.meta.dither.x_offset
        y_offset = input_model.meta.dither.y_offset

    dithered_v2, dithered_v3 = idltov23(x_offset, y_offset)

    # v23toworld requires a wavelength along with v2, v3, but value does not affect return
    v23toworld = input_model.meta.wcs.get_transform('v2v3', 'world')
    dithered_ra, dithered_dec, _ = v23toworld(dithered_v2, dithered_v3, 0.0)

    x, y = input_model.meta.wcs.backward_transform(dithered_ra, dithered_dec, middle_wl)

    if dispaxis == HORIZONTAL:
        return y
    else:
        return x


def psf_profile(input_model, psf_ref_name, specwcs_ref_name, middle_wl, location,
                wl_array, optimize_shifts=True, model_nod_pair=True):
    """Create a spatial profile from a PSF reference.

    Currently only works on MIRI LRS-FIXEDSLIT exposures.
    Input data must be point source.

    The extraction routine can support multiple sources for
    simultaneous extraction, but for this first version, we will assume
    one source only, located at the planned position (dither RA/Dec), and
    return a single profile.

    Parameters
    ----------
    input_model : data model
        This can be either the input science file or one SlitModel out of
        a list of slits.
    psf_ref_name : str
        PSF reference filename.
    specwcs_ref_name : str
        Reference file containing information on the spectral trace.
    middle_wl : float or None
        Wavelength value to use as the center of the trace. If not provided,
        the wavelength at the center of the bounding box will be used.
    location : float or None
        Spatial index to use as the center of the trace.  If not provided,
        the location at the center of the bounding box will be used.
    wl_array : ndarray
        Array of wavelength values, matching the input model data shape, for
        each pixel in the array.
    optimize_shifts : bool, optional
        If True, the spatial location of the trace will be optimized via
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
    # Check if data is resampled
    if hasattr(input_model.meta.cal_step, 'resample'):
        resampled = str(input_model.meta.cal_step.resample) == 'COMPLETE'
    else:
        # guess from the wcs
        resampled = 'v2v3' not in input_model.meta.wcs.available_frames

    # Check input exposure type
    exp_type = input_model.meta.exposure.type
    if exp_type == 'MIR_LRS-FIXEDSLIT' and not resampled and specwcs_ref_name != 'N/A':
        trace, wave_trace, wavetab = open_specwcs(specwcs_ref_name, exp_type)
    else:
        trace, wave_trace, wavetab = None, None, None

    # Read in reference files
    psf_model = open_psf(psf_ref_name, exp_type)

    dispaxis = input_model.meta.wcsinfo.dispersion_direction
    wcs = input_model.meta.wcs
    bbox = wcs.bounding_box
    center_x = np.mean(bbox[0])
    center_y = np.mean(bbox[1])

    # Determine the location using the WCS
    if middle_wl is None or np.isnan(middle_wl):
        _, _, middle_wl = wcs(center_x, center_y)
    if dispaxis == HORIZONTAL:
        nominal = center_y
    else:
        nominal = center_x
    if location is None or np.isnan(location):
        location = nominal

    y0 = int(np.ceil(bbox[1][0]))
    y1 = int(np.ceil(bbox[1][1]))
    x0 = int(np.round(bbox[0][0]))
    x1 = int(np.round(bbox[0][1]))
    if input_model.data.ndim == 3:
        cutout = input_model.data[0, y0:y1, x0:x1]
        cutout_var = input_model.var_rnoise[0, y0:y1, x0:x1]
    else:
        cutout = input_model.data[y0:y1, x0:x1]
        cutout_var = input_model.var_rnoise[y0:y1, x0:x1]
    cutout_wl = wl_array[y0:y1, x0:x1]

    # Get the effective index into the 1D PSF wavelengths from the data wavelengths
    psf_wave = psf_model.wave
    sort_idx = np.argsort(psf_wave)
    valid_wave = np.isfinite(psf_wave[sort_idx])
    wave_idx = np.interp(cutout_wl, psf_wave[sort_idx][valid_wave], sort_idx[valid_wave],
                         left=np.nan, right=np.nan)

    # Get the spatial shift for the output wavelengths
    psf_subpix = psf_model.meta.psf.subpix
    if trace is not None:
        # Perform fit of reference trace and corresponding wavelength
        # The wavelength for the reference trace does not exactly line up
        # exactly with the PSF data
        cs = CubicSpline(wavetab, trace)
        cen_shift = cs(middle_wl)
        shift = location - cen_shift
        log.debug(f'For this wavelength, the reference trace location is at {cen_shift:.2f}')
        log.debug(f'Shift applied to reference trace: {shift:.2f}')

        # adjust the trace to the slit region
        trace_cutout = trace - bbox[0][0]
        trace_shift = trace_cutout + shift

        psf_interp = interpolate.interp1d(wavetab, trace_shift, fill_value="extrapolate")
        psf_location = psf_interp(cutout_wl)

    else:
        psf_location = location

    psf_shift = psf_model.meta.psf.center_col - (psf_location * psf_subpix)

    # Check if we need to add a negative nod pair trace
    if model_nod_pair:
        nod_subtracted = str(input_model.meta.cal_step.back_sub) == 'COMPLETE'
        if not nod_subtracted:
            log.info('Input data was not nod-subtracted. '
                     'A negative trace will not be modeled.')
            nod_offset = None
        else:
            nod_center = nod_pair_location(input_model, middle_wl, dispaxis)
            if np.isnan(nod_center):
                log.warning('Nod center could not be estimated from the WCS.')
                log.warning('The negative nod will not be modeled.')
                nod_offset = None
            else:
                nod_offset = location - nod_center
    else:
        nod_offset = None

    # Get an index grid for the data cutout
    data_shape = cutout.shape
    _y, _x = np.mgrid[:data_shape[0], :data_shape[1]]

    # Add the wavelength and spatial shifts to the coordinates to map to
    if dispaxis == HORIZONTAL:
        xidx = wave_idx
        yidx = _y * psf_subpix + psf_shift
    else:
        xidx = _x * psf_subpix + psf_shift
        yidx = wave_idx

    # If desired, add additional spatial shifts to the starting locations of
    # the primary trace (and negative nod pair trace if necessary)
    if optimize_shifts:
        log.info('Optimizing trace locations')
        extra_shift, nod_offset = optimize.minimize(
            _profile_residual, [0.0, nod_offset],
            (cutout, cutout_var, xidx, yidx, psf_subpix, psf_model.data, dispaxis)).x
        location += extra_shift
    else:
        extra_shift = 0.0

    log.info(f'Centering profile on spectrum at {location:.2f}, wavelength {middle_wl:.2f}')
    if nod_offset is not None:
        log.info(f'Also modeling a negative trace at {location - nod_offset:.2f} '
                 f'(offset: {nod_offset:.2f})')

    # Make a spatial profile from the shifted PSF data
    sprofiles = _make_cutout_profile(xidx, yidx, psf_subpix, psf_model.data,
                                     dispaxis, extra_shift=extra_shift,
                                     nod_offset=nod_offset)

    # Make the output profile, matching the input data
    data_shape = input_model.data.shape
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
