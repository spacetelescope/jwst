import logging
import numpy as np

from scipy.interpolate import CubicSpline
from scipy import interpolate
from scipy import ndimage
from astropy.io import fits

from stdatamodels.jwst.datamodels import MiriLrsPsfModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""


def open_specwcs(specwcs_ref_name: str, exp_type: str):
    """Open the specwcs reference file.

    Currently only works on MIRI LRS-FIXEDSLIT exposures.

    Parameters
    ----------
    specwcs_ref_name : str
        The name of the specwcs reference file. This file contains
        information of the trace location. For MIRI LRS-FIXEDSlIT it
        is a FITS file containing the x,y center of the trace.
    ext_type : str
        The exposure type of the data.

    Returns
    -------
    trace, wave_trace, wavetab
        Center of the trace in x and y for a given wavelength.

    """
    if exp_type == 'MIR_LRS-FIXEDSLIT':
        # use fits to read file (the datamodel does not have all that is needed)
        ref = fits.open(specwcs_ref_name)

        with ref:
            lrsdata = np.array([d for d in ref[1].data])
            # Get the zero point from the reference data.
            # The zero_point is X, Y  (which should be COLUMN, ROW)
            # These are 1-indexed in CDP-7 (i.e., SIAF convention) so must be converted to 0-indexed
            # for lrs_fixedslit
            zero_point = ref[0].header['imx'] - 1, ref[0].header['imy'] - 1

        # In the lrsdata reference table, X_center,Y_center, wavelength  relative to zero_point

        xcen = lrsdata[:, 0]
        ycen = lrsdata[:, 1]
        wavetab = lrsdata[:, 2]
        trace = xcen + zero_point[0]
        wave_trace = ycen + zero_point[1]

    else:
        raise NotImplementedError(f'Specwcs files for EXP_TYPE {exp_type} '
                                  f'are not supported.')

    return trace, wave_trace, wavetab
    

def open_psf(psf_refname: str, exp_type: str):
    """Open the PSF reference file.

    Parameters
    ----------
    psf_ref_name : str
        The name of the psf reference file. 
    ext_type : str
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
        raise NotImplementedError(f'PSF files for EXP_TYPE {exp_type} '
                                  f'are not supported.')
    return psf_model 


def psf_profile(input_model, psf_ref_name, specwcs_ref_name, middle_wl, location):
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
    # Check input exposure type
    exp_type = input_model.meta.exposure.type
    if exp_type != 'MIR_LRS-FIXEDSLIT':
        raise NotImplementedError(f'PSF extraction is not supported for '
                                  f'EXP_TYPE {exp_type}')

    # Read in reference files
    trace, wave_trace, wavetab = open_specwcs(specwcs_ref_name, exp_type)
    psf_model = open_psf(psf_ref_name, exp_type)

    dispaxis = input_model.meta.wcsinfo.dispersion_direction
    wcs = input_model.meta.wcs
    bbox = wcs.bounding_box
    center_x = np.mean(bbox[0])
    center_y = np.mean(bbox[1])

    # Determine the location using the WCS
    if middle_wl is None:
        _, _, middle_wl = wcs(center_x, center_y)
    if location is None:
        if dispaxis == HORIZONTAL:
            location = center_y
        else:
            location = center_x

    y0 = int(np.ceil(bbox[1][0]))
    y1 = int(np.ceil(bbox[1][1]))
    x0 = int(np.round(bbox[0][0]))
    x1 = int(np.round(bbox[0][1]))
    cutout = input_model.data[y0:y1, x0:x1]

    # Perform fit of reference trace and corresponding wavelength
    # The wavelength for the reference trace does not exactly line up exactly with the data
    cs = CubicSpline(wavetab, trace)
    cen_shift = cs(middle_wl)
    shift = location - cen_shift
    log.info(f'Centering profile on spectrum at {location}, wavelength {middle_wl}')
    log.info(f'For this wavelength, the reference trace location is at {cen_shift}')
    log.info(f'Shift to apply to ref trace: {shift}')

    # todo - if possible, fix this for s2d -
    #  cen_shift is wrong, wavelengths don't match PSF

    # adjust the trace to the slit region
    trace_cutout = trace - bbox[0][0]
    trace_shift = trace_cutout + shift
    psf_wave = psf_model.wave

    # trace_shift: for each wavelength in the PSF, this is the shift in x to apply
    #   to the PSF image to shift it to fall on the source.
    # wavetab : this is the wavelength corresponding to the trace.
    #   This wavelength may not match exactly to the PSF.

    # Determine what the shifts per row are for the wavelengths
    # given by the model PSF
    psf_subpix = psf_model.meta.psf.subpix

    psf_interp = interpolate.interp1d(wavetab, trace_shift, fill_value="extrapolate")
    psf_shift = psf_interp(psf_wave)
    psf_shift = psf_model.meta.psf.center_col - (psf_shift * psf_subpix)

    # Note: this assumes that data wavelengths are identical to PSF wavelengths
    data_shape = cutout.shape
    _y, _x = np.mgrid[:data_shape[0], :data_shape[1]]

    # Scale cross-dispersion coordinates by subpixel value and shift by trace
    if dispaxis == HORIZONTAL:
        if data_shape[1] != psf_shift.size:
            log.error('Data shape does not match PSF reference.')
            log.error('Optimal extraction must be performed on cal files.')
            raise NotImplementedError('Optimal extraction not implemented for resampled data.')

        _y_sh = _y * psf_subpix + psf_shift
        sprofile = ndimage.map_coordinates(psf_model.data, [_y_sh, _x], order=1)
    else:
        if data_shape[0] != psf_shift.size:
            log.error('Data shape does not match PSF reference.')
            log.error('Optimal extraction must be performed on cal files.')
            raise NotImplementedError('Optimal extraction not implemented for resampled data.')

        _x_sh = _x * psf_subpix + psf_shift[:, np.newaxis]
        sprofile = ndimage.map_coordinates(psf_model.data, [_y, _x_sh], order=1)

    # Normalize the spatial profile at each dispersion element
    if dispaxis == HORIZONTAL:
        psum = np.sum(sprofile, axis=0)
        sprofile[:, psum > 0] = sprofile[:, psum > 0] / psum[psum > 0]
        sprofile[:, psum <= 0] = 0.0
    else:
        psum = np.sum(sprofile, axis=1)
        sprofile[psum > 0, :] = sprofile[psum > 0, :] / psum[psum > 0, None]
        sprofile[psum <= 0, :] = 0.0
    sprofile[~np.isfinite(sprofile)] = 0.0
    sprofile[sprofile < 0] = 0.0

    # Make the output profile, matching the input data
    data_shape = input_model.data.shape
    profile = np.full(data_shape, 0.0)
    output_y = _y + y0
    output_x = _x + x0
    valid = (output_y >= 0) & (output_y < y1) & (output_x >= 0) & (output_x < x1)
    profile[output_y[valid], output_x[valid]] = sprofile[valid]

    if dispaxis == HORIZONTAL:
        limits = (y0, y1)
    else:
        limits = (x0, x1)
    return profile, *limits
