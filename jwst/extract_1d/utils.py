import abc
import logging
import copy
import json
import math
import numpy as np
from stdatamodels import DataModel
from stdatamodels.jwst import datamodels
from typing import Union, Tuple, NamedTuple, List
from astropy.modeling import polynomial
from astropy.io import fits
from gwcs import WCS
from ..lib.wcs_utils import get_wavelengths

# Collection of useful routines

HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""

def locn_from_wcs(
         dispaxis,
         wcs,
         input_model: DataModel,
         slit: Union[DataModel, None],
         targ_ra: Union[float, None],
         targ_dec: Union[float, None],) -> Union[Tuple[int, float, float], None]:

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
    tuple (middle, middle_wl, locn) or None
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

    # WFSS data are not currently supported by this function
    if input_model.meta.exposure.type in ['NIS_WFSS', 'NRC_WFSS', 'NRC_GRISM']:
        log.warning("Can't use target coordinates to get location of spectrum "
                    f"for exp type {input_model.meta.exposure.type}")
        return

    bb = wcs.bounding_box  # ((x0, x1), (y0, y1))
    print('****** ',bb)
    if bb is None:
        if slit is None:
            shape = input_model.data.shape
        else:
            shape = slit.data.shape
        bb = wcs_bbox_from_shape(shape)

    print('dispaxis', dispaxis) 
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

    print('middle wl', middle_wl)
    if input_model.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_MSASPEC',
                                              'NRS_BRIGHTOBJ']:
        if slit is None:
            xpos = input_model.source_xpos
            ypos = input_model.source_ypos
        else:
            xpos = slit.source_xpos
            ypos = slit.source_ypos
        slit2det = wcs.get_transform('slit_frame', 'detector')
        x_y = slit2det(xpos, ypos, middle_wl)
        log.info("Using source_xpos and source_ypos to center extraction.")

    elif input_model.meta.exposure.type == 'MIR_LRS-FIXEDSLIT':
        try:
            if slit is None:
                dithra = input_model.meta.dither.dithered_ra
                dithdec = input_model.meta.dither.dithered_dec
            else:
                dithra = slit.meta.dither.dithered_ra
                dithdec = slit.meta.dither.dithered_dec
                
            x_y = wcs.backward_transform(dithra, dithdec, middle_wl)
        except AttributeError:
            log.warning("Dithered pointing location not found in wcsinfo. "
                        "Defaulting to TARG_RA / TARG_DEC for centering.")
            return


    # locn is the XD location of the spectrum:
    if dispaxis == HORIZONTAL:
        locn = x_y[1]
    else:
        locn = x_y[0]

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

    # If the target is at the edge of the image or at the edge of the non-NaN area, we can't use the WCS to find the
    # location of the target spectrum.
    if locn < lower or locn > upper:
        log.warning(f"WCS implies the target is at {locn:.2f}, which is outside the bounding box,")
        log.warning("so we can't get spectrum location using the WCS")
        locn = None

    return middle, middle_wl, locn



def get_target_coordinates(
            input_model: DataModel, slit: Union[DataModel, None]
    ) -> Tuple[Union[float, None], Union[float, None]]:
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


def replace_bad_values(
        data: np.ndarray,
        input_dq: Union[np.ndarray, None],
        wl_array: np.ndarray
) -> np.ndarray:
    """Replace values flagged with DO_NOT_USE or that have NaN wavelengths.

    Parameters
    ----------
    data : ndarray
        The science data array.

    input_dq : ndarray or None
        If not None, this will be checked for flag value DO_NOT_USE.  The
        science data will be set to NaN for every pixel that is flagged
        with DO_NOT_USE in `input_dq`.

    wl_array : ndarray, 2-D
        Wavelengths corresponding to `data`.  For any element of this
        array that is NaN, the corresponding element in `data` will be
        set to NaN.

    Returns
    -------
    ndarray
        A possibly modified copy of `data`.  If no change was made, this
        will be a view rather than a copy.  Values that are set to NaN
        should not be included when doing the 1-D spectral extraction.

    """
    mask = np.isnan(wl_array)

    if input_dq is not None:
        bad_mask = np.bitwise_and(input_dq, dqflags.pixel['DO_NOT_USE']) > 0
        mask = np.logical_or(mask, bad_mask)

    if np.any(mask):
        mod_data = data.copy()
        mod_data[mask] = np.nan
        return mod_data

    return data


def nans_at_endpoints(
        wavelength: np.ndarray,
        dq: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, slice]:
    """Flag NaNs in the wavelength array.

    Extended summary
    ----------------
    Both input arrays should be 1-D and have the same shape.
    If NaNs are present at endpoints of `wavelength`, the arrays will be
    trimmed to remove the NaNs.  NaNs at interior elements of `wavelength`
    will be left in place, but they will be flagged with DO_NOT_USE in the
    `dq` array.

    Parameters
    ----------
    wavelength : ndarray
        Array of wavelengths, possibly containing NaNs.

    dq : ndarray
        Data quality array.

    Returns
    -------
    wavelength, dq : ndarray
        The returned `dq` array may have NaNs flagged with DO_NOT_USE,
        and both arrays may have been trimmed at either or both ends.

    slc : slice
        The slice to be applied to other output arrays to match the modified
        shape of the wavelength array.
    """
    # The input arrays will not be modified in-place.
    new_wl = wavelength.copy()
    new_dq = dq.copy()
    nelem = wavelength.shape[0]
    slc = slice(nelem)

    nan_mask = np.isnan(wavelength)
    new_dq[nan_mask] = np.bitwise_or(new_dq[nan_mask], dqflags.pixel['DO_NOT_USE'])
    not_nan = np.logical_not(nan_mask)
    flag = np.where(not_nan)

    if len(flag[0]) > 0:
        n_trimmed = flag[0][0] + nelem - (flag[0][-1] + 1)

        if n_trimmed > 0:
            log.debug(f"Output arrays have been trimmed by {n_trimmed} elements")

            slc = slice(flag[0][0], flag[0][-1] + 1)
            new_wl = new_wl[slc]
            new_dq = new_dq[slc]
    else:
        new_dq |= dqflags.pixel['DO_NOT_USE']

    return new_wl, new_dq, slc


def setup_data(
        input_model: DataModel,
#        integ: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray]:
    """Pull out the data from the slit to perform extraction

    Parameters
    ----------
    input_model : data model
        The input science model.

    slit : one slit from a MultiSlitModel (or similar), or None
        If slit is None, the data array is input_model.data; otherwise,
        the data array is slit.data.
        In the former case, if `integ` is zero or larger, the spectrum
        will be extracted from the 2-D slice input_model.data[integ].

    integ : int
        For the case that input_model is a SlitModel or a CubeModel,
        `integ` is the integration number.  If the integration number is
        not relevant (i.e. the data array is 2-D), `integ` should be -1.



    """

    # hook for lrs_slitless data 
    #if integ > -1:
    #    log.info(f"Extracting integration {integ + 1}")
    #    data = input_model.data[integ]
    #    var_poisson = input_model.var_poisson[integ]
    #    var_rnoise = input_model.var_rnoise[integ]
    #    var_flat = input_model.var_flat[integ]
    #    input_dq = input_model.dq[integ]

    bbox = input_model.meta,wcs.bounding_box
    x0, x1 = bbox[0]
    y0, y1 = bbox[1]
    i1, i2, j1, j2 = (int(np.ceil(y0)), int(np.ceil(y1)), int(np.round(x0)), int(np.round(x1)))
    data = model.data[i1:i2, j1:j2]
    var_poisson = model.var_poisson[i1:i2, j1:j2]
    var_rnoise = model.var_rnoise[i1:i2, j1:j2]
    dq = model.dq[i1:i2, j1:j2]
    
    #  Ensure variance arrays have been populated. If not, zero fill.
    if np.shape(var_poisson) != np.shape(data):
        var_poisson = np.zeros_like(data)
        var_rnoise = np.zeros_like(data)
        
    if np.shape(var_flat) != np.shape(data):
        var_flat = np.zeros_like(data)
        
    if dq.size == 0:
        dq = None

    # not sure we need the wl_array - might do this differently 
    wl_array = get_wavelengths(input_model)
    wl_array = wl_array[i1:i2, j1:j2]
    data = replace_bad_values(data, dq, wl_array)

    # Use uniform weights for now, setting invalid values to zero.
    weights = np.ones(data.shape)
    weights[dq%2 == 1] = 0
    weights[np.isnan(var_rnoise)] = 0

    var_rnoise[np.isnan(var_rnoise)] = 1e30
    var_rnoise[var_rnoise <= 0] = 1e30

    var_poisson[np.isnan(var_poisson)] = 1e30
    var_poisson[var_poisson <= 0] = 1e30
    data[np.isnan(data)] = 0

    return data, dq, var_poisson, var_rnoise, weights
