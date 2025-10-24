import logging
import multiprocessing as mp
import warnings

import numpy as np
from astropy.modeling.mappings import Mapping
from scipy import sparse

from jwst.lib.winclip import get_clipped_pixels
from jwst.wfss_contam.sens1d import create_1d_sens

log = logging.getLogger(__name__)


__all__ = ["disperse"]


def _determine_native_wl_spacing(
    x0_sky,
    y0_sky,
    sky_to_imgxy,
    imgxy_to_grismxy,
    order,
    wmin,
    wmax,
    oversample_factor=2,
):
    """
    Determine the wavelength spacing necessary to adequately sample the dispersed frame.

    Parameters
    ----------
    x0_sky : float or ndarray
        RA of the input pixel position in direct image and segmentation map
    y0_sky : float or ndarray
        Dec of the input pixel position in direct image and segmentation map
    sky_to_imgxy : astropy model
        Transform from sky to image coordinates
    imgxy_to_grismxy : astropy model
        Transform from image to grism coordinates
    order : int
        Spectral order number
    wmin : float
        Minimum wavelength for dispersed spectra
    wmax : float
        Maximum wavelength for dispersed spectra
    oversample_factor : int, optional
        Factor by which to oversample the wavelength grid

    Returns
    -------
    lambdas : ndarray
        Wavelengths at which to compute dispersed pixel values

    Notes
    -----
    It was found that the native wavelength spacing varies by a few percent or less
    across the detector for both NIRCam and NIRISS. This function has the capability to
    take in many x0, y0 at once and take the median to get the wavelengths,
    but typically it's okay to just put in any x0, y0 pair.
    """
    # Get x/y positions in the grism image corresponding to wmin and wmax:
    # Convert to x/y in the direct image frame
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, 1, order)
    # then convert to x/y in the grism image frame.
    xwmin, ywmin = imgxy_to_grismxy(x0_xy, y0_xy, wmin, order)
    xwmax, ywmax = imgxy_to_grismxy(x0_xy, y0_xy, wmax, order)
    dxw = xwmax - xwmin
    dyw = ywmax - ywmin

    # Create list of wavelengths on which to compute dispersed pixels
    dw = np.abs((wmax - wmin) / (dyw - dxw))
    dlam = np.median(dw / oversample_factor)
    lambdas = np.arange(wmin, wmax + dlam, dlam)
    return lambdas


def _disperse_onto_grism(x0_sky, y0_sky, sky_to_imgxy, imgxy_to_grismxy, lambdas, order):
    """
    Compute x/y positions in the grism image for the set of desired wavelengths.

    Parameters
    ----------
    x0_sky : ndarray
        RA of the input pixel position in direct image and segmentation map
    y0_sky : ndarray
        Dec of the input pixel position in direct image and segmentation map
    sky_to_imgxy : astropy model
        Transform from sky to image coordinates
    imgxy_to_grismxy : astropy model
        Transform from image to grism coordinates
    lambdas : ndarray
        Wavelengths at which to compute dispersed pixel values
    order : int
        Spectral order number

    Returns
    -------
    x0s : ndarray
        X coordinates of dispersed pixels in the grism image
    y0s : ndarray
        Y coordinates of dispersed pixels in the grism image
    lambdas : ndarray
        Wavelengths corresponding to each dispersed pixel
    """
    # x/y in image frame of grism image is the same for all wavelengths
    x0_sky = np.repeat(x0_sky[np.newaxis, :], len(lambdas), axis=0)
    y0_sky = np.repeat(y0_sky[np.newaxis, :], len(lambdas), axis=0)

    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, lambdas, order)
    del x0_sky, y0_sky

    # Convert to x/y in grism frame.
    lambdas = np.repeat(lambdas[:, np.newaxis], x0_xy.shape[1], axis=1)
    x0s, y0s = imgxy_to_grismxy(x0_xy, y0_xy, lambdas, order)
    # x0s, y0s now have shape (n_lam, n_pixels)
    return x0s, y0s, lambdas


def _collect_outputs_by_source(xs, ys, counts, source_ids_per_pixel):
    """
    Collect the dispersed pixel values into separate images for each source.

    Parameters
    ----------
    xs : ndarray
        X coordinates of dispersed pixels
    ys : ndarray
        Y coordinates of dispersed pixels
    counts : ndarray
        Count rates of dispersed pixels
    source_ids_per_pixel : int array
        Source IDs of the dispersed pixels

    Returns
    -------
    outputs_by_source : dict
        Dictionary containing dispersed images and bounds for each source ID
    """
    # Collect the outputs source-by-source
    outputs_by_source = {}
    source_ids = np.unique(source_ids_per_pixel)
    for this_sid in source_ids:
        this_sid_idx = source_ids_per_pixel == this_sid
        this_xs = xs[this_sid_idx]
        this_ys = ys[this_sid_idx]
        this_flxs = counts[this_sid_idx]
        if len(this_xs) == 0:
            continue

        img, bounds = _build_dispersed_image_of_source(this_xs, this_ys, this_flxs)
        outputs_by_source[this_sid] = {
            "bounds": bounds,
            "image": img,
        }
    return outputs_by_source


def _build_dispersed_image_of_source(x, y, flux):
    """
    Convert a flattened list of pixels to a 2-D grism image of that source.

    Parameters
    ----------
    x : ndarray
        X coordinates of pixels in the grism image
    y : ndarray
        Y coordinates of pixels in the grism image
    flux : ndarray
        Fluxes of pixels in the grism image

    Returns
    -------
    _type_
        _description_
    """
    minx = int(min(x))
    maxx = int(max(x))
    miny = int(min(y))
    maxy = int(max(y))
    a = sparse.coo_matrix(
        (flux, (y - miny, x - minx)), shape=(maxy - miny + 1, maxx - minx + 1)
    ).toarray()
    bounds = [minx, maxx, miny, maxy]
    return a, bounds


def disperse(
    xs,
    ys,
    fluxes,
    source_ids_per_pixel,
    order,
    wmin,
    wmax,
    sens_waves,
    sens_resp,
    direct_image_wcs,
    grism_wcs,
    naxis,
    oversample_factor=2,
    phot_per_lam=True,
):
    """
    Compute the dispersed image pixel values from the direct image.

    Parameters
    ----------
    xs : ndarray
        Flat array of X coordinates of pixels in the direct image
    ys : ndarray
        Flat array of Y coordinates of pixels in the direct image
    fluxes : ndarray
        Fluxes of the pixels in the direct image corresponding to xs, ys
    source_ids_per_pixel : int array
        Source IDs of the input pixels in the segmentation map
    order : int
        Spectral order number
    wmin : float
        Minimum wavelength for dispersed spectra
    wmax : float
        Maximum wavelength for dispersed spectra
    sens_waves : float array
        Wavelength array from photom reference file
    sens_resp : float array
        Response (flux calibration) array from photom reference file
    direct_image_wcs : WCS object
        WCS object for the direct image and segmentation map
    grism_wcs : WCS object
        WCS object for the grism image
    naxis : tuple
        Dimensions of the grism image (naxis[0], naxis[1])
    oversample_factor : int, optional
        Factor by which to oversample the wavelength grid
    phot_per_lam : bool
        If True, it is assumed that the photometric response is calibrated per wavelength,
        as is done for NIRCam, and therefore we need to multiply by dlam to un-do the calibration
        in this step.
        If False, it is assumed that the response is calibrated per pixel, and there is no need
        to multiply by dlam. This is what's done for NIRISS.

    Returns
    -------
    outputs_by_source : dict
        Dictionary containing dispersed images and bounds for each source ID
        in the specified spectral order.
    """
    n_input_sources = np.unique(source_ids_per_pixel).size
    log.debug(
        f"{mp.current_process()} dispersing {n_input_sources} "
        f"sources in order {order} with total number of pixels: {len(xs)}"
    )
    width = 1.0
    height = 1.0
    x0 = xs + 0.5 * width
    y0 = ys + 0.5 * height
    del xs, ys

    # Set up the transforms we need from the input WCS objects
    sky_to_imgxy = grism_wcs.get_transform("world", "detector")
    imgxy_to_grismxy = grism_wcs.get_transform("detector", "grism_detector")

    # We only need the x,y outputs of imgxy_to_grismxy
    # Making the number of outputs dynamic handles legacy WCS objects that did not pass
    # the x0, y0, and order through the transform unmodified like the current version does.
    n_outputs = len(imgxy_to_grismxy.outputs)
    imgxy_to_grismxy = imgxy_to_grismxy | Mapping((0, 1), n_inputs=n_outputs)

    # Find RA/Dec of the input pixel position in direct image
    x0_sky, y0_sky = direct_image_wcs(x0, y0, with_bounding_box=False)
    del x0, y0

    # native spacing does not change much over the detector, so just put in one x0, y0
    lambdas = _determine_native_wl_spacing(
        x0_sky[0],
        y0_sky[0],
        sky_to_imgxy,
        imgxy_to_grismxy,
        order,
        wmin,
        wmax,
        oversample_factor=oversample_factor,
    )
    nlam = len(lambdas)
    dlam = lambdas[1] - lambdas[0]

    x0s, y0s, lambdas = _disperse_onto_grism(
        x0_sky,
        y0_sky,
        sky_to_imgxy,
        imgxy_to_grismxy,
        lambdas,
        order,
    )
    del x0_sky, y0_sky

    # If none of the dispersed pixel indexes are within the image frame,
    # return a null result without wasting time doing other computations
    if x0s.min() >= naxis[0] or x0s.max() < 0 or y0s.min() >= naxis[1] or y0s.max() < 0:
        return

    source_ids_per_pixel = np.repeat(source_ids_per_pixel[np.newaxis, :], nlam, axis=0)
    fluxes = np.repeat(fluxes[np.newaxis, :], nlam, axis=0)

    # Discretize x and y coordinates to integer pixel values, keeping track of the fractional area
    # that each pixel contributes to the final grism image.
    # The resulting x, y coordinate pairs are non-unique: there are multiple wavelengths
    # that contribute to each pixel.
    padding = 1
    xs, ys, areas, index = get_clipped_pixels(x0s, y0s, padding, naxis[0], naxis[1], width, height)
    del x0s, y0s

    lambdas = np.take(lambdas, index)
    fluxes = np.take(fluxes, index)
    source_ids_per_pixel = np.take(source_ids_per_pixel, index)
    del index

    # compute 1D sensitivity array corresponding to list of wavelengths
    sens, no_cal = create_1d_sens(lambdas, sens_waves, sens_resp)
    del lambdas

    # Compute countrates for dispersed pixels.
    # The input direct image data is already photometrically calibrated,
    # so we need to basically apply a reverse flux calibration here.
    # Divide out the response values to convert from Mjy/sr to DN/s.
    # Note that the photom reference files are constructed differently for NIRCam and NIRISS
    # in the way they handle the dispersion: NIRCam uses per/wavelength, NIRISS uses per-pixel.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero")
        if phot_per_lam:
            # NIRCam case. Oversampling is accounted for by the spacing of dlam.
            counts = fluxes * areas * dlam * 1e4 / sens
        else:
            # NIRISS case. Oversampling must be handled directly to avoid double-counting.
            counts = fluxes * areas / (sens * oversample_factor)
    counts[no_cal] = 0.0  # set to zero where no flux cal info available
    del fluxes, areas, sens, dlam, no_cal

    outputs_by_source = _collect_outputs_by_source(xs, ys, counts, source_ids_per_pixel)
    del xs, ys, counts, source_ids_per_pixel
    n_out = len(outputs_by_source)
    log.debug(
        f"{mp.current_process()} finished order {order} with {n_out} "
        f"sources that overlap with the output frame "
        f"(out of {n_input_sources} input sources)"
    )
    return outputs_by_source
