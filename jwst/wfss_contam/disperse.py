import numpy as np
import warnings
import multiprocessing as mp

from scipy import sparse

from jwst.lib.winclip import get_clipped_pixels
from jwst.wfss_contam.sens1d import create_1d_sens

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def disperse(
    xs,
    ys,
    fluxes,
    source_ids_per_pixel,
    order,
    pivlam,
    wmin,
    wmax,
    sens_waves,
    sens_resp,
    seg_wcs,
    grism_wcs,
    naxis,
    oversample_factor=2,
    xoffset=0,
    yoffset=0,
):
    """
    Compute dispersed pixel values for sources identified in the segmentation map.

    Parameters
    ----------
    xs : float array
        X coordinates of pixels in the segmentation map
    ys : float array
        Y coordinates of pixels in the segmentation map
    fluxes : float array
        Fluxes of sources in the segmentation map
    source_ids_per_pixel : int array
        Source IDs corresponding to each pixel in the segmentation map
    order : int
        Spectral order number to process
    pivlam : float
        Pivot wavelength for the spectral order
    wmin : float
        Minimum wavelength for dispersed spectra
    wmax : float
        Maximum wavelength for dispersed spectra
    sens_waves : float array
        Wavelength array from photom reference file
    sens_resp : float array
        Response (flux calibration) array from photom reference file
    seg_wcs : WCS object
        WCS object for the segmentation map
    grism_wcs : WCS object
        WCS object for the grism image
    naxis : tuple
        Dimensions of the grism image (naxis[0], naxis[1])
    oversample_factor : int, optional
        Factor by which to oversample the wavelength grid (default is 2)
    xoffset : float, optional
        X offset to apply to the dispersed pixel positions (default is 0)
    yoffset : float, optional
        Y offset to apply to the dispersed pixel positions (default is 0)

    Returns
    -------
    outputs_by_source : dict
        Dictionary containing dispersed images and bounds for each source ID
        in the specified spectral order.
    """
    log.debug(
        f"{mp.current_process()} dispersing {np.unique(source_ids_per_pixel).size} "
        f"sources in order {order}"
    )
    log.debug(f"total number of pixels: {len(xs)}")
    width = 1.0
    height = 1.0
    x0 = xs + 0.5 * width
    y0 = ys + 0.5 * height

    # Compute the WCS transforms
    # Setup the transforms we need from the input WCS objects
    sky_to_imgxy = grism_wcs.get_transform("world", "detector")
    imgxy_to_grismxy = grism_wcs.get_transform("detector", "grism_detector")

    # Get x/y positions in the grism image corresponding to wmin and wmax:
    # Start with RA/Dec of the input pixel position in segmentation map,
    x0_sky, y0_sky = seg_wcs(x0, y0)
    # then convert to x/y in the direct image frame corresponding
    # to the grism image,
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, 1, order)
    # then finally convert to x/y in the grism image frame
    xwmin, ywmin = imgxy_to_grismxy(x0_xy + xoffset, y0_xy + yoffset, wmin, order)
    xwmax, ywmax = imgxy_to_grismxy(x0_xy + xoffset, y0_xy + yoffset, wmax, order)
    dxw = xwmax - xwmin
    dyw = ywmax - ywmin

    # Create list of wavelengths on which to compute dispersed pixels
    lams = np.array([pivlam] * len(fluxes))
    dw = np.abs((wmax - wmin) / (dyw - dxw))
    dlam = np.median(dw / oversample_factor)  # TODO: validate that just taking median is ok
    lambdas = np.arange(wmin, wmax + dlam, dlam)
    n_lam = len(lambdas)

    # Compute lists of x/y positions in the grism image for
    # the set of desired wavelengths:
    # x/y in image frame of grism image is the same for all wavelengths
    x0_sky = np.repeat(x0_sky[np.newaxis, :], n_lam, axis=0)
    y0_sky = np.repeat(y0_sky[np.newaxis, :], n_lam, axis=0)
    source_ids_per_pixel = np.repeat(source_ids_per_pixel[np.newaxis, :], n_lam, axis=0)
    fluxes = np.repeat(fluxes[np.newaxis, :], n_lam, axis=0)
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, lambdas, order)

    # Convert to x/y in grism frame.
    # lambdas needs same shape as x0_xy to be indexed by np.take below
    lambdas = np.repeat(lambdas[:, np.newaxis], x0_xy.shape[1], axis=1)
    x0s, y0s = imgxy_to_grismxy(x0_xy + xoffset, y0_xy + yoffset, lambdas, order)
    # x0s, y0s now have shape (n_lam, n_pixels)

    # Compute arrays of dispersed pixel locations and areas
    padding = 1
    # If none of the dispersed pixel indexes are within the image frame,
    # return a null result without wasting time doing other computations
    if x0s.min() >= naxis[0] or x0s.max() < 0 or y0s.min() >= naxis[1] or y0s.max() < 0:
        # log.info(f"No dispersed pixels within image frame for order {order}.")
        return

    xs, ys, areas, index = get_clipped_pixels(x0s, y0s, padding, naxis[0], naxis[1], width, height)
    lams = np.take(lambdas, index)
    fluxes = np.take(fluxes, index)

    # compute 1D sensitivity array corresponding to list of wavelengths
    # TODO: what wavelength unit does this expect?
    sens, no_cal = create_1d_sens(lams, sens_waves, sens_resp)

    # Compute countrates for dispersed pixels. Note that dispersed pixel
    # values are naturally in units of physical fluxes, so we divide out
    # the sensitivity (flux calibration) values to convert to units of
    # countrate (DN/s).
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero")
        counts = fluxes * lams * areas / (sens * oversample_factor)
    counts[no_cal] = 0.0  # set to zero where no flux cal info available

    # keep track of source IDs for each dispersed pixel
    dispersed_source_ids = np.take(source_ids_per_pixel, index)

    # Build the dispersed image and make a slit model for each source
    outputs_by_source = {}
    source_ids = np.unique(dispersed_source_ids)
    for this_sid in source_ids:
        this_sid_idx = dispersed_source_ids == this_sid
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
    log.debug(
        f"{mp.current_process()} finished order {order} with {len(outputs_by_source)} "
        "sources in the output frame"
    )
    return outputs_by_source


def _build_dispersed_image_of_source(x, y, flux):
    minx = int(min(x))
    maxx = int(max(x))
    miny = int(min(y))
    maxy = int(max(y))
    a = sparse.coo_matrix(
        (flux, (y - miny, x - minx)), shape=(maxy - miny + 1, maxx - minx + 1)
    ).toarray()
    bounds = [minx, maxx, miny, maxy]
    return a, bounds
