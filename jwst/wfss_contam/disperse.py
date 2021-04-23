from scipy.interpolate import interp1d
import numpy as np
from ..lib.winclip import get_clipped_pixels
import sys, time


def dispersed_pixel(x0, y0, width, height, lams, flxs, order, wmin, wmax,
                    seg_wcs, grism_wcs, ID, naxis, oversample_factor=2,
                    extrapolate_SED=False, xoffset=0, yoffset=0):
    """This function take a list of pixels and disperses them using the information contained
    in the grism image WCS object and returns a list of pixels and fluxes.

    Parameters
    ----------
    x0: list
        A list of n x-coordinates of the centers of the pixels.
    y0: list
        A list of n y-coordinates of the centers of the pixels.
    width: int
        Width of the pixels to be dispersed.
    height: int
        Width of the pixels to be dispersed.
    f0: list
        A list of n flux (flam) for each of the pixels contained in x0s,y0s.
        The entries in the list are wavelength and flux.
    order: int
        The spectral order to disperse.
    seg_wcs: WCS object
        The WCS object of the segmentation map.
    grism_wcs: WCS object
        The WCS object of the grism image.
    ID: int
        The ID of the object this is for.
    oversample_factor: int
        The amount of oversampling required above that of the input spectra or natural dispersion,
        whichever is smaller. Default=2.
    extrapolate_SED: bol
        Whether to allow for the SED of the object to be extrapolated when it does not fully cover the
        needed wavelength range. Default if False.
    xoffset int
        A pixel offset to apply when computing the dispersion (for padded images for example)
    yoffset int
        A pixel offset to apply when computing the dispersion (for padded images for example)

    Output
    ------
    xs: array
        A list of x-coordinates dispersed pixel
    ys: array
        A list of y-coordinates dispersed pixel
    areas: array
        A list of the areas of the incident pixel thatm when dispersed falls on the dispersed pixel
    lams: array
        A list of the wavelength of each of the dispersed pixels
    counts: array
        A list of counts for each of the dispersed pixels
    ID: array
        A list containing the ID. Returned for bookkeeping convenience.
    """

    # Setup the transforms we need from the input WCS objects
    # img_to_grism = grism_wcs.get_transform('detector', 'grism_detector')
    sky_to_imgxy = grism_wcs.get_transform('world', 'detector')
    imgxy_to_grismxy = grism_wcs.get_transform('detector', 'grism_detector')
    offset_x = -imgxy_to_grismxy.offset_1
    offset_y = -imgxy_to_grismxy.offset_2

    if len(lams) > 1:
        if extrapolate_SED is False:
            f = interp1d(lams, flxs, fill_value=0., bounds_error=False)
        else:
            f = interp1d(lams, flxs, fill_value="extrapolate", bounds_error=False)
    else:
        def f(x):
            return flxs[0]

    # Get the x/y positions corresponding to wmin and wmax
    # xwmin, ywmin, _, _, _ = img_to_grism(x0, y0, wmin, order)
    # xwmax, ywmax, _, _, _ = img_to_grism(x0, y0, wmax, order)
    # xwmin, ywmin = img_to_grism(x0+offset_x, y0+offset_y, wmin, order)
    # xwmax, ywmax = img_to_grism(x0+offset_x, y0+offset_y, wmax, order)
    # RA/Dec of pixel position in segmentation map
    x0_sky, y0_sky = seg_wcs(x0, y0)
    # Pixel position in direct image frame of grism image
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, 1, order)
    # Now pixel positions in grism image
    xwmin, ywmin = imgxy_to_grismxy(x0_xy + offset_x, y0_xy + offset_y, wmin, order)
    xwmax, ywmax = imgxy_to_grismxy(x0_xy + offset_x, y0_xy + offset_y, wmax, order)
    # print("xwmin, ywmin:", xwmin, ywmin)
    # print("xwmax, ywmax:", xwmax, ywmax)
    dxw = xwmax - xwmin
    dyw = ywmax - ywmin

    # Compute the delta-wave per pixel
    dw = np.abs((wmax - wmin) / (dyw - dxw))
    # print("dw:", dw)

    # Use a natural wavelength scale or the wavelength scale of the input SED/spectrum,
    # whichever is smaller, divided by oversampling requested
    input_dlam = np.median(lams[1:] - lams[:-1])
    if input_dlam < dw:
        # print("input_dlam is < dw")
        dlam = input_dlam / oversample_factor
    else:
        # print("use computed dw")
        dlam = dw / oversample_factor

    lambdas = np.arange(wmin, wmax + dlam, dlam)
    n_lam = len(lambdas)

    # Compute lists of x/y positions in the grism image for
    # the set of desired wavelengths
    # x0s, y0s, _, _, _ = img_to_grism([x0]*n_lam, [y0]*n_lam, lambdas, [order]*n_lam)
    # x0s, y0s = img_to_grism([x0+offset_x]*n_lam, [y0+offset_y]*n_lam, lambdas, [order]*n_lam)
    # RA/Dec of segmentation map pixel positions
    x0_sky, y0_sky = seg_wcs([x0] * n_lam, [y0] * n_lam)
    # Pixel positions in direct image frame of grism image
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, lambdas, [order] * n_lam)
    x0s, y0s = imgxy_to_grismxy(x0_xy + offset_x, y0_xy + offset_y, lambdas, [order] * n_lam)

    padding = 1
    xs, ys, areas, index = get_clipped_pixels(
        x0s, y0s,
        padding,
        naxis[0], naxis[1],
        width, height
    )
    lams = np.take(lambdas, index)
    counts = f(lams) * areas * np.abs(dlam)

    if xs.size <= 1:
        return None

    return xs, ys, areas, lams, counts, ID
