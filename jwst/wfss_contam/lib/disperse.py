from scipy.interpolate import interp1d
import numpy as np
from .polyclip import polyclip


def dispersed_pixel(x0s, y0s, f0, order, wmin, wmax, seg_wcs, grism_wcs, ID, oversample_factor=2,
                    extrapolate_SED=False, xoffset=0, yoffset=0):
    """This function take a list of pixels and disperses them using the information contained
    in the grism image WCS object and returns a list of pixels and fluxes.

    Parameters
    ----------
    x0s: list
        A list of n x-coordinates.
    y0s: list
        A list of n y-coordinates.
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
    #img_to_grism = grism_wcs.get_transform('detector', 'grism_detector')
    sky_to_imgxy = grism_wcs.get_transform('world', 'detector')
    imgxy_to_grismxy = grism_wcs.get_transform('detector', 'grism_detector')
    offset_x = -imgxy_to_grismxy.offset_1
    offset_y = -imgxy_to_grismxy.offset_2

    if len(f0[0]) > 1:
        #print("f0:", f0, len(f0[0]), len(f0[1]))
        if extrapolate_SED is False:
            f = interp1d(f0[0], f0[1], fill_value=0., bounds_error=False)
        else:
            f = interp1d(f0[0], f0[1], fill_value="extrapolate", bounds_error=False)
    else:
        #print("f0: is not > 1")
        #f = lambda x: f0[1][0]
        def f(x):
            return f0[1][0]

    #print("x0s:", x0s)
    #print("y0s:", y0s)
    #print("f0:", f0)

    # Mean x/y of input pixel coords
    x0 = np.mean(x0s)
    y0 = np.mean(y0s)
    #print("x0, y0:", x0, y0)

    # deltas relative to mean x/y coords
    # typcially results in a list like -0.5, +0.5, +0.5, -0.5 when the input
    # simply brackets a single pixel
    dx0s = [t-x0 for t in x0s]
    dy0s = [t-y0 for t in y0s]
    #print("dx0s:", dx0s)
    #print("dy0s:", dy0s)

    # Get the x/y positions corresponding to wmin and wmax
    #xwmin, ywmin, _, _, _ = img_to_grism(x0, y0, wmin, order)
    #xwmax, ywmax, _, _, _ = img_to_grism(x0, y0, wmax, order)
    #xwmin, ywmin = img_to_grism(x0+offset_x, y0+offset_y, wmin, order)
    #xwmax, ywmax = img_to_grism(x0+offset_x, y0+offset_y, wmax, order)
    # RA/Dec of pixel position in segmentation map
    x0_sky, y0_sky = seg_wcs(x0, y0)
    # Pixel position in direct image frame of grism image
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, 1, order)
    # Now pixel positions in grism image
    xwmin, ywmin = imgxy_to_grismxy(x0_xy+offset_x, y0_xy+offset_y, wmin, order)
    xwmax, ywmax = imgxy_to_grismxy(x0_xy+offset_x, y0_xy+offset_y, wmax, order)
    #print("xwmin, ywmin:", xwmin, ywmin)
    #print("xwmax, ywmax:", xwmax, ywmax)
    dxw = xwmax - xwmin
    dyw = ywmax - ywmin

    # Compute the delta-wave per pixel
    dw = np.abs((wmax - wmin) / (dyw - dxw))
    #print("dw:", dw)

    # Use a natural wavelength scale or the wavelength scale of the input SED/spectrum,
    # whichever is smaller, divided by oversampling requested
    input_dlam = np.median(f0[0][1:] - f0[0][:-1])
    if input_dlam < dw:
        #print("input_dlam is < dw")
        dlam = input_dlam / oversample_factor
    else:
        #print("use computed dw")
        dlam = dw / oversample_factor

    lambdas = np.arange(wmin, wmax + dlam, dlam)
    n_lam = len(lambdas)
    #print("n_lam:", n_lam)

    # Compute lists of x/y positions in the grism image for
    # the set of desired wavelengths
    #x0s, y0s, _, _, _ = img_to_grism([x0]*n_lam, [y0]*n_lam, lambdas, [order]*n_lam)
    #x0s, y0s = img_to_grism([x0+offset_x]*n_lam, [y0+offset_y]*n_lam, lambdas, [order]*n_lam)
    # RA/Dec of segmentation map pixel positions
    x0_sky, y0_sky = seg_wcs([x0]*n_lam, [y0]*n_lam)
    # Pixel positions in direct image frame of grism image
    x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, lambdas, [order]*n_lam)
    x0s, y0s = imgxy_to_grismxy(x0_xy+offset_x, y0_xy+offset_y, lambdas, [order]*n_lam)
    #print("x0s:", x0s)
    #print("y0s:", y0s)

    padding = 1
    left = x0s.astype(np.int32) - padding
    right = x0s.astype(np.int32) + padding
    bottom = y0s.astype(np.int32) - padding
    top = y0s.astype(np.int32) + padding

    px = np.array([x0s+dx0s[0], x0s+dx0s[1], x0s+dx0s[2], x0s+dx0s[3]], dtype=np.float32).transpose().ravel()
    py = np.array([y0s+dy0s[0], y0s+dy0s[1], y0s+dy0s[2], y0s+dy0s[3]], dtype=np.float32).transpose().ravel()

    lams = np.array([[ll, 0, 0, 0] for ll in lambdas], dtype=np.float32).transpose().ravel()

    poly_inds = np.arange(0, (n_lam+1)*4, 4, dtype=np.int32)
    n_poly = len(x0s)
    n = len(lams)  # number of pixels we are "dropping", e.g. number of wav bins
    n *= 2

    index = np.zeros(n, dtype=np.int32)
    x = np.zeros(n, dtype=np.int32)
    y = np.zeros(n, dtype=np.int32)
    areas = np.zeros(n, dtype=np.float32)
    nclip_poly = np.array([0], np.int32)
    polyclip.polyclip_multi4(left, right, bottom, top, px, py, n_poly, poly_inds,
                             x, y, nclip_poly, areas, index)

    xs = x[0:nclip_poly[0]]
    ys = y[0:nclip_poly[0]]
    areas = areas[0:nclip_poly[0]]
    lams = np.take(lambdas, index)[0:len(xs)]
    # factor of 10000 because dlam is in micron and we want Angstrom with to apply f(lams)
    #counts = f(lams) * areas * sens(lams) * np.abs(dlam) * 10000.
    counts = f(lams) * areas * np.abs(dlam) * 10000.
    vg = (xs >= 0) & (ys >= 0)

    if len(xs[vg]) == 0:
        return np.array([0]), np.array([0]), np.array([0]), np.array([0]), np.array([0]), 0

    #print("xs:", xs[vg])
    #print("ys:", ys[vg])
    #print("areas:", areas[vg])
    return xs[vg], ys[vg], areas[vg], lams[vg], counts[vg], ID
