from scipy.interpolate import interp1d
import numpy as np
from .polyclip import polyclip


def dispersed_pixel(x0s, y0s, f0, order, C, ID, oversample_factor=2, extrapolate_SED=False, xoffset=0, yoffset=0):
    """This function take a list of pixels and disperses them using the information contained
    in a GRISMCONF file, and returns a list of pixel and fluxes.

    Parameters
    ----------
    x0s: list
        A list of n x-coordinates.
    y0s: list
        A list of n y-coordinates.
    f0: list
        A list of n flux (flam) for each of the pixels contained in x0s,y0s.
    order: str
        The name of the spectral order to disperse.
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
    if len(f0[0]) > 1:
        #print("f0:",f0,len(f0[0]),len(f0[1]))
        if extrapolate_SED is False:
            f = interp1d(f0[0], f0[1], fill_value=0., bounds_error=False)
        else:
            f = interp1d(f0[0], f0[1], fill_value="extrapolate", bounds_error=False)
    else:
        #f = lambda x: f0[1][0]
        def f(x): return f0[1][0]

    # Retrieve sensitivity calibration reference data for the spectral order
    # Old aXe style Config ref file. Needs updating
    # to use JWST photom spectral flux calibration curves
    sens = C.SENS[order]

    # Mean x/y of input pixel list
    x0 = np.mean(x0s)
    y0 = np.mean(y0s)

    dx0s = [t-x0 for t in x0s]
    dy0s = [t-y0 for t in y0s]

    # Figuring out a few things about size of order, dispersion and wavelengths to use
    # Old aXe style Config info. Needs updating to use WCS waverange info.
    wmin = C.WRANGE[order][0]
    wmax = C.WRANGE[order][1]

    # Old aXe style Config info. Need to update.
    t0 = C.INVDISPL(order, x0+xoffset, y0+yoffset, wmin)
    t1 = C.INVDISPL(order, x0+xoffset, y0+yoffset, wmax)

    # Old aXe style Config info. Need to update.
    dx0 = C.DISPX(order, x0+xoffset, y0+yoffset, t0) - C.DISPX(order, x0+xoffset, y0+yoffset, t1)
    dx1 = C.DISPY(order, x0+xoffset, y0+yoffset, t0) - C.DISPY(order, x0+xoffset, y0+yoffset, t1)

    dw = np.abs((wmax - wmin) / (dx1 - dx0))

    # Use a natural wavelength scale or the wavelength scale of the input SED/spectrum,
    # whichever is smaller, divided by oversampling requested
    input_dlam = np.median(f0[0][1:] - f0[0][:-1])
    if input_dlam < dw:
        dlam = input_dlam / oversample_factor
    else:
        dlam = dw / oversample_factor

    lambdas = np.arange(wmin, wmax + dlam, dlam)

    # Old aXe style Config info. Need to update.
    dS = C.INVDISPL(order, x0 + xoffset, y0 + yoffset, lambdas)

    m = len(lambdas)

    # Old aXe style Config info. Need to update.
    dXs = C.DISPX(order, x0 + xoffset, y0 + yoffset, dS)
    dYs = C.DISPY(order, x0 + xoffset, y0 + yoffset, dS)

    x0s = x0 + dXs
    y0s = y0 + dYs

    padding = 1
    left = x0s.astype(np.int32) - padding
    right = x0s.astype(np.int32) + padding
    bottom = y0s.astype(np.int32) - padding
    top = y0s.astype(np.int32) + padding

    px = np.array([x0s+dx0s[0], x0s+dx0s[1], x0s+dx0s[2], x0s+dx0s[3]], dtype=np.float32).transpose().ravel()
    py = np.array([y0s+dy0s[0], y0s+dy0s[1], y0s+dy0s[2], y0s+dy0s[3]], dtype=np.float32).transpose().ravel()

    lams = np.array([[ll, 0, 0, 0] for ll in lambdas], dtype=np.float32).transpose().ravel()

    poly_inds = np.arange(0, (m+1)*4, 4, dtype=np.int32)

    n_poly = len(x0s)

    n = len(lams)  # number of pixels we are "dropping", e.g. number of wav bins

    n = 2 * n

    index = np.zeros(n, dtype=np.int32)
    x = np.zeros(n, dtype=np.int32)
    y = np.zeros(n, dtype=np.int32)
    areas = np.zeros(n, dtype=np.float32)
    nclip_poly = np.array([0], np.int32)
    polyclip.polyclip_multi4(left, right, bottom, top, px, py, n_poly, poly_inds, x, y, nclip_poly, areas, index)

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

    return xs[vg], ys[vg], areas[vg], lams[vg], counts[vg], ID
