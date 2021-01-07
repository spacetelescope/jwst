# Python implementation: anand@stsci.edu 6 Mar 2013
# Algorithm: eelliott@stsci.edu -  Applied Optics, Vol 44, No. 8 10 March 2005
# Sabatke et al.
# Erin Elliott's analytical hexagon-aperture PSF, page 1361 equation 5
# Coordinate center at center of symmetry, and flat edge along xi axis
#     ---  eta
#   /     \ ^
#   \     / |
#     ---   -> xi
# hex(xi,eta) = g(xi,eta) + g(-xi,eta)
#
import logging
import numpy as np


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def g_eeAG(xi, eta, **kwargs):
    """
    Short Summary
    -------------
    Calculate the Fourier transform of one half of a hexagon that is bisected
    from one corner to its diametrically opposite corner.
    { DG: how does this compare to  g_eeGEN() ? }

    Parameters
    ----------
    xi: 2D float array
        hexagon's coordinate center at center of symmetry, along flat edge

    eta: 2D float array
        hexagon's coordinate center at center of symmetry, normal to xi

    c (optional, via **kwargs): tuple(float, float)
        coordinates of center

    pixel (optional, via **kwargs): float
        pixel scale

    d (optional, via **kwargs): float
        flat-to-flat distance across hexagon

    lambda (optional, via **kwargs): float
        wavelength

    minus: (optional, via **kwargs) boolean
        if set, use flipped sign of xi in calculation

    Returns
    -------
    g: 2D complex array
        Fourier transform of one half of a hexagon.
    """
    c = kwargs['c']
    pixel = kwargs['pixel']
    d = kwargs['d']
    lam = kwargs['lam']
    xi = (d / lam) * pixel * (xi - c[0])
    eta = (d / lam) * pixel * (eta - c[1])

    if kwargs['minus'] is True:
        xi = -1 * xi
    i = 1j
    Pi = np.pi

    g1 = np.exp(-i * Pi * (2 * eta / np.sqrt(3) + xi))
    g2 = (np.sqrt(3) * eta - 3 * xi)
    g3 = (np.exp(i * Pi * np.sqrt(3) * eta) - np.exp(i * Pi * \
            (4 * eta / np.sqrt(3) + xi)))
    g4 = (np.sqrt(3) * eta + 3 * xi)
    g5 = (np.exp(i * Pi * eta / np.sqrt(3)) - np.exp(i * Pi * xi))
    g6 = (4 * Pi * Pi * (eta * eta * eta - 3 * eta * xi * xi))
    g = g1 * (g2 * g3 + g4 * g5) / g6

    return g


def glimit(xi, eta, **kwargs):
    """
    Short Summary
    -------------
    Calculate the analytic limit of the Fourier transform of one half of the
        hexagon along eta=0.

    Parameters
    ----------
    xi: 2D float array
        hexagon's coordinate center at center of symmetry, along flat edge

    eta: 2D float array
        hexagon's coordinate center at center of symmetry, normal to xi;
        not currently used

    c (optional, via **kwargs): tuple(float, float)
        coordinates of center

    pixel (optional, via **kwargs): float
        pixel scale

    d (optional, via **kwargs): float
        flat-to-flat distance across hexagon

    lam: (optional, via **kwargs): float
        wavelength

    minus: (optional, via **kwargs) boolean
        if set, use flipped sign of xi in calculation

    Returns
    -------
    g: complex
        analytic limit of the Fourier transform of one half of the hexagon
        along eta=0
    """
    c = kwargs['c']
    pixel = kwargs['pixel']
    d = kwargs['d']
    lam = kwargs['lam']
    xi = (d / lam) * pixel * (xi - c[0])

    if kwargs['minus'] is True:
        xi = -1 * xi

    Pi = np.pi

    g1 = (np.exp(-1j * Pi * xi) / (2 * np.sqrt(3) * Pi * Pi * xi * xi))
    g2 = (-1 + 1j * Pi * xi + np.exp(1j * Pi * xi) - 2j * Pi * xi *
            np.exp(1J * Pi * xi))
    g = g1 * g2

    return g


def centralpix_limit():
    """
    Short Summary
    -------------
    Calculate the analytic limit of the Fourier transform of one half of the
    hexagon at the origin.

    Parameters
    ----------
    None

    Returns
    -------
    g: float
        analytic limit of the Fourier transform of one half of the hexagon
        at the origin.
    """
    g = np.sqrt(3) / 4.0

    return g


def mas2rad(mas):
    """
    Short Summary
    -------------
    Convert angle in milli arc-sec to radians

    Parameters
    ----------
    mas: float
        angle in milli arc-sec

    Returns
    -------
    rad: float
        angle in radians
    """
    rad = mas * (10**(-3)) / (3600 * 180 / np.pi)
    return rad


def hex_eeAG(s=(121, 121), c=None, d=0.80, lam=4.3e-6,
             pitch=mas2rad(65)):
    """
    Short Summary
    -------------
    Calculate the hexagonal hole Fourier transform by adding the transforms
    of the 2 symmetric parts.

    Parameters
    ----------
    s: (int,int) tuple
        size of hexagonal primary beam

    c: (float,float) tuple
        location of center of hexagonal primary beam

    d: float
        flat-to-flat distance across hexagon

    lam: float
        wavelength

    pitch: float
        sampling pitch in radians in image plane

    Returns
    -------
    np.abs(hex_complex): 2D float array
        hexagonal hole Fourier transform by adding the transforms
    """
    if c is None:
        c = float(s[0]) / 2.0 - 0.5, float(s[1]) / 2.0 - 0.5

    log.debug('hex_eeAG: center: %s, s: %s', c, s)

    h1 = np.fromfunction(g_eeAG, s, d=d, c=c, lam=lam, pixel=pitch, minus=False)
    h2 = np.fromfunction(g_eeAG, s, d=d, c=c, lam=lam, pixel=pitch, minus=True)
    hex_complex = h1 + h2

    # There will be a strip of NaNs down the middle (eta-axis)
    (xnan, ynan) = np.where(np.isnan(hex_complex))

    # The "yval" will be the same for all points;
    # loop over the xi values to replace NaN strip with limiting behavior.
    for index, val in enumerate(xnan):
        h1 = glimit(xnan[index], ynan[index], d=d, c=c, lam=lam, pixel=pitch,
                    minus=False)
        h2 = glimit(xnan[index], ynan[index], d=d, c=c, lam=lam, pixel=pitch,
                    minus=True)
        hex_complex[xnan[index], ynan[index]] = h1 + h2

    (xnan, ynan) = np.where(np.isnan(hex_complex))

    # Replace NaN strip with limiting behavior; the same for both halves
    hex_complex[xnan[:], ynan[:]] = 2.0 * centralpix_limit()

    if (log.getEffectiveLevel() <= logging.DEBUG):
        hr = hex_complex.real
        hi = hex_complex.imag
        log.debug('hex_eeAG: hr.min: %s, hr.mean: %s, hr.max: %s',
                   hr.min(), hr.mean(), hr.max())
        log.debug('hex_eeAG: hi.min: %s, hi.mean: %s, hi.max: %s',
                   hi.min(), hi.mean(), hi.max())

    return np.abs(hex_complex)
