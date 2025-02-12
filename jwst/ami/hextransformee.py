#  Module for calculation of the hexagon-aperture PSFs
import numpy as np


def gfunction(xi, eta, **kwargs):
    """
    Fourier transform a half-hexagon.

    By half-hexagon, it is meant that
    the hexagon is bisected from one corner to its diametrically opposite corner.

    Parameters
    ----------
    xi : 2D float array
        Hexagon's coordinate center at center of symmetry, along flat edge

    eta : 2D float array
        Hexagon's coordinate center at center of symmetry, normal to xi

    **kwargs : dict
        Keyword arguments
        c (optional, via **kwargs): tuple(float, float)
            coordinates of center

        pixel (optional, via **kwargs): float
            pixel scale

        d (optional, via **kwargs): float
            flat-to-flat distance across hexagon

        lam (optional, via **kwargs): float
            wavelength

        minus: (optional, via **kwargs) boolean
            if set, use flipped sign of xi in calculation

    Returns
    -------
    g*affine2d.distortphase(xi,eta): 2D complex array
        Fourier transform of one half of a hexagon.
    """
    c = kwargs["c"]
    pixel = kwargs["pixel"]
    d = kwargs["d"]
    lam = kwargs["lam"]
    xi = (d / lam) * pixel * (xi - c[0])
    eta = (d / lam) * pixel * (eta - c[1])
    affine2d = kwargs["affine2d"]

    i = 1j
    pi = np.pi

    xip, etap = affine2d.distort_f_args(xi, eta)

    if kwargs["minus"] is True:
        xip = -1 * xip

    g = (
        np.exp(-i * pi * (2 * etap / np.sqrt(3) + xip))
        * (
            (np.sqrt(3) * etap - 3 * xip)
            * (np.exp(i * pi * np.sqrt(3) * etap) - np.exp(i * pi * (4 * etap / np.sqrt(3) + xip)))
            + (np.sqrt(3) * etap + 3 * xip)
            * (np.exp(i * pi * etap / np.sqrt(3)) - np.exp(i * pi * xip))
        )
        / (4 * pi * pi * (etap * etap * etap - 3 * etap * xip * xip))
    )

    return g * affine2d.distortphase(xi, eta)


def hextransform(s=None, c=None, d=None, lam=None, pitch=None, affine2d=None):
    """
    Calculate the complex array analytical transform of a (distorted if necessary) hexagon.

    Parameters
    ----------
    s : (int,int) tuple
        Size of hexagonal primary beam

    c : (float,float) tuple
        Location of center of hexagonal primary beam

    d : float
        Flat-to-flat distance across hexagon

    lam : float
        Vavelength

    pitch : float
        Sampling pitch in radians in image plane

    affine2d : Affine2d object
        Distortion object

    Returns
    -------
    hex_complex : 2D complex array
        Complex array analytical transform of a hexagon
    """
    if c is None:
        c = (float(s[0]) / 2.0 - 0.5, float(s[1]) / 2.0 - 0.5)

    # deal with central pixel singularity:
    c_adjust = c

    eps_offset = 1.0e-8

    # The fractional part of the offsets:
    d0, d1 = (c[0] - int(c[0]), c[1] - int(c[1]))

    # Are they very small (i.e. 'almost exactly' centered on 0)?
    if abs(d0) < 0.5 * eps_offset:  # might have the singular central pixel here
        c_adjust[0] = c[0] + eps_offset

    if abs(d1) < 0.5 * eps_offset:  # might have the singular central pixel here
        c_adjust[1] = c[1] + eps_offset

    hex_complex = np.fromfunction(
        gfunction,
        s,
        d=d,
        c=c_adjust,
        lam=lam,
        pixel=pitch,
        affine2d=affine2d,
        minus=False,
    ) + np.fromfunction(
        gfunction,
        s,
        d=d,
        c=c_adjust,
        lam=lam,
        pixel=pitch,
        affine2d=affine2d,
        minus=True,
    )
    hex_complex[int(c[0]), int(c[1])] = np.sqrt(3) / 2.0

    return hex_complex
