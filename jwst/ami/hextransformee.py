#  Module for calculation of the hexagon-aperture PSFs
import numpy as np


def gfunction(xi, eta, **kwargs):
    """
    Fourier transform a half-hexagon.

    By half-hexagon, it is meant that
    the hexagon is bisected from one corner to its diametrically opposite corner.
    TODO: The kwargs are not optional, so there's no reason they should be kwargs at all.

    Parameters
    ----------
    xi : 2D float array
        Hexagon's coordinate center at center of symmetry, along flat edge
    eta : 2D float array
        Hexagon's coordinate center at center of symmetry, normal to xi
    **kwargs : dict
        Keyword arguments, which at present are NOT optional:
        c : tuple(float, float), required
            Coordinates of center
        pixel : float, required
            Pixel scale
        d : float, required
            Flat-to-flat distance across hexagon
        lam : float, required
            Wavelength
        affine2d : Affine2d object, required
        minus : bool, required
            If True, use flipped sign of xi in calculation

    Returns
    -------
    2D complex array
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

    TODO: Many of the inputs are not optional, and None is not a valid input.
    So they should not have defaults. The function signature
    will need to be re-ordered to have the required inputs first.
    TODO: This function modifies input parameter c (offsets it by a small amount)
    and there is no real reason for this. Could lead to bugs in the future
    and should be fixed.
    TODO: Why does this function force the center pixel to have value sqrt(3)/2
    but does not normalize the rest of the PSF accordingly?
    TODO: Functionality was written to handle the case that c is None, but it fails
    with TypeError: 'tuple' object does not support item assignment at the line
    c_adjust[0] = c[0] + eps_offset. This should either be fixed, or the option
    to have c=None should be removed.

    Parameters
    ----------
    s : (int,int) tuple, required
        Size of output hexagonal beam in pixels
    c : (float,float) tuple, required
        Location of center of hexagonal primary beam in pixels.
        If None, the center is assumed to be at the center of the array.
    d : float, required
        Flat-to-flat distance across hexagon in meters
    lam : float, required
        Wavelength of the observation in meters
    pitch : float, required
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
