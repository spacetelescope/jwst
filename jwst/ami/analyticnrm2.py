#! /usr/bin/env  python
# Heritage mathematia nb from Alex & Laurent
# Python by Alex Greenbaum & Anand Sivaramakrishnan Jan 2013
# updated May 2013 to include hexagonal envelope

from . import hexee

import logging
import numpy as np
import scipy.special
from . import leastsqnrm

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def Jinc(x, y):
    """
    Short Summary
    -------------
    Compute 2d Jinc for given coordinates

    Parameters
    ----------
    x,y: floats
        input coordinates

    Returns
    -------
    jinc_2d: float array
        2d Jinc at the given coordinates, with NaNs replaced by pi/4.
    """
    R = (Jinc.d / Jinc.lam) * Jinc.pitch *  \
           np.sqrt((x - Jinc.offx)*(x - Jinc.offx) + \
           (y - Jinc.offy)*(y - Jinc.offy))

    jinc_2d = leastsqnrm.replacenan(scipy.special.jv(1, np.pi * R)/(2.0 * R))

    return jinc_2d


def phasor(kx, ky, hx, hy, lam, phi, pitch):
    """
    Short Summary
    -------------
    Calculate wavefront for a single hole ??

    Parameters
    ----------
    kx, ky: float
        image plane coords in units of sampling pitch (oversampled, or not)

    hx, hy: float
        hole centers in meters

    lam: float
        wavelength

    phi: float
        distance of fringe from hole center in units of waves

    pitch: float
        sampling pitch in radians in image plane

    Returns
    -------
    phasor: complex
        Calculate wavefront for a single hole
    """
    return np.exp(-2 * np.pi * 1j * ((pitch * kx * hx + pitch * ky * hy)
               / lam + (phi / lam)))


def interf(kx, ky):
    """
    Short Summary
    -------------
    Calculate interference for all holes.

    Parameters
    ----------
    kx, ky: float, float
        x-component and y-component of image plane (spatial frequency) vector

    Returns
    -------
    interference: 2D complex array
        interference for all holes
    """
    interference = 0j
    for hole, ctr in enumerate(interf.ctrs):
        interference += phasor((kx - interf.offx), (ky - interf.offy),
                               ctr[0], ctr[1], interf.lam,
                               interf.phi[hole], interf.pitch)

    return interference


def ASF(pixel, fov, oversample, ctrs, d, lam, phi, centering=(0.5, 0.5)):
    """
    Short Summary
    -------------
    Calculate the Amplitude Spread Function (a.k.a. image plane complex
    amplitude) for a circular aperture

    Parameters
    ----------
    pixel: float
        pixel scale

    fov: integer
        number of detector pixels on a side

    oversample: integer
        oversampling factor

    ctrs: float, float
        coordinates of hole center

    d: float
        hole diameter

    lam: float
        wavelength

    phi: float
        distance of fringe from hole center in units of waves

    centering: string
        if set to 'PIXELCENTERED' or unspecified, the offsets will be set to
        (0.5,0.5); if set to 'PIXELCORNER', the offsets will be set to
        (0.0,0.0).

    Returns
    -------
    asf: 2D complex array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a circular aperture
    """
    if centering == 'PIXELCENTERED':
        off_x = 0.5
        off_y = 0.5
    elif centering == 'PIXELCORNER':
        off_x = 0.0
        off_y = 0.0
    else:
        off_x, off_y = centering

    # log.debug('ASF centering %s:', centering)
    # log.debug('ASF offsets %s %s:', off_x, off_y)

    # Jinc parameters
    Jinc.lam = lam
    Jinc.offx = oversample * fov / 2.0 - off_x # in pixels
    Jinc.offy = oversample * fov / 2.0 - off_y
    Jinc.pitch = pixel / float(oversample)
    Jinc.d = d

    primarybeam = np.fromfunction(Jinc, (int((oversample * fov)),
                                         int((oversample * fov))))
    primarybeam = primarybeam.transpose()

    # interference terms' parameters
    interf.lam = lam
    interf.offx = oversample * fov / 2.0 - off_x # in pixels
    interf.offy = oversample * fov / 2.0 - off_y
    interf.pitch = pixel / float(oversample)
    interf.ctrs = ctrs
    interf.phi = phi

    fringing = np.fromfunction(interf, (int((oversample * fov)),
                                        int((oversample * fov))))
    fringing = fringing.transpose()

    asf = primarybeam * fringing

    return asf


def ASFfringe(pixel, fov, oversample, ctrs, d, lam, phi, centering=(0.5, 0.5)):
    """
    Short Summary
    -------------
    Amplitude Spread Function (a.k.a. image plane complex amplitude)
    for a fringe

    Parameters
    ----------
    pixel: float
        pixel scale

    fov: integer
        number of detector pixels on a side

    oversample: integer
        oversampling factor

    ctrs: 2D float array
        centers of holes

    d: float
        hole diameter

    lam: float
        wavelength

    phi: float
        distance of fringe from hole center in units of waves

    centering: string
        if set to 'PIXELCENTERED' or unspecified, the offsets will be set to
        (0.5,0.5); if set to 'PIXELCORNER', the offsets will be set to
        (0.0,0.0).

    Returns
    -------
    fringing: 2D complex array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a fringe
    """
    if centering == 'PIXELCENTERED':
        off_x = 0.5
        off_y = 0.5
    elif centering == 'PIXELCORNER':
        off_x = 0.0
        off_y = 0.0
    else:
        off_x, off_y = centering

    # log.debug('ASFfringe centering %s:', centering)
    # log.debug('ASFfringe offsets %s %s:', off_x, off_y)

    # Jinc parameters
    Jinc.lam = lam
    Jinc.offx = oversample * fov / 2.0 - off_x # in pixels
    Jinc.offy = oversample * fov / 2.0 - off_y
    Jinc.pitch = pixel / float(oversample)
    Jinc.d = d

    # interference terms' parameters
    interf.lam = lam
    interf.offx = oversample * fov / 2.0 - off_x # in pixels
    interf.offy = oversample * fov / 2.0 - off_y
    interf.pitch = pixel / float(oversample)
    interf.ctrs = ctrs
    interf.phi = phi

    fringing = np.fromfunction(interf, (int((oversample * fov)),
                                        int((oversample * fov))))
    fringing = fringing.transpose()

    return fringing


def ASFhex(pixel, fov, oversample, ctrs, d, lam, phi, centering='PIXELCENTERED'):
    """
    Short Summary
    -------------
    Amplitude Spread Function (a.k.a. image plane complex amplitude)
    for a hexagonal aperture

    Parameters
    ----------
    pixel: float
        pixel scale

    fov: integer
        number of detector pixels on a side

    oversample: integer
        oversampling factor

    ctrs: 2D float array
        centers of holes

    d: float
        flat-to-flat distance across hexagon

    lam: float
        wavelength

    phi: float
        distance of fringe from hole center in units of waves

    centering: string
        type of centering

    Returns
    -------
    asf: 2D complex array
        Amplitude Spread Function (a.k.a. image plane complex amplitude) for
        a hexagonal aperture
    """
    # log.debug('centering: %s', centering)

    if centering == 'PIXELCENTERED':
        off_x = 0.5
        off_y = 0.5
    elif centering == 'PIXELCORNER':
        off_x = 0.0
        off_y = 0.0
    else:
        off_x, off_y = centering

    #Hex kwargs
    offx = (float(oversample * fov) / 2.0) - off_x # in pixels
    offy = (float(oversample * fov) / 2.0) - off_y

    # log.debug('ASF offsets for x and y in pixels: %s %s', offx, offy)
    # log.debug('ASF centering:%s', centering)

    pitch = pixel / float(oversample)

    # interference terms' parameters
    interf.lam = lam
    interf.offx = (oversample * fov) / 2.0 - off_x # in pixels
    interf.offy = (oversample * fov) / 2.0 - off_y
    interf.pitch = pixel / float(oversample)
    interf.ctrs = ctrs
    interf.phi = phi

    primarybeam = hexee.hex_eeAG(s=(oversample * fov, oversample * fov),
                                 c=(offx, offy), d=d, lam=lam, pitch=pitch)

    fringing = np.fromfunction(interf, (int((oversample * fov)),
                                        int((oversample * fov))))
    fringing = fringing.transpose()

    asf = primarybeam * fringing

    return asf


def PSF(pixel, fov, oversample, ctrs, d, lam, phi, centering='PIXELCENTERED',
        shape='circ'):
    """
    Short Summary
    -------------
    Calculate the PSF for the requested shape

    Parameters
    ----------
    pixel: float
        pixel scale

    fov: integer
        number of detector pixels on a side

    oversample: integer
        oversampling factor

    ctrs: 2D float array
        centers of holes

    d: float
        hole diameter for 'circ'; flat-to-flat distance across for 'hex'

    lam: float
        wavelength

    phi: float
        distance of fringe from hole center in units of waves

    centering: string
        type of centering

    shape: string
        shape of hole; possible values are 'circ', 'hex', and 'fringe'

    Returns
    -------
    PSF - 2D float array
    """
    if shape == 'circ':
        asf = ASF(pixel, fov, oversample, ctrs, d, lam, phi, centering)
    elif shape == 'hex':
        asf = ASFhex(pixel, fov, oversample, ctrs, d, lam, phi, centering)
    elif shape == 'fringe': # Alex: "not needed,only used for visualization"
        asf = ASFfringe(pixel, fov, oversample, ctrs, d, lam, phi, centering)
    else:
        log.critical('Pupil shape %s not supported', shape)

    # log.debug('-----------------')
    # log.debug(' PSF Parameters: ')
    # log.debug('-----------------')
    # log.debug('pixel: %s, fov: %s, oversampling: %s', pixel, fov, oversample)
    # log.debug('d: %s, wavelength: %s, pistons: %s, shape: %s', d, lam, phi,
    #           shape)

    PSF_ = asf * asf.conj()

    return PSF_.real
