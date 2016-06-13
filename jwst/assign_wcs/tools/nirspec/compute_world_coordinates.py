"""
A simple tool to read in the output of extract2d and
apply the WCS transforms to all pixels in a slit.
For each slit it writes the results as a cube with three planes
(wavelength, ra, dec) in a separate fits extension.
The file is saved with an suffix "world_coordinates".

Requested by the NIRSPEC team.

"""
from __future__ import absolute_import, division, unicode_literals, print_function

import numpy as np
from astropy.io import fits
from jwst import datamodels


def compute_world_coordinates(fname, output=None):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    fixed slit observation. For each slit the output is a cube
    of three planes (wavelengh, ra, dec), saved as a separate
    FITS extension.

    Parameters
    ----------
    fname : str
        The name of a file with extracted slits, i.e. the output
        of extract2d.
    output : str
        The name of the output file. If None the root of the input
        file is used with an extension world_coordinates.

    Examples
    --------
    >>> compute_world_coordinates('nrs1_fixed_assign_wcs_extract_2d.fits')

    """
    model = models.MultiSlitModel(fname)
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    for slit in model.slits:
        ysize, xsize = slit.data.shape
        y, x = np.mgrid[: ysize, : xsize]
        ra, dec, lam = slit.meta.wcs(x, y)
        world_coordinates = np.array([lam, ra, dec])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = 'ote_x, arcsec'
        imhdu.header['PLANE3'] = 'ote_y, arcsec'
        imhdu.header['SLIT'] = slit.name
        hdulist.append(imhdu)
    if output is not None:
        base, ext = os.path.splitext(output)
        if ext != "fits":
            ext = "fits"
        if not base.endswith('world_coordinates'):
            "".join([base, '_world_coordinates'])
        "".join([base, ext])
    else:
        root = model.meta.filename.split('_')
        output = "".join([root[0], '_world_coordinates', '.fits'])
    hdulist.writeto(output)
    del hdulist
    model.close()


def compute_msa_coordinates(fname, output=None):
    """
    Computes wavelengths, and relative MSA coordinates of a NIRSPEC
    fixed slit observation. For each slit the output is a cube
    of three planes (wavelengh, MSA_X, MSA_Y), saved as a separate
    FITS extension.

    Parameters
    ----------
    fname : str
        The name of a file with extracted slits, i.e. the output
        of extract2d.
    output : str
        The name of the output file. If None the root of the input
        file is used with a suffix "msa".

    Examples
    --------
    >>> compute_msa_coordinates('nrs1_fixed_assign_wcs_extract_2d.fits')

    """
    model = models.MultiSlitModel(fname)
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'msa'
    hdulist.append(phdu)
    for slit in model.slits:
        ysize, xsize = slit.data.shape
        y, x = np.mgrid[: ysize, : xsize]
        det2msa = slit.meta.wcs.get_transform('detector', 'msa')
        x, y, lam = det2msa(x, y)
        msa_coordinates = np.array([lam, x, y])
        imhdu = fits.ImageHDU(data=msa_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = 'msa_x, relative to center of slit (0, 0)'
        imhdu.header['PLANE3'] = 'msa_y, relative to center of slit (0, 0)'
        imhdu.header['SLIT'] = slit.name
        hdulist.append(imhdu)
    if output is not None:
        base, ext = os.path.splitext(output)
        if ext != "fits":
            ext = "fits"
        if not base.endswith('msa'):
            "".join([base, '_msa'])
        "".join([base, ext])
    else:
        root = model.meta.filename.split('_')
        output = "".join([root[0], '_msa', '.fits'])
    hdulist.writeto(output)
    del hdulist
    model.close()
