#! /usr/bin/env python
"""
A simple tool to read in the output of extract2d (FS and MOS) or assign_wcs (IFU) and
apply the WCS transforms to all pixels in a slit. For each slit it writes the results
as a cube with three planes (wavelength, ra, dec) in a separate fits extension.
The file is saved with an suffix "world_coordinates".

Requested by the NIRSPEC team.

Build 6 testing.
"""
from __future__ import absolute_import, division, unicode_literals, print_function
import os.path
import numpy as np
from astropy.io import fits
from gwcs import wcstools
from ... import datamodels
from .. import nirspec


imaging_modes = supported_modes = ['nrs_taconfirm', 'nrs_brightobj', 'nrs_bota', 'nrs_tacq', 'nrs_focus',
                                   'nrs_lamp', 'nrs_mimf', 'nrs_image', 'nrs_confirm', 'nrs_taslit']


def ifu_coords(fname, output=None):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    IFU exposure after running assign_wcs.

    Parameters
    ----------
    fname : str
        The name of a file after assign_wcs was run.
    output : str
        The name of the output file. If None the root of the input
        file is used with an extension world_coordinates.

    Examples
    --------
    >>> compute_world_coordinates('nrs1_fixed_assign_wcs_extract_2d.fits')

    """
    model = datamodels.ImageModel(fname)
    if model.meta.exposure.type.lower() != 'nrs_ifu':
        raise ValueError("Expected an IFU observation,"
                         "(EXP_TYPE=NRS_IFU), got {0}".format(model.meta.exposure.type))
    ifu_slits = nirspec.nrs_ifu_wcs(model)

    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    output_frame = ifu_slits[0].available_frames[-1]
    for i, slit in enumerate(ifu_slits):
        x, y = wcstools.grid_from_domain(slit.domain)
        ra, dec, lam = slit(x, y)
        detector2slit = slit.get_transform('detector', 'slit_frame')
        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([np.rot(lam), np.rot90(ra), np.rot90(dec), np.rot90(sy)])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = '{0}_x, arcsec'.format(output_frame)
        imhdu.header['PLANE3'] = '{0}_y, arcsec'.format(output_frame)
        imhdu.header['PLANE4'] = 'slit_y, relative to center (0, 0)'
        imhdu.header['SLIT'] = "SLIT_{0}".format(i)
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


def compute_world_coordinates(fname, output=None):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    FS or MOS observation after running extract_2d.

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
    model = datamodels.MultiSlitModel(fname)
    if model.meta.exposure.type.lower() not in ['nrs_fixedslit', 'nrs_msaspec']:
        raise ValueError("Expected a FS or MOS observation,"
                         "(EXP_TYPE=NRS_FIXEDSLIT, NRS_MSASPEC),"
                         "got {0}".format(model.meta.exposure.type))
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    output_frame = model.slits[0].meta.wcs.available_frames[-1]
    for slit in model.slits:
        #x, y = wcstools.grid_from_domain(slit.meta.wcs.domain)
        xstart, xend = slit.xstart, slit.xstart + slit.xsize
        ystart, yend = slit.ystart, slit.ystart + slit.ysize
        x, y = np.mgrid[xstart: xend, ystart: yend]
        ra, dec, lam = slit.meta.wcs(x, y)
        detector2slit = slit.meta.wcs.get_transform('detector', 'slit_frame')

        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([np.rot90(lam), np.rot90(ra), np.rot90(dec), np.rot90(sy)])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = '{0}_x, arcsec'.format(output_frame)
        imhdu.header['PLANE3'] = '{0}_y, arcsec'.format(output_frame)
        imhdu.header['PLANE4'] = 'slit_y, relative to center (0, 0)'
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


def imaging_coords(fname, output=None):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    imaging exposure after running assign_wcs.

    Parameters
    ----------
    fname : str
        The name of a file after assign_wcs was run.
    output : str
        The name of the output file. If None the root of the input
        file is used with an extension world_coordinates.

    Examples
    --------
    >>> compute_world_coordinates('nrs1_fixed_assign_wcs_extract_2d.fits')

    """
    model = datamodels.ImageModel(fname)
    if model.meta.exposure.type.lower() not in imaging_modes:
        raise ValueError("Observation mode {0} is not supported.".format(model.meta.exposure.type))
    ifu_slits = nirspec.nrs_ifu_wcs(model)

    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    output_frame = ifu_slits[0].available_frames[-1]
    for i, slit in enumerate(ifu_slits):
        x, y = wcstools.grid_from_domain(slit.domain)
        ra, dec, lam = slit(x, y)
        world_coordinates = np.array([lam, ra, dec])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = '{0}_x, arcsec'.format(output_frame)
        imhdu.header['PLANE3'] = '{0}_y, arcsec'.format(output_frame)
        imhdu.header['SLIT'] = "SLIT_{0}".format(i)
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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Create NIRSPEC world coordinates file.")
    parser.add_argument('mode', type=str, help='Observing mode: one of "ifu", "fs", "mos", "imaging"')
    parser.add_argument('filename', help='Name of file after assign_wcs (IFU) or extract_2d was run (FS and MOS)')
    res = parser.parse_args()
    if res.mode.lower() == "ifu":
        ifu_coords(res.filename)
    elif res.mode.lower() in ["fs", "mos", "msa"]:
        compute_world_coordinates(res.filename)
    elif res.mode.lower() == "imaging":
        imaging_coords(res.filename)
    else:
        print("Invalid mode: {0} ".format(res.mode))
