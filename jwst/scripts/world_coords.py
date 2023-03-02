#!/usr/bin/env python

# Copyright (C) 2018 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.


"""
A tool to read in the output of extract2d (FS and MOS) or assign_wcs (IFU) and
apply the WCS transforms to all pixels in a slit. For each slit it writes the results
as a cube with four planes (wavelength, ra, dec, y_slit) in a separate fits extension.
The file is saved with a suffix 'world_coordinates'.

Requested by the NIRSPEC team.

Build 7.1 testing.
"""

import argparse
import os
import os.path
import sys
import warnings

import numpy as np
from astropy.io import fits
from gwcs import wcstools
from gwcs.utils import _toindex

from stdatamodels.jwst.datamodels.util import open as dmopen

from jwst.assign_wcs import nirspec


def main(filenames, mode):
    """
    Loop over the input files and apply the proper method to
    create the output file.

    Parameters
    ----------
    filenames : list of strings
        A list of input and output file names
    mode : str or None
        The exposure mode to use if exposure type
        is not in the image header
    """
    methods = {
        'nrs_ifu': ifu_coords,
        'nrs_fixedslit': compute_world_coordinates,
        'nrs_msaspec': compute_world_coordinates,
        'nrs_taconfirm': imaging_coords,
        'nrs_brightobj': imaging_coords,
        'nrs_bota': imaging_coords,
        'nrs_tacq': imaging_coords,
        'nrs_focus': imaging_coords,
        'nrs_lamp': imaging_coords,
        'nrs_mimf': imaging_coords,
        'nrs_image': imaging_coords,
        'nrs_confirm': imaging_coords,
        'nrs_taslit': imaging_coords,
    }

    ok = True
    outputs = build_output_files(filenames)
    for (input_file, output_file) in zip(filenames, outputs):
        if not os.path.exists(input_file):
            warn_user(input_file, 'is not found')
            ok = False
            continue

        try:
            model = dmopen(input_file, pass_invalid_values=True)
        except IOError:
            warn_user(input_file, 'is not valid fits')
            ok = False
            continue

        method = get_update_method(methods, model, mode)
        if method is None:
            ok = False
        else:
            hdulist = method(model)
            hdulist.writeto(output_file, overwrite=True)
            del hdulist
        model.close()

    if not ok:
        msg = 'Error in input file'
        if len(filenames) > 1:
            msg += 's'
        raise ValueError(msg)


def build_output_files(filenames):
    """
    Create a list of output files of the same length
    as the list of input files

    Parameters
    ----------
    filenames : list of strings
        A list of imput and output filenames

    returns a list of output file names
    removes output file from filenames
    """

    # Check to see if the last file is an output file
    output = filenames[-1]
    if os.path.exists(output) and not os.path.isdir(output):
        output = None
    else:
        filenames.pop()

    # Handle the case of a single input and output
    # Punt the other cases to the code below
    if output is not None and not os.path.isdir(output):
        if len(filenames) == 1:
            return [output]
        else:
            # Let other code handle the missing file error
            filenames.append(output)
            outputs = []
    else:
        outputs = []

    # Track if output file was a directory
    if output is not None and os.path.isdir(output):
        output_dir = output
    else:
        output_dir = None

    # Build an output filename for each input
    for filename in filenames:
        # Build output base and extensions
        dirname, root = os.path.split(filename)
        base, ext = os.path.splitext(root)
        if ext != '.fits':
            ext = '.fits'
        if not base.endswith('world_coordinates'):
            base = '_'.join([base, 'world_coordinates'])

        # Substitute output directory if present
        if output_dir is not None:
            dirname = output_dir

        # Combine directory, base, and extension
        outputs.append(os.path.join(dirname, base + ext))

    return outputs


def compute_world_coordinates(model):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    FS or MOS observation after running extract_2d.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract2d.

    returns a fits hdulist
    """
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    output_frame = model.slits[0].meta.wcs.available_frames[-1]

    for slit in model.slits:
        # slit.x(y)start are 1-based, turn them to 0-based for extraction
        # xstart, xend = slit.xstart - 1, slit.xstart -1 + slit.xsize
        # ystart, yend = slit.ystart - 1, slit.ystart -1 + slit.ysize
        # y, x = np.mgrid[ystart: yend, xstart: xend]
        x, y = wcstools.grid_from_bounding_box(
            slit.meta.wcs.bounding_box, step=(1, 1), center=True
        )

        ra, dec, lam = slit.meta.wcs(x, y)
        detector2slit = slit.meta.wcs.get_transform('detector', 'slit_frame')

        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([lam, ra, dec, sy])  # , x, y])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = '{0}_x, arcsec'.format(output_frame)
        imhdu.header['PLANE3'] = '{0}_y, arcsec'.format(output_frame)
        imhdu.header['PLANE4'] = 'slit_y, relative to center (0, 0)'
        imhdu.header['SLIT'] = slit.name

        # add the overall subarray offset
        imhdu.header['CRVAL1'] = slit.xstart - 1 + model.meta.subarray.xstart
        imhdu.header['CRVAL2'] = slit.ystart - 1 + model.meta.subarray.ystart

        imhdu.header['CRPIX1'] = 1
        imhdu.header['CRPIX2'] = 1
        imhdu.header['CTYPE1'] = 'pixel'
        imhdu.header['CTYPE2'] = 'pixel'
        hdulist.append(imhdu)

    return hdulist


def get_update_method(methods, model, mode):
    """
    Find the method to convert the input file to output

    Parameters
    ----------
    methods : dictionary
        A dictionary whose keys are exposure type strings and
        whose values are the update functions associated with them
    model : ImageModel
        A model whose metadata possibly contains the exposure type
    mode : str
        A mode used to select the method if the exposure type
        is not found in the model

    returns the update function
    """
    exp_type = model.meta.exposure.type
    if exp_type is None and mode is not None:
        mode = 'nrs_' + mode.lower()
        for candidate_type in methods.keys():
            if candidate_type.startswith(mode):
                if exp_type is None:
                    exp_type = candidate_type
                else:
                    exp_type = None
                    break

        if exp_type is None:
            warn_user(mode, 'does not match an exposure type')
        else:
            model.meta.exposure.type = exp_type.upper()

    exp_type = exp_type.lower()
    method = methods.get(exp_type)
    if method is None:
        warnings.warn(exp_type, 'is not a supported exposure type')

    return method


def ifu_coords(model):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    IFU exposure after running assign_wcs.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract2d.


    returns a fits hdulist
    """
    ifu_slits = nirspec.nrs_ifu_wcs(model)

    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)
    output_frame = ifu_slits[0].available_frames[-1]

    for i, slit in enumerate(ifu_slits):
        x, y = wcstools.grid_from_bounding_box(slit.bounding_box, (1, 1), center=True)
        ra, dec, lam = slit(x, y)
        detector2slit = slit.get_transform('detector', 'slit_frame')
        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([lam, ra, dec, sy])

        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header['PLANE1'] = 'lambda, microns'
        imhdu.header['PLANE2'] = '{0}_x, arcsec'.format(output_frame)
        imhdu.header['PLANE3'] = '{0}_y, arcsec'.format(output_frame)
        imhdu.header['PLANE4'] = 'slit_y, relative to center (0, 0)'
        imhdu.header['SLIT'] = 'SLIT_{0}'.format(i)

        # -1 and +1 are to express this in 1-based coordinates
        imhdu.header['CRVAL1'] = (
            model.meta.subarray.xstart - 1 + int(_toindex(slit.bounding_box[0][0])) + 1
        )
        imhdu.header['CRVAL2'] = (
            model.meta.subarray.xstart - 1 + int(_toindex(slit.bounding_box[1][0])) + 1
        )

        # Input coordinates will be 1-based.
        imhdu.header['CRPIX1'] = 1
        imhdu.header['CRPIX2'] = 1
        imhdu.header['CTYPE1'] = 'pixel'
        imhdu.header['CTYPE2'] = 'pixel'
        hdulist.append(imhdu)
    return hdulist


def imaging_coords(model):
    """
    Computes wavelengths, and space coordinates of a NIRSPEC
    imaging exposure after running assign_wcs.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract2d.


    returns a fits hdulist
    """
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['filename'] = model.meta.filename
    phdu.header['data'] = 'world coordinates'
    hdulist.append(phdu)

    output_frame = model.available_frames[-1]
    bb = model.meta.wcs.bounding_box
    x, y = wcstools.grid_from_bounding_box(bb, step=(1, 1), center=True)
    ra, dec, lam = slit(x + 1, y + 1)  # noqa
    world_coordinates = np.array([lam, ra, dec])

    imhdu = fits.ImageHDU(data=world_coordinates)
    imhdu.header['PLANE1'] = 'lambda, microns'
    imhdu.header['PLANE2'] = '{0}_x, deg'.format(output_frame)
    imhdu.header['PLANE3'] = '{0}_y, deg'.format(output_frame)
    imhdu.header['SLIT'] = 'SLIT_{0}'.format(i)  # noqa

    # add the overall subarray offset
    imhdu.header['CRVAL1'] = model.meta.subarray.xstart - 1 + int(_toindex(bb[0][0]))
    imhdu.header['CRVAL2'] = model.meta.subarray.ystart - 1 + int(_toindex(bb[1][0]))

    # Input coordinates will be 1-based.
    imhdu.header['CRPIX1'] = 1
    imhdu.header['CRPIX2'] = 1
    imhdu.header['CTYPE1'] = 'pixel'
    imhdu.header['CTYPE2'] = 'pixel'
    hdulist.append(imhdu)
    return hdulist


def warn_user(*argv):
    """
    Send a warning message to stderr

    Parameters
    ----------
    argv : strings
        One or more strings to be printed
    """
    print(' '.join(argv), file=sys.stderr)


if __name__ == '__main__':
    short_description = 'Create NIRSPEC world coordinates file'
    long_description = """

A tool to read in the output of extract2d (FS and MOS) or assign_wcs
(IFU) and apply the WCS transforms to all pixels in a slit. For each
slit it writes the results as a cube with four planes (wavelength, ra,
dec, y_slit) in a separate fits extension. The file is saved with a
suffix 'world_coordinates'.

This script can be run in one of three ways:

With a single input and output filename.

> world_coord input output

With one or more input filenames. The output names are generated from
the inputs.

> world_coord input_1 input_2 ... input_n

With an output directory name.

> world_coord input_1 input_2 ... output_directory

The output filename cannot be the name of an existing file. The output
directory must exist.

Optionally you can choose the exposure type with the -mode flag. This
flag can be abbreviated to -m. The -mode is only used if exposure type
is not found in the primary header. The mode is the exposure type
without the nrs_ prefix or any unique prefix of it.
"""

    parser = argparse.ArgumentParser(
        description=short_description,
        epilog=long_description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-mode', help='set exposure mode')

    parser.add_argument(
        'filenames', help='Name of input and/or output files', nargs='*'
    )

    res = parser.parse_args()
    main(res.filenames, res.mode)
