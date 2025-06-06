#!/usr/bin/env python

"""
Read in the output of extract_2d and apply the WCS transforms to all pixels in a slit.

The extract_2d file can be (FS and MOS) or assign_wcs (IFU).

For each slit it writes the results as a cube with four planes (wavelength, ra, dec, y_slit)
in a separate fits extension.

The file is saved with a suffix 'world_coordinates'.
"""

# Licensed under a 3-clause BSD style license - see LICENSE

import argparse
import os
import warnings
import logging

import numpy as np
from astropy.io import fits
from gwcs import wcstools
from gwcs.utils import _toindex

from stdatamodels.jwst.datamodels.util import open as dmopen

from jwst.assign_wcs import nirspec


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def main():
    """Loop over the input files and apply the proper method to create the output file."""
    short_description = "Create NIRSPEC world coordinates file"
    long_description = """

    A tool to read in the output of extract_2d (FS and MOS) or assign_wcs
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

    parser.add_argument("filenames", help="Name of input and/or output files", nargs="*")

    parser.add_argument("-m", "--mode", default=None, help="set exposure mode")

    res = parser.parse_args()
    filenames = res.filenames
    mode = res.mode

    methods = {
        "nrs_ifu": ifu_coords,
        "nrs_fixedslit": compute_world_coordinates,
        "nrs_msaspec": compute_world_coordinates,
        "nrs_taconfirm": imaging_coords,
        "nrs_brightobj": imaging_coords,
        "nrs_bota": imaging_coords,
        "nrs_tacq": imaging_coords,
        "nrs_focus": imaging_coords,
        "nrs_lamp": imaging_coords,
        "nrs_mimf": imaging_coords,
        "nrs_image": imaging_coords,
        "nrs_confirm": imaging_coords,
        "nrs_taslit": imaging_coords,
    }

    ok = True
    outputs = build_output_files(filenames)
    for input_file, output_file in zip(filenames, outputs, strict=True):
        if not os.path.exists(input_file):  # noqa: PTH110
            warn_user(input_file, "is not found")
            ok = False
            continue

        try:
            model = dmopen(input_file, pass_invalid_values=True)
            if mode is None:
                mode = model.meta.exposure.type
        except OSError:
            warn_user(input_file, "is not valid fits")
            ok = False
            continue

        method = get_update_method(methods, model, mode)
        if method is None:
            ok = False
        else:
            hdulist = method(model)
            hdulist.writeto(output_file, overwrite=True)
            logging.info(f"File written: {output_file}")
            del hdulist
        model.close()

    if not ok:
        msg = "Error in input file"
        if len(filenames) > 1:
            msg += "s"
        raise ValueError(msg)


def build_output_files(filenames):
    """
    Create a list of output files of the same length as the list of input files.

    Parameters
    ----------
    filenames : list of strings
        A list of input and output filenames

    Returns
    -------
    outputs : list
        Output file names
    """
    # Check to see if the last file is an output file
    output = filenames[-1]
    if os.path.exists(output) and not os.path.isdir(output):  # noqa: PTH110, PTH112
        output = None
    else:
        filenames.pop()

    # Handle the case of a single input and output
    # Punt the other cases to the code below
    if output is not None and not os.path.isdir(output):  # noqa: PTH112
        if len(filenames) == 1:
            return [output]
        else:
            # Let other code handle the missing file error
            filenames.append(output)
            outputs = []
    else:
        outputs = []

    # Track if output file was a directory
    if output is not None and os.path.isdir(output):  # noqa: PTH112
        output_dir = output
    else:
        output_dir = None

    # Build an output filename for each input
    for filename in filenames:
        # Build output base and extensions
        dirname, root = os.path.split(filename)
        base, ext = os.path.splitext(root)  # noqa: PTH122
        if ext != ".fits":
            ext = ".fits"
        if not base.endswith("world_coordinates"):
            base = "_".join([base, "world_coordinates"])

        # Substitute output directory if present
        if output_dir is not None:
            dirname = output_dir

        # Combine directory, base, and extension
        outputs.append(os.path.join(dirname, base + ext))  # noqa: PTH118

    return outputs


def compute_world_coordinates(model):
    """
    Compute wavelengths and space coordinates of an NRS FS or MOS observation after extract_2d.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract_2d.

    Returns
    -------
    hdulist : fits.HDUList
        The list of HDUs object.
    """
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["filename"] = model.meta.filename
    phdu.header["data"] = "world coordinates"
    hdulist.append(phdu)
    output_frame = model.slits[0].meta.wcs.available_frames[-1]

    for slit in model.slits:
        # slit.x(y)start are 1-based, turn them to 0-based for extraction
        # xstart, xend = slit.xstart - 1, slit.xstart -1 + slit.xsize
        # ystart, yend = slit.ystart - 1, slit.ystart -1 + slit.ysize
        # y, x = np.mgrid[ystart: yend, xstart: xend]
        x, y = wcstools.grid_from_bounding_box(slit.meta.wcs.bounding_box, step=(1, 1), center=True)

        ra, dec, lam = slit.meta.wcs(x, y)
        detector2slit = slit.meta.wcs.get_transform("detector", "slit_frame")

        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([lam, ra, dec, sy])  # , x, y])
        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header["PLANE1"] = "lambda, microns"
        imhdu.header["PLANE2"] = f"{output_frame}_x, arcsec"
        imhdu.header["PLANE3"] = f"{output_frame}_y, arcsec"
        imhdu.header["PLANE4"] = "slit_y, relative to center (0, 0)"
        imhdu.header["SLIT"] = slit.name

        # add the overall subarray offset
        imhdu.header["CRVAL1"] = slit.xstart - 1 + model.meta.subarray.xstart
        imhdu.header["CRVAL2"] = slit.ystart - 1 + model.meta.subarray.ystart

        imhdu.header["CRPIX1"] = 1
        imhdu.header["CRPIX2"] = 1
        imhdu.header["CTYPE1"] = "pixel"
        imhdu.header["CTYPE2"] = "pixel"
        hdulist.append(imhdu)

    return hdulist


def get_update_method(methods, model, mode):
    """
    Find the method to convert the input file to output.

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

    Returns
    -------
    method : function
        Correct update function to call
    """
    exp_type = model.meta.exposure.type
    if exp_type is None and mode is not None:
        mode = "nrs_" + mode.lower()
        for candidate_type in methods.keys():
            if candidate_type.startswith(mode):
                if exp_type is None:
                    exp_type = candidate_type
                else:
                    exp_type = None
                    break

        if exp_type is None:
            warn_user(mode, "does not match an exposure type")
        else:
            model.meta.exposure.type = exp_type.upper()

    exp_type = exp_type.lower()
    method = methods.get(exp_type)
    if method is None:
        warnings.warn(exp_type, "is not a supported exposure type", stacklevel=2)

    return method


def ifu_coords(model):
    """
    Compute wavelengths, and space coordinates of a NIRSPEC IFU exposure after running assign_wcs.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract_2d.

    Returns
    -------
    hdulist : fits.HDUList
        Fits HDU list object
    """
    ifu_slits = nirspec.nrs_ifu_wcs(model)

    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["filename"] = model.meta.filename
    phdu.header["data"] = "world coordinates"
    hdulist.append(phdu)
    output_frame = ifu_slits[0].available_frames[-1]

    for i, slit in enumerate(ifu_slits):
        x, y = wcstools.grid_from_bounding_box(slit.bounding_box, (1, 1), center=True)
        ra, dec, lam = slit(x, y)
        detector2slit = slit.get_transform("detector", "slit_frame")
        sx, sy, ls = detector2slit(x, y)
        world_coordinates = np.array([lam, ra, dec, sy])

        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header["PLANE1"] = "lambda, microns"
        imhdu.header["PLANE2"] = f"{output_frame}_x, arcsec"
        imhdu.header["PLANE3"] = f"{output_frame}_y, arcsec"
        imhdu.header["PLANE4"] = "slit_y, relative to center (0, 0)"
        imhdu.header["SLIT"] = f"SLIT_{i}"

        # -1 and +1 are to express this in 1-based coordinates
        imhdu.header["CRVAL1"] = (
            model.meta.subarray.xstart - 1 + int(_toindex(slit.bounding_box[0][0])) + 1
        )
        imhdu.header["CRVAL2"] = (
            model.meta.subarray.xstart - 1 + int(_toindex(slit.bounding_box[1][0])) + 1
        )

        # Input coordinates will be 1-based.
        imhdu.header["CRPIX1"] = 1
        imhdu.header["CRPIX2"] = 1
        imhdu.header["CTYPE1"] = "pixel"
        imhdu.header["CTYPE2"] = "pixel"
        hdulist.append(imhdu)
    return hdulist


def imaging_coords(model):
    """
    Compute wavelengths and space coordinates of a NIRSPEC imaging exposure after assign_wcs.

    Parameters
    ----------
    model : ImageModel
        An image with extracted slits, i.e. the output of extract_2d.

    Returns
    -------
    hdulist : fits.HDUList
        Fits HDU list object
    """
    hdulist = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["filename"] = model.meta.filename
    phdu.header["data"] = "world coordinates"
    hdulist.append(phdu)

    output_frame = model.available_frames[-1]
    for i, slit in enumerate(output_frame):
        bb = model.meta.wcs.bounding_box
        x, y = wcstools.grid_from_bounding_box(bb, step=(1, 1), center=True)
        ra, dec, lam = slit(x + 1, y + 1)
        world_coordinates = np.array([lam, ra, dec])

        imhdu = fits.ImageHDU(data=world_coordinates)
        imhdu.header["PLANE1"] = "lambda, microns"
        imhdu.header["PLANE2"] = f"{output_frame}_x, deg"
        imhdu.header["PLANE3"] = f"{output_frame}_y, deg"
        imhdu.header["SLIT"] = f"SLIT_{i}"

        # add the overall subarray offset
        imhdu.header["CRVAL1"] = model.meta.subarray.xstart - 1 + int(_toindex(bb[0][0]))
        imhdu.header["CRVAL2"] = model.meta.subarray.ystart - 1 + int(_toindex(bb[1][0]))

        # Input coordinates will be 1-based.
        imhdu.header["CRPIX1"] = 1
        imhdu.header["CRPIX2"] = 1
        imhdu.header["CTYPE1"] = "pixel"
        imhdu.header["CTYPE2"] = "pixel"
        hdulist.append(imhdu)
    return hdulist


def warn_user(*argv):
    """Send a warning message to stderr."""
    logging.warning(*argv)


if __name__ == "__main__":
    main()
