#!/usr/bin/env python
# Copyright (C) 2023 Association of Universities for Research in Astronomy (AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
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
#
"""
Apply adjustments to an imaging JWST GWCS object. Input files can be
``DataModel``s serialized to ASDF or ASDF-in-FITS files, or they could be
"simple" ASDF files storing just the GWCS model in the root of the ASDF file
under the key ``'wcs'``, i.e., ``asdf_file.tree['wcs']``.

Examples
--------

From command line::

    % adjust_wcs [-h] (-u | --suffix SUFFIX | -f FILE) [--overwrite]
                 [-r RA_DELTA [units]] [-d DEC_DELTA [units]]
                 [-o ROLL_DELTA [units]] [-s SCALE_FACTOR]
                 [-v] [input ...]

    % adjust_wcs data_model_*_cal.fits -u -s 1.002

    % adjust_wcs data_model_*_cal.fits --suffix wcsadj -s 1.002 -r 0.2 arcsec -o 0.0023 deg
"""
import argparse
import glob
import logging
import os
import sys

import jwst
from jwst.tweakreg.utils import adjust_wcs
from jwst.assign_wcs.util import update_fits_wcsinfo
from astropy import units


_ANGLE_PARS = ["-r", "--ra_delta", "-d", "--dec_delta", "-o", "--roll_delta"]


# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.NullHandler())


def _replace_suffix(file, new_suffix):
    sep = '_'
    directory, base = os.path.split(file)
    root, ext = os.path.splitext(base)
    if sep in root:
        name_parts = root.rpartition(sep)
        root = ''.join(name_parts[:-1] + (new_suffix, ))
    else:
        root = root + sep + new_suffix
    base = root + ext

    new_file_name = os.path.join(directory, base)

    return new_file_name


def _is_float(arg):
    try:
        float(arg)
        return True
    except ValueError:
        pass
    return False


def _is_unit(arg):
    if arg.startswith("-"):
        return False
    try:
        units.Unit(arg)
        return True
    except ValueError:
        pass
    return False


def angle(arg):
    args = arg.strip().split(' ')
    if len(args) == 1:
        return float(args[0])
    elif len(args) == 2:
        return units.Quantity(float(args[0]), units.Unit(args[1]), dtype=float)
    else:
        raise argparse.ArgumentTypeError()


def main():
    if len(sys.argv) <= 1:
        raise ValueError("Missing required arguments.")

    # Parse input parameters
    parser = argparse.ArgumentParser(
        prog="adjust_wcs",
        description="Apply adjustments to data model's WCS"
    )

    parser.add_argument(
        "arg0",
        nargs='*',
        metavar="input",
        type=str,
        help="Input data model file(s)."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-u",
        "--update",
        action="store_true",
        help="Update input data models."
    )
    group.add_argument(
        "--suffix",
        type=str,
        help="Suffix for output files."
    )
    group.add_argument(
        "-f",
        "--file",
        type=str,
        help="Output file name."
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files (ignored with -u/--update)."
    )

    group_wcs = parser.add_argument_group(
        title='WCS adjustment parameters'
    )
    group_wcs.add_argument(
        "-r",
        "--ra_delta",
        type=angle,
        default=0.0,
        help="Delta RA (longitude), degrees."
    )
    group_wcs.add_argument(
        "-d",
        "--dec_delta",
        type=angle,
        default=0.0,
        help="Delta DEC (latitude), degrees."
    )
    group_wcs.add_argument(
        "-o",
        "--roll_delta",
        type=angle,
        default=0.0,
        help="Delta roll angle, degrees."
    )
    group_wcs.add_argument(
        "-s",
        "--scale_factor",
        type=float,
        default=1.0,
        help="Scale factor."
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="v{0}".format(jwst.__version__)
    )

    # pre-process argv for units and negative floats:
    argv = sys.argv[1:]
    argv_new = [os.path.basename(sys.argv[0])]
    while argv:
        argv_new.append(argv.pop(0))

        # check whether next argument is a float:
        if (argv_new[-1] in _ANGLE_PARS and len(argv) >= 1 and
                _is_float(argv[0])):
            angle_value = argv.pop(0)
            # check whether next argument is a unit:
            if len(argv) >= 1 and _is_unit(argv[0]):
                angle_unit = argv.pop(0)
                argv_new.append(f" {angle_value} {angle_unit}")
            else:
                # assume angle units are degrees
                argv_new.append(f" {angle_value}")

    options = parser.parse_args(argv_new)

    files = []
    for f in options.arg0:
        files.extend(glob.glob(f))

    if options.file and len(files) > 1 :
        parser.error(
            "argument -f/--file: not allowed with multiple input files"
        )

    wcs_pars = [
        options.ra_delta,
        options.dec_delta,
        options.roll_delta,
        options.scale_factor
    ]

    if wcs_pars == [0.0, 0.0, 0.0, 1.0]:
        logger.info(
            "All WCS adjustment parameters ('ra_delta', 'dec_delta',"
            "roll_delta, and scale_factor) have default values"
            "(identical transformation)."
        )
        logger.info("Nothing to do. Quitting.")
        return

    for fname in files:
        data_model = jwst.datamodels.open(fname)

        # Compute adjusted WCS:
        out_wcs = adjust_wcs(
            wcs=data_model.meta.wcs,
            delta_ra=options.ra_delta,
            delta_dec=options.dec_delta,
            delta_roll=options.roll_delta,
            scale_factor=options.scale_factor
        )

        # Save results:
        data_model.meta.wcs = out_wcs
        # Also update FITS representation in input exposures for
        # subsequent reprocessing by the end-user.
        try:
            update_fits_wcsinfo(
                data_model,
                max_pix_error=0.001,
                npoints=64
            )
        except (ValueError, RuntimeError) as e:
            logger.warning(
                "Failed to update 'meta.wcsinfo' with FITS SIP "
                f'approximation. Reported error is:\n"{e.args[0]}"'
            )

        if options.update:
            data_model.save(fname, overwrite=True)
            logger.info(f"Updated data model '{fname}'")
        else:
            if options.file is None:
                output_file = _replace_suffix(fname, options.suffix)
            else:
                output_file = options.file

            data_model.save(output_file, overwrite=options.overwrite)
            logger.info(f"Saved data model to '{output_file}'")


if __name__ == '__main__':
    main()
