#!/usr/bin/env python

"""Add velocity aberration correction information to the FITS files provided."""

import logging
import sys
import warnings
import argparse
from pathlib import Path

from jwst.lib.set_velocity_aberration import add_dva

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def parse_args(args):
    """
    Parse given arguments.

    Parameters
    ----------
    args : command-line input
        The arguments to parse

    Returns
    -------
    parser : object
        Parsed arguments
    """
    description_text = """
    Add velocity aberration correction information to FITS files.
    Input via command line can be one or more FITS files.
    It assumes the following keywords are present in the file header:
    JWST_DX (km/sec),
    JWST_DY (km/sec),
    JWST_DZ (km/sec),
    RA_REF (deg),
    DEC_REF (deg).

    The keywords added are:
    VA_SCALE (dimensionless scale factor)

    It does not currently place the new keywords in any particular location
    in the header other than what is required by the standard.
    """
    parser = argparse.ArgumentParser(
        prog="set_velocity_aberration",
        description=description_text,
    )
    parser.add_argument(
        "filename",
        nargs="+",
        help="FITS file(s) to which to add the velocity aberration correction information.",
    )
    parser.add_argument(
        "-f",
        "--force_level1bmodel",
        choices=[1, 0],
        default=1,
        type=int,
        help="Force the input file to be treated as a Level1bModel. Options 0 (do not force)"
        "or 1 (yes, force). Default is 1.",
    )
    return parser.parse_args(args)


def main():
    """Parse arguments and add velocity aberration correction information to the files provided."""
    args = parse_args(sys.argv[1:])
    for filename in args.filename:
        add_dva(filename)


def deprecated_name():
    """Raise warning if filename.* is no longer used, and provide correct one."""
    filename = Path(__file__)
    warnings.warn(
        f"usage of `{filename.name}` is deprecated; use `{filename.stem}` instead", stacklevel=2
    )

    main()


if __name__ == "__main__":
    main()
