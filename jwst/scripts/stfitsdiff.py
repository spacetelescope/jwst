#!/usr/bin/env python

"""Allow command stfitsdiff be used from terminal."""

import ast
import logging
import sys
from argparse import ArgumentParser

from jwst.regtest.st_fitsdiff import STFITSDiffBeta as STFITSDiff


logging.basicConfig(level=logging.INFO, format="", datefmt="", stream=sys.stdout)


def _is_number(s):
    is_number = True
    try:
        int(s)
    except ValueError:
        is_number = False
    return is_number


def _check_arglist(arg_list, list_values):
    if arg_list is not None:
        if "," in arg_list:
            list_values.extend(arg_list.split(sep=","))
        else:
            list_values.append(arg_list)
    return list_values


def main():
    """Find the differences between the files and report."""
    # Parse command line arguments.
    parser = ArgumentParser(
        description="Get the differences between two fits files and report,",
        epilog="e.g. $ stfitsdiff jw000_rate_a.fits jw000_rate_b.fits",
    )

    # Required arguments

    parser.add_argument("file_a", type=str, help="The file to compare.")

    parser.add_argument("file_b", type=str, help="The file to compare from.")

    # Optional arguments (with corresponding default values)

    parser.add_argument(
        "-ih",
        "--ignore_hdus",
        dest="ignore_hdus",
        action="store",
        default=None,
        help="HDU names to ignore when comparing two FITS files.",
    )

    parser.add_argument(
        "-ik",
        "--ignore_keywords",
        dest="ignore_keywords",
        action="store",
        default=None,
        help="Header keywords to ignore when comparing two headers.",
    )

    parser.add_argument(
        "-ic",
        "--ignore_comments",
        dest="ignore_comments",
        action="store",
        default=None,
        help="A list of header keywords whose comments should be ignored in the comparison.",
    )

    parser.add_argument(
        "-if",
        "--ignore_fields",
        dest="ignore_fields",
        action="store",
        default=None,
        help="""The (case-insensitive) names of any table columns to ignore, if any table data
                             is to be compared.""",
    )

    parser.add_argument(
        "-n",
        "--numdiffs",
        dest="numdiffs",
        action="store",
        default=10,
        type=int,
        help="""The number of pixel/table values to output when reporting HDU data differences.
                             Though the count of differences is the same either way, this allows
                             controlling the number of different values that are kept in memory
                             or output. If a negative value is given, then numdiffs is treated
                             as unlimited.""",
    )

    parser.add_argument(
        "-r",
        "--rtol",
        dest="rtol",
        action="store",
        default=1e-5,
        type=float,
        help="The relative difference to allow when comparing two float values either in "
        "header values, image arrays, or table columns.",
    )

    parser.add_argument(
        "-a",
        "--atol",
        dest="atol",
        action="store",
        default=1e-7,
        type=float,
        help="The allowed absolute difference.",
    )

    parser.add_argument(
        "-ib",
        "--ignore_blanks",
        dest="ignore_blanks",
        action="store_false",
        default=True,
        help="Ignore extra whitespace at the end of string values either in headers or data.",
    )

    parser.add_argument(
        "-ibc",
        "--ignore_blank_cards",
        dest="ignore_blank_cards",
        action="store_false",
        default=True,
        help="Ignore all cards that are blank, i.e. they only contain whitespace.",
    )

    parser.add_argument(
        "-rpl",
        "--report_pixel_loc_diffs",
        dest="report_pixel_loc_diffs",
        action="store_true",
        default=False,
        help="Report the numdiffs pixel locations where differences are found.",
    )

    parser.add_argument(
        "-et",
        "--extension_tolerances",
        dest="extension_tolerances",
        action="store",
        default=None,
        help="""Provide a different relative and absolute tolerance for the given extensions
                             (use no spaces and double quotes encasing the whole dictionary), e.g.
                             --extension_tolerances="{'sci':{'rtol':1e-3,'atol':1e-2}}}"
                             It does not matter if the keys in the dictionary are upper or lower
                             case. The key for only providing main header tolerances is 'primary'.
                             The key 'headers' can be used to provide a special tolerance for all
                             extension headers in the HDU list.
                             The key 'default' is optional, i.e. if it is not provided then the
                             default values will be used, otherwise the default value will be the
                             one in the dictionary.""",
    )

    # Get the arguments

    args = parser.parse_args()
    file_a = args.file_a
    file_b = args.file_b

    stfitsdiff_default_kwargs = {
        "ignore_hdus": [],
        "ignore_keywords": [],
        "ignore_comments": [],
        "ignore_fields": [],
        "numdiffs": args.numdiffs,
        "rtol": args.rtol,
        "atol": args.atol,
        "ignore_blanks": args.ignore_blanks,
        "ignore_blank_cards": args.ignore_blank_cards,
        "report_pixel_loc_diffs": args.report_pixel_loc_diffs,
    }

    ignore_hdus = ["ASDF"]
    if args.ignore_hdus is not None:
        ignore_hdus = _check_arglist(args.ignore_hdus, ignore_hdus)
        # Check if there are any number extensions in ignore_hdus
        for exti, ext in enumerate(ignore_hdus):
            if _is_number(ext):
                ignore_hdus[exti] = int(ext)
    stfitsdiff_default_kwargs["ignore_hdus"] = ignore_hdus

    ignore_keywords = ["DATE", "CAL_VER", "CAL_VCS", "CRDS_VER", "CRDS_CTX", "NAXIS1", "TFORM*"]
    if args.ignore_keywords is not None:
        ignore_keywords = _check_arglist(args.ignore_keywords, ignore_keywords)
    stfitsdiff_default_kwargs["ignore_keywords"] = ignore_keywords

    ignore_comments = ["DATE", "CAL_VER", "CAL_VCS", "CRDS_VER", "CRDS_CTX", "NAXIS1", "TFORM*"]
    if args.ignore_comments is not None:
        ignore_comments = _check_arglist(args.ignore_comments, ignore_comments)
    stfitsdiff_default_kwargs["ignore_comments"] = ignore_comments

    if args.ignore_fields is not None:
        stfitsdiff_default_kwargs["ignore_fields"] = _check_arglist(
            args.ignore_fields, stfitsdiff_default_kwargs["ignore_fields"]
        )

    # If provided, make sure the extension_tolerances is a dictionary and not a string
    err_msg = """
    Dictionary format error. Use no spaces and double quotes encasing the whole dictionary, e.g.
    --extension_tolerances="{'sci':{'rtol':1e-3,'atol':1e-2},'err':{'rtol':1e-1,'atol':1e-2}}"
    """
    if args.extension_tolerances is not None:
        try:
            stfitsdiff_default_kwargs["extension_tolerances"] = ast.literal_eval(
                args.extension_tolerances
            )
        except (NameError, TypeError, ValueError, SyntaxError):
            logging.error(err_msg)
            exit()

    # Find the differences
    logging.info("\n* STScI Custom FITSDiff")
    try:
        diff = STFITSDiff(file_a, file_b, **stfitsdiff_default_kwargs)
        logging.info(diff.report())
    except (NameError, TypeError, ValueError, SyntaxError):
        logging.error(err_msg)


if __name__ == "__main__":
    main()
