#!/usr/bin/env python

"""Pointing verification."""

import argparse
import logging
from sys import stdout

from jwst.lib import pointing_summary

__all__ = []  # type: ignore[var-annotated]


# Begin execution
def main():
    """Run pointing verification."""
    parser = argparse.ArgumentParser(
        description="""Summarize various pointing information in a table. Compare
                    the calculated V1 and REFPOINT pointing with the proposed TARGET pointing.
                    E.g.
                    $ pointing_summary exp1.fits
                    $ pointing_summary *.fits
                    """
    )

    parser.add_argument(
        "exposures", type=str, nargs="+", help="List of JWST data files to examine."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=stdout,
        help="File to write summary table to. Default is standard output.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity. Specifying multiple times adds more output.",
    )
    parser.add_argument(
        "--extra-meta",
        type=str,
        nargs="+",
        help="Extra meta information from the exposures to add to the result table",
    )

    args = parser.parse_args()

    # Configure logging
    log_handler = logging.StreamHandler()
    logger = logging.getLogger("jwst")
    logger.addHandler(log_handler)
    log_levels = [logging.WARNING, logging.INFO, logging.DEBUG]

    # Set output detail.
    level = log_levels[min(len(log_levels) - 1, args.verbose)]
    logger.setLevel(level)

    # Process the file list.
    logger.info("Starting pointing summary.")
    deltas = pointing_summary.calc_deltas(args.exposures, extra_meta=args.extra_meta)
    deltas.write(args.output, format="ascii.ecsv")
    logger.info("........Summary completed.")


if __name__ == "__main__":
    main()
