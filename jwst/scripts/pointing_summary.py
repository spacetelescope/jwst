#!/usr/bin/env python

"""Pointing verification

Review contents of a set of given models for pointing information.
Compare the calculated V1 and REFPOINT pointing with the proposed TARGET pointing.

Examples
--------
>>> pointing_summary exp1.fits

>>> pointing_summary *.fits
"""
import argparse
import logging
from sys import stdout

from jwst.lib import pointing_summary

log_handler = logging.StreamHandler()
logger = logging.getLogger('jwst')
logger.addHandler(log_handler)
LogLevels = [logging.WARNING, logging.INFO, logging.DEBUG]


# Begin execution
def main():
    parser = argparse.ArgumentParser(
        description='Summarize various pointing information in a table.'
    )

    parser.add_argument(
        'exposures', type=str, nargs='+',
        help='List of JWST data files to examine.'
    )
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=stdout,
        help='File to write summary table to. Default is standard output.'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='Increase verbosity. Specifying multiple times adds more output.'
    )
    parser.add_argument(
        '--extra-meta', type=str, nargs='+',
        help='Extra meta information from the exposures to add to the result table'
    )

    args = parser.parse_args()

    # Set output detail.
    level = LogLevels[min(len(LogLevels) - 1, args.verbose)]
    logger.setLevel(level)

    # Process the file list.
    logger.info('Starting pointing summary.')
    deltas = pointing_summary.calc_deltas(args.exposures, extra_meta=args.extra_meta)
    deltas.write(args.output, format='ascii.ecsv')
    logger.info('........Summary completed.')


if __name__ == '__main__':
    main()
