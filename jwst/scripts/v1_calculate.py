#!/usr/bin/env python

"""V1 Calculation based on time and engineering database info
"""
import argparse
import logging
from sys import stdout

from astropy.time import Time

import jwst.lib.set_telescope_pointing as stp
from jwst.lib import v1_calculate

# Configure logging
logger = logging.getLogger('jwst')
logger.propagate = False
logger_handler = logging.StreamHandler()
logger.addHandler(logger_handler)
logger_format_debug = logging.Formatter('%(levelname)s:%(filename)s::%(funcName)s: %(message)s')

# Available reduce functions
REDUCE_FUNCS_MAPPING = {
    'all': stp.all_pointings,
    'first': stp.first_pointing,
    'average': stp.pointing_from_average
}
REDUCE_FUNCS = list(REDUCE_FUNCS_MAPPING.keys())


# Begin execution
def main():
    parser = argparse.ArgumentParser(
        description='Calculate V1 over a time period.'
    )

    parser.add_argument(
        'time_sources', type=str, nargs='+',
        help=('Either a list of JWST data files to retrieve the timing from'
              ' or a start and end time to retrieve pointing information for.'
              )
    )
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=stdout,
        help='File to write V1 calculation table to. Default is standard output.'
    )
    parser.add_argument(
        '--pointing', type=str, choices=REDUCE_FUNCS, default='average',
        help=('Which pointing(s) to use within the specified time interval.'
              f' Choices: {REDUCE_FUNCS}'
              ' (default: %(default)s)'
              )
    )
    parser.add_argument(
        '--method',
        type=stp.Methods, choices=list(stp.Methods), default=stp.Methods.default,
        help='Algorithmic method to use. Default: %(default)s'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='Increase verbosity. Specifying multiple times adds more output.'
    )
    parser.add_argument(
        '--override-transforms', type=str, default=None,
        help='Transform matrices to use instead of calculated'
    )
    parser.add_argument(
        '--engdb-url', type=str, default=None,
        help=('The engineering database to use. If unspecified, the value'
              ' of the environment ENG_BASE_URL is used. Otherwise, the production'
              ' system is used.'
              )
    )

    args = parser.parse_args()

    # Set output detail.
    level = stp.LOGLEVELS[min(len(stp.LOGLEVELS) - 1, args.verbose)]
    logger.setLevel(level)
    if level <= logging.DEBUG:
        logger_handler.setFormatter(logger_format_debug)

    override_transforms = args.override_transforms
    if override_transforms:
        override_transforms = stp.Transforms.from_asdf(override_transforms)

    # Determine whether the sources are time specifications or a file list.
    if len(args.time_sources) == 2:
        try:
            obsstart = Time(args.time_sources[0])
            obsend = Time(args.time_sources[1])
        except ValueError:
            input_as_files = True
        else:
            logger.info(f'Retrieving V1 over the time span {obsstart.isot} - {obsend.isot}')
            input_as_files = False
            if args.pointing != 'all':
                logger.warning(
                    'V1 pointings have been requested over a time range.'
                    ' However, the \'pointing\' option is not \'all\'.'
                    '\nThere will only be a single result returned. Is this what was desired?'
                    '\nSuggestion: Use \'--pointing=all\''
                )
    else:
        input_as_files = True

    # Process the file list.
    logger.info('Starting V1 calculation...')
    if input_as_files:
        v1s = v1_calculate.v1_calculate_from_models(
            args.time_sources, engdb_url=args.engdb_url,
            reduce_func=REDUCE_FUNCS_MAPPING[args.pointing],
            method=args.method, override_transforms=override_transforms
        )
    else:
        v1s = v1_calculate.v1_calculate_over_time(
            obsstart.mjd, obsend.mjd, engdb_url=args.engdb_url,
            reduce_func=REDUCE_FUNCS_MAPPING[args.pointing],
            method=args.method, override_transforms=override_transforms
        )

    formatted = v1_calculate.simplify_table(v1s)
    formatted.write(args.output, format='ascii.ecsv')
    logger.info('...V1 calculation completed.')


if __name__ == '__main__':
    main()
