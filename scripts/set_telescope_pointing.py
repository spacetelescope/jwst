#!/usr/bin/env python
"""
Set the initial world coordinate system for JWST exposures

The JWST engineering database is queried for the JWST observatory
orientation parameters, and converts that orientation to a WCS
for a list of exposures.

Copyright
---------
Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
"""

import argparse
import logging
import sys

from jwst.lib.set_telescope_pointing import add_wcs

logger = logging.getLogger('jwst')
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Update basic WCS information in JWST exposures from the engineering database.'
    )
    parser.add_argument(
        'exposures', type=str, nargs='+',
        help='List of JWST exposures to update.'
    )
    parser.add_argument(
        '--strict-time', action='store_true',
        help='Force pointings to be within the observation time.'
    )
    parser.add_argument(
        '--allow-default', action='store_true',
        help='If pointing information cannot be determine, use header information.'
    )
    parser.add_argument(
        '--siaf', type=str, default=None,
        help='SIAF PRD sqlite database file. If not specified, default is to use `$XML_DATA/prd.db`'
    )

    args = parser.parse_args()

    for filename in args.exposures:
        logger.info(
            '\n------'
            'Setting pointing for {}'.format(filename)
        )
        try:
            add_wcs(
                filename,
                siaf_path=args.siaf,
                strict_time=args.strict_time,
                strict_pointing=not args.allow_default
            )
        except ValueError as exception:
            logger.info('Cannot determine pointing information: ' + str(exception))
