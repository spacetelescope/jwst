#!/usr/bin/env python

# Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

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

# STDLIB
import argparse
import sys
import warnings

# LOCAL
from jwst import fits_generator


def main():
    if '--version' in sys.argv:
        sys.stdout.write('%s\n' % fits_generator.__version__)
        sys.exit(0)

    parser = argparse.ArgumentParser(
        """Adds data to a raw FITS file to correspond to a given file type."""
    )
    parser.add_argument(
        'files',
        metavar='file',
        nargs='+',
        help="""A list of files, containing 1 input fits file, and 0 or
                more data files to pull data from.  Most output file types
                require at least a target and program data file.""",
    )
    parser.add_argument(
        '--type',
        '-t',
        nargs=1,
        default=None,
        help="""Output file type.  If not provided, the result will be a level
        1B FITS file of the appropriate instrument type, determined
        automatically from the input FITS file.""",
    )
    parser.add_argument(
        '--output',
        '-o',
        nargs=1,
        default=None,
        help="""Output filename.  If not provided, a default one will be
                generated based on the content.""",
    )
    parser.add_argument(
        '--clobber',
        '-c',
        action='store_true',
        help='Overwrite the output file if it exists.',
    )
    parser.add_argument(
        '--verify', '-v', action='store_true', help='Verify the output after generation'
    )

    args = parser.parse_args()

    if args.type is not None:
        filetype = args.type[0]
    else:
        filetype = None

    hdulist = fits_generator.generate(args.files, filetype=filetype, verify=args.verify)
    filename = args.output
    if filename is None:
        try:
            filename = fits_generator.guess_filename(hdulist)
        except ValueError as e:
            warnings.warn(str(e))
            filename = 'output.fits'
    hdulist.writeto(filename, clobber=args.clobber, output_verify='silentfix')


if __name__ == '__main__':
    main()
