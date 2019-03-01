#!/usr/bin/env python

# Copyright (C) 2015-2016 Association of Universities for Research in Astronomy (AURA)

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

'''
This is a script to update JWST data files to include
barycentric and heliocentric times as computed from the spacecraft
times. It also updates the extension that contains the times
for each group by adding the necessary columns to that table.
'''

import argparse
import logging

import astropy.coordinates
import astropy.io.fits as fits
import astropy.time
import astropy.units as u
import numpy as np

from jwst import timeconversion

# Setup the logging.
logging.basicConfig(level=logging.INFO)


def set_bary_helio_times(filename, jwstpos=None):
    '''
    Compute the barycentric and heliocentric times for the
    given file and update the contents of the file to contain
    these values.
    '''
    logging.info('Starting set_bary_helio_times task')
    hdul = fits.open(filename, mode='update')
    pheader = hdul[0].header
    # Obtain the necessary info from the header
    targcoord = astropy.coordinates.SkyCoord(
        ra=pheader['TARG_RA'],
        dec=pheader['TARG_DEC'],
        frame='fk5',
        unit=(u.hourangle, u.deg)
    )
    starttime = pheader['EXPSTART']
    middletime = pheader['EXPMID']
    endtime = pheader['EXPEND']
    times = astropy.time.Time(
        np.array([starttime, middletime, endtime]),
        format='mjd', scale='utc').tt.mjd
    ((bstrtime, bmidtime, bendtime), (hstrtime, hmidtime, hendtime)) = \
        timeconversion.compute_bary_helio_time(
            targetcoord=(targcoord.ra.degree, targcoord.dec.degree),
            times=times
        )
    pheader['BSTRTIME'] = bstrtime
    pheader['BMIDTIME'] = bmidtime
    pheader['BENDTIME'] = bendtime
    pheader['BARTDELT'] = (bstrtime - starttime) * 86400.
    pheader['HSTRTIME'] = hstrtime
    pheader['HMIDTIME'] = hmidtime
    pheader['HENDTIME'] = hendtime
    pheader['HELIDELT'] = (hstrtime - starttime) * 86400.

    # Now modify the table
    try:
        tabhdu = hdul['GROUP']
    except KeyError:
        logging.info('No GROUP extension found. Ignoring GROUP calculations')
        pass
    else:
        logging.info('Calculating GROUP times')
        tendtimes = tabhdu.data['group_end_time']

        # replace colon as separator between date and time to be consistent
        tendtimeslist = [item[:10] + 'T' + item[11:] for item in tendtimes]
        astropy_endtimes = astropy.time.Time(
            tendtimeslist, format='isot', scale='utc'
        )
        mjdtimes = astropy_endtimes.tt.mjd
        try:
            btimes, htimes = timeconversion.compute_bary_helio_time(
                targetcoord=(targcoord.ra.degree, targcoord.dec.degree),
                times=mjdtimes
            )
        except Exception as exception:
            # Cleanup before exiting
            hdul.flush()
            hdul.close()
            raise Exception('Error in calculating times.') from exception
        else:
            bcol = fits.Column(
                name='bary_end_time', format='D', unit='MJD', array=btimes
            )
            hcol = fits.Column(
                name='helio_end_time', format='D', unit='MJD', array=htimes
            )

        binhdu = fits.BinTableHDU.from_columns(
            tabhdu.columns + fits.ColDefs([bcol, hcol])
        )
        binhdu.header['EXTNAME'] = 'GROUP'
        hdul['GROUP'] = binhdu

    hdul.flush()
    hdul.close()
    logging.info('Completed set_bary_helio_times task')


def main(args):

    if args.jwstpos:
        jwstposstr = args.jwstpos.split(',')
        jwstpos = [float(item) for item in jwstposstr]
    else:
        jwstpos = None
    for filename in args.filenames:
        set_bary_helio_times(filename, jwstpos=jwstpos)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add corresponding barycentric and heliocentric times to JWST data files"
    )
    parser.add_argument('filenames', metavar='Filename', nargs='+')
    parser.add_argument(
        '--jwstpos',
        help="x,y,z JWST position relative to Earth's center in km as comma separated string of values"
    )
    args = parser.parse_args()
    main(args)
