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

"""
This contains the functions that actually compute the time travel difference
between two frames of reference.
"""

import os
import os.path
import numpy as np
from jplephem.spk import SPK
import astropy.time as atime
import astropy.coordinates as acoord
import astropy.constants
import pymssql
import scipy.interpolate as sciint
import logging

# Configure logging
logger = logging.getLogger(__name__)

JDOFFSET = 2400000.5

# Find path to ephemeris from environmental variable.
try:
    EPHEMERIS_PATH = os.environ['JPL_EPHEMERIS']
except KeyError:
    raise KeyError('Environmental variable JPL_EPHEMERIS not found')

# The following three functions are only needed if reading from a file
# Currently it is obtained from a DMS database


def read_jwst_ephemeris(times):
    """Compute the x,y,z positions of JWST.

    Parameters
    ----------
    times : ndarray, 1-D, float
        The times (Julian Day Number) at which the positions of JWST will
        be calculated.

    Returns
    -------
    outpos : ndarray, 2-D, float
        This array has shape (nt, 3), where `nt` is the length of the input
        `times`.  For each index `i` in `times`, `outpos[i, :]` is the
        position (rectangular coordinates) of JWST in kilometers with
        respect to the center of the Earth.
    """

    f = open(os.path.join(EPHEMERIS_PATH, 'ephem_10days.txt'), 'r')
    # First get header info.
    mdict = get_jwst_ephem_metadata(f)
    if ((mdict['CENTER_NAME'] != 'Earth') or
        (mdict['OBJECT_NAME'] != 'JWST') or
        (mdict['REF_FRAME'] != 'EME2000') or
        (mdict['TIME_SYSTEM'] != 'UTC')):
        raise RuntimeError("JWST ephemeris file metadata doesn't match expectations")
    # Determine size of time data; used to estimate location of desired times
    start = f.tell()
    time1, position1 = parse_jwst_ephem_line(f.readline())
    time2, position2 = parse_jwst_ephem_line(f.readline())
    f.seek(1, os.SEEK_END)
    end = f.tell()
    starttime = atime.Time(mdict['START_TIME']).jd
    endtime = atime.Time(mdict['STOP_TIME']).jd
    interval = (time2 - time1) * 1440 # in minutes
    if interval > 10:
        raise ValueError("time interval in JWST ephemeris file is too large")
    # Determine min, max of times sought
    mintime = times.min()
    maxtime = times.max()
    if (mintime < starttime) or (maxtime > endtime):
        raise ValueError(
            'Some of times provided are out of the range of times in the JWST ephemeris file')
    # Determine fractional locations
    minpos = int((mintime - starttime) / (endtime - starttime) * (end - start) + start - 5000)
    maxpos = int((maxtime - starttime) / (endtime - starttime) * (end - start) + start + 5000)
    if minpos < start:
        minpos = start
    if maxpos > end:
        maxpos = end
    f.seek(minpos)
    text = f.read(maxpos - minpos)
    f.close()

    # Convert this to lines
    lines = text.split('\n')
    # Throw away truncated lines
    if not lines[0].startswith('20') or len(lines[0] < 125):
        lines = lines[1:]
    if len(lines[-1]) < 80: # Can be looser here since
                            # we don't care about the last 3 numbers.
        lines = lines[:-1]
    # Turn into numpy arrays
    nrecords = len(lines)
    etimes = np.zeros(nrecords)
    epos = np.zeros((nrecords, 3))
    for i, line in enumerate(lines):
        etimes[i], epos[i, :] = parse_jwst_ephem_line(line)
    outpos = np.zeros((len(times), 3))
    for i, t in enumerate(times):
        tdiff = np.abs(etimes - t)
        tdiffmin = np.where(tdiff == tdiff.min())[0][0]
        if tdiff[tdiffmin] * 1440 > 6:
            raise RuntimeError(
                "didn't find time close enough; check the code for sufficient margin")
        outpos[i, :] = epos[tdiffmin, :]
    return outpos


def get_jwst_ephem_metadata(fh):
    """Parse the metadata from the ephemeris file.

    Parameters
    ----------
    fh: file handle
        File handle for the JWST 10-day ephemeris file.  The lines in this
        file have the form 'keyword = value'.

    Returns
    -------
    mdict: dictionary
        Each key is the word before '=' in a line read from the input file,
        and the value is the string that follows '='.
    """

    mdict = {}
    while True:
        line = fh.readline()
        if line.startswith('META_STOP'):
            return mdict
        elif line.startswith('COMMENT') or line.startswith('META_START'):
            pass
        else:
            lparts = [item.strip() for item in line.split('=')]
            mdict[lparts[0]] = lparts[1]


def parse_jwst_ephem_line(line):
    """Parse the relevant elements of a JWST ephemeris file line.

    Parameters
    ----------
    line : str
        One line that was read from a JWST ephemeris file.

    Returns
    -------
    float
        The time (Julian Day Number) at which JWST had the specified
        coordinates.

    ndarray, 2-D, float
        The coordinates (km) of JWST at the specified time.
    """

    lparts = line.split()
    time = atime.Time(lparts[0])
    position = np.array([float(item) for item in lparts[1:4]])
    return time.jd, position

#-------end of file based jwst ephem routines

def jwst_ephem_interp(t, padding=3):
    """Interpolate within an array of JWST positions.

    Given the values of time (etab[:, 0]), x, y, z (etab[:, 1], etab[:, 2],
    etab[:, 3]) obtained from an ephemeris, apply cubic interpolation to
    obtain x, y, z values for the requested time(s).

    Parameters
    ----------
    t: ndarray, 1-D, float
        An array of times in MJD.

    padding : int
        The number of entries required before and after the time range of
        interest.  The default is 3.

    Returns
    -------
    tuple of three floats
        The coordinates (km) of JWST at the specified time.
    """

    logger.debug('times="{}"'.format(t))
    etab = get_jwst_ephemeris()

    # Select only the portion of the table with relevant times.
    mask_min = etab[:, 0] >= t.min()
    mask_max = etab[:, 0] <= t.max()
    mask = np.logical_and(mask_min, mask_max)
    try:
        first_idx = np.nonzero(mask_min)[0][0]
    except IndexError:
        first_idx = len(mask_min) - 1
    try:
        last_idx = np.nonzero(mask_max)[0][-1]
    except IndexError:
        last_idx = 0
    first_idx -= padding
    last_idx += padding
    fit_type = 'cubic'
    if first_idx < 0 or last_idx >= len(mask):
        logger.warning('Times extend outside range of ephemeris.'
                       ' Extrapolation will be used.')
        fit_type = 'linear'
        first_idx = max(first_idx, 0)
        last_idx = min(last_idx, len(mask) - 1)
    logger.debug('idxs = {} {}'.format(first_idx, last_idx))
    for diff in range(padding):
        mask[first_idx + diff] = True
        mask[last_idx - diff] = True

    roi = etab[mask]
    logger.debug('roi = {}'.format(roi.shape))

    extrapolate = 'extrapolate' if fit_type == 'linear' else None
    fx = sciint.interp1d(
        roi[:, 0], roi[:, 1], kind=fit_type,
        assume_sorted=True, fill_value=extrapolate
    )
    fy = sciint.interp1d(
        roi[:, 0], roi[:, 2], kind=fit_type,
        assume_sorted=True, fill_value=extrapolate
    )
    fz = sciint.interp1d(
        roi[:, 0], roi[:, 3], kind=fit_type,
        assume_sorted=True, fill_value=extrapolate
    )
    return fx(t), fy(t), fz(t)


# xxx Is this comment correct?
# the following is only needed if obtaining from a file instead of DB, not used yet
def get_jwst_ephemeris():
    """Read the contents of the JWST database.

    Extracts predicted JWST ephemeris from DMS database (the whole thing
    for now) and returns a numpy 2d array.

    Returns
    -------
    etab : ndarray, 2-D, float
        The contents of the DMS database for JWST times and positions.
        The shape is (n, 7), where `n` is the number of entries in the
        database.  For each entry `i`, the values are as follows:
        etab[i, 0] is the time
        etab[i, 1] is the X component of JWST position
        etab[i, 2] is the Y component of JWST position
        etab[i, 3] is the Z component of JWST position
        etab[i, 4] is the X component of the velocity of JWST
        etab[i, 5] is the Y component of the velocity of JWST
        etab[i, 6] is the Z component of the velocity of JWST
    """

    if 'METRICS_SERVER' in os.environ:
        eserver = os.environ['METRICS_SERVER']
    else:
        eserver = 'JWDMSDEVDBVM1'
    if 'METRICS_DB' in os.environ:
        edb = os.environ['METRICS_DB']
    else:
        edb = 'jwdpmetrics5'
    logger.info(
        'Ephemeris connect info:'
        ' eserver={}'
        ' edb={}'.format(eserver, edb)
    )
    logger.debug(
        'Ephemeris connect info:'
        ' eserver={}'
        ' edb={}'.format(eserver, edb)
    )
    conn = pymssql.connect(server=eserver, database=edb)
    cur = conn.cursor()
    cur.execute('select ephem_time,jwst_x,jwst_y,jwst_z,jwst_dx,jwst_dy,jwst_dz from predictephemeris')
    etab = np.array(cur.fetchall())
    return etab


def get_jwst_position(times, jwstpos, debug=False):
    """Get the position of JWST.

    This returns the pair of relative positions to
    the barycenter and heliocenter in that order
    as a tuple of two arrays with 3 components

    Parameters
    ----------
    times : float, or 1-D ndarray of float
        A time or an array of times, in MJD (TT).

    jwstpos : 1-D array, or None
        `jwstpos` is a 3 element vector (list, tuple, whatever) in km if it
        is provided.

    debug : bool
        This is for testing.  If `debug` is `True`, the values returned
        will be the barycenter of the Earth and the center of the Sun (km),
        both with respect to the solar-system barycenter.

    Returns
    -------
    jwst_bary_vectors : ndarray, float
        Equatorial (ICRS) rectangular coordinates of JWST with respect to
        the solar-system barycenter, in km.  If `times` is a float,
        `jwst_bary_vectors` will have shape (3,); otherwise, it will have
        shape (3, nt), where `nt` is the length of the input `times`.

    jwst_sun_vectors : ndarray, float
        Equatorial (ICRS) rectangular coordinates of JWST with respect to
        the center of the Sun, in km.  If `times` is a float,
        `jwst_sun_vectors` will have shape (3,); otherwise, it will have
        shape (3, nt), where `nt` is the length of the input `times`.
    """

    ekernel = SPK.open(os.path.join(EPHEMERIS_PATH, 'de430.bsp'))
    # JPL ephemeris uses JD_tt based times
    barysun_baryearth_pos = ekernel[0, 3].compute(times + JDOFFSET)
    barysun_centersun_pos = ekernel[0, 10].compute(times + JDOFFSET)
    baryearth_centerearth_pos = ekernel[3, 399].compute(times + JDOFFSET)
    ekernel.close()
    if debug:
        return barysun_baryearth_pos, barysun_centersun_pos
    barysun_centerearth_pos = barysun_baryearth_pos + baryearth_centerearth_pos
    centersun_centerearth_pos = barysun_centerearth_pos - barysun_centersun_pos
    if jwstpos is not None:
        centerearth_jwst = jwstpos
    else:
        centerearth_jwst = jwst_ephem_interp(times)
    return (barysun_centerearth_pos + centerearth_jwst), \
            (centersun_centerearth_pos + centerearth_jwst)


def get_target_vector(targetcoord):
    """Convert spherical coordinates to rectangular.

    Returns a unit vector given ra and dec that astropy coordinates can handle.

    Parameters
    ----------
    targetcoord : two-element ndarray (or tuple, list), float
        The right ascension and declination of the target, in degrees.

    Returns
    -------
    ndarray, float
        A three-element ndarray, a unit vector pointing toward the target.
    """

    ra, dec = targetcoord
    coord = acoord.SkyCoord(ra, dec, distance=1, frame='icrs', unit='deg')
    cartcoord = coord.represent_as(acoord.CartesianRepresentation)
    x = cartcoord.x.value
    y = cartcoord.y.value
    z = cartcoord.z.value
    vector = np.array([x, y, z])
    return vector / np.sqrt((vector**2).sum())


def compute_bary_helio_time(targetcoord, times, jwstpos=None):
    """Correct the times of observation by subtracting light travel time.

    Extended summary
    ----------------
    The end point computational routine to compute the distance of JWST
    to the barycenter (or the Sun) projected onto the unit vector to the
    target and determine the relative light travel time that results.

    Parameters
    ----------
    targetcoord : two-element ndarray (or tuple, list), float
        The right ascension and declination of the target, in degrees.

    times : float, or 1-D ndarray of float
        A time or an array of times, in MJD (TT).

    jwstpos : list, tuple, ndarray, or None
        This will normally be None (which is the default).
        If not None, this is a three-element vector, giving the rectangular
        coordinates of JWST in km with respect to the center of the Earth.
        When jwstpos is provided, it overrides what would have been obtained
        from the JWST ephemeris. This is useful for regression testing.

    Returns
    -------
    barycentric times : float, or 1-D ndarray of float
        The time or times, corrected to the solar-system barycenter.

    heliocentric times : float, or 1-D ndarray of float
        The time or times, corrected to the center of the Sun.
    """

    tvector = get_target_vector(targetcoord)
    jwst_bary_vectors, jwst_sun_vectors = get_jwst_position(times, jwstpos)
    proj_bary_dist = (tvector * jwst_bary_vectors.transpose()).sum(axis=-1)
    proj_sun_dist = (tvector * jwst_sun_vectors.transpose()).sum(axis=-1)
    cspeed = astropy.constants.c.value / 1000.
    return times + proj_bary_dist / cspeed / 86400., times + proj_sun_dist / cspeed / 86400.
