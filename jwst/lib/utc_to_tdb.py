import logging

from astropy.io import fits
from astropy.time import Time

"""
If we can't import timeconversion, functions in this module will be called
instead, using astropy.coordinates for the location of the solar-system
barycenter, and a linear function based on header keywords for the location
of JWST with respect to the Earth.
"""
try:
    from jwst import timeconversion
    USE_TIMECONVERSION = True
except Exception:
    import numpy as np
    import astropy.coordinates as acoord
    import astropy.constants
    USE_TIMECONVERSION = False

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# This is an interface to call
#    timeconversion.compute_bary_helio_time(targetcoord, times)

def utc_tdb(filename, update_tdb=True, update_velosys=True,
            delta_t=1000., ndecimals=2):
    """Convert start, mid, end times from UTC to TDB.

    Parameters
    ----------
    filename: str
        Name of an input FITS file containing an INT_TIMES table.

    update_tdb: bool
        If True, the UTC times in the INT_TIMES table will be converted
        to TDB, and the three TDB columns in the INT_TIMES table will be
        updated with the computed values.

    update_velosys: bool
        If True, the radial velocity of JWST will be computed at the
        middle of the exposure, and the value in meters / second will
        be assigned to keyword VELOSYS in the first SCI extension.

    delta_t: float
        This is used only if `update_velosys` is True.  The TDB times
        corresponding to two UTC times (that are centered on the middle
        of the exposure and separated by `delta_t` seconds) will be
        computed, and the radial velocity will be computed from the
        difference between the two TDB times.

    ndecimals: int
        This is used only if `update_velosys` is True.  The radial
        velocity will be rounded to `ndecimals` decimal places.

    Returns
    -------
    tuple of three numpy arrays
        If `update_tdb` is True, the returned value will be the TDB time
        or times, expressed as MJD and with time scale TDB; otherwise,
        a tuple of three zeros (floats) will be returned.
    """

    if not update_tdb and not update_velosys:
        log.warning("Both update_tdb and update_velosys are False; "
                    "there's nothing to do.")
        return (0., 0., 0.)

    log.info("Processing file %s", filename)
    fd = fits.open(filename, mode="update")

    targetcoord = (fd[0].header["targ_ra"], fd[0].header["targ_dec"])

    if update_velosys:
        # Compute VELOSYS at the middle of the exposure.
        try:
            expmid = fd[0].header["expmid"]
            update_keyword = True
        except KeyError as e:
            log.warning(str(e))
            log.warning("Can't update VELOSYS without the time.")
            update_keyword = False

        if update_keyword:
            half_dt = delta_t / (86400. * 2.)       # and convert to days
            tt_rv_times = np.array([to_tt(expmid - half_dt),
                                    to_tt(expmid + half_dt)])
            if USE_TIMECONVERSION:
                tdb_rv_times = timeconversion.compute_bary_helio_time(
                                    targetcoord, tt_rv_times)[0]
            else:
                (eph_time, jwst_pos, jwst_vel) = get_jwst_keywords(fd)
                jwstpos_rv = linear_pos(tt_rv_times, eph_time,
                                        jwst_pos, jwst_vel)
                tdb_rv_times = compute_bary_helio_time2(
                                    targetcoord, tt_rv_times, jwstpos_rv)[0]

            # This is almost exactly the same as delta_t, but in units of days.
            delta_tt = tt_rv_times[1] - tt_rv_times[0]
            # This will be close to delta_t (in days), but it can differ due to
            # the motion of the earth along the line of sight to or from the
            # target.
            delta_tdb = tdb_rv_times[1] - tdb_rv_times[0]
            time_diff = (delta_tt - delta_tdb) * 86400.     # seconds
            radial_velocity = time_diff * astropy.constants.c.value / delta_t
            radial_velocity = round(radial_velocity, ndecimals)
            log.info("radial velocity = {}".format(radial_velocity))

            scihdu = find_hdu(fd, "sci")
            fd[scihdu].header["velosys"] = \
                    (radial_velocity, "Radial velocity wrt barycenter [m / s]")

    if update_tdb:
        try:
            hdunum = find_hdu(fd, "int_times")
        except RuntimeError as e:
            log.warning(str(e))
            fd.close()
            return (0., 0., 0.)

        # TT, MJD
        tt_start_times = to_tt(fd[hdunum].data.field("int_start_MJD_UTC"))
        tt_mid_times = to_tt(fd[hdunum].data.field("int_mid_MJD_UTC"))
        tt_end_times = to_tt(fd[hdunum].data.field("int_end_MJD_UTC"))

        # Function compute_bary_helio_time returns both barycentric and
        # heliocentric times; the "[0]" extracts the former.
        if USE_TIMECONVERSION:
            log.debug("Using the timeconversion module.")
            tdb_start_times = timeconversion.compute_bary_helio_time(
                                    targetcoord, tt_start_times)[0]
            tdb_mid_times = timeconversion.compute_bary_helio_time(
                                    targetcoord, tt_mid_times)[0]
            tdb_end_times = timeconversion.compute_bary_helio_time(
                                    targetcoord, tt_end_times)[0]
        else:
            log.warning("Couldn't import the timeconversion module;")
            log.warning("using astropy.coordinates, and "
                        "JWST position and velocity keywords.")
            (eph_time, jwst_pos, jwst_vel) = get_jwst_keywords(fd)
            jwstpos = linear_pos(tt_start_times, eph_time, jwst_pos, jwst_vel)
            tdb_start_times = compute_bary_helio_time2(
                                    targetcoord, tt_start_times, jwstpos)[0]
            jwstpos = linear_pos(tt_mid_times, eph_time, jwst_pos, jwst_vel)
            tdb_mid_times = compute_bary_helio_time2(
                                    targetcoord, tt_mid_times, jwstpos)[0]
            jwstpos = linear_pos(tt_end_times, eph_time, jwst_pos, jwst_vel)
            tdb_end_times = compute_bary_helio_time2(
                                    targetcoord, tt_end_times, jwstpos)[0]

        try:
            # TDB, MJD
            fd[hdunum].data.field("int_start_BJD_TDB")[:] = \
                        tdb_start_times.copy()
            fd[hdunum].data.field("int_mid_BJD_TDB")[:] = tdb_mid_times.copy()
            fd[hdunum].data.field("int_end_BJD_TDB")[:] = tdb_end_times.copy()
        except KeyError:
            log.warning("One or more of the *BJD_TDB columns do not exist,")
            log.warning("so the INT_TIMES table was not updated.")
    else:
        (tdb_start_times, tdb_mid_times, tdb_end_times) = (0., 0., 0.)

    fd.close()

    return (tdb_start_times, tdb_mid_times, tdb_end_times)


def find_hdu(fd, extname):
    """Find the HDU with name extname.

    Parameters
    ----------
    fd: fits HDUList object
        List of HDUs in the input file.

    extname: str
        The extension name to be found in fd.

    Returns
    -------
    hdunum: int
        The HDU number (0 is the primary HDU) with EXTNAME = `extname`.
    """

    extname_uc = extname.upper()
    nhdu = len(fd)
    hdunum = None
    for i in range(1, nhdu):
        hdr = fd[i].header
        if "EXTNAME" in hdr and hdr["EXTNAME"] == extname_uc:
            hdunum = i
            break
    if hdunum is None or len(fd[hdunum].data) < 1:
        fd.close()
        raise RuntimeError("An extension with name {} is required."
                           .format(extname))

    return hdunum


def to_tt(utc_time):
    """Convert UTC to TT.

    Parameters
    ----------
    utc_time: float or numpy array
        Time or array of times, expressed as MJD with time scale UTC.

    Returns
    -------
    float or numpy array
        The time or times, expressed as MJD but with time scale TT.
    """

    temp = Time(utc_time, format="mjd", scale="utc")

    return temp.tt.value


""" ###
The following functions are only used if the timeconversion module could not
be imported, i.e. if USE_TIMECONVERSION is False.
"""

def get_jwst_position2(times, jwstpos, use_jpl_ephemeris=False):
    '''
    This returns the pair of relative positions from
    the barycenter and heliocenter to JWST in that order
    as a tuple of two arrays, each of shape (len(times), 3).
    '''

    t = Time(times, format="mjd", scale="tt")

    if use_jpl_ephemeris:
        from astropy.coordinates import solar_system_ephemeris
        solar_system_ephemeris.set('jpl')

    # Vectors from the solar-system barycenter to the center of the Earth.
    bary_earth = acoord.get_body_barycentric("earth", t)

    # Vectors from the solar-system barycenter to the center of the Sun.
    bary_sun = acoord.get_body_barycentric("sun", t)

    # Vectors from the center of the Sun to the center of the Earth.
    sun_earth = bary_earth - bary_sun

    # Convert to ordinary numpy arrays of 3-element vectors, in km.

    barysun_centerearth_pos = np.empty((len(t), 3), dtype=np.float64)
    barysun_centerearth_pos[:, 0] = bary_earth.x.si.value / 1000.
    barysun_centerearth_pos[:, 1] = bary_earth.y.si.value / 1000.
    barysun_centerearth_pos[:, 2] = bary_earth.z.si.value / 1000.

    centersun_centerearth_pos = np.empty((len(t), 3), dtype=np.float64)
    centersun_centerearth_pos[:, 0] = sun_earth.x.si.value / 1000.
    centersun_centerearth_pos[:, 1] = sun_earth.y.si.value / 1000.
    centersun_centerearth_pos[:, 2] = sun_earth.z.si.value / 1000.

    centerearth_jwst = jwstpos

    return (barysun_centerearth_pos + centerearth_jwst), \
            (centersun_centerearth_pos + centerearth_jwst)

def get_target_vector2(targetcoord):
    '''
    returns a unit vector given ra and dec that astropy coordinates can handle
    '''
    ra, dec = targetcoord
    coord = acoord.SkyCoord(ra, dec, distance=1, frame='icrs', unit='deg')
    cartcoord = coord.represent_as(acoord.CartesianRepresentation)
    x = cartcoord.x.value
    y = cartcoord.y.value
    z = cartcoord.z.value
    vector = np.array([x, y, z])
    return vector / np.sqrt((vector**2).sum())

def compute_bary_helio_time2(targetcoord, times, jwstpos,
                            use_jpl_ephemeris=False):
    '''
    The end point computational routine to compute the distance of JWST
    to the sun (or barycenter) projected onto the unit vector to the
    target and determine the relative light travel time that results.

    times is assumed to be MJD_TT, and it can be either a scalar value
    or an array.

    jwstpos should be a 3-element vector (list, tuple, ndarray), or an array
    of such vectors, in km.  The shape of jwstpos
    should be (3,) or (1, 3) or (len(times), 3).
    jwstpos overrides what would have been obtained from the JWST ephemeris.
    This is useful for regression testing.
    '''
    tvector = get_target_vector2(targetcoord)
    jwst_bary_vectors, jwst_sun_vectors = get_jwst_position2(
                        times, jwstpos, use_jpl_ephemeris)
    proj_bary_dist = (tvector * jwst_bary_vectors).sum(axis=-1)
    proj_sun_dist = (tvector * jwst_sun_vectors).sum(axis=-1)
    cspeed = astropy.constants.c.value / 1000.
    return times + proj_bary_dist / cspeed / 86400., times + proj_sun_dist / cspeed / 86400.


def get_jwst_keywords(fd):

    # TT MJD time at which JWST has position jwst_pos and velocity jwst_vel.
    eph_time = to_tt(fd[0].header['eph_time'])

    jwst_pos = np.array((fd[0].header['jwst_x'],
                         fd[0].header['jwst_y'],
                         fd[0].header['jwst_z']), dtype=np.float64)

    jwst_vel = np.array((fd[0].header['jwst_dx'],
                         fd[0].header['jwst_dy'],
                         fd[0].header['jwst_dz']), dtype=np.float64)

    return(eph_time, jwst_pos, jwst_vel)


def linear_pos(tt_times, eph_time, jwst_pos, jwst_vel):
    """Compute JWST position as a linear function of time.

    Parameters
    ----------
    tt_times: float or ndarray
        Time or array of times, expressed as MJD with time scale TT.

    eph_time: float
        The time (MJD, TT) at which JWST had the position and velocity
        given by `jwst_pos` and `jwst_vel` respectively.

    jwst_pos: ndarray
        A three-element vector in rectangular coordinates, giving the
        position of JWST in km at time `eph_time`.

    jwst_vel: ndarray
        A three-element vector in rectangular coordinates, giving the
        velocity of JWST in km/s at time `eph_time`.

    Returns
    -------
    ndarray
        An array of shape (n_times, 3), giving the position of JWST at
        each of the times in `tt_times`.  `n_times` is the length of
        `tt_times`, or 1 if `tt_times` is a float.
    """

    jwst_pos_2d = jwst_pos.reshape((3, 1))
    jwst_vel_2d = jwst_vel.reshape((3, 1))
    dt = tt_times - eph_time
    jwstpos = (jwst_pos_2d + dt * jwst_vel_2d).transpose()

    return jwstpos
