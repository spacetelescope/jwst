"""Set Telescope Pointing from quaternions"""
import logging
import numpy as np
from numpy import cos, sin

from namedlist import namedlist

from jwst.lib.engdb_tools import (
    ENGDB_BASE_URL,
    ENGDB_Service,
)

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# Degree, radian, angle transformations
R2D = 180./np.pi
D2R = np.pi/180.
A2R = D2R/3600.
R2A = 3600.*R2D

# SIAF container
SIAF = namedlist(
    'SIAF',
    ['v2ref', 'v3ref', 'v3idlyang', 'vparity'],
    default=None
)

# Pointing container
Pointing = namedlist(
    'Pointing',
    ['q', 'j2fgs_matrix', 'fsmcorr', 'obstime'],
    default=None
)

# WCS reference container
WCSRef = namedlist(
    'WCSRef',
    ['ra', 'dec', 'pa'],
    default=None
)


def update_wcs(hdu, default_pa_v3=0.):
    """Update WCS pointing information

    Given a FITS Header Data Unit, determine the simple WCS parameters
    from the SIAF keywords in the HDU and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    hdu: astropy.io.fits HDU-like
        The Header Data Unit to update, in-place

    default_pa_v3: float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.
    """

    # Get the SIAF and observation parameters
    header = hdu.header
    obsstart = float(header['EXPSTART'])
    obsend = float(header['EXPEND'])
    siaf = SIAF(
        v2ref=float(header['V2_REF']),
        v3ref=float(header['V3_REF']),
        v3idlyang=float(header['V3I_YANG']),
        vparity=int(header['VPARITY'])
    )

    # Get the pointing information
    try:
        pointing = get_pointing(obsstart, obsend)
    except ValueError as exception:
        ra = header['TARG_RA']
        dec = header['TARG_DEC']
        roll = default_pa_v3

        logger.warning(
            'Cannot retrieve telescope pointing.'
            '\n{}'
            '\nUsing TARG_RA={}, TARG_DEC={} and PA_V3={} '
            'to set pointing.'.format(exception, ra, dec, roll)
        )

        wcsinfo = (ra, dec, roll)
        vinfo = (ra, dec, roll)

    else:
        # compute relevant WCS information
        logger.info('Successful read of engineering quaternions:')
        logger.info('\tPointing = {}'.format(pointing))
        wcsinfo, vinfo = calc_wcs(siaf, pointing)

    logger.info('Aperture WCS info: {}'.format(wcsinfo))
    logger.info('V1 WCS info: {}'.format(vinfo))

    # Update V1 pointing
    header['RA_V1'] = vinfo.ra
    header['DEC_V1'] = vinfo.dec
    header['PA_V3'] = vinfo.pa

    # Update Aperture pointing
    header['CRVAL1'] = wcsinfo.ra
    header['CRVAL2'] = wcsinfo.dec
    header['PC1_1'] = -np.cos(wcsinfo.pa * D2R)
    header['PC1_2'] = np.sin(wcsinfo.pa * D2R)
    header['PC2_1'] = np.sin(wcsinfo.pa * D2R)
    header['PC2_2'] = np.cos(wcsinfo.pa * D2R)
    header['RA_REF'] = wcsinfo.ra
    header['DEC_REF'] = wcsinfo.dec
    header['ROLL_REF'] = compute_local_roll(
        vinfo.pa, wcsinfo.ra, wcsinfo.dec, siaf.v2ref, siaf.v3ref
    )
    header['WCSAXES'] = len(header['CTYPE*'])

    logger.info('Header successfully updated')


def calc_wcs(siaf, pointing):
    """Transform from the given SIAF information and Pointing
    the aperture and V1 wcs

    Parameters
    ----------
    siaf: SIAF
        The SIAF transformation. See ref:`Notes` for further details
    pointing: Pointing
        The telescope pointing. See ref:`Notes` for further details

    Returns
    -------
    (wcsinfo, vinfo): (WCSRef, WCSRef)
        A 2-tuple is returned with the WCS pointing for
        the aperture and the V1 axis

    Notes
    -----

    The SIAF information is as follows:

    v2ref (arcsec), v3ref (arcsec), v3idlyang (deg), vidlparity
    (+1 or -1), are the relevant siaf parameters. The assumed
    units are shown in parentheses.

    It is assumed that the siaf ref position is the corresponding WCS
    reference position.

    The `Pointing` information is as follows:

    Parameter q is the SA_ZATTEST<n> engineering parameters where
    n ranges from 1 to 4.

    Parameter j2fgs_matrix is the transformation matrix specified by
    engineering parameters SA_ZRFGS2J<n><m> where both n and m range
    from 1 to 3. This is to be provided as a 1d list using this order:
    11, 21, 31, 12, 22, 32, 13, 23, 33

    Parameter fsmcorr are two values provided as a list consisting of:
    [SA_ZADUCMDX, SA_ZADUCMDY]
    """

    # Determine the ECI to J-frame matrix
    m_eci2j = calc_eci2j_matrix(pointing.q)

    # Calculate the J-frame to FGS! ICS matrix
    m_j2fgs1 = calc_j2fgs1_matrix(pointing.j2fgs_matrix)

    # Calculate the FSM corrections to the SI_FOV frame
    m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(pointing.fsmcorr)

    # Calculate the FGS1 ICS to SI-FOV matrix
    m_fgs1_to_sifov = calc_fgs1_to_sifov_mastrix()

    # Calculate ECI to SI FOV
    m_eci2sifov = np.dot(
        m_sifov_fsm_delta,
        np.dot(
            m_fgs1_to_sifov,
            np.dot(
                m_j2fgs1,
                m_eci2j
            )
        )
    )

    # Calculate SI FOV to V1 matrix
    m_sifov2v = calc_sifov2v_matrix(PITCH_V2)

    # Calculate the complete transform to the V1 reference
    m_eci2v = np.dot(m_sifov2v, m_eci2sifov)

    # Calculate the V1 wcs information
    vinfo = WCSRef()

    # V1 RA/Dec is the first row of the transform
    vinfo.ra, vinfo.dec = vector_to_ra_dec(m_eci2v[0])

    # V3 is the third row of the transformation
    v3info = WCSRef()
    v3info.ra, v3info.dec = vector_to_ra_dec(m_eci2v[2])

    # Calculate the V3 position angle
    vinfo.pa = calc_position_angle(vinfo, v3info)

    # Calculate the SIAF transform matrix
    m_v2siaf = calc_v2siaf_matrix(siaf)

    # Calculate the full ECI to SIAF transform matrix
    m_eci2siaf = np.dot(m_v2siaf, m_sifov2v)

    # Calculate the Aperture WCS
    wcsinfo = calc_aperture_wcs(siaf, m_eci2siaf)

    # Convert all WCS to degrees
    wcsinfo = WCSRef(
        ra=wcsinfo.ra * R2D,
        dec=wcsinfo.dec * R2D,
        pa=wcsinfo.pa * R2D
    )
    vinfo = WCSRef(
        ra=vinfo.ra * R2D,
        dec=vinfo.dec * R2D,
        pa=vinfo.pa * R2D
    )

    # That's all folks
    return (wcsinfo, vinfo)


def get_pointing(obsstart, obsend, result_type='first'):
    """
    Get telescope pointing engineering data.

    Parameters
    ----------
    obsstart, obsend: float
        MJD observation start/end times

    result_type: str
        What to return. Possible options are:
            `first`: Return the first non-zero matricies
            `all`: Return all non-zero matricies within
                   the given range.

    Returns
    -------
    pointing: Pointing or [Pointing(, ...)]
        The engineering pointing parameters.
        If the `result_type` is `all`, a list
        of pointings will be returned

    Raises
    ------
    ValueError
        Cannot retrieve engineering information

    Notes
    -----
    For the moment, the first found values will be used.
    This will need be re-examined when more information is
    available.
    """
    logger.info(
        'Determining pointing between observations times (mjd):'
        '\n\tobsstart = {}'
        '\n\tobsend = {}'.format(obsstart, obsend)
    )
    logger.info(
        'Querying engineering DB: {}'.format(ENGDB_BASE_URL)
    )
    try:
        engdb = ENGDB_Service()
    except Exception as exception:
        raise ValueError(
            'Cannot open engineering DB connection'
            '\nException: {}'.format(
                exception
            )
        )
    params = {
        'SA_ZATTEST1':  None,
        'SA_ZATTEST2':  None,
        'SA_ZATTEST3':  None,
        'SA_ZATTEST4':  None,
        'SA_ZRFGS2J11': None,
        'SA_ZRFGS2J21': None,
        'SA_ZRFGS2J31': None,
        'SA_ZRFGS2J12': None,
        'SA_ZRFGS2J22': None,
        'SA_ZRFGS2J32': None,
        'SA_ZRFGS2J13': None,
        'SA_ZRFGS2J23': None,
        'SA_ZRFGS2J33': None,
        'SA_ZADUCMDX':  None,
        'SA_ZADUCMDY':  None,
    }
    for param in params:
        try:
            params[param] = engdb.get_values(
                param, obsstart, obsend, time_format='mjd', include_obstime=True
            )
        except Exception as exception:
            raise ValueError(
                'Cannot retrive {} from engineering.'
                '\nFailure was {}'.format(
                    param,
                    exception
                )
            )

    # Find the first set of non-zero values
    results = []
    for idx in range(len(params['SA_ZATTEST1'])):
        values = [
            params[param][idx].value
            for param in params
        ]
        if any(values):
            pointing = Pointing()

            # The tagged obstime will come from the SA_ZATTEST1 mneunonic
            pointing.obstime = params['SA_ZATTEST1'][idx].obstime

            # Fill out the matricies
            pointing.q = np.array([
                params['SA_ZATTEST1'][idx].value,
                params['SA_ZATTEST2'][idx].value,
                params['SA_ZATTEST3'][idx].value,
                params['SA_ZATTEST4'][idx].value,
            ])

            pointing.j2fgs_matrix = np.array([
                params['SA_ZRFGS2J11'][idx].value,
                params['SA_ZRFGS2J21'][idx].value,
                params['SA_ZRFGS2J31'][idx].value,
                params['SA_ZRFGS2J12'][idx].value,
                params['SA_ZRFGS2J22'][idx].value,
                params['SA_ZRFGS2J32'][idx].value,
                params['SA_ZRFGS2J13'][idx].value,
                params['SA_ZRFGS2J23'][idx].value,
                params['SA_ZRFGS2J33'][idx].value,
            ])

            pointing.fsmcorr = np.array([
                params['SA_ZADUCMDX'][idx].value,
                params['SA_ZADUCMDY'][idx].value,

            ])

            results.append(pointing)

            # Short circuit if all we're looking for is the first.
            if result_type == 'first':
                break

    if not len(results):
        raise ValueError(
                'No non-zero quanternion found '
                'in the DB between MJD {} and {}'.format(obsstart, obsend)
            )

    if result_type == 'first':
        return results[0]
    else:
        return results
