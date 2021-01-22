"""Set Telescope Pointing from quaternions"""
from copy import copy
import logging
from math import (cos, sin)
import os.path
import sqlite3
from collections import defaultdict, namedtuple

from astropy.time import Time
from scipy.interpolate import interp1d
import numpy as np

from ..assign_wcs.util import update_s_region_keyword, calc_rotation_matrix
from ..assign_wcs.pointing import v23tosky
from ..datamodels import Level1bModel
from ..lib.engdb_tools import ENGDB_Service
from ..lib.pipe_utils import is_tso
from .exposure_types import IMAGING_TYPES, FGS_GUIDE_EXP_TYPES

TYPES_TO_UPDATE = set(list(IMAGING_TYPES) + FGS_GUIDE_EXP_TYPES)

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Conversion from seconds to MJD
SECONDS2MJD = 1 / 24 / 60 / 60

# Default transformation matricies
FGS12SIFOV_DEFAULT = np.array(
    [[0.9999994955442, 0.0000000000000, 0.0010044457459],
     [0.0000011174826, 0.9999993811310, -0.0011125359826],
     [-0.0010044451243, 0.0011125365439, 0.9999988766756]]
)

J2FGS_MATRIX_DEFAULT = np.array(
    [[0.999997425983907, 0, -0.002268926080840],
     [0., 1., 0.],
     [0.002268926080840, 0., 0.999997425983907]]
)

SIFOV2V_DEFAULT = np.array(
    [[0.99999742598, 0., 0.00226892608],
     [0., 1., 0.],
     [-0.00226892608, 0., 0.99999742598]]
)


# Degree, radian, angle transformations
R2D = 180./np.pi
D2R = np.pi/180.
A2R = D2R/3600.
R2A = 3600.*R2D

# SIAF container
# The names should correspond to the names in the ``wcsinfo`` schema.
# It is populated by the SIAF values in the PRD database based
# on APERNAME and UseAfterDate and used to populate the keywords
# in Level1bModel data models.
SIAF = namedtuple("SIAF", ["v2_ref", "v3_ref", "v3yangle", "vparity",
                           "crpix1", "crpix2", "cdelt1", "cdelt2",
                           "vertices_idl"])
# Set default values for the SIAF.
# Values which are needed by the pipeline are set to None which
# triggers a ValueError if missing in the SIAF database.
# Quantities not used by the pipeline get a default value -
# FITS keywords and aperture vertices.
SIAF.__new__.__defaults__ = (None, None, None, None, 0, 0, 3600, 3600,
                             (0, 1, 1, 0, 0, 0, 1, 1))

# Pointing container
Pointing = namedtuple("Pointing", ["q", "j2fgs_matrix", "fsmcorr", "obstime"])
Pointing.__new__.__defaults__ = (None, None, None, None)

# Transforms
Transforms = namedtuple("Transforms",
                        [
                            'm_eci2j',            # ECI to J-Frame
                            'm_j2fgs1',           # J-Frame to FGS1
                            'm_sifov_fsm_delta',  # FSM correction
                            'm_fgs12sifov',       # FGS1 to SIFOV
                            'm_eci2sifov',        # ECI to SIFOV
                            'm_sifov2v',          # SIFOV to V1
                            'm_eci2v',            # ECI to V
                            'm_v2siaf',           # V to SIAF
                            'm_eci2siaf'          # ECI to SIAF
                        ])
Transforms.__new__.__defaults__ = (None, None, None, None, None,
                                   None, None, None, None)
# WCS reference container
WCSRef = namedtuple('WCSRef', ['ra', 'dec', 'pa'])
WCSRef.__new__.__defaults__ = (None, None, None)


def add_wcs(filename, default_pa_v3=0., siaf_path=None, engdb_url=None,
            tolerance=60, allow_default=False, reduce_func=None,
            dry_run=False, **transform_kwargs):
    """Add WCS information to a FITS file.

    Telescope orientation is attempted to be obtained from
    the engineering database. Failing that, a default pointing
    is used based on proposal target.

    The FITS file is updated in-place.

    Parameters
    ----------
    filename: str
        The path to a data file.

    default_pa_v3: float
        The V3 position angle to use if the pointing information
        is not found.

    siaf_path: str or file-like object or None
        The path to the SIAF database.

    engdb_url: str or None
        URL of the engineering telemetry database REST interface.

    tolerance: int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default: bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func: func or None
        Reduction function to use on values.

    dry_run: bool
        Do not write out the modified file.

    transform_kwargs: dict
        Keyword arguments used by matrix calculation routines

    Notes
    -----
    This function adds absolute pointing information to the FITS files provided
    to it on the command line (one or more).

    Currently it only uses a constant value for the engineering keywords
    since the Engineering Database does not yet contain them.

    It starts by populating the headers with values from the SIAF database.
    It adds the following keywords to all files:

    V2_REF (arcseconds)
    V3_REF (arcseconds)
    VPARITY (+1 or -1)
    V3I_YANG (decimal degrees)

    In addition for ``IMAGING_MODES`` it adds these keywords:

    CRPIX1 (pixels)
    CRPIX2 (pixels)
    CDELT1 (deg/pix)
    CDELT2 (deg/pix)
    CUNIT1 (str)
    CUNIT2
    CTYPE1
    CTYPE2
    WCSAXES

    The keywords computed and added to all files are:

    RA_V1
    DEC_V1
    PA_V3
    RA_REF
    DEC_REF
    ROLL_REF
    S_REGION

    In addition the following keywords are computed and added to IMAGING_MODES only:

    CRVAL1
    CRVAL2
    PC1_1
    PC1_2
    PC2_1
    PC2_2

    It does not currently place the new keywords in any particular location
    in the header other than what is required by the standard.
    """
    logger.info('Updating WCS info for file {}'.format(filename))
    with Level1bModel(filename) as model:
        update_wcs(
            model,
            default_pa_v3=default_pa_v3,
            siaf_path=siaf_path,
            engdb_url=engdb_url,
            tolerance=tolerance,
            allow_default=allow_default,
            reduce_func=reduce_func,
            **transform_kwargs
        )

        try:
            if model.meta.target.type.lower() == 'moving':
                update_mt_kwds(model)
        except AttributeError:
            pass

        model.meta.model_type = None

        if dry_run:
            logger.info('Dry run requested; results are not saved.')
        else:
            model.save(filename)
    logger.info('...update completed')

def update_mt_kwds(model):
    """Add/update the Moving target header keywords

    If the target type is "moving_target" check for the moving target position
    table. If this is available calculate the moving target position keywords
    and insert or update MT_RA & MT_DEC.
    """

    if model.hasattr('moving_target'):
        time_mt = model.moving_target.time
        exp_midpt_mjd = model.meta.exposure.mid_time
        # check to see if the midpoint of the observation is contained within
        # the timerange of the MT table
        if time_mt[0] <= exp_midpt_mjd <= time_mt[-1]:
            ra = model.moving_target.moving_target_RA
            dec = model.moving_target.moving_target_Dec
            f_ra = interp1d(time_mt, ra)
            f_dec = interp1d(time_mt, dec)
            model.meta.wcsinfo.mt_ra = f_ra(exp_midpt_mjd).item(0)
            model.meta.wcsinfo.mt_dec = f_dec(exp_midpt_mjd).item(0)
        else:
            logger.info('Exposure midpoint {} is not in the moving_target '
                        'table range of {} to {}'.format(exp_midpt_mjd, time_mt[0], time_mt[-1]))
            return
    else:
        logger.info("Moving target position table not found in the file")
        return

    logger.info("Moving target RA and Dec updated.")
    return model


def update_wcs(model, default_pa_v3=0., default_roll_ref=0., siaf_path=None, engdb_url=None,
               tolerance=60, allow_default=False,
               reduce_func=None, **transform_kwargs):
    """Update WCS pointing information

    Given a `jwst.datamodels.DataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the roll ref angle.

    siaf_path : str
        The path to the SIAF file, i.e. ``XML_DATA`` env variable.

    engdb_url: str or None
        URL of the engineering telemetry database REST interface.

    tolerance: int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default: bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func: func or None
        Reduction function to use on values.

    transform_kwargs: dict
        Keyword arguments used by matrix calculation routines.
    """

    # If the type of exposure is not FGS, then attempt to get pointing
    # from telemetry.
    try:
        exp_type = model.meta.exposure.type.lower()
    except AttributeError:
        exp_type = None
    aperture_name = model.meta.aperture.name.upper()
    if aperture_name != "UNKNOWN":
        logger.info("Updating WCS for aperture {}".format(aperture_name))
        useafter = model.meta.observation.date
        siaf = _get_wcs_values_from_siaf(aperture_name, useafter, siaf_path)
        populate_model_from_siaf(model, siaf)
    else:
        logger.warning("Aperture name is set to 'UNKNOWN'. "
                       "WCS keywords will not be populated from SIAF.")
        siaf = SIAF()

    if exp_type in FGS_GUIDE_EXP_TYPES:
        update_wcs_from_fgs_guiding(
            model, default_roll_ref=default_roll_ref
        )
    else:
        update_wcs_from_telem(
            model, default_pa_v3=default_pa_v3, siaf=siaf, engdb_url=engdb_url,
            tolerance=tolerance, allow_default=allow_default,
            reduce_func=reduce_func, **transform_kwargs
        )


def update_wcs_from_fgs_guiding(model, default_roll_ref=0.0, default_vparity=1, default_v3yangle=0.0):
    """ Update WCS pointing from header information

    For Fine Guidance guiding observations, nearly everything
    in the `wcsinfo` meta information is already populated,
    except for the PC matrix. This function updates the PC
    matrix based on the rest of the `wcsinfo`.

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update.

    default_pa_v3 : float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.

    default_vparity: int
        The default `VIdlParity` to use and should
        be either "1" or "-1". "1" is the
        default since FGS guiding will be using the
        OSS aperture.
    """

    logger.info('Updating WCS for Fine Guidance.')

    # Get position angle
    try:
        roll_ref = model.meta.wcsinfo.roll_ref if model.meta.wcsinfo.roll_ref is not None else default_roll_ref
    except AttributeError:
        logger.warning(
            'Keyword `ROLL_REF` not found. Using {} as default value'.format(
                default_roll_ref
            )
        )
        roll_ref = default_roll_ref

    roll_ref = np.deg2rad(roll_ref)

    # Get VIdlParity
    try:
        vparity = model.meta.wcsinfo.vparity
    except AttributeError:
        logger.warning(
            'Keyword "VPARITY" not found. Using {} as default value'.format(
                default_vparity
            )
        )
        vparity = default_vparity

    try:
        v3i_yang = model.meta.wcsinfo.v3yangle
    except AttributeError:
        logger.warning(f'Keyword "V3I_YANG" not found. Using {default_v3yangle} as default value.')

        v3i_yang = default_v3yangle

    (
        model.meta.wcsinfo.pc1_1,
        model.meta.wcsinfo.pc1_2,
        model.meta.wcsinfo.pc2_1,
        model.meta.wcsinfo.pc2_2
    ) = calc_rotation_matrix(roll_ref, np.deg2rad(v3i_yang), vparity=vparity)


def update_wcs_from_telem(
    model, default_pa_v3=0., siaf=None, engdb_url=None,
    tolerance=0, allow_default=False,
    reduce_func=None, **transform_kwargs
):
    """Update WCS pointing information

    Given a `jwst.datamodels.DataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update.

    default_pa_v3 : float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.

    siaf_path : str
        The path to the SIAF file, i.e. ``XML_DATA`` env variable.

    engdb_url: str or None
        URL of the engineering telemetry database REST interface.

    tolerance: int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default: bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func: func or None
        Reduction function to use on values.

    transform_kwargs: dict
        Keyword arguments used by matrix calculation routines.
    """

    logger.info('Updating wcs from telemetry.')

    # Get the SIAF and observation parameters
    obsstart = model.meta.exposure.start_time
    obsend = model.meta.exposure.end_time
    if None in siaf:
        # Check if any of "v2_ref", "v3_ref", "v3yangle", "vparity" is None
        # and raise an error. The other fields have default values.
        raise ValueError('Insufficient SIAF information found in header.')

    # Setup default WCS info if actual pointing and calculations fail.
    wcsinfo = WCSRef(
        model.meta.target.ra,
        model.meta.target.dec,
        default_pa_v3
    )
    vinfo = wcsinfo

    # Get the pointing information
    try:
        pointing = get_pointing(obsstart, obsend, engdb_url=engdb_url,
                                tolerance=tolerance, reduce_func=reduce_func)
    except ValueError as exception:
        if not allow_default:
            raise
        else:
            logger.warning(
                'Cannot retrieve valid telescope pointing.'
                ' Default pointing parameters will be used.'
                '\nException is {}'.format(exception)
            )
            logger.info("Setting ENGQLPTG keyword to PLANNED")
            model.meta.visit.pointing_engdb_quality = "PLANNED"
    else:
        # compute relevant WCS information
        logger.info('Successful read of engineering quaternions:')
        logger.info('\tPointing = {}'.format(pointing))
        model.meta.visit.pointing_engdb_quality = "CALCULATED"
        try:
            wcsinfo, vinfo = calc_wcs(pointing, siaf, **transform_kwargs)
            logger.info("Setting ENGQLPTG keyword to CALCULATED")
        except Exception as e:
            logger.warning(
                'WCS calculation has failed and will be skipped.'
                'Default pointing parameters will be used.'
                '\nException is {}'.format(e)
            )
            if not allow_default:
                raise
            else:
                logger.info("Setting ENGQLPTG keyword to PLANNED")
                model.meta.visit.pointing_engdb_quality = "PLANNED"
    logger.info('Aperture WCS info: {}'.format(wcsinfo))
    logger.info('V1 WCS info: {}'.format(vinfo))

    # Update V1 pointing
    model.meta.pointing.ra_v1 = vinfo.ra
    model.meta.pointing.dec_v1 = vinfo.dec
    model.meta.pointing.pa_v3 = vinfo.pa

    # Update Aperture pointing
    model.meta.aperture.position_angle = wcsinfo.pa
    model.meta.wcsinfo.ra_ref = wcsinfo.ra
    model.meta.wcsinfo.dec_ref = wcsinfo.dec
    model.meta.wcsinfo.roll_ref = compute_local_roll(
        vinfo.pa, wcsinfo.ra, wcsinfo.dec, siaf.v2_ref, siaf.v3_ref
    )
    if model.meta.exposure.type.lower() in TYPES_TO_UPDATE:
        model.meta.wcsinfo.crval1 = wcsinfo.ra
        model.meta.wcsinfo.crval2 = wcsinfo.dec
        (
            model.meta.wcsinfo.pc1_1,
            model.meta.wcsinfo.pc1_2,
            model.meta.wcsinfo.pc2_1,
            model.meta.wcsinfo.pc2_2
        ) = calc_rotation_matrix(
            np.deg2rad(model.meta.wcsinfo.roll_ref),
            np.deg2rad(model.meta.wcsinfo.v3yangle),
            vparity=siaf.vparity
        )

    # Calculate S_REGION with the footprint
    # information
    try:
        update_s_region(model, siaf)
    except Exception as e:
        logger.warning(
            'Calculation of S_REGION failed and will be skipped.'
            '\nException is {}'.format(e)
        )


def update_s_region(model, siaf):
    """Update ``S_REGION`` sky footprint information.

    The ``S_REGION`` keyword is intended to store the spatial footprint of
    an observation using the VO standard STCS representation.

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update in-place.
    siaf : namedtuple
        The ``SIAF`` tuple withg values populated from the PRD database.
    """
    vertices = siaf.vertices_idl
    xvert = vertices[:4]
    yvert = vertices[4:]
    logger.info(
        "Vertices for aperture {0}: {1}".format(model.meta.aperture.name, vertices)
    )
    # Execute IdealToV2V3, followed by V23ToSky
    from ..transforms.models import IdealToV2V3
    vparity = model.meta.wcsinfo.vparity
    v3yangle = model.meta.wcsinfo.v3yangle

    # V2_ref and v3_ref should be in arcsec
    idltov23 = IdealToV2V3(
        v3yangle,
        model.meta.wcsinfo.v2_ref, model.meta.wcsinfo.v3_ref,
        vparity
    )
    v2, v3 = idltov23(xvert, yvert)  # in arcsec

    # hardcode wrapping angles for V2 and RA here. Could be made more
    # flexible if needed.
    v23tosky_tr = v23tosky(model, wrap_v2_at=180, wrap_lon_at=360)
    ra_vert, dec_vert = v23tosky_tr(v2, v3)
    # Do not do any sorting, use the vertices in the SIAF order.
    footprint = np.array([ra_vert, dec_vert]).T
    update_s_region_keyword(model, footprint)


def calc_wcs_over_time(obsstart, obsend, engdb_url=None, tolerance=60, reduce_func=None,
                       siaf=None, **transform_kwargs):
    """Calculate V1 and WCS over a time period

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    reduce_func : func or None
        Reduction function to use on values.
        If None, the average pointing is returned.

    siaf : SIAF or None
        The SIAF transformation. If `None`, a unit transformation is used.

    transform_kwargs : dict
        Keyword arguments used by matrix calculation routines

    Returns
    -------
    obstimes, wcsinfos, vinfos : [], [WCSRef[,...]], [WCSRef[,...]]
        A 3-tuple is returned with the WCS pointings for
        the aperture and the V1 axis
    """
    # Setup structures
    obstimes = list()
    wcsinfos = list()
    vinfos = list()

    # Calculate WCS
    try:
        pointings = get_pointing(obsstart, obsend, engdb_url=engdb_url,
                             tolerance=tolerance, reduce_func=reduce_func)
    except ValueError:
        logger.warning("Cannot get valid engineering mnemonics from engineering database")
        raise
    if not isinstance(pointings, list):
        pointings = [pointings]
    for pointing in pointings:
        wcsinfo, vinfo = calc_wcs(pointing, siaf=siaf, **transform_kwargs)
        obstimes.append(pointing.obstime)
        wcsinfos.append(wcsinfo)
        vinfos.append(vinfo)

    return obstimes, wcsinfos, vinfos


def calc_wcs(pointing, siaf=None, **transform_kwargs):
    """Transform from the given SIAF information and Pointing
    the aperture and V1 wcs

    Parameters
    ----------
    pointing : Pointing
        The telescope pointing. See ref:`Notes` for further details

    siaf : SIAF or None
        The SIAF transformation. See ref:`Notes` for further details.
        If `None`, unit transformation is used.

    transform_kwargs : dict
        Keyword arguments used by matrix calculation routines

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
    if siaf is None:
        siaf = SIAF()

    # Calculate transforms
    tforms = calc_transforms(pointing, siaf, **transform_kwargs)

    # Calculate the V1 WCS information
    vinfo = calc_v1_wcs(tforms.m_eci2v)

    # Calculate the Aperture WCS
    wcsinfo = calc_aperture_wcs(tforms.m_eci2siaf)

    # That's all folks
    return (wcsinfo, vinfo)


def calc_transforms(pointing, siaf, fsmcorr_version='latest', fsmcorr_units='arcsec', j2fgs_transpose=True):
    """Calculate transforms from pointing to SIAF

    Given the spacecraft pointing parameters and the
    aperture-specific SIAF, calculate all the transforms
    necessary to produce WCS information.

    Parameters
    ----------
    pointing : Pointing
        Observatory pointing information

    siaf : SIAF
        Aperture information

    fsmcorr_version : str
        The version of the FSM correction calculation to use.
        See :ref:`calc_sifov_fsm_delta_matrix`

    fsmcorr_units : str
        Units of the FSM correction values. Default is 'arcsec'.
        See :ref:`calc_sifov_fsm_delta_matrix`

    j2fgs_transpose : bool
        Transpose the `j2fgs1` matrix.

    Returns
    -------
    transforms: Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    The matrix transform pipeline to convert from ECI J2000 observatory
    qauternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =           # The complete transformation
            M_v1_to_siaf      *   # V1 to SIAF
            M_sifov_to_v1     *   # Science Instruments Aperture to V1
            M_sifov_fsm_delta *   # Fine Steering Mirror correction
            M_fgs1_to_sifov   *   # FGS1 to Science Instruments Aperture
            M_j_to_fgs1       *   # J-Frame to FGS1
            M_eci_to_j        *   # ECI to J-Frame
    """

    # Determine the ECI to J-frame matrix
    m_eci2j = calc_eci2j_matrix(pointing.q)

    # Calculate the J-frame to FGS! ICS matrix
    m_j2fgs1 = calc_j2fgs1_matrix(pointing.j2fgs_matrix, transpose=j2fgs_transpose)

    # Calculate the FSM corrections to the SI_FOV frame
    m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        pointing.fsmcorr, fsmcorr_version=fsmcorr_version, fsmcorr_units=fsmcorr_units
    )

    # Calculate the FGS1 ICS to SI-FOV matrix
    m_fgs12sifov = calc_fgs1_to_sifov_mastrix()

    # Calculate SI FOV to V1 matrix
    m_sifov2v = calc_sifov2v_matrix()

    # Calculate ECI to SI FOV
    m_eci2sifov = np.dot(
        m_sifov_fsm_delta,
        np.dot(
            m_fgs12sifov,
            np.dot(
                m_j2fgs1,
                m_eci2j
            )
        )
    )

    # Calculate the complete transform to the V1 reference
    m_eci2v = np.dot(
        m_sifov2v,
        m_eci2sifov
    )

    # Calculate the SIAF transform matrix
    m_v2siaf = calc_v2siaf_matrix(siaf)

    # Calculate the full ECI to SIAF transform matrix
    m_eci2siaf = np.dot(
        m_v2siaf,
        m_eci2v
    )

    tforms = Transforms(m_eci2j=m_eci2j, m_j2fgs1=m_j2fgs1, m_sifov_fsm_delta=m_sifov_fsm_delta,
                        m_fgs12sifov=m_fgs12sifov, m_eci2sifov=m_eci2sifov, m_sifov2v=m_sifov2v,
                        m_eci2v=m_eci2v, m_v2siaf=m_v2siaf, m_eci2siaf=m_eci2siaf
                        )
    return tforms


def calc_v1_wcs(m_eci2v):
    """Calculate the V1 wcs information

    Parameters
    ----------
    m_eci2v: np.array((3, 3))
        The ECI to V transformation matrix

    Returns
    -------
    vinfo: WCSRef
        The V1 wcs pointing
    """
    # V1 RA/Dec is the first row of the transform
    v1_ra, v1_dec = vector_to_ra_dec(m_eci2v[0])
    vinfo = WCSRef(v1_ra, v1_dec, None)

    # V3 is the third row of the transformation
    v3_ra, v3_dec = vector_to_ra_dec(m_eci2v[2])
    v3info = WCSRef(v3_ra, v3_dec, None)

    # Calculate the V3 position angle
    v1_pa = calc_position_angle(vinfo, v3info)

    # Convert to degrees
    vinfo = WCSRef(
        ra=vinfo.ra * R2D,
        dec=vinfo.dec * R2D,
        pa=v1_pa * R2D
    )

    return vinfo


def calc_aperture_wcs(m_eci2siaf):
    """Calculate the aperture WCS

    Parameters
    ----------
    m_eci2siaf: np.array((3, 3))
        The ECI to SIAF transformation matrix

    Returns
    -------
    wcsinfo: WCSRef
        The aperturn wcs information
    """

    # Calculate point on sky.
    # Note, the SIAF referenct point is hardcoded to
    # (0, 0). The calculation is left here in case
    # this is not desired.

    siaf_x = 0. * A2R
    siaf_y = 0. * A2R
    refpos = np.array(
        [siaf_x,
         siaf_y,
         np.sqrt(1.-siaf_x * siaf_x - siaf_y * siaf_y)]
    )
    msky = np.dot(m_eci2siaf.transpose(), refpos)
    wcs_ra, wcs_dec = vector_to_ra_dec(msky)

    # Calculate the position angle
    vysiaf = np.array([0., 1., 0.])
    myeci = np.dot(m_eci2siaf.transpose(), vysiaf)

    # The Y axis of the aperture is given by
    vy_ra, vy_dec = vector_to_ra_dec(myeci)

    # The VyPA @ xref,yref is given by
    y = cos(vy_dec) * sin(vy_ra-wcs_ra)
    x = sin(vy_dec) * cos(wcs_dec) - \
        cos(vy_dec) * sin(wcs_dec) * cos((vy_ra - wcs_ra))
    wcs_pa = np.arctan2(y, x)

    # Convert all WCS to degrees
    wcsinfo = WCSRef(
        ra=wcs_ra * R2D,
        dec=wcs_dec * R2D,
        pa=wcs_pa * R2D
    )

    return wcsinfo


def calc_eci2j_matrix(q):
    """Calculate ECI to J-frame matrix from quaternions

    Parameters
    ----------
    q: np.array(q1, q2, q3, q4)
        Array of quaternions from the engineering database

    Returns
    -------
    transform: np.array((3, 3))
        The transform matrix representing the transformation
        from observatory orientation to J-Frame
    """
    q1, q2, q3, q4 = q
    transform = np.array(
        [[1. - 2.*q2*q2 - 2.*q3*q3,
          2.*(q1*q2 + q3*q4),
          2.*(q3*q1 - q2*q4)],
         [2.*(q1*q2 - q3*q4),
          1. - 2.*q3*q3 - 2.*q1*q1,
          2.*(q2*q3 + q1*q4)],
         [2.*(q3*q1 + q2*q4),
          2.*(q2*q3 - q1*q4),
          1. - 2.*q1*q1 - 2.*q2*q2]]
    )

    return transform


def calc_j2fgs1_matrix(j2fgs_matrix, transpose=False):
    """Calculate the J-frame to FGS1 transformation

    Parameters
    ----------
    j2fgs_matrix: n.array((9,))
        Matrix parameters from the engineering database.
        If all zeros, a predefined matrix is used.

    transpose: bool
        Transpose the resulting matrix.

    Returns
    -------
    transform: np.array((3, 3))
        The transformation matrix
    """
    if np.isclose(j2fgs_matrix, 0.).all():
        logger.warning(
            'J-Frame to FGS1 engineering parameters are all zero.'
            '\nUsing default matrix'
        )
        m_partial = np.asarray(
            [
                [0., 1., 0.],
                [0., 0., 1.],
                [1., 0., 0.]
            ]
        )
        transform = np.dot(
            m_partial,
            J2FGS_MATRIX_DEFAULT
        )

    else:
        logger.info(
            'Using J-Frame to FGS1 engineering parameters'
            ' for the J-Frame to FGS1 transformation.'
        )
        transform = np.array(j2fgs_matrix).reshape((3, 3))
        if transpose:
            transform = transform.transpose()

    return transform


def calc_sifov_fsm_delta_matrix(fsmcorr, fsmcorr_version='latest', fsmcorr_units='arcsec'):
    """Calculate Fine Steering Mirror correction matrix

    Parameters
    ----------
    fsmcorr : np.array((2,))
        The FSM correction parameters:
            0: SA_ZADUCMDX
            1: SA_ZADUCMDY

    fsmcorr_version : str
        The version of the FSM correction calculation to use.
        Versions available:
            latest: The state-of-art. Currently `v2`
            v2: Update 201708 to use actual spherical calculations
            v1: Original linear approximation

    fsmcorr_units : str
        The units of the FSM correction values. Default is `arcsec`.

    Returns
    -------
    transform: np.array((3, 3))
        The transformation matrix
    """
    version = fsmcorr_version.lower()
    units = fsmcorr_units.lower()
    logger.debug('Using version {}'.format(version))
    logger.debug('Using units {}'.format(units))

    x = fsmcorr[0]  # SA_ZADUCMDX
    y = fsmcorr[1]  # SA_ZADUCMDY

    # If FSMCORR values are in arcsec, convert to radians
    if units == 'arcsec':
        x *= D2R / 3600.
        y *= D2R / 3600.

    # `V1`: Linear approximation calculation
    if version == 'v1':
        transform = np.array(
            [
                [1.,       x/22.01, y/21.68],
                [-x/22.01, 1.,      0.],
                [-y/21.68, 0.,      1.]
            ]
        )

    # Default or `V2`: Direct spherical calculation
    # Note: With the "0.0" in the lower middle Y transform
    else:
        if version not in ('latest', 'v2'):
            logger.warning(
                'Unknown version "{}" specified.'
                ' Using the latest (spherical) calculation.'
            )
        m_x_partial = np.array(
            [
                [1., 0.,      0.],
                [0., cos(x),  sin(x)],
                [0., -sin(x), cos(x)]
            ]
        )
        m_y_partial = np.array(
            [
                [cos(y), 0., -sin(y)],
                [0.,     1., 0.],
                [sin(y), 0., cos(y)]
            ]
        )
        transform = np.dot(m_x_partial, m_y_partial)

    return transform


def calc_fgs1_to_sifov_mastrix():
    """
    Calculate the FGS! to SI-FOV matrix

    Currently, this is a defined matrix
    """
    m_partial = np.array(
        [[0, 0, 1],
         [1, 0, 0],
         [0, 1, 0]]
    )

    transform = np.dot(m_partial, FGS12SIFOV_DEFAULT)
    return transform


def calc_sifov2v_matrix():
    """Calculate the SI-FOV to V-Frame matrix

    This is currently defined as the inverse Euler rotation
    about an angle of 7.8 arcmin. Here returns the pre-calculate
    matrix.
    """
    return SIFOV2V_DEFAULT


def calc_v2siaf_matrix(siaf):
    """Calculate the SIAF transformation matrix

    Parameters
    ----------
    siaf: SIAF
        The SIAF parameters

    Returns
    -------
    transform: np.array((3, 3))
        The V1 to SIAF transformation matrix
    """
    v2, v3, v3idlyang, vparity = (siaf.v2_ref, siaf.v3_ref,
                                  siaf.v3yangle, siaf.vparity)
    v2 *= A2R
    v3 *= A2R
    v3idlyang *= D2R

    mat = np.array(
        [[cos(v3)*cos(v2),
          cos(v3)*sin(v2),
          sin(v3)],
         [-cos(v3idlyang)*sin(v2)+sin(v3idlyang)*sin(v3)*cos(v2),
          cos(v3idlyang)*cos(v2)+sin(v3idlyang)*sin(v3)*sin(v2),
          -sin(v3idlyang)*cos(v3)],
         [-sin(v3idlyang)*sin(v2)-cos(v3idlyang)*sin(v3)*cos(v2),
          sin(v3idlyang)*cos(v2)-cos(v3idlyang)*sin(v3)*sin(v2),
          cos(v3idlyang)*cos(v3)]])
    pmat = np.array([[0., vparity, 0.],
                     [0., 0., 1.],
                     [1., 0., 0.]])

    transform = np.dot(pmat, mat)
    return transform


def calc_position_angle(v1, v3):
    """Calculate V3 position angle @V1

    Parameters
    ----------
    v1: WCSRef
        The V1 wcs parameters

    v3: WCSRef
        The V3 wcs parameters

    Returns
    -------
    v3_pa: float
      The V3 position angle, in radians
    """
    y = cos(v3.dec) * sin(v3.ra-v1.ra)
    x = sin(v3.dec) * cos(v1.dec) - \
        cos(v3.dec) * sin(v1.dec) * cos((v3.ra - v1.ra))
    v3_pa = np.arctan2(y, x)

    return v3_pa


def get_pointing(obsstart, obsend, engdb_url=None,
                 tolerance=60, reduce_func=None):
    """
    Get telescope pointing engineering data.

    Parameters
    ----------
    obsstart, obsend: float
        MJD observation start/end times

    engdb_url: str or None
        URL of the engineering telemetry database REST interface.

    tolerance: int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    reduce_func: func or None
        Reduction function to use on values.
        If None, the average pointing is returned.

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
    if reduce_func is None:
        reduce_func = pointing_from_average

    logger.info(
        'Determining pointing between observations times (mjd):'
        'obsstart = {obsstart} obsend = {obsend}'
        '\nTelemetry search tolerance = {tolerance}'
        '\nReduction function = {reduce_func}'
        ''.format(
            obsstart=obsstart, obsend=obsend, tolerance=tolerance, reduce_func=reduce_func
        )
    )

    mnemonics = get_mnemonics(obsstart, obsend, tolerance, engdb_url=engdb_url)
    reduced = reduce_func(mnemonics)

    logger.debug(f'Memonics found:\n{mnemonics}')
    logger.info(f'Reduced set of pointings:\n{reduced}')

    return reduced


def vector_to_ra_dec(v):
    """Returns tuple of spherical angles from unit direction Vector

    Parameters
    ----------
    v: [v0, v1, v2]

    Returns
    -------
    ra, dec: float, float
        The spherical angles, in radians
    """
    ra = np.arctan2(v[1], v[0])
    dec = np.arcsin(v[2])
    if ra < 0.:
        ra += 2. * np.pi
    return(ra, dec)


def compute_local_roll(pa_v3, ra_ref, dec_ref, v2_ref, v3_ref):
    """
    Computes the position angle of V3 (measured N to E)
    at the center af an aperture.

    Parameters
    ----------
    pa_v3 : float
        Position angle of V3 at (V2, V3) = (0, 0) [in deg]
    v2_ref, v3_ref : float
        Reference point in the V2, V3 frame [in arcsec]
    ra_ref, dec_ref : float
        RA and DEC corresponding to V2_REF and V3_REF, [in deg]

    Returns
    -------
    new_roll : float
        The value of ROLL_REF (in deg)

    """
    v2 = np.deg2rad(v2_ref / 3600)
    v3 = np.deg2rad(v3_ref / 3600)
    ra_ref = np.deg2rad(ra_ref)
    dec_ref = np.deg2rad(dec_ref)
    pa_v3 = np.deg2rad(pa_v3)

    M = np.array(
        [[cos(ra_ref) * cos(dec_ref),
          -sin(ra_ref) * cos(pa_v3) + cos(ra_ref) * sin(dec_ref) * sin(pa_v3),
          -sin(ra_ref) * sin(pa_v3) - cos(ra_ref) * sin(dec_ref) * cos(pa_v3)],
         [sin(ra_ref) * cos(dec_ref),
          cos(ra_ref) * cos(pa_v3) + sin(ra_ref) * sin(dec_ref) * sin(pa_v3),
          cos(ra_ref) * sin(pa_v3) - sin(ra_ref) * sin(dec_ref) * cos(pa_v3)],
         [sin(dec_ref),
          -cos(dec_ref) * sin(pa_v3),
          cos(dec_ref) * cos(pa_v3)]]
    )

    return _roll_angle_from_matrix(M, v2, v3)


def _roll_angle_from_matrix(matrix, v2, v3):
    X = -(matrix[2, 0] * np.cos(v2) + matrix[2, 1] * np.sin(v2)) * np.sin(v3) + matrix[2, 2] * np.cos(v3)
    Y = (matrix[0, 0] * matrix[1, 2] - matrix[1, 0] * matrix[0, 2]) * np.cos(v2) + \
        (matrix[0, 1] * matrix[1, 2] - matrix[1, 1] * matrix[0, 2]) * np.sin(v2)
    new_roll = np.rad2deg(np.arctan2(Y, X))
    if new_roll < 0:
        new_roll += 360
    return new_roll


def _get_wcs_values_from_siaf(aperture_name, useafter, prd_db_filepath=None):
    """
    Query the SIAF database file and get WCS values.

    Given an ``APERTURE_NAME`` and a ``USEAFTER`` date query the SIAF database
    and extract the following keywords:
    ``V2Ref``, ``V3Ref``, ``V3IdlYAngle``, ``VIdlParity``,
    ``XSciRef``, ``YSciRef``, ``XSciScale``, ``YSciScale``,
    ``XIdlVert1``, ``XIdlVert2``, ``XIdlVert3``, ``XIdlVert4``,
    ``YIdlVert1``, ``YIdlVert2``, ``YIdlVert3``, ``YIdlVert4``

    Parameters
    ----------
    aperture_name : str
        The name of the aperture in the data file (``model.meta.aperture.name``).
    useafter : str
        The date of observation (``model.meta.date``)
    prd_db_filepath : str
        The path to the SIAF (PRD database) file.
        If None, attempt to get the path from the ``XML_DATA`` environment variable.

    Returns
    -------
    siaf : namedtuple
        The SIAF namedtuple with values from the PRD database.
    """
    if prd_db_filepath is None:
        try:
            prd_db_filepath = os.path.join(os.environ['XML_DATA'], "prd.db")
        except KeyError:
            message = "Unknown path to PRD DB file or missing env variable ``XML_DATA``."
            logger.info(message)
            raise KeyError(message)
    if not os.path.exists(prd_db_filepath):
        message = "Invalid path to PRD DB file: {0}".format(prd_db_filepath)
        logger.info(message)
        raise OSError(message)
    prd_db_filepath = "file:{0}?mode=ro".format(prd_db_filepath)
    logger.info("Using SIAF database from {}".format(prd_db_filepath))
    logger.info("Quering SIAF for aperture "
                "{0} with USEAFTER {1}".format(aperture_name, useafter))
    aperture = (aperture_name, useafter)

    RESULT = {}
    PRD_DB = False
    try:
        PRD_DB = sqlite3.connect(prd_db_filepath, uri=True)

        cursor = PRD_DB.cursor()
        cursor.execute("SELECT Apername, V2Ref, V3Ref, V3IdlYAngle, VIdlParity, "
                       "XSciRef, YSciRef, XSciScale, YSciScale, "
                       "XIdlVert1, XIdlVert2, XIdlVert3, XIdlVert4, "
                       "YIdlVert1, YIdlVert2, YIdlVert3, YIdlVert4 "
                       "FROM Aperture WHERE Apername = ? and UseAfterDate <= ? ORDER BY UseAfterDate LIMIT 1",
                       aperture)
        for row in cursor:
            RESULT[row[0]] = tuple(row[1:17])
        PRD_DB.commit()
    except (sqlite3.Error, sqlite3.OperationalError) as err:
        print("Error: " + err.args[0])
        raise
    finally:
        if PRD_DB:
            PRD_DB.close()
    logger.info("loaded {0} table rows from {1}".format(len(RESULT), prd_db_filepath))
    default_siaf = SIAF()
    if RESULT:
        # This populates the SIAF tuple with the values from the database.
        # The last 8 values returned from the database are the vertices.
        # They are wrapped in a list and assigned to SIAF.vertices_idl.
        values = list(RESULT.values())[0]
        vert = values[-8:]
        values = list(values[: - 8])
        values.append(vert)
        # If any of "crpix1", "crpix2", "cdelt1", "cdelt2", "vertices_idl" is None
        # reset ot to the default value.
        for i in range(4, 8):
            if values[i] is None:
                values[i] = default_siaf[i]
        siaf = SIAF(*values)
        return siaf
    else:
        return default_siaf


def get_mnemonics(obsstart, obsend, tolerance, engdb_url=None):
    """Retrieve pointing mnemonics from the engineering database

    Parameters
    ----------
    obsstart, obsend: float
        MJD observation start/end times

    tolerance: int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    engdb_url: str or None
        URL of the engineering telemetry database REST interface.

    Returns
    -------
    mnemonics: {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Raises
    ------
    ValueError
        Cannot retrieve engineering information

    """
    try:
        engdb = ENGDB_Service(base_url=engdb_url)
    except Exception as exception:
        raise ValueError(
            'Cannot open engineering DB connection'
            '\nException: {}'.format(
                exception
            )
        )
    logger.info(
        'Querying engineering DB: {}'.format(engdb.base_url)
    )

    mnemonics = {
        'SA_ZATTEST1':  None,
        'SA_ZATTEST2':  None,
        'SA_ZATTEST3':  None,
        'SA_ZATTEST4':  None,
        'SA_ZRFGS2J11': None,
        'SA_ZRFGS2J12': None,
        'SA_ZRFGS2J13': None,
        'SA_ZRFGS2J21': None,
        'SA_ZRFGS2J22': None,
        'SA_ZRFGS2J23': None,
        'SA_ZRFGS2J31': None,
        'SA_ZRFGS2J32': None,
        'SA_ZRFGS2J33': None,
        'SA_ZADUCMDX':  None,
        'SA_ZADUCMDY':  None,
    }

    # Retrieve the mnemonics from the engineering database.
    # Check for whether the bracket values are used and
    # within tolerance.
    for mnemonic in mnemonics:
        try:
            mnemonics[mnemonic] = engdb.get_values(
                mnemonic, obsstart, obsend,
                time_format='mjd', include_obstime=True,
                include_bracket_values=True
            )
        except Exception as exception:
            raise ValueError(
                'Cannot retrive {} from engineering.'
                '\nFailure was {}'.format(
                    mnemonic,
                    exception
                )
            )

        # If more than two points exist, throw off the bracket values.
        # Else, ensure the bracket values are within the allowed time.
        if len(mnemonics[mnemonic]) > 2:
            mnemonics[mnemonic] = mnemonics[mnemonic][1:-1]
        else:
            logger.warning(
                'Mnemonic {} has no telemetry within the observation time.'
                '\nAttempting to use bracket values within {} seconds'.format(
                    mnemonic, tolerance
                )
            )
            tolerance_mjd = tolerance * SECONDS2MJD
            allowed_start = obsstart - tolerance_mjd
            allowed_end = obsend + tolerance_mjd
            allowed = [
                value
                for value in mnemonics[mnemonic]
                if allowed_start <= value.obstime.mjd <= allowed_end
            ]
            if not len(allowed):
                raise ValueError(
                    'No telemetry exists for mnemonic {} within {} and {}'.format(
                        mnemonic,
                        Time(allowed_start, format='mjd').isot,
                        Time(allowed_end, format='mjd').isot
                    )
                )
            mnemonics[mnemonic] = allowed

    # All mnemonics must have some values.
    if not all([len(mnemonic) for mnemonic in mnemonics.values()]):
        raise ValueError('Incomplete set of pointing mnemonics')

    return mnemonics


def all_pointings(mnemonics):
    """V1 of making pointings

    Parameters
    ==========
    mnemonics: {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointings: [Pointing[,...]]
        List of pointings.
    """
    pointings = []
    filled = fill_mnemonics_chronologically(mnemonics)
    for obstime, mnemonics_at_time in filled.items():

        # Fill out the matricies
        q = np.array([
            mnemonics_at_time['SA_ZATTEST1'].value,
            mnemonics_at_time['SA_ZATTEST2'].value,
            mnemonics_at_time['SA_ZATTEST3'].value,
            mnemonics_at_time['SA_ZATTEST4'].value,
        ])

        j2fgs_matrix = np.array([
            mnemonics_at_time['SA_ZRFGS2J11'].value,
            mnemonics_at_time['SA_ZRFGS2J12'].value,
            mnemonics_at_time['SA_ZRFGS2J13'].value,
            mnemonics_at_time['SA_ZRFGS2J21'].value,
            mnemonics_at_time['SA_ZRFGS2J22'].value,
            mnemonics_at_time['SA_ZRFGS2J23'].value,
            mnemonics_at_time['SA_ZRFGS2J31'].value,
            mnemonics_at_time['SA_ZRFGS2J32'].value,
            mnemonics_at_time['SA_ZRFGS2J33'].value,
        ])

        fsmcorr = np.array([
            mnemonics_at_time['SA_ZADUCMDX'].value,
            mnemonics_at_time['SA_ZADUCMDY'].value,

        ])
        pointing = Pointing(q=q, obstime=obstime, j2fgs_matrix=j2fgs_matrix,
                            fsmcorr=fsmcorr)
        pointings.append(pointing)

    if not len(pointings):
        raise ValueError('No non-zero quanternion found.')

    return pointings


def populate_model_from_siaf(model, siaf):
    """
    Populate the WCS keywords of a Level1bModel from the SIAF.

    Parameters
    ----------
    model : `~jwst.datamodels.Level1bModel`
        Input data as Level1bModel.
    siaf : namedtuple
        The WCS keywords read in from the SIAF.
    """
    # Update values from the SIAF for all exposures.
    model.meta.wcsinfo.v2_ref = siaf.v2_ref
    model.meta.wcsinfo.v3_ref = siaf.v3_ref
    model.meta.wcsinfo.v3yangle = siaf.v3yangle
    model.meta.wcsinfo.vparity = siaf.vparity

    # For imaging modes, also update the basic FITS WCS keywords
    if model.meta.exposure.type.lower() in TYPES_TO_UPDATE:
        logger.info('Setting basic FITS WCS keywords for imaging')
        model.meta.wcsinfo.ctype1 = 'RA---TAN'
        model.meta.wcsinfo.ctype2 = 'DEC--TAN'
        model.meta.wcsinfo.wcsaxes = 2
        model.meta.wcsinfo.cunit1 = "deg"
        model.meta.wcsinfo.cunit2 = "deg"
        model.meta.wcsinfo.crpix1 = siaf.crpix1
        model.meta.wcsinfo.crpix2 = siaf.crpix2
        model.meta.wcsinfo.cdelt1 = siaf.cdelt1 / 3600  # in deg
        model.meta.wcsinfo.cdelt2 = siaf.cdelt2 / 3600  # in deg
        model.meta.coordinates.reference_frame = "ICRS"

    # For TSO exposures, also populate XREF_SCI/YREF_SCI keywords,
    # which are used by the Cal pipeline to determine the
    # location of the source.
    # Note that we use a combination of the is_tso function and
    # a check on EXP_TYPE, because there are rare corner cases
    # where EXP_TIME=NRC_TSGRISM, TSOVISIT=False, NINTS=1, which
    # normally return False, but we want to treat it as TSO anyway.
    if is_tso(model) or model.meta.exposure.type.lower() in ['nrc_tsimage', 'nrc_tsgrism']:
        logger.info('TSO exposure:')
        logger.info(' setting xref_sci to {}'.format(siaf.crpix1))
        logger.info(' setting yref_sci to {}'.format(siaf.crpix2))
        model.meta.wcsinfo.siaf_xref_sci = siaf.crpix1
        model.meta.wcsinfo.siaf_yref_sci = siaf.crpix2


def first_pointing(mnemonics):
    """Return first pointing

    Parameters
    ==========
    mnemonics: {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing: Pointing
        First pointing.

    """
    pointings = all_pointings(mnemonics)
    return pointings[0]


def pointing_from_average(mnemonics):
    """Determine single pointing from average of available pointings

    Parameters
    ==========
    mnemonics: {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing: Pointing
        Pointing from average.

    """
    # Get average observation time. This is keyed off the q0 quaternion term, SA_ZATTEST1
    times = [
        eng_param.obstime.unix
        for eng_param in mnemonics['SA_ZATTEST1']
    ]
    goodtimes = []
    for this_time in times:
        if this_time != 0.0:
            goodtimes.append(this_time)
    if len(goodtimes) > 0:
        obstime = Time(np.average(goodtimes), format='unix')
    else:
        raise ValueError("No valid times in range")
    # Get averages for all the mnemonics.
    mnemonic_averages = {}
    zero_mnemonics = []
    for mnemonic in mnemonics:
        values = [
            eng_param.value
            for eng_param in mnemonics[mnemonic]
        ]
        # Weed out mnemonic entries that are zero
        # 'SA_ZADUCMDX' and 'SA_ZADUCMDY' can be zero
        if mnemonic not in ['SA_ZADUCMDX', 'SA_ZADUCMDY']:
            good_mnemonic = []
            for this_value in values:
                if this_value != 0.0:
                    good_mnemonic.append(this_value)
            if len(good_mnemonic) > 0:
                mnemonic_averages[mnemonic] = np.average(good_mnemonic)
            else:
                zero_mnemonics.append(mnemonic)
        else:
            mnemonic_averages[mnemonic] = np.average(values)
    # Raise exception if there are mnemonics with only zeros in the time range
    if len(zero_mnemonics):
        logger.warning("The following engineering mnemonics only contained zeros in the requested time interval:")
        badmnemonicsstring = ' '.join(zero_mnemonics)
        logger.info(badmnemonicsstring)
        raise ValueError("Bad telemetry values")
    # Fill out the pointing matrices.
    q = np.array([
        mnemonic_averages['SA_ZATTEST1'],
        mnemonic_averages['SA_ZATTEST2'],
        mnemonic_averages['SA_ZATTEST3'],
        mnemonic_averages['SA_ZATTEST4']
    ])

    j2fgs_matrix = np.array([
        mnemonic_averages['SA_ZRFGS2J11'],
        mnemonic_averages['SA_ZRFGS2J12'],
        mnemonic_averages['SA_ZRFGS2J13'],
        mnemonic_averages['SA_ZRFGS2J21'],
        mnemonic_averages['SA_ZRFGS2J22'],
        mnemonic_averages['SA_ZRFGS2J23'],
        mnemonic_averages['SA_ZRFGS2J31'],
        mnemonic_averages['SA_ZRFGS2J32'],
        mnemonic_averages['SA_ZRFGS2J33']
    ])

    fsmcorr = np.array([
        mnemonic_averages['SA_ZADUCMDX'],
        mnemonic_averages['SA_ZADUCMDY']

    ])
    pointing = Pointing(obstime=obstime, q=q, j2fgs_matrix=j2fgs_matrix,
                        fsmcorr=fsmcorr)
    # That's all folks
    return pointing


def fill_mnemonics_chronologically(mnemonics, filled_only=True):
    """Return time-ordered mnemonic list with progressive values

    The different set of mnemonics used for observatory orientation
    appear at different cadences. This routine creates a time-ordered dictionary
    with all the mnemonics for each time found in the engineering. For mnemonics
    missing for a particular time, the last previous value is used.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]]}

    filled_only : bool
        Only return a matrix where observation times have all the mnemonics defined.

    Returns
    -------
    filled_by_time : {obstime: {mnemonic: value}}
    """
    # Collect all information by observation time and order.
    by_obstime = defaultdict(dict)
    n_mnemonics = len(mnemonics)
    for mnemonic, values in mnemonics.items():
        for value in values:
            by_obstime[value.obstime][mnemonic] = value
    by_obstime = sorted(by_obstime.items())

    # Created the filled matrix
    filled = dict()
    last_obstime = dict()
    for obstime, mnemonics_at_time in by_obstime:
        last_obstime.update(mnemonics_at_time)
        if len(last_obstime) >= n_mnemonics or not filled_only:

            # Engineering data may be present, but all zeros.
            # Filter out this situation also.
            if filled_only:
                values = [
                    value.value
                    for value in last_obstime.values()
                ]
                if not any(values):
                    continue

            filled[obstime] = copy(last_obstime)

    return filled
