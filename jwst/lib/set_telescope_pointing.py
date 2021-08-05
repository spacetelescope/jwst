"""Set Telescope Pointing from quaternions"""
import sys

import asdf
from collections import defaultdict, namedtuple
from copy import copy, deepcopy
import dataclasses
from enum import Enum
import logging
from math import (asin, atan2, cos, sin, sqrt)
import typing

from astropy.time import Time
import numpy as np
from scipy.interpolate import interp1d

from .exposure_types import IMAGING_TYPES, FGS_GUIDE_EXP_TYPES
from .set_velocity_aberration import compute_va_effects_vector
from .siafdb import SIAF, SiafDb
from ..assign_wcs.util import update_s_region_keyword, calc_rotation_matrix
from ..assign_wcs.pointing import v23tosky
from .. import datamodels
from ..lib.engdb_tools import ENGDB_Service
from ..lib.pipe_utils import is_tso

TYPES_TO_UPDATE = set(list(IMAGING_TYPES) + FGS_GUIDE_EXP_TYPES)

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]


# The available methods for transformation
class Methods(Enum):
    COARSE_TR_202107 = ('coarse_tr_202107', 'calc_transforms_coarse_tr_202107')
    GSCMD_J3PAGS = ('gscmd', 'calc_transforms_gscmd_j3pags')
    GSCMD_V3PAGS = ('gscmd_v3pags', 'calc_transforms_gscmd_v3pags')
    OPS_TR_202107 = ('ops_tr_202107', 'calc_transforms_ops_tr_202107')
    ORIGINAL = ('original', 'calc_transforms_original')
    TR_202105 = ('tr_202105', 'calc_transforms_tr202105')  # Was FULL
    TR_202105_VA = ('tr_202105_va', 'calc_transforms_velocity_abberation_tr202105')  # Was FULLVA
    TRACK_TR_202107 = ('track_tr_202107', 'calc_transforms_track_tr_202107')

    # Aliases
    default = OPS_TR_202107  # Algorithm to use by default. Used by Operations.
    COARSE = COARSE_TR_202107  # Default algorithm under PCS_MODE COARSE.
    OPS = OPS_TR_202107  # Default algorithm for use by Operations.
    TRACK = TRACK_TR_202107  # Default algorithm under PCS_MODE TRACK/FINEGUIDE.

    def __new__(cls: object, value: str, func_name: str):
        obj = object.__new__(cls)
        obj._value_ = value
        obj._func_name = func_name
        return obj

    @property
    def func(self):
        return globals()[self._func_name]

    def __str__(self):
        return self.value


# FGS id to aperture name
FGSId2Aper = {1: 'FGS1_FULL_OSS', 2: 'FGS2_FULL_OSS'}

# Definition of th J3 Ideal Y-Angle
J3IDLYANGLE = -1.25  # Degrees

# Conversion from seconds to MJD
SECONDS2MJD = 1 / 24 / 60 / 60

# Default transformation matricies
FGS12SIFOV_DEFAULT = np.array(
    [[0.9999994955442, 0.0000000000000, 0.0010044457459],
     [0.0000011174826, 0.9999993811310, -0.0011125359826],
     [-0.0010044451243, 0.0011125365439, 0.9999988766756]]
)

J2FGS_MATRIX_DEFAULT = np.array([
    [-0.0010044400033, 0.9999994955442, 0.0000033964915],
    [0.0033814583568, 0.0000000000000, 0.9999942828533],
    [0.9999937784005, 0.0010044457459, -0.0033814566510]
])

SIFOV2V_DEFAULT = np.array(
    [[0.99999742598, 0., 0.00226892608],
     [0., 1., 0.],
     [-0.00226892608, 0., 0.99999742598]]
)

MX2Z = np.array(
    [[0, 1, 0],
     [0, 0, 1],
     [1, 0, 0]]
)
MZ2X = np.array(
    [[0, 0, 1],
     [1, 0, 0],
     [0, 1, 0]]
)

# Degree, radian, angle transformations
R2D = 180. / np.pi
D2R = np.pi / 180.
A2R = D2R / 3600.
R2A = 3600. * R2D
PI2 = np.pi * 2.

# Pointing container
Pointing = namedtuple('Pointing', ['q', 'j2fgs_matrix', 'fsmcorr', 'obstime', 'gs_commanded', 'fgsid', 'gs_position'])
Pointing.__new__.__defaults__ = ((None,) * 5)


# Transforms
@dataclasses.dataclass
class Transforms:
    m_eci2fgs1: np.array = None         # ECI to FGS1
    m_eci2gs: np.array = None           # ECI to Guide Star
    m_eci2j: np.array = None            # ECI to J-Frame
    m_eci2siaf: np.array = None         # ECI to SIAF
    m_eci2sifov: np.array = None        # ECI to SIFOV
    m_eci2v: np.array = None            # ECI to V
    m_fgs12fgsx: np.array = None        # FGS1 to FGSx transformation
    m_fgs12sifov: np.array = None       # FGS1 to SIFOV
    m_gs2gsapp: np.array = None         # Velocity abberation
    m_j2fgs1: np.array = None           # J-Frame to FGS1
    m_sifov_fsm_delta: np.array = None  # FSM correction
    m_sifov2v: np.array = None          # SIFOV to V1
    m_v2siaf: np.array = None           # V to SIAF
    override: object = None             # Override values. Either another Transforms or dict-like object
    """Transformation matrices"""

    @classmethod
    def from_asdf(cls, asdf_file):
        """Create Transforms from AsdfFile

        Parameters
        ----------
        asdf_file : Stream-like or `asdf.AsdfFile`
            The asdf to create from.

        Returns
        -------
        transforms : Transforms
            The Transforms instance.
        """
        if isinstance(asdf_file, asdf.AsdfFile):
            transforms = asdf_file.tree['transforms']
        else:
            with asdf.open(asdf_file, copy_arrays=True, lazy_load=False) as af:
                transforms = af.tree['transforms']

        return cls(**transforms)

    def to_asdf(self):
        """Serialize to AsdfFile

        Returns
        -------
        asdf_file : asdf.AsdfFile
            The ASDF serialization.

        Notes
        -----
        The `override` transforms are not serialized, since the values of this transform
        automatically represent what is in the override.
        """
        self_dict = dataclasses.asdict(self)
        del self_dict['override']  # Do not serialize the override transforms
        asdf_file = asdf.AsdfFile({'transforms': self_dict})
        return asdf_file

    def write_to_asdf(self, path):
        """Serialize to a file path

        Parameters
        ----------
        path : Stream-like
        """
        asdf_file = self.to_asdf()
        asdf_file.write_to(path, all_array_storage='inline')

    def __getattribute__(self, name):
        """If an override has been specified, return that value regardless

        Notes
        -----
        This dunder method is called for ALL attributes. Tread carefully.
        """
        # If the attribute is not a field, just return its value. Like NOW.
        if name.startswith('_') or name not in self._fields or name == 'override':
            return object.__getattribute__(self, name)

        override = self.override
        override_value = getattr(override, name) if override else None
        return override_value if override_value is not None else object.__getattribute__(self, name)

    def __post_init__(self):
        """Post-initialization of a DataClass"""

        # Create a simple list of fields to check against.
        self._fields = [field.name for field in dataclasses.fields(self)]


# WCS reference container
WCSRef = namedtuple('WCSRef', ['ra', 'dec', 'pa'])
WCSRef.__new__.__defaults__ = (None, None, None)


@dataclasses.dataclass
class TransformParameters:
    """Parameters required the calculations

    Parameters
    ----------
    allow_default : bool
        If telemetry cannot be determined, use existing
        information in the observation's header.

    default_pa_v3 : float
        The V3 position angle to use if the pointing information
        is not found.

    dry_run : bool
        Do not write out the modified file.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    fsmcorr_version : str
        The version of the FSM correction calculation to use.
        See :ref:`calc_sifov_fsm_delta_matrix`

    fsmcorr_units : str
        Units of the FSM correction values. Default is 'arcsec'.
        See :ref:`calc_sifov_fsm_delta_matrix`

    guide_star_wcs : WCSRef
        Guide star WCS info, typically from the input model.

    j2fgs_transpose : bool
        Transpose the `j2fgs1` matrix.

    jwst_velocity : numpy.array
        The [DX, DY, DZ] barycentri velocity vector

    method : `Methods`
        The method, or algorithm, to use in calculating the transform.
        If not specified, the default method is used.

    override_transforms : `Transforms`
        If set, matrices that should be used instead of the calculated one.

    pcs_mode : str
        The tracking mode in use.

    pointing : `Pointing`
        The quaternion.

    reduce_func : func or None
        Reduction function to use on values.

    siaf : `SIAF`
        The SIAF information for the input model

    siaf_db : `SiafDb`
        The SIAF database

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    useafter : str
        The date of observation (``model.meta.date``)
    """
    allow_default: bool = False
    default_pa_v3: float = 0.
    dry_run: bool = False
    engdb_url: str = None
    fsmcorr_version: str = 'latest'
    fsmcorr_units: str = 'arcsec'
    guide_star_wcs: WCSRef = WCSRef()
    j2fgs_transpose: bool = True
    jwst_velocity: np.array = None
    method: Methods = None
    override_transforms: Transforms = None
    pcs_mode: str = None
    pointing: Pointing = None
    reduce_func: typing.Callable = None
    siaf: SIAF = None
    siaf_db: SiafDb = None
    tolerance: float = 60.
    useafter: str = None
    v3pa_at_gs: float = None

    def as_reprdict(self):
        """Return a dict where all values are REPR of their values"""
        d = dataclasses.asdict(self)
        r = {key: repr(value) for key, value in d.items()}
        return r


def add_wcs(filename, default_pa_v3=0., siaf_path=None, engdb_url=None,
            tolerance=60, allow_default=False, reduce_func=None,
            dry_run=False, save_transforms=None, **transform_kwargs):
    """Add WCS information to a FITS file.

    Telescope orientation is attempted to be obtained from
    the engineering database. Failing that, a default pointing
    is used based on proposal target.

    The FITS file is updated in-place.

    Parameters
    ----------
    filename : str
        The path to a data file.

    default_pa_v3 : float
        The V3 position angle to use if the pointing information
        is not found.

    siaf_path : str or file-like object or None
        The path to the SIAF database. See `SiafDb` for more information.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default : bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func : func or None
        Reduction function to use on values.

    dry_run : bool
        Do not write out the modified file.

    save_transforms : Path-like or None
        File to save the calculated transforms to.

    transform_kwargs : dict
        Keyword arguments used by matrix calculation routines.

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
    logger.info('Updating WCS info for file %s', filename)
    with datamodels.Level1bModel(filename) as model:
        t_pars, transforms = update_wcs(
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
            logger.info('Saving updated model %s', filename)
            model.save(filename)
            if transforms and save_transforms:
                logger.info('Saving transform matrices to %s', save_transforms)
                transforms.write_to_asdf(save_transforms)

    logger.info('...update completed')


def update_mt_kwds(model):
    """Add/update the Moving target header keywords

    If the target type is "moving_target" check for the moving target position
    table. If this is available calculate the moving target position keywords
    and insert or update MT_RA & MT_DEC.
    """

    if model.hasattr('moving_target'):
        time_mt = Time(model.moving_target.time, format='isot')
        time_mt = [t.mjd for t in time_mt]
        exp_midpt_mjd = model.meta.exposure.mid_time
        # check to see if the midpoint of the observation is contained within
        # the timerange of the MT table
        if time_mt[0] <= exp_midpt_mjd <= time_mt[-1]:
            ra = model.moving_target.mt_apparent_RA
            dec = model.moving_target.mt_apparent_Dec
            f_ra = interp1d(time_mt, ra)
            f_dec = interp1d(time_mt, dec)
            model.meta.wcsinfo.mt_ra = f_ra(exp_midpt_mjd).item(0)
            model.meta.wcsinfo.mt_dec = f_dec(exp_midpt_mjd).item(0)
            model.meta.target.ra = f_ra(exp_midpt_mjd).item(0)
            model.meta.target.dec = f_dec(exp_midpt_mjd).item(0)
        else:
            logger.info('Exposure midpoint %s is not in the moving_target '
                        'table range of %s to %s', exp_midpt_mjd, time_mt[0], time_mt[-1])
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
        The path to the SIAF database. See `SiafDb` for more information.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default : bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func : func or None
        Reduction function to use on values.

    transform_kwargs : dict
        Keyword arguments used by matrix calculation routines.

    Returns
    -------
    t_pars, transforms : TransformParameters, Transforms
        The parameters and transforms calculated. May be
        None for either if telemetry calculations were not
        performed.
    """
    t_pars = transforms = None  # Assume telemetry is not used.

    # Open the SIAF dabase
    with SiafDb(siaf_path) as siaf_db:

        # If the type of exposure is not FGS, then attempt to get pointing
        # from telemetry.
        try:
            exp_type = model.meta.exposure.type.lower()
        except AttributeError:
            exp_type = None
        aperture_name = model.meta.aperture.name.upper()
        useafter = model.meta.observation.date
        if aperture_name != "UNKNOWN":
            logger.info("Updating WCS for aperture %s", aperture_name)
            siaf = siaf_db.get_wcs(aperture_name, useafter)
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
            t_pars = TransformParameters(
                default_pa_v3=default_pa_v3, siaf=siaf, engdb_url=engdb_url,
                tolerance=tolerance, allow_default=allow_default,
                reduce_func=reduce_func, siaf_db=siaf_db, useafter=useafter,
                **transform_kwargs
            )
            transforms = update_wcs_from_telem(model, t_pars)

    return t_pars, transforms


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

    default_vparity : int
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
        logger.warning('Keyword `ROLL_REF` not found. Using %s as default value', default_roll_ref)
        roll_ref = default_roll_ref

    roll_ref = np.deg2rad(roll_ref)

    # Get VIdlParity
    try:
        vparity = model.meta.wcsinfo.vparity
    except AttributeError:
        logger.warning('Keyword "VPARITY" not found. Using %s as default value', default_vparity)
        vparity = default_vparity

    try:
        v3i_yang = model.meta.wcsinfo.v3yangle
    except AttributeError:
        logger.warning('Keyword "V3I_YANG" not found. Using %s as default value.', default_v3yangle)
        v3i_yang = default_v3yangle

    (
        model.meta.wcsinfo.pc1_1,
        model.meta.wcsinfo.pc1_2,
        model.meta.wcsinfo.pc2_1,
        model.meta.wcsinfo.pc2_2
    ) = calc_rotation_matrix(roll_ref, np.deg2rad(v3i_yang), vparity=vparity)


def update_wcs_from_telem(model, t_pars: TransformParameters):
    """Update WCS pointing information

    Given a `jwst.datamodels.DataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update. The update is done in-place.

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms or None
        If available, the transformation matrices.
    """
    logger.info('Updating wcs from telemetry.')
    transforms = None  # Assume no transforms are calculated.

    # observation parameters
    obsstart = model.meta.exposure.start_time
    obsend = model.meta.exposure.end_time

    # Check that the V-frame to SIAF refpoint has been defined.
    v2siaf_pars = [getattr(t_pars.siaf, attr) for attr in ['v2_ref', 'v3_ref', 'v3yangle', 'vparity']]
    if None in v2siaf_pars:
        raise ValueError('Insufficient SIAF information found in header.')

    # Setup default WCS info if actual pointing and calculations fail.
    wcsinfo = WCSRef(
        model.meta.target.ra,
        model.meta.target.dec,
        t_pars.default_pa_v3
    )
    vinfo = wcsinfo

    # Get Guide Star information
    t_pars.guide_star_wcs = WCSRef(
        model.meta.guidestar.gs_ra,
        model.meta.guidestar.gs_dec,
        model.meta.guidestar.gs_v3_pa_science
    )
    t_pars.pcs_mode = model.meta.guidestar.gs_pcs_mode
    logger.debug('guide_star_wcs from model: %s', t_pars.guide_star_wcs)
    logger.debug('PCS_MODE: %s', t_pars.pcs_mode)

    # Get jwst velocity
    t_pars.jwst_velocity = np.array([
        model.meta.ephemeris.velocity_x_bary,
        model.meta.ephemeris.velocity_y_bary,
        model.meta.ephemeris.velocity_z_bary,
    ])
    logger.debug('JWST Velocity: %s', t_pars.jwst_velocity)

    # Retrieve tracking mode.

    # Get the pointing information
    try:
        pointing = get_pointing(obsstart, obsend, engdb_url=t_pars.engdb_url,
                                tolerance=t_pars.tolerance, reduce_func=t_pars.reduce_func)
    except ValueError as exception:
        if not t_pars.allow_default:
            raise
        else:
            logger.warning(
                'Cannot retrieve valid telescope pointing.'
                ' Default pointing parameters will be used.'
            )
            logger.warning('Exception is %s', exception)
            logger.info("Setting ENGQLPTG keyword to PLANNED")
            model.meta.visit.engdb_pointing_quality = "PLANNED"
    else:
        # compute relevant WCS information
        logger.info('Successful read of engineering quaternions:')
        logger.info('\tPointing: %s', pointing)
        try:
            t_pars.pointing = pointing
            wcsinfo, vinfo, transforms = calc_wcs(t_pars)
            pointing_engdb_quality = f'CALCULATED_{t_pars.method.value.upper()}'
            logger.info('Setting ENGQLPTG keyword to %s', pointing_engdb_quality)
            model.meta.visit.engdb_pointing_quality = pointing_engdb_quality
        except Exception as e:
            logger.warning(
                'WCS calculation has failed and will be skipped.'
                'Default pointing parameters will be used.'
            )
            logger.warning('Exception is %s', e)
            if not t_pars.allow_default:
                raise
            else:
                logger.info("Setting ENGQLPTG keyword to PLANNED")
                model.meta.visit.engdb_pointing_quality = "PLANNED"
    logger.info('Aperture WCS info: %s', wcsinfo)
    logger.info('V1 WCS info: %s', vinfo)

    # Update V1 pointing
    model.meta.pointing.ra_v1 = vinfo.ra
    model.meta.pointing.dec_v1 = vinfo.dec
    model.meta.pointing.pa_v3 = vinfo.pa

    # Update Aperture pointing
    model.meta.aperture.position_angle = wcsinfo.pa
    model.meta.wcsinfo.ra_ref = wcsinfo.ra
    model.meta.wcsinfo.dec_ref = wcsinfo.dec
    model.meta.wcsinfo.roll_ref = compute_local_roll(
        vinfo.pa, wcsinfo.ra, wcsinfo.dec, t_pars.siaf.v2_ref, t_pars.siaf.v3_ref
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
            vparity=t_pars.siaf.vparity
        )

    # Calculate S_REGION with the footprint
    # information
    try:
        update_s_region(model, t_pars.siaf)
    except Exception as e:
        logger.warning('Calculation of S_REGION failed and will be skipped.')
        logger.warning('Exception is %s', e)

    return transforms


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
    logger.info("Vertices for aperture %s: %s", model.meta.aperture.name, vertices)

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


def calc_wcs_over_time(obsstart, obsend, t_pars: TransformParameters):
    """Calculate V1 and WCS over a time period

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

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
        pointings = get_pointing(obsstart, obsend, engdb_url=t_pars.engdb_url,
                                 tolerance=t_pars.tolerance, reduce_func=t_pars.reduce_func)
    except ValueError:
        logger.warning("Cannot get valid engineering mnemonics from engineering database")
        raise
    if not isinstance(pointings, list):
        pointings = [pointings]
    for pointing in pointings:
        t_pars.pointing = pointing
        wcsinfo, vinfo, transforms = calc_wcs(t_pars)
        obstimes.append(pointing.obstime)
        wcsinfos.append(wcsinfo)
        vinfos.append(vinfo)

    return obstimes, wcsinfos, vinfos


def calc_wcs(t_pars: TransformParameters):
    """Transform from the given SIAF information and Pointing
    the aperture and V1 wcs

    Parameters
    ----------
    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    wcsinfo, vinfo, transforms : WCSRef, WCSRef, Transforms
        A 3-tuple is returned with the WCS pointing for
        the aperture and the V1 axis, and the transformation matrices.

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
    if t_pars.siaf is None:
        t_pars.siaf = SIAF()

    # Calculate transforms
    tforms = calc_transforms(t_pars)

    # Calculate the V1 WCS information
    vinfo = calc_wcs_from_matrix(tforms.m_eci2v)

    # Calculate the Aperture WCS
    wcsinfo = calc_aperture_wcs(tforms.m_eci2siaf)

    # That's all folks
    return wcsinfo, vinfo, tforms


def calc_transforms(t_pars: TransformParameters):
    """Calculate transforms  which determine reference point celestial WCS

    Given the spacecraft pointing parameters and the
    aperture-specific SIAF, calculate all the transforms
    necessary to produce WCS information.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : `Transforms`
        The list of coordinate matrix transformations
    """
    t_pars.method = t_pars.method if t_pars.method else Methods.default
    transforms = t_pars.method.func(t_pars)
    return transforms


def calc_transforms_tr202105(t_pars: TransformParameters):
    """Calculate transforms which determine reference point celestial WCS from the original, pre-JSOCINT-555 algorithm

    Given the spacecraft pointing parameters and the aperture-specific SIAF,
    calculate all the transforms necessary to produce the celestial World
    Coordinate system information.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    The matrix transform pipeline to convert from ECI J2000 observatory
    qauternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =       # Complete transformation
            M_v_to_siaf    *  # V to SIAF
            M_eci_to_v        # ECI to V

        where

        M_eci_to_v =
            M_sifov_to_v   *  # SIFOV to V
            M_z_to_j       *  # ICS-frame to SIAF-frame
            M_eci_to_sifov    # ECI to SIFOV

        M_eci_to_sifov =
            M_z_to_x          *  # Rotate from Z out to X out
            M_sifov_fsm_delta *  # SIFOV correstion due to Fast Steering Mirror offsets
            M_fgs1_to_sifov   *  # FGS1 to SIFOV
            M_j_to_fgs1       *  # J-frame to FGS1
            M_eci_to_j           # ECI to J-Frame

    """
    logger.info('Calculating transforms using FULL quaternion method...')
    t_pars.method = Methods.TR_202105
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # ---
    # Calculate base matrices from telemetry, siaf database, and/or constants
    # ___

    # Determine the ECI to J-frame matrix
    t.m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)

    # Calculate the J-frame to FGS1 ICS matrix
    t.m_j2fgs1 = calc_j2fgs1_matrix(t_pars.pointing.j2fgs_matrix, transpose=t_pars.j2fgs_transpose)

    # Calculate the FGS1 ICS to SI-FOV matrix
    t.m_fgs12sifov = calc_fgs1_to_sifov_matrix(t_pars.siaf_db, t_pars.useafter)

    # Calculate the FSM corrections to the SI_FOV frame
    t.m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        t_pars.pointing.fsmcorr, fsmcorr_version=t_pars.fsmcorr_version, fsmcorr_units=t_pars.fsmcorr_units
    )

    # Calculate SI FOV to V1 matrix
    t.m_sifov2v = calc_sifov2v_matrix()

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # ---
    # Calculate matrices based on previously computed matrices
    # ---

    # Calculate ECI to SIFOV complete transformation
    t.m_eci2sifov = np.linalg.multi_dot([
        MZ2X, t.m_sifov_fsm_delta, t.m_fgs12sifov, t.m_j2fgs1, t.m_eci2j
    ])
    logger.debug('m_eci2sifov: %s', t.m_eci2sifov)

    # Calculate the complete transform to the V1 reference
    t.m_eci2v = np.dot(t.m_sifov2v, t.m_eci2sifov)
    logger.debug('m_eci2v: %s', t.m_eci2v)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_coarse_tr_202107(t_pars: TransformParameters):
    """Calculate transforms for  guiding as per TR presented in 2021-07

    This implements equation 42 from Technical Report JWST-STScI-003222, SM-12. 2021-07
    From Section 5:

    In COARSE mode the measured attitude of the J-frame of the spacecraft is
    determined by the star tracker and inertial gyroscopes attitude
    measurements and is converted to an estimated guide star inertial attitude
    using the equations in section 3.2. The V-frame attitude then is
    determined.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    The matrix transform pipeline to convert from ECI J2000 observatory
    qauternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =
            transpose(M_v_to_fgsx)  *
            transpose(M_fgsx_to_gs) *
            M_x_to_z                *
            M_eci_to_gs

        where

            M_fgsx_to_v = FGSx to V-frame
            M_gs_to_fgsx = Guide star to FGSx
            M_eci_to_gs = ECI to Guide star
    """
    logger.info('Calculating transforms using TR 202107 COARSE Tracking method...')
    t_pars.method = Methods.COARSE_TR_202107

    # Determine the M_eci_to_gs matrix. Since this is a full train, the matrix
    # is returned as part of the full Transforms object. Many of the required
    # matrices are already determined as part of this calculation.
    t = calc_m_eci2gs(t_pars)

    # Determine the M_fgsx_to_v matrix
    siaf = t_pars.siaf_db.get_wcs(FGSId2Aper[t_pars.pointing.fgsid])
    t.m_v2fgsx = calc_v2siaf_matrix(siaf)

    # Determine M_eci_to_v frame.
    t.m_eci2v = np.linalg.multi_dot([np.transpose(t.m_v2fgsx), np.transpose(t.m_fgsx2gs), MX2Z, t.m_eci2gs])
    logger.debug('M_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate full transformation
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_track_tr_202107(t_pars: TransformParameters):
    """Calculate transforms for TRACK/FINEGUIDE guiding as per TR presented in 2021-07

    This implements equation 43 from Technical Report JWST-STScI-003222, SM-12. 2021-07
    From Section 5:

    Under guide star control the guide star position is measured relative to
    the V-frame. The V3 position angle at the guide star is derived from the
    measured J-frame attitude. Then using the corrected guide star catalog
    position the data be used to determine the inertial V-frame attitude on the
    sky.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    The matrix transform pipeline to convert from ECI J2000 observatory
    qauternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =       # Complete transformation
            M_v_to_siaf    *  # V to SIAF
            M_eci_to_v        # ECI to V

        where

            M_eci_to_v = Conversion of the attitude to a DCM
    """
    logger.info('Calculating transforms using TR 202107 TRACK/FINEGUIDE Tracking method...')
    t_pars.method = Methods.TRACK_TR_202107
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Determine V3PA@GS
    v3pags = calc_v3pags(t_pars)
    t_pars.guide_star_wcs = WCSRef(t_pars.guide_star_wcs.ra, t_pars.guide_star_wcs.dec, v3pags)

    # Transform the guide star location in ideal detector coordinates to the telescope/V23 frame.
    gs_pos_v23 = trans_fgs2v(t_pars.pointing.fgsid, t_pars.pointing.gs_position, t_pars.siaf_db)

    # Calculate the M_eci2v matrix. This is the attitude matrix of the observatory
    # relative to the guide star.
    t.m_eci2v = calc_attitude_matrix(t_pars.guide_star_wcs, v3pags, gs_pos_v23)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_ops_tr_202107(t_pars: TransformParameters):
    """Calculate transforms as per TR presented in 2021-07

    This implements the ECI-to-SIAF transformation from Technical Report JWST-STScI-003222, SM-12. 2021-07
    The actual implementation depends on the guide star mode, represented by the header keyword PCS_MODE.
    For COARSE or NONE, the method COARSE_TR_202107 is used.
    For TRACK or FINEGUIDE, the method TRACK_TR_202107 is used.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations
    """
    if t_pars.pcs_mode is None or t_pars.pcs_mode in ['NONE', 'COARSE']:
        return calc_transforms_coarse_tr_202107(t_pars)
    elif t_pars.pcs_mode in ['TRACK', 'FINEGUIDE']:
        return calc_transforms_track_tr_202107(t_pars)
    else:
        raise ValueError(f'Invalid PCS_MODE: {t_pars.pcs_mode}. Should be in ["NONE", "COARSE", "TRACK", "FINEGUIDE"]')


def calc_transforms_velocity_abberation_tr202105(t_pars: TransformParameters):
    """Calculate transforms which determine reference point celestial WCS from the original, pre-JSOCINT-555 algorithm

    Given the spacecraft pointing parameters and the aperture-specific SIAF,
    calculate all the transforms necessary to produce the celestial World
    Coordinate system information.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
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
            M_z_to_x          *   # Transposition
            M_fgs1_to_sifov   *   # FGS1 to Science Instruments Aperture
            M_gs2gsapp        *   # Apply Velocity Aberration
            M_j_to_fgs1       *   # J-Frame to FGS1
            M_eci_to_j        *   # ECI to J-Frame

    """
    logger.info('Calculating transforms using FULLVA quaternion method with velocity aberration...')
    t_pars.method = Methods.TR_202105_VA
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Determine the ECI to J-frame matrix
    t.m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)

    # Calculate the J-frame to FGS1 ICS matrix
    t.m_j2fgs1 = calc_j2fgs1_matrix(t_pars.pointing.j2fgs_matrix, transpose=t_pars.j2fgs_transpose)

    t.m_eci2fgs1 = np.dot(t.m_j2fgs1, t.m_eci2j)
    logger.debug('m_eci2fgs1: %s', t.m_eci2fgs1)

    # Calculate Velocity Aberration corrections
    t.m_gs2gsapp = calc_gs2gsapp_tr_202105(t.m_eci2fgs1, t_pars.jwst_velocity)

    # Calculate the FSM corrections to the SI_FOV frame
    t.m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        t_pars.pointing.fsmcorr, fsmcorr_version=t_pars.fsmcorr_version, fsmcorr_units=t_pars.fsmcorr_units
    )

    # Calculate the FGS1 ICS to SI-FOV matrix
    t.m_fgs12sifov = calc_fgs1_to_sifov_matrix(t_pars.siaf_db, t_pars.useafter)

    # Calculate SI FOV to V1 matrix
    t.m_sifov2v = calc_sifov2v_matrix()

    # Calculate ECI to SI FOV
    t.m_eci2sifov = np.linalg.multi_dot(
        [MZ2X, t.m_sifov_fsm_delta, t.m_fgs12sifov, t.m_gs2gsapp, t.m_j2fgs1, t.m_eci2j]
    )
    logger.debug('m_eci2sifov: %s', t.m_eci2sifov)

    # Calculate the complete transform to the V1 reference
    t.m_eci2v = np.dot(t.m_sifov2v, t.m_eci2sifov)
    logger.debug('m_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_original(t_pars: TransformParameters):
    """Calculate transforms which determine reference point celestial WCS from the original, pre-JSOCINT-555 algorithm

    Given the spacecraft pointing parameters and the aperture-specific SIAF,
    calculate all the transforms necessary to produce the celestial World
    Coordinate system information.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
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
            M_z_to_x          *   # Transposition
            M_fgs1_to_sifov   *   # FGS1 to Science Instruments Aperture
            M_j_to_fgs1       *   # J-Frame to FGS1
            M_eci_to_j        *   # ECI to J-Frame

    """
    logger.info('Calculating transforms using ORIGINAL method...')
    t_pars.method = Methods.ORIGINAL
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Determine the ECI to J-frame matrix
    t.m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)

    # Calculate the J-frame to FGS1 ICS matrix
    t.m_j2fgs1 = calc_j2fgs1_matrix(t_pars.pointing.j2fgs_matrix, transpose=t_pars.j2fgs_transpose)

    # Calculate the FSM corrections to the SI_FOV frame
    t.m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        t_pars.pointing.fsmcorr, fsmcorr_version=t_pars.fsmcorr_version, fsmcorr_units=t_pars.fsmcorr_units
    )

    # Calculate the FGS1 ICS to SI-FOV matrix
    # The SIAF is explicitly not used because this was not accounted for in
    # the original specification.
    t.m_fgs12sifov = calc_fgs1_to_sifov_matrix(siaf_db=None, useafter=None)

    # Calculate SI FOV to V1 matrix
    t.m_sifov2v = calc_sifov2v_matrix()

    # Calculate ECI to SI FOV
    t.m_eci2sifov = np.linalg.multi_dot(
        [t.m_sifov_fsm_delta, MZ2X, t.m_fgs12sifov, t.m_j2fgs1, t.m_eci2j]
    )
    logger.debug('m_eci2sifov: %s', t.m_eci2sifov)

    # Calculate the complete transform to the V1 reference
    t.m_eci2v = np.dot(t.m_sifov2v, t.m_eci2sifov)
    logger.debug('m_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_gscmd_j3pags(t_pars: TransformParameters):
    """Calculate transforms which determine reference point celestial WCS using the JSOCINT-555 fix and J3PA@GS

    Given the spacecraft pointing parameters and the
    aperture-specific SIAF, calculate all the transforms
    necessary to produce WCS information. This is using a modified
    algorithm to fix a documentation error and to override a shift
    being introduced by the OTB simulator.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    This algorithm is in response to an apparent misbehavior of the OTB simulator.
    See JSOCINT-555 for full details. The workaround is to ignore the quaternion
    and FGS1 J-frame telemetry. Instead use the commanded guide star telemetry.

    The matrix transform pipeline to convert from ECI J2000 observatory
    qauternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =                    # The complete transformation
            M_v1_to_siaf               *   # V1 to SIAF
            M_sifov_to_v1              *   # Science Instruments Aperture to V1
            M_z_to_x                   *   # Transposition
            M_sifov_fsm_delta          *   # Fine Steering Mirror correction
            M_fgs1_to_sifov_fgs1siaf   *   # FGS1 to Science Instruments Aperture
            M_eci_to_fgs1                  # ECI to FGS1 using commanded information
    """
    logger.info('Calculating transforms using GSCMD with J3PA@GS method...')
    t_pars.method = Methods.GSCMD_J3PAGS
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Determine the ECI to FGS1 ICS using guide star telemetry.
    t.m_eci2fgs1 = calc_eci2fgs1_j3pags(t_pars)

    # Calculate the FGS1 ICS to SI-FOV matrix
    t.m_fgs12sifov = calc_fgs1_to_sifov_matrix(siaf_db=t_pars.siaf_db, useafter=t_pars.useafter)

    # Calculate the FSM corrections to the SI_FOV frame
    t.m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        t_pars.pointing.fsmcorr, fsmcorr_version=t_pars.fsmcorr_version, fsmcorr_units=t_pars.fsmcorr_units
    )

    # Calculate SI FOV to V1 matrix
    t.m_sifov2v = calc_sifov2v_matrix()

    # Calculate ECI to SI FOV
    t.m_eci2sifov = np.linalg.multi_dot(
        [MZ2X, t.m_sifov_fsm_delta, t.m_fgs12sifov, t.m_eci2fgs1]
    )
    logger.debug('m_eci2sifov: %s', t.m_eci2sifov)

    # Calculate the complete transform to the V1 reference
    t.m_eci2v = np.dot(t.m_sifov2v, t.m_eci2sifov)
    logger.debug('m_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_gscmd_v3pags(t_pars: TransformParameters):
    """Calculate transforms which determine reference point celestial WCS using the JSOCINT-555 fix

    Given the spacecraft pointing parameters and the
    aperture-specific SIAF, calculate all the transforms
    necessary to produce WCS information. This is using a modified
    algorithm to fix a documentation error and to override a shift
    being introduced by the OTB simulator.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    This algorithm is in response to an apparent misbehavior of the OTB simulator.
    See JSOCINT-555 for full details. The workaround is to ignore the quaternion
    and FGS1 J-frame telemetry. Instead use the commanded guide star telemetry.

    The matrix transform pipeline to convert from ECI J2000 observatory
    quaternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =                    # The complete transformation
            M_v1_to_siaf               *   # V1 to SIAF
            M_sifov_to_v1              *   # Science Instruments Aperture to V1
            M_z_to_x                   *   # Transposition
            M_sifov_fsm_delta          *   # Fine Steering Mirror correction
            M_fgs1_to_sifov_fgs1siaf   *   # FGS1 to Science Instruments Aperture
            M_eci_to_fgs1                  # ECI to FGS1 using commanded information
    """
    logger.info('Calculating transforms using GSCMD using V3PA@GS method...')
    t_pars.method = Methods.GSCMD_V3PAGS
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Determine the ECI to FGS1 ICS using guide star telemetry.
    t.m_eci2fgs1 = calc_eci2fgs1_v3pags(t_pars)

    # Calculate the FSM corrections to the SI_FOV frame
    t.m_sifov_fsm_delta = calc_sifov_fsm_delta_matrix(
        t_pars.pointing.fsmcorr, fsmcorr_version=t_pars.fsmcorr_version, fsmcorr_units=t_pars.fsmcorr_units
    )

    # Calculate the FGS1 ICS to SI-FOV matrix
    t.m_fgs12sifov = calc_fgs1_to_sifov_matrix(siaf_db=t_pars.siaf_db, useafter=t_pars.useafter)

    # Calculate SI FOV to V1 matrix
    t.m_sifov2v = calc_sifov2v_matrix()

    # Calculate ECI to SI FOV
    t.m_eci2sifov = np.linalg.multi_dot(
        [MZ2X, t.m_sifov_fsm_delta, t.m_fgs12sifov, t.m_eci2fgs1]
    )
    logger.debug('m_eci2sifov: %s', t.m_eci2sifov)

    # Calculate the complete transform to the V1 reference
    t.m_eci2v = np.dot(t.m_sifov2v, t.m_eci2sifov)
    logger.debug('m_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.dot(t.m_v2siaf, t.m_eci2v)
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_eci2fgs1_j3pags(t_pars: TransformParameters):
    """Calculate full ECI to FGS1 matrix based on commanded guidestar information and J3PA@GS

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters

    Return
    ------
    m_eci2fgs1 : numpy.array
        The ECI to FGS1ics matrix

    Modifies
    --------
    t_pars.guide_star_wcs
        Update PA to calculated if originally None.
    """
    # Calculate J3 position.
    m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)
    j3_ra, j3_dec = vector_to_angle(m_eci2j[2])
    j3wcs = WCSRef(j3_ra, j3_dec, None)
    logger.debug('j3wcs: %s', j3wcs)

    # Calculate J3PA@GS
    logger.debug('Incoming guide_star_wcs: %s', t_pars.guide_star_wcs)
    if not t_pars.guide_star_wcs.pa:
        gs_wcs_ra = WCSRef(t_pars.guide_star_wcs.ra * D2R,
                           t_pars.guide_star_wcs.dec * D2R, None)
        pa = calc_position_angle(gs_wcs_ra, j3wcs)
        pa *= R2D
        t_pars.guide_star_wcs = WCSRef(t_pars.guide_star_wcs.ra, t_pars.guide_star_wcs.dec, pa)

    logger.debug('guide_star_wcs: %s', t_pars.guide_star_wcs)
    logger.debug('gs_commanded: %s', t_pars.pointing.gs_commanded)

    m_gs_commanded = calc_m_gs_commanded(t_pars.guide_star_wcs, J3IDLYANGLE, t_pars.pointing.gs_commanded)

    # Need to invert the martrix
    m_gs_commanded = np.linalg.inv(m_gs_commanded)

    logger.debug('m_gs_commanded: %s', m_gs_commanded)

    m_eci2fgs1 = np.dot(MX2Z, m_gs_commanded)

    logger.debug('m_eci2fgs1: %s', m_eci2fgs1)
    return m_eci2fgs1


def calc_gs2gsapp(m_eci2gsics, jwst_velocity):
    """Calculate the Velocity Aberration correction

    This implements equation 37 from Technical Report JWST-STScI-003222, SM-12. 2021-07
    From Section 3.2.5:

    The velocity aberration correction is applied in the direction of the of
    the guide star. The matrix that translates from ECI to the guide star ICS
    frame is M_(ECI→GSICS), where the GS position vector is along the z-axis in
    the guide star ICS frame.

    Parameters
    ----------
    m_eci2gsics : numpy.array(3, 3)
        The the ECI to Guide Star transformation matrix, in the ICS frame.

    jwst_velocity : numpy.array([dx, dy, dz])
        The barycentric velocity of JWST.

    Returns
    -------
    m_gs2gsapp : numpy.array(3, 3)
        The velocity aberration correction matrix.
    """
    # Check velocity. If present, negate the velocity since
    # the desire is to remove the correction.
    if jwst_velocity is None or any(jwst_velocity == None):  # noqa Syntax needed for numpy arrays.
        logger.warning('Velocity: %s contains None. Cannot calculate aberration. Returning identity matrix', jwst_velocity)
        return np.identity(3)
    velocity = -1 * jwst_velocity

    # Eq. 31: Guide star position vector
    uz = np.array([0., 0., 1.])
    u_gseci = np.dot(np.transpose(m_eci2gsics), uz)

    # Eq. 33: Compute the apparent shift due to velocity aberration.
    try:
        scale_factor, u_gseci_app = compute_va_effects_vector(*velocity, u_gseci)
    except TypeError:
        logger.warning('Failure in computing velocity aberration. Returning identity matrix.')
        logger.warning('Exception: %s', sys.exc_info())
        return np.identity(3)

    # Eq. 36: Rotate from ICS into the guide star frame.
    u_gs_app = np.dot(m_eci2gsics, u_gseci_app)

    # Eq. 37: Compute the M_gs2gsapp matrix
    u_prod = np.cross(uz, u_gs_app)
    u_prod_mag = np.linalg.norm(u_prod)
    a_hat = u_prod / u_prod_mag
    m_a_hat = np.array([[0., -a_hat[2], a_hat[1]],
                        [a_hat[2], 0., -a_hat[0]],
                        [-a_hat[1], a_hat[0], 0.]])
    theta = np.arcsin(u_prod_mag)

    m_gs2gsapp = np.identity(3) \
        - (m_a_hat * np.sin(theta)) \
        + (2 * m_a_hat**2 * np.sin(theta / 2.)**2)

    logger.debug('m_gs2gsapp: %s', m_gs2gsapp)
    return m_gs2gsapp


def calc_gs2gsapp_tr_202105(m_eci2fgs1, jwst_velocity):
    """Calculate the Velocity Aberration correction

    Parameters
    ----------
    m_eci2fgs : numpy.array(3, 3)
        The matrix defining the FGS1, guide star, position

    jwst_velocity : numpy.array([dx, dy, dz])
        The barycentric velocity of JWST

    Returns
    -------
    m_gs2gsapp : numpy.array(3, 3)
        The corection matrix

    Notes
    -----
    The algorithm is from an NGAS memo describing the procedure. The steps referred
    to in the commentary are in direct relation to the memo's text.

    Step 1 is the calculation of the Meci2fgs1 matrix, which is the input.
    """
    ux = np.array([1., 0., 0.])

    # Step 2: Compute guide star unit vector
    u = np.dot(np.transpose(m_eci2fgs1), ux)

    # Step 3: Apply VA
    try:
        scale_factor, u_j2000 = compute_va_effects_vector(*jwst_velocity, u)
    except TypeError:
        logger.warning('Failure in computing velocity aberration. Returning identity matrix.')
        logger.warning('Exception: %s', sys.exc_info())
        return np.identity(3)

    # Step 4: Transform to guide star attitude frame
    u_gs = np.dot(m_eci2fgs1, u_j2000)

    # Step 5: Compute apparent attitude correction
    u_prod = np.cross(ux, u_gs)
    a_hat = u_prod / np.linalg.norm(u_prod)
    a_hat_sym = np.array([[0., -a_hat[2], a_hat[1]],
                          [a_hat[2], 0., -a_hat[0]],
                          [-a_hat[1], a_hat[0], 0.]])
    theta = np.arccos(np.dot(ux, u_gs) / np.linalg.norm(u_gs))

    m_gs2gsapp = np.identity(3) - a_hat_sym * np.sin(theta) + 2 * a_hat_sym**2 * np.sin(theta / 2.)**2

    logger.debug('m_gs2gsapp: %s', m_gs2gsapp)
    return m_gs2gsapp


def calc_eci2fgs1_v3pags(t_pars: TransformParameters):
    """Calculate full ECI to FGS1 matrix based on commanded guidestar information andd V3PA@GS

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.
    """
    fgs1_siaf = t_pars.siaf_db.get_wcs('FGS1_FULL_OSS', useafter=t_pars.useafter)
    logger.debug('fgs1_siaf: %s', fgs1_siaf)

    logger.debug('Incoming guide_star_wcs: %s', t_pars.guide_star_wcs)
    if not t_pars.guide_star_wcs.pa:
        pa = calc_v3pa_at_gs_from_original(t_pars)
        t_pars.guide_star_wcs = WCSRef(t_pars.guide_star_wcs.ra, t_pars.guide_star_wcs.dec, pa)

    logger.debug('guide_star_wcs: %s', t_pars.guide_star_wcs)
    logger.debug('gs_commanded: %s', t_pars.pointing.gs_commanded)

    m_gs_commanded = calc_m_gs_commanded(t_pars.guide_star_wcs, fgs1_siaf.v3yangle, t_pars.pointing.gs_commanded)

    # Need to invert the martrix
    m_gs_commanded = np.linalg.inv(m_gs_commanded)

    logger.debug('m_gs_commanded: %s', m_gs_commanded)

    m_eci2fgs1 = np.dot(MX2Z, m_gs_commanded)

    logger.debug('m_eci2fgs1: %s', m_eci2fgs1)
    return m_eci2fgs1


def calc_attitude_matrix(wcs, yangle, position):
    """Calculate the DCM attitude from known positions and roll angles.

    This implements Appendix A from Technical Report JWST-STScI-003222, SM-12. 2021-07

    Parameters
    ----------
    wcs : WCSRef
        The guide star position

    yangle : float
        The IdlYangle of the point in question.

    position : numpy.array(2)
        The position in Ideal frame.

    Returns
    -------
    m : np.array(3,3)
        The transformation matrix
    """
    # Convert to radians
    ra = wcs.ra * D2R
    dec = wcs.dec * D2R
    yangle_ra = yangle * D2R
    pos_rads = position * A2R
    v2 = pos_rads[0]
    v3 = pos_rads[1]

    # Create the matrices
    r1 = np.array([
        [cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)],
        [(-sin(ra) * cos(yangle_ra)) + (cos(ra) * sin(dec) * sin(yangle_ra)),
         (cos(ra) * cos(yangle_ra)) + (sin(ra) * sin(dec) * sin(yangle_ra)),
         -cos(dec) * sin(yangle_ra)],
        [(-sin(ra) * sin(yangle_ra)) - (cos(ra) * sin(dec) * cos(yangle_ra)),
         (cos(ra) * sin(yangle_ra)) - (sin(ra) * sin(dec) * cos(yangle_ra)),
         cos(dec) * cos(yangle_ra)]
    ])

    r2 = np.array([
        [cos(v2) * cos(v3), -sin(v2), -cos(v2) * sin(v3)],
        [sin(v2) * cos(v3), cos(v2), -sin(v2) * sin(v3)],
        [sin(v3), 0., cos(v3)]
    ])

    # Final transformation
    m = np.dot(r2, r1)

    logger.debug('attitude DCM: %s', m)
    return m


def calc_m_gs_commanded(guide_star_wcs, yangle, gs_commanded):
    """Calculate the guide star matrix

    Parameters
    ----------
    guide_star_wcs : WCSRef
        The guide star position

    yangle : float
        The IdlYangle of the point in question.

    gs_commanded : numpy.array(2)
        The commanded position from telemetry

    Returns
    -------
    m_gs_commanded : np.array(3,3)
        The transformation matrix
    """

    # Define the individal rotations
    def r1(a):
        r = np.array([
            [1, 0, 0],
            [0, cos(a), -sin(a)],
            [0, sin(a), cos(a)]
        ])
        return r

    def r2(a):
        r = np.array([
            [cos(a), 0, sin(a)],
            [0, 1, 0],
            [-sin(a), 0, cos(a)]
        ])
        return r

    def r3(a):
        r = np.array([
            [cos(a), -sin(a), 0],
            [sin(a), cos(a), 0],
            [0, 0, 1]
        ])
        return r

    # Convert to radians
    ra = guide_star_wcs.ra * D2R
    dec = guide_star_wcs.dec * D2R
    pa = guide_star_wcs.pa * D2R
    yangle_ra = yangle * D2R
    gs_commanded_rads = gs_commanded * A2R

    # Calculate
    m = np.linalg.multi_dot([r3(ra), r2(-dec), r1(-(pa + yangle_ra)), r2(gs_commanded_rads[1]), r3(-gs_commanded_rads[0])])

    logger.debug('m: %s', m)
    return m


def calc_v3pa_at_gs_from_original(t_pars: TransformParameters) -> float:
    """Calculate V3PA@GS by precomputing V1 using original method

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    v3pa_at_gs : float
    """

    # Calculate V1 using ORIGINAL method
    tforms = calc_transforms_original(deepcopy(t_pars))
    vinfo = calc_wcs_from_matrix(tforms.m_eci2v)

    # Convert to radians
    ra_v1 = vinfo.ra * D2R
    dec_v1 = vinfo.dec * D2R
    pa_v3 = vinfo.pa * D2R
    ra_gs = t_pars.guide_star_wcs.ra * D2R
    dec_gs = t_pars.guide_star_wcs.dec * D2R

    # Calculate
    x = -sin(ra_v1) * sin(pa_v3) - cos(ra_v1) * sin(dec_v1) * cos(pa_v3)
    y = cos(ra_v1) * sin(pa_v3) - sin(ra_v1) * sin(dec_v1) * cos(pa_v3)

    ra_v3 = atan2(y, x)
    dec_v3 = asin(cos(dec_v1) * cos(pa_v3))

    v3pa_at_gs = calc_position_angle(WCSRef(ra_gs, dec_gs, None), WCSRef(ra_v3, dec_v3, None))
    v3pa_at_gs *= R2D

    logger.debug('v3pa_at_gs: %s', v3pa_at_gs)
    return v3pa_at_gs


def calc_wcs_from_matrix(m):
    """Calculate the WCS information from a DCM.

    Parameters
    ----------
    m: np.array((3, 3))
        The DCM matrix to extract WCS information from

    Returns
    -------
    wcs : WCSRef
        The WCS.
    """
    # V1 RA/Dec is the first row of the transform
    v1_ra, v1_dec = vector_to_angle(m[0])
    wcs = WCSRef(v1_ra, v1_dec, None)

    # V3 is the third row of the transformation
    v3_ra, v3_dec = vector_to_angle(m[2])
    v3wcs = WCSRef(v3_ra, v3_dec, None)

    # Calculate the V3 position angle
    v1_pa = calc_position_angle(wcs, v3wcs)

    # Convert to degrees
    wcs = WCSRef(
        ra=wcs.ra * R2D,
        dec=wcs.dec * R2D,
        pa=v1_pa * R2D
    )

    logger.debug('wcs: %s', wcs)
    return wcs


def calc_aperture_wcs(m_eci2siaf):
    """Calculate the aperture WCS

    Parameters
    ----------
    m_eci2siaf : np.array((3, 3))
        The ECI to SIAF transformation matrix

    Returns
    -------
    wcsinfo : WCSRef
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
         np.sqrt(1. - siaf_x * siaf_x - siaf_y * siaf_y)]
    )
    msky = np.dot(m_eci2siaf.transpose(), refpos)
    wcs_ra, wcs_dec = vector_to_angle(msky)

    # Calculate the position angle
    vysiaf = np.array([0., 1., 0.])
    myeci = np.dot(m_eci2siaf.transpose(), vysiaf)

    # The Y axis of the aperture is given by
    vy_ra, vy_dec = vector_to_angle(myeci)

    # The VyPA @ xref,yref is given by
    y = cos(vy_dec) * sin(vy_ra - wcs_ra)
    x = sin(vy_dec) * cos(wcs_dec) - \
        cos(vy_dec) * sin(wcs_dec) * cos((vy_ra - wcs_ra))
    wcs_pa = np.arctan2(y, x)
    if wcs_pa < 0.:
        wcs_pa += PI2

    # Convert all WCS to degrees
    wcsinfo = WCSRef(
        ra=wcs_ra * R2D,
        dec=wcs_dec * R2D,
        pa=wcs_pa * R2D
    )

    logger.debug('wcsinfo: %s', wcsinfo)
    return wcsinfo


def calc_eci2j_matrix(q):
    """Calculate ECI to J-frame matrix from quaternions

    Parameters
    ----------
    q : np.array(q1, q2, q3, q4)
        Array of quaternions from the engineering database

    Returns
    -------
    transform : np.array((3, 3))
        The transform matrix representing the transformation
        from observatory orientation to J-Frame
    """
    q1, q2, q3, q4 = q
    transform = np.array(
        [[1. - 2. * q2 * q2 - 2. * q3 * q3,
          2. * (q1 * q2 + q3 * q4),
          2. * (q3 * q1 - q2 * q4)],
         [2. * (q1 * q2 - q3 * q4),
          1. - 2. * q3 * q3 - 2. * q1 * q1,
          2. * (q2 * q3 + q1 * q4)],
         [2. * (q3 * q1 + q2 * q4),
          2. * (q2 * q3 - q1 * q4),
          1. - 2. * q1 * q1 - 2. * q2 * q2]]
    )

    logger.debug('quaternion: %s', transform)
    return transform


def calc_j2fgs1_matrix(j2fgs_matrix, transpose=True):
    """Calculate the J-frame to FGS1 transformation

    Parameters
    ----------
    j2fgs_matrix : n.array((9,))
        Matrix parameters from the engineering database.
        If all zeros, a predefined matrix is used.

    transpose : bool
        Transpose the resulting matrix.

    Returns
    -------
    transform : np.array((3, 3))
        The transformation matrix

    Notes
    -----
    The parameter `transpose` is defaulted to `True` because the
    matrix, as defined in the engineering telemetry, is actually for
    FGS1-to-J-frame. However, all documentation has always
    referred to this J-to-FGS1.
    """
    if np.isclose(j2fgs_matrix, 0.).all():
        logger.warning('J-Frame to FGS1 engineering parameters are all zero.')
        logger.warning('Using default matrix')
        transform = J2FGS_MATRIX_DEFAULT

    else:
        logger.info(
            'Using J-Frame to FGS1 engineering parameters'
            ' for the J-Frame to FGS1 transformation.'
        )
        transform = np.array(j2fgs_matrix).reshape((3, 3))

    if transpose:
        logger.info('Transposing the J-Frame to FGS matrix.')
        transform = transform.transpose()

    logger.debug('j2fgs1: %s', transform)
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
    transform : np.array((3, 3))
        The transformation matrix
    """
    version = fsmcorr_version.lower()
    units = fsmcorr_units.lower()
    logger.debug('Using version %s', version)
    logger.debug('Using units %s', units)

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
                [1., x / 22.01, y / 21.68],
                [-x / 22.01, 1., 0.],
                [-y / 21.68, 0., 1.]
            ]
        )

    # Default or `V2`: Direct spherical calculation
    # Note: With the "0.0" in the lower middle Y transform
    else:
        if version not in ('latest', 'v2'):
            logger.warning(
                'Unknown version "%s" specified. Using the latest (spherical) calculation.', version
            )
        m_x_partial = np.array(
            [
                [1., 0., 0.],
                [0., cos(x), sin(x)],
                [0., -sin(x), cos(x)]
            ]
        )
        m_y_partial = np.array(
            [
                [cos(y), 0., -sin(y)],
                [0., 1., 0.],
                [sin(y), 0., cos(y)]
            ]
        )
        transform = np.dot(m_x_partial, m_y_partial)

    logger.debug('fsm_delta_matrix: %s', transform)
    return transform


def calc_fgs1_to_sifov_matrix(siaf_db=None, useafter=None):
    """
    Calculate the FGS1 to SI-FOV matrix

    Based on algorithm determined in JSOCINT-555 for a more
    accurate determination of the matrix using the SIAF.
    """
    if siaf_db is None:
        logger.warning('No SIAF database specified. Using the default FGS1_to_SIFOV matrix')
        return FGS12SIFOV_DEFAULT

    try:
        fgs1_siaf = siaf_db.get_wcs('FGS1_FULL_OSS', useafter)
    except (KeyError, OSError, TypeError, RuntimeError):
        logger.warning('Cannot read a SIAF database. Using the default FGS1_to_SIFOV matrix')
        return FGS12SIFOV_DEFAULT

    v2 = fgs1_siaf.v2_ref * A2R
    v3 = (fgs1_siaf.v3_ref + 7.8 * 60.0) * A2R
    y = fgs1_siaf.v3yangle * D2R

    m = np.array([
        [cos(v2) * cos(y) + sin(v2) * sin(v3) * sin(y), cos(v2) * sin(y) - sin(v2) * sin(v3) * cos(y), sin(v2) * cos(v3)],
        [-cos(v3) * sin(y), cos(v3) * cos(y), sin(v3)],
        [-sin(v2) * cos(y) + cos(v2) * sin(v3) * sin(y), -sin(v2) * sin(y) - cos(v2) * sin(v3) * cos(y), cos(v2) * cos(v3)]
    ])

    logger.debug('FGS1_to_SIFOV from siaf is %s', m)
    return m


def calc_sifov2v_matrix():
    """Calculate the SI-FOV to V-Frame matrix

    This is currently defined as the inverse Euler rotation
    about an angle of 7.8 arcmin. Here returns the pre-calculate
    matrix.
    """
    return SIFOV2V_DEFAULT


def calc_v2siaf_matrix(siaf):
    """Calculate the SIAF transformation matrix

    This implements equation 10 from Technical Report JWST-STScI-003222, SM-12. 2021-07
    From Section 3.1:

    The V to SIAF parameters V3IdlYang, V2Ref, V3Ref, and VIdlParity are
    defined and their usage explained in SIAF2017. The parameter values for
    each aperture are specified in the Project Reference Database (PRD).

    Parameters
    ----------
    siaf : SIAF
        The SIAF parameters, where angles are in arcseconds/degrees

    Returns
    -------
    transform : np.array((3, 3))
        The V1 to SIAF transformation matrix

    """
    v2, v3, v3idlyang, vparity = (siaf.v2_ref, siaf.v3_ref,
                                  siaf.v3yangle, siaf.vparity)
    v2 *= A2R
    v3 *= A2R
    v3idlyang *= D2R

    mat = np.array(
        [[cos(v3) * cos(v2),
          cos(v3) * sin(v2),
          sin(v3)],
         [-cos(v3idlyang) * sin(v2) + sin(v3idlyang) * sin(v3) * cos(v2),
          cos(v3idlyang) * cos(v2) + sin(v3idlyang) * sin(v3) * sin(v2),
          -sin(v3idlyang) * cos(v3)],
         [-sin(v3idlyang) * sin(v2) - cos(v3idlyang) * sin(v3) * cos(v2),
          sin(v3idlyang) * cos(v2) - cos(v3idlyang) * sin(v3) * sin(v2),
          cos(v3idlyang) * cos(v3)]])
    pmat = np.array([[0., vparity, 0.],
                     [0., 0., 1.],
                     [1., 0., 0.]])

    transform = np.dot(pmat, mat)

    logger.debug('transform: %s', transform)
    return transform


def calc_position_angle(target, point):
    """Calculate point position angle @target

    Parameters
    ----------
    target : WCSRef
        The TARGET wcs parameters

    point : WCSRef
        The POINT wcs parameters

    Returns
    -------
    point_pa : float
      The POINT position angle, in radians
    """
    y = cos(point.dec) * sin(point.ra - target.ra)
    x = sin(point.dec) * cos(target.dec) - \
        cos(point.dec) * sin(target.dec) * cos((point.ra - target.ra))
    point_pa = np.arctan2(y, x)
    if point_pa < 0:
        point_pa += PI2
    if point_pa >= PI2:
        point_pa -= PI2

    logger.debug('Given target: %s, point: %s, then PA: %s', target, point, point_pa)
    return point_pa


def get_pointing(obsstart, obsend, engdb_url=None,
                 tolerance=60, reduce_func=None):
    """
    Get telescope pointing engineering data.

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

    Returns
    -------
    pointing : Pointing or [Pointing(, ...)]
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

    logger.info('Determining pointing between observations times (mjd):')
    logger.info('obsstart: %s obsend: %s', obsstart, obsend)
    logger.info('Telemetry search tolerance: %s', tolerance)
    logger.info('Reduction function: %s', reduce_func)

    mnemonics = get_mnemonics(obsstart, obsend, tolerance, engdb_url=engdb_url)
    reduced = reduce_func(mnemonics)

    logger.log(DEBUG_FULL, 'Memonics found:')
    logger.log(DEBUG_FULL, '%s', mnemonics)
    logger.info('Reduced set of pointings:')
    logger.info('%s', reduced)

    return reduced


def vector_to_angle(v):
    """Returns tuple of spherical angles from unit direction Vector

    This implements equations 8 & 9 from Technical Report JWST-STScI-003222, SM-12. 2021-07

    Parameters
    ----------
    v : [v0, v1, v2]

    Returns
    -------
    alpha, delta : float, float
        The spherical angles, in radians
    """
    alpha = np.arctan2(v[1], v[0])
    delta = np.arcsin(v[2])
    if alpha < 0.:
        alpha += 2. * np.pi
    return(alpha, delta)


def angle_to_vector(alpha, delta):
    """Convert spherical angles to unit vector

    This implements equation 7 from Technical Report JWST-STScI-003222, SM-12. 2021-07

    Parameters
    ----------
    alpha, delta : float
        Spherical angles in radians

    Returns
    -------
    v : [float, float, float]
        Unit vector
    """
    v0 = sin(alpha) * cos(delta)
    v1 = sin(alpha) * sin(delta)
    v2 = cos(delta)

    return [v0, v1, v2]


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


def get_mnemonics(obsstart, obsend, tolerance, engdb_url=None):
    """Retrieve pointing mnemonics from the engineering database

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    Returns
    -------
    mnemonics : {mnemonic: [value[,...]][,...]}
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
            '\nException: {}'.format(exception)
        )
    logger.info(
        'Querying engineering DB: %s', engdb.base_url
    )

    mnemonics = {
        'SA_ZATTEST1': None,
        'SA_ZATTEST2': None,
        'SA_ZATTEST3': None,
        'SA_ZATTEST4': None,
        'SA_ZRFGS2J11': None,
        'SA_ZRFGS2J12': None,
        'SA_ZRFGS2J13': None,
        'SA_ZRFGS2J21': None,
        'SA_ZRFGS2J22': None,
        'SA_ZRFGS2J23': None,
        'SA_ZRFGS2J31': None,
        'SA_ZRFGS2J32': None,
        'SA_ZRFGS2J33': None,
        'SA_ZADUCMDX': None,
        'SA_ZADUCMDY': None,
        'SA_ZFGGSCMDX': None,
        'SA_ZFGGSCMDY': None,
        'SA_ZFGGSPOSX': None,
        'SA_ZFGGSPOSY': None,
        'SA_ZFGDETID': None
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
            logger.warning('Mnemonic %s has no telemetry within the observation time.', mnemonic)
            logger.warning('Attempting to use bracket values within %s seconds', tolerance)

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
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointings : [Pointing[,...]]
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

        gs_commanded = np.array([
            mnemonics_at_time['SA_ZFGGSCMDX'].value,
            mnemonics_at_time['SA_ZFGGSCMDY'].value

        ])

        gs_position = np.array([
            mnemonics_at_time['SA_ZFGGSPOSX'].value,
            mnemonics_at_time['SA_ZFGGSPOSY'].value

        ])

        fgsid = mnemonics_at_time['SA_ZFGDETID'].value

        pointing = Pointing(q=q, obstime=obstime, j2fgs_matrix=j2fgs_matrix,
                            fsmcorr=fsmcorr, gs_commanded=gs_commanded,
                            fgsid=fgsid, gs_position=gs_position)
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
        logger.info(' setting xref_sci to %s', siaf.crpix1)
        logger.info(' setting yref_sci to %s', siaf.crpix2)
        model.meta.wcsinfo.siaf_xref_sci = siaf.crpix1
        model.meta.wcsinfo.siaf_yref_sci = siaf.crpix2


def first_pointing(mnemonics):
    """Return first pointing

    Parameters
    ==========
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing : Pointing
        First pointing.

    """
    pointings = all_pointings(mnemonics)
    return pointings[0]


def pointing_from_average(mnemonics):
    """Determine single pointing from average of available pointings

    Parameters
    ==========
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing : Pointing
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
        # Weed out mnemonic entries that are zero, though some are OK to be zero.
        if mnemonic not in ['SA_ZADUCMDX', 'SA_ZADUCMDY', 'SA_ZFGGSCMDX', 'SA_ZFGGSCMDY',
                            'SA_ZFGGSPOSX', 'SA_ZFGGSPOSY']:
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

    gs_commanded = np.array([
        mnemonic_averages['SA_ZFGGSCMDX'],
        mnemonic_averages['SA_ZFGGSCMDY']

    ])

    gs_position = np.array([
        mnemonic_averages['SA_ZFGGSPOSX'],
        mnemonic_averages['SA_ZFGGSPOSY']
    ])

    # For FGS ID, just take the first one.
    fgsid = mnemonics['SA_ZFGDETID'][0].value

    pointing = Pointing(obstime=obstime, q=q, j2fgs_matrix=j2fgs_matrix,
                        fsmcorr=fsmcorr, gs_commanded=gs_commanded,
                        fgsid=fgsid, gs_position=gs_position)
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


def calc_estimated_gs_wcs(t_pars: TransformParameters):
    """Calculate the estimated guide star RA/DEC/Y-angle

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    gs_wcs : WCSRef
        Estimated RA, Dec, and Y-angle. All in degrees.

    Notes
    -----

    This implements equation 15, 16, 17 from Technical Report JWST-STScI-003222, SM-12.
    """

    # Determine the ECI to Guide star transformation
    t = calc_m_eci2gs(t_pars)
    m_eci2gs = t.m_eci2gs

    # Determine the wcs
    wcs = calc_wcs_from_matrix(m_eci2gs)

    logger.debug('wcs: %s', wcs)
    return wcs


def calc_v3pags(t_pars: TransformParameters):
    """Calculate the V3 Position Angle at the Guide Star

    This implements equation 18 from Technical Report JWST-STScI-003222, SM-12. 2021-07

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    v3pags : float
        The V3 Position Angle at the Guide Star, in degrees
    """

    # Determine Guides Star estimated WCS information.
    gs_wcs = calc_estimated_gs_wcs(t_pars)

    # Retrieve the Ideal Y-angle for the desired FGS
    fgs_siaf = t_pars.siaf_db.get_wcs(FGSId2Aper[t_pars.pointing.fgsid], useafter=t_pars.useafter)

    # Calculate V3PAGS
    v3pags = gs_wcs.pa - fgs_siaf.v3yangle

    logger.debug('v3pags: %s', v3pags)
    return v3pags


def calc_m_eci2gs(t_pars: TransformParameters):
    """Calculate the M_eci2gs matrix as per TR presented in 2021-07

    This implements equation 14 from Technical Report JWST-STScI-003222, SM-12.
    From Section 3.2:

    The equation is formed by inverting the equation in Section 5.9.1.2 of
    SE-20 which converts from the attitude specified at the Guide Star into the
    commanded spacecraft J-frame attitude. With the inversion, the commanded
    J-frame attitude quaternion is replaced in this equation by the matrix
    derived from the measured J-frame attitude quaternion.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The calculated transforms. The target transform is
        `transforms.m_eci2gs`. See the notes for other transforms
        used and calculated.

    Notes
    -----
    The transform train needed to calculate M_eci_to_gs is

        M_eci_to_gs =
            M_z_to_x               *
            M_gsics_to_gsappics    *
            M_fgsics_to_gsics      *
            M_fgs1ics_to_M_fgsics  *
            M_j_to_fgs1ics         *
            M_eci_to_j

        where

            M_eci_to_gs = ECI to Guide Star Ideal Frame
            M_gsics_to_gsappics = Velocity Aberration correction
            M_fgsics_to_gsics = Convert from the FGS ICS frame to Guide Star ICS frame
            M_fgs1ics_to_fgsics = Convert from the FGS1 ICS frame to the in-use FGS ICS frame
            M_j_to_fgs1ics = Convert from J frame to FGS1 ICS frame
            M_eci_to_j = ECI (quaternion) to J-frame
    """

    # Initial state of the transforms
    t = Transforms(override=t_pars.override_transforms)

    t.m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)
    t.m_j2fgs1 = calc_j2fgs1_matrix(t_pars.pointing.j2fgs_matrix, t_pars.j2fgs_transpose)
    t.m_fgs12fgsx = calc_m_fgs12fgsx(t_pars.pointing.fgsid, t_pars.siaf_db)
    t.m_fgsx2gs = calc_m_fgsx2gs(t_pars.pointing.gs_commanded)

    # Apply the Velocity Aberration. To do so, the M_eci2gsics matrix must be created. This
    # is used to calculate the aberration matrix.
    # Also, since the aberration is to be removed, the velocity is negated.
    m_eci2gsics = np.linalg.multi_dot([t.m_fgsx2gs, t.m_fgs12fgsx, t.m_j2fgs1, t.m_eci2j])
    logger.debug('m_eci2gsics: %s', m_eci2gsics)
    t.m_gs2gsapp = calc_gs2gsapp(m_eci2gsics, t_pars.jwst_velocity)

    # Put it all together
    t.m_eci2gs = np.linalg.multi_dot([MZ2X, t.m_gs2gsapp, m_eci2gsics])
    logger.debug('m_eci2gs: %s', t.m_eci2gs)

    # That's all folks
    return t


def calc_m_fgs12fgsx(fgsid, siaf_db):
    """Calculate the FGS1 to FGSx matrix

    This implements equation 23 from Technical Report JWST-STScI-003222, SM-12, 2021-07
    From Section 3.2.3:

    A selected guide star being used, could be in FGS 1 or FGS 2. The JWST ACS
    always uses the FGS 1 ICS frame to calculate the commanded spacecraft
    J-frame attitude and is used in the attitude control loop. If the specified
    guide star is in FGS 2, its position will be converted to the FGS 1 ICS
    using an on board FGS2 to FGS1 k-constant matrix. Here we are creating the
    FGS1 ICS to FGSj ICS DCM which converts from the FGS1 ICS frame to the FGSj
    ICS frame using SIAF parameters for the FGSs.

    Parameters
    ----------
    fgsid : int
        The id of the FGS in actually in use. 1 or 2

    siaf_db : SiafDb
        The SIAF database.

    Returns
    -------
    m_fgs12fgsx : numpy.array(3, 3)
        The DCM to transform from FGS1 ICS frame to the desired FGS frame

    """
    # If the in-use FGS is FGS1, no transformation is necessary.
    # Simply return the identity matrix.
    if fgsid == 1:
        m_fgs12fgsx = np.identity(3)
        logger.debug('FGS1 is in use, the identity matrix is returned: %s', m_fgs12fgsx)
        return m_fgs12fgsx

    if fgsid != 2:
        raise ValueError(f'fgsid == {fgsid} is invalid. Must be 1 or 2')

    # FGS2 is in use. Calculate the transform from FGS1 to FGS2
    fgs1_siaf = siaf_db.get_wcs(FGSId2Aper[1])
    fgs2_siaf = siaf_db.get_wcs(FGSId2Aper[2])
    m_fgs1 = calc_v2siaf_matrix(fgs1_siaf)
    m_fgs2 = calc_v2siaf_matrix(fgs2_siaf)

    m_fgs12fgsx = np.dot(m_fgs2, m_fgs1.transpose())

    logger.debug('m_fgs12fgsx: %s', m_fgs12fgsx)
    return m_fgs12fgsx


def calc_m_fgsx2gs(gs_commanded):
    """Calculate the FGS1 to commanded Guide Star frame

    This implements equation 25 from Technical Report JWST-STScI-003222, SM-12, 2021-07

    Parameters
    ----------
    gs_commanded : numpy.array(2)
        The Guide Star commanded position, in arcseconds

    Returns
    -------
    m_fgsx2gs : numpy.array(3, 3)
        The DCM transform from FGSx (1 or 2) to Guide Star ICS frame
    """
    m_gs2fgsx = calc_m_gs2fgsx(gs_commanded)
    m_fgsx2gs = m_gs2fgsx.transpose()

    logger.debug('m_fgsx2gs: %s', m_fgsx2gs)
    return m_fgsx2gs


def calc_m_gs2fgsx(gs_commanded):
    """Calculate the Guides Star frame to FGSx ICS frame

    This implements equation 26 from Technical Report JWST-STScI-003222, SM-12, 2021-07

    Parameters
    ----------
    gs_commanded : numpy.array(2)
        The commanded position of the guide stars, in arcseconds

    Returns
    -------
    m_gs2fgsx : numpy.array(3, 3)
        The guide star to FGSx transformation
    """
    in_rads = gs_commanded * A2R
    x, y = in_rads
    m_x = np.array([
        [cos(-x), 0., -sin(-x)],
        [0., 1., 0.],
        [sin(-x), 0., cos(-x)]
    ])
    m_y = np.array([
        [1.0, 0.0, 0.0],
        [0., cos(y), sin(y)],
        [0., -sin(y), cos(y)]
    ])
    m_gs2fgsx = np.dot(m_y, m_x)

    logger.debug('m_gs2fgsx: %s', m_gs2fgsx)
    return m_gs2fgsx


def trans_fgs2v(fgsid, ideal, siaf_db):
    """Transform an Ideal coordinate to V coordinates

    Parameters
    ----------
    fgsid : int
        The FGS in use.

    ideal : numpy.array(2)
        The Ideal coordinates in arcseconds

    siaf_db : SiafDb
        The SIAF database.

    Returns
    -------
    v : numpy.array(2)
        The V-frame coordinates in arcseconds
    """
    ideal_rads = ideal * A2R
    ideal_vec = cart_to_vector(ideal_rads)
    siaf = siaf_db.get_wcs(FGSId2Aper[fgsid])
    m_v2fgs = calc_v2siaf_matrix(siaf)
    v_vec = np.dot(m_v2fgs.transpose(), ideal_vec)
    v_rads = np.array(vector_to_angle(v_vec))
    v = v_rads * R2A

    logger.debug('FGS%s %s -> V %s', fgsid, ideal, v)
    return v


def cart_to_vector(coord):
    """Convert Cartesian to a unit vector

    This implements equation 4 from Technical Report JWST-STScI-003222, SM-12. 2021-07

    Parameters
    ----------
    coord : numpy.array(2)
        The Cartesian coordinate.

    Returns
    -------
    vector : numpy.array(3)
        The vector version
    """
    vector = np.array([
        coord[0],
        coord[1],
        sqrt(1 - coord[0]**2 - coord[1]**2)
    ])

    return vector
