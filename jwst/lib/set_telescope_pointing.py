"""Set Telescope Pointing from Observatory Engineering Telemetry

Calculate and update the pointing-related and world coordinate system-related
keywords. Given a time period, usually defined by an exposure, the engineering
mnemonic database is queried for observatory orientation. The orientation
defines the sky coordinates a particular point on the observatory is pointed to.
Then, using a set of matrix transformations, the sky coordinates of the
reference pixel of a desired aperture is calculated.

The transformations are defined by the Technical Reference JWST-STScI-003222,
SM-12. This document has undergone a number of revisions. The current version
implemented is based on an internal email version Rev. C, produced 2021-11.

There are a number of algorithms, or *methods*, that have been implemented.
Most represent the historical refinement of the algorithm. Until the technical
reference is finalized, all methods will remain in the code. The default,
state-of-the art algorithm is represented by method ``OPS_TR_202111``,
implemented by
`~jwst.lib.set_telescope_pointing.calc_transforms_ops_tr_202111`.

Interface
=========

The primary usage is through the command line interface
`set_telescope_pointing.py`. Operating on a list of JWST Level 1b exposures,
this command updates the world coordinate system keywords with the values
necessary to translate from aperture pixel to sky coordinates.

Access to the JWST Engineering Mnemonic database is required. See the
:ref:`Engineering Database Interface<engdb>` for more information.

Programmatically, the command line is implemented by the function
`~jwst.lib.set_telescope_pointing.add_wcs`, which calls the basic function
`~jwst.lib.set_telescope_pointing.calc_wcs`. The available methods are defined
by `~jwst.lib.set_telescope_pointing.Methods`.

There are two data structures used to maintain the state of the transformation.
`~jwst.lib.set_telescope_pointing.TransformParameters` contains the parameters
needed to perform the transformations.
`~jwst.lib.set_telescope_pointing.Transforms` contains the calculated
transformation matrices.

Transformation Matrices
-----------------------

All the transformation matrices, as defined by
`~jwst.lib.set_telescope_pointing.Transforms`, are Direction Cosine Matrices
(DCM). A DCM contains the Euler rotation angles that represent the sky
coordinates for a particular frame-of-reference. The initial DCM is provided
through the engineering telemetry and represents where in the sky either the
Fine Guidance Sensor (FGS) or star tracker is pointed to. Then, through a set
of transformations, the DCM for the reference point of the target aperture
is calculated.

"""
import functools
import sys

import asdf
from collections import defaultdict, namedtuple
from copy import copy
import dataclasses
from enum import Enum
import logging
from math import (cos, sin, sqrt)
import typing

from astropy import units as U
from astropy.table import Table
from astropy.time import Time
import numpy as np
from scipy.interpolate import interp1d

from stdatamodels.jwst import datamodels

from .exposure_types import IMAGING_TYPES, FGS_GUIDE_EXP_TYPES
from .set_velocity_aberration import compute_va_effects_vector
from .siafdb import SIAF, SiafDb
from ..assign_wcs.util import update_s_region_keyword, calc_rotation_matrix
from ..assign_wcs.pointing import v23tosky
from ..lib.engdb_tools import ENGDB_Service
from ..lib.pipe_utils import is_tso

__all__ = [
    'Methods',
    'TransformParameters',
    'Transforms',
    'WCSRef',
    'add_wcs',
    'calc_transforms',
    'calc_transforms_ops_tr_202111',
    'calc_wcs',
    'calc_wcs_over_time',
    'update_wcs',
]

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]

# Datamodels that can be updated, normally
EXPECTED_MODELS = (datamodels.Level1bModel, datamodels.ImageModel, datamodels.CubeModel)

# Exposure types that can be updated, normally
TYPES_TO_UPDATE = set(list(IMAGING_TYPES) + FGS_GUIDE_EXP_TYPES)

# Mnemonics for each transformation method.
# dict where value indicates whether the mnemonic is required or not.
COURSE_TR_202111_MNEMONICS = {
    'SA_ZATTEST1': True,
    'SA_ZATTEST2': True,
    'SA_ZATTEST3': True,
    'SA_ZATTEST4': True,
    'SA_ZRFGS2J11': True,
    'SA_ZRFGS2J12': True,
    'SA_ZRFGS2J13': True,
    'SA_ZRFGS2J21': True,
    'SA_ZRFGS2J22': True,
    'SA_ZRFGS2J23': True,
    'SA_ZRFGS2J31': True,
    'SA_ZRFGS2J32': True,
    'SA_ZRFGS2J33': True,
    'SA_ZADUCMDX': False,
    'SA_ZADUCMDY': False,
    'SA_ZFGGSCMDX': False,
    'SA_ZFGGSCMDY': False,
    'SA_ZFGDETID': False,
}

TRACK_TR_202111_MNEMONICS = {
    **COURSE_TR_202111_MNEMONICS,
    'SA_ZFGGSPOSX': False,
    'SA_ZFGGSPOSY': False,
}

FGS_ACQ_EXP_TYPES = ['fgs_acq1', 'fgs_acq2']
FGS_ACQ_MNEMONICS = {
    'IFGS_ACQ_DETXCOR': True,
    'IFGS_ACQ_DETYCOR': True,
    'IFGS_ACQ_DETXSIZ': True,
    'IFGS_ACQ_DETYSIZ': True,
    'IFGS_ACQ_XPOSG': True,
    'IFGS_ACQ_YPOSG': True,
}

FGS_GUIDED_EXP_TYPES = ['fgs_fineguide', 'fgs_track']
FGS_GUIDED_MNEMONICS = {
    'IFGS_TFGGS_X',
    'IFGS_TFGGS_Y',
    'IFGS_TFGDET_XCOR',
    'IFGS_TFGDET_YCOR',
    'IFGS_TFGDET_XSIZ',
    'IFGS_TFGDET_YSIZ',
}

FGS_ID_EXP_TYPES = ['fgs_id-image', 'fgs_id-stack']
FGS_ID_MNEMONICS = {
    'IFGS_ID_XPOSG',
    'IFGS_ID_YPOSG',
    'IFGS_ID_DETXCOR',
    'IFGS_ID_DETYCOR',
    'IFGS_ID_DETXSIZ',
    'IFGS_ID_DETYSIZ',
}

# FGS ACQ1/ACQ2 modes, a dedicated range of mmenonics need to be present.
# These define those ranges. Key is the ACQ exposure type.
FGS_ACQ_MINVALUES = {
    'fgs_acq1': 1,
    'fgs_acq2': 4
}
FGS_ACQ_SLICES = {
    'fgs_acq1': slice(0, 3),
    'fgs_acq2': slice(3, 8)
}
FGS_ACQ_WINDOW_INDEX = {
    'fgs_acq1': 0,
    'fgs_acq2': -1
}


# The available methods for transformation
class Methods(Enum):
    """Available methods to calculate V1 and aperture WCS information

    Current state-of-art is OPS_TR_202111. This method chooses either COARSE_TR_202111 or
    TRACK_TR_202111 depending on the guidance mode, as specified by header keyword PCS_MODE.
    """
    #: COARSE tracking mode algorithm, TR version 2021-11.
    COARSE_TR_202111 = ('coarse_tr_202111', 'calc_transforms_coarse_tr_202111', 'calc_wcs_tr_202111', COURSE_TR_202111_MNEMONICS)
    #: Method to use in OPS to use TR version 2021-11
    OPS_TR_202111 = ('ops_tr_202111', 'calc_transforms_ops_tr_202111', 'calc_wcs_tr_202111', TRACK_TR_202111_MNEMONICS)
    #: TRACK and FINEGUIDE mode algorithm, TR version 2021-11
    TRACK_TR_202111 = ('track_tr_202111', 'calc_transforms_track_tr_202111', 'calc_wcs_tr_202111', TRACK_TR_202111_MNEMONICS)

    # Aliases
    #: Algorithm to use by default. Used by Operations.
    default = OPS_TR_202111
    #: Default algorithm under PCS_MODE COARSE.
    COARSE = COARSE_TR_202111
    #: Default algorithm for use by Operations.
    OPS = OPS_TR_202111
    #: Default algorithm under PCS_MODE TRACK/FINEGUIDE/MOVING.
    TRACK = TRACK_TR_202111

    def __new__(cls: object, value: str, func_name: str, calc_func: str, mnemonics: dict):
        obj = object.__new__(cls)
        obj._value_ = value
        obj._func_name = func_name
        obj._calc_func = calc_func
        obj._mnemonics = mnemonics
        return obj

    @property
    def calc_func(self):
        """Function associated with the method"""
        return globals()[self._calc_func]

    @property
    def func(self):
        """Function associated with the method"""
        return globals()[self._func_name]

    @property
    def mnemonics(self):
        return self._mnemonics

    def __str__(self):
        return self.value


# FGS id to aperture name
FGSId2Aper = {1: 'FGS1_FULL_OSS', 2: 'FGS2_FULL_OSS'}

# FGS Ids
FGSIDS = [1, 2]

# Definition of th J3 Ideal Y-Angle
J3IDLYANGLE = -1.25  # Degrees

# Conversion from seconds to MJD
SECONDS2MJD = 1 / 24 / 60 / 60

# Default transformation matrices
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

# Define the transformation matrices to move between the Idealized Coordinate System (ICS)
# and the Idealized Coordinate System (Idl). ICS is the spacecraft-centric system used by
# all frames up through the V-frame. Idl is used by the instruments.
# Reference: Eqs. 1 & 2 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
M_idl2ics = MX2Z = np.array(
    [[0, 1, 0],
     [0, 0, 1],
     [1, 0, 0]]
)
M_ics2idl = MZ2X = np.array(
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
# Attributes are as follows. Except for the observation time, all values
# are retrieved from the engineering data.
#    q            : Quaternion of the FGS.
#    j2fgs_matrix : J-frame to FGS transformation.
#    fsmcorr      : Fine Steering Mirror position.
#    obstime      : Mid-point time of the observation at which all other values have been calculated.
#    gs_commanded : Guide star position as originally commanded.
#    fgsid        : FGS in use, 1 or 2.
#    gs_position  : X/Y guide star position in the FGS.
Pointing = namedtuple('Pointing', ['q', 'j2fgs_matrix', 'fsmcorr', 'obstime', 'gs_commanded', 'fgsid', 'gs_position'])
Pointing.__new__.__defaults__ = ((None,) * 5)

# Guide Star ACQ pointing container
# Attributes are as follows. All values are retrieved from the engineering.
#    position : X/Y position of the guide star within the acquisition window of the FGS.
#    corner   : X/Y corner of the acquisition window within the FGS.
#    size     : X/Y size of the acquisition window.
GuideStarPosition = namedtuple('GuideStarPosition', ['position', 'corner', 'size'])
GuideStarPosition.__new__.__defaults__ = ((None,) * 3)


# Transforms
@dataclasses.dataclass
class Transforms:
    """The matrices used in calculation of the M_eci2siaf transformation
    """
    #: ECI to FGS1
    m_eci2fgs1: np.array = None
    #: ECI to Guide Star
    m_eci2gs: np.array = None
    #: ECI to J-Frame
    m_eci2j: np.array = None
    #: ECI to SIAF
    m_eci2siaf: np.array = None
    #: ECI to SIFOV
    m_eci2sifov: np.array = None
    #: ECI to V
    m_eci2v: np.array = None
    #: FGSX to Guide Stars transformation
    m_fgsx2gs: np.array = None
    #: FGS1 to SIFOV
    m_fgs12sifov: np.array = None
    #: Velocity aberration
    m_gs2gsapp: np.array = None
    #: J-Frame to FGS1
    m_j2fgs1: np.array = None
    #: FSM correction
    m_sifov_fsm_delta: np.array = None
    #: SIFOV to V1
    m_sifov2v: np.array = None
    #: V to SIAF
    m_v2siaf: np.array = None
    #: Override values. Either another Transforms or dict-like object
    override: object = None

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
    """
    #: If telemetry cannot be determined, use existing information in the observation's header.
    allow_default: bool = False
    #: The V3 position angle to use if the pointing information is not found.
    default_pa_v3: float = 0.
    #: Detector in use.
    detector: str = None
    #: Do not write out the modified file.
    dry_run: bool = False
    #: URL of the engineering telemetry database REST interface.
    engdb_url: str = None
    #: Exposure type
    exp_type: str = None
    #: FGS to use as the guiding FGS. If None, will be set to what telemetry provides.
    fgsid: int = None
    #: The version of the FSM correction calculation to use. See `calc_sifov_fsm_delta_matrix`
    fsmcorr_version: str = 'latest'
    #: Units of the FSM correction values. Default is 'arcsec'. See `calc_sifov_fsm_delta_matrix`
    fsmcorr_units: str = 'arcsec'
    #: Guide star WCS info, typically from the input model.
    guide_star_wcs: WCSRef = WCSRef()
    #: Transpose the `j2fgs1` matrix.
    j2fgs_transpose: bool = True
    #: The [DX, DY, DZ] barycentri velocity vector
    jwst_velocity: np.array = None
    #: The method, or algorithm, to use in calculating the transform. If not specified, the default method is used.
    method: Methods = Methods.default
    #: Observation end time
    obsend: float = None
    #: Observation start time
    obsstart: float = None
    #: If set, matrices that should be used instead of the calculated one.
    override_transforms: Transforms = None
    #: The tracking mode in use.
    pcs_mode: str = None
    #: The observatory orientation, represented by the ECI quaternion, and other engineering mnemonics
    pointing: Pointing = None
    #: Reduction function to use on values.
    reduce_func: typing.Callable = None
    #: The SIAF information for the input model
    siaf: SIAF = None
    #: The SIAF database
    siaf_db: SiafDb = None
    #: If no telemetry can be found during the observation,
    #: the time, in seconds, beyond the observation time to search for telemetry.
    tolerance: float = 60.
    #: The date of observation (`jwst.datamodels.JwstDataModel.meta.date`)
    useafter: str = None
    #: V3 position angle at Guide Star (`jwst.datamodels.JwstDataModel.meta.guide_star.gs_v3_pa_science`)
    v3pa_at_gs: float = None

    def as_reprdict(self):
        """Return a dict where all values are REPR of their values"""
        r = dict((field.name, repr(getattr(self, field.name))) for field in dataclasses.fields(self))
        return r

    def update_pointing(self):
        """Update pointing information"""
        self.pointing = get_pointing(self.obsstart, self.obsend,
                                     mnemonics_to_read=self.method.mnemonics,
                                     engdb_url=self.engdb_url,
                                     tolerance=self.tolerance, reduce_func=self.reduce_func)


def add_wcs(filename, allow_any_file=False, force_level1bmodel=False,
            default_pa_v3=0., siaf_path=None, prd=None, engdb_url=None,
            fgsid=None, tolerance=60, allow_default=False, reduce_func=None,
            dry_run=False, save_transforms=None, **transform_kwargs):
    """Add WCS information to a JWST DataModel.

    Telescope orientation is attempted to be obtained from
    the engineering database. Failing that, a default pointing
    is used based on proposal target.

    The file is updated in-place.

    Parameters
    ----------
    filename : str
        The path to a data file.

    allow_any_file : bool
        Attempt to add the WCS information to any type of file.
        The default, `False`, only allows modifications of files that contain
        known datamodels of `Level1bmodel`, `ImageModel`, or `CubeModel`.

    force_level1bmodel : bool
        If not `allow_any_file`, and the input file model is unknown,
        open the input file as a Level1bModel regardless.

    default_pa_v3 : float
        The V3 position angle to use if the pointing information
        is not found.

    siaf_path : str or file-like object or None
        The path to the SIAF database. See `SiafDb` for more information.

    prd : str
        The PRD version from the `pysiaf` to use.
        `siaf_path` overrides this value.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    fgsid : int or None
        When in COARSE mode, the FGS to use as the guider reference.
        If None, use what is provided in telemetry.

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

    This function adds absolute pointing information to the JWST
    datamodels provided. By default, only Stage 1 and Stage 2a exposures are
    allowed to be updated. These have the suffixes of "uncal", "rate", and
    "rateints" representing datamodels Level1bModel, ImageModel, and CubeModel.
    Any higher level product, from Stage 2b and beyond, that has had the
    `assign_wcs` step applied, have improved WCS information. Running
    this task on such files will potentially corrupt the WCS.

    It starts by populating the headers with values from the SIAF database.
    It adds the following keywords to all files:

    V2_REF (arcseconds)
    V3_REF (arcseconds)
    VPARITY (+1 or -1)
    V3I_YANG (decimal degrees)

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
    try:
        model = datamodels.open(filename, guess=allow_any_file)
    except TypeError:
        if force_level1bmodel:
            logger.warning(f'Input {filename} is an unknown model, opening as a Level1bModel.')
            model = datamodels.Level1bModel(filename)
        else:
            raise

    try:
        if type(model) not in EXPECTED_MODELS:
            logger.warning(f'Input {model} is not of an expected type (uncal, rate, rateints)'
                           '\n    Updating pointing may have no effect or detrimental effects on the WCS information,'
                           '\n    especially if the input is the result of Level2b or higher calibration.')
            if not allow_any_file:
                raise TypeError(f'Input model {model} is not one of {EXPECTED_MODELS} and `allow_any_file` is `False`.'
                                '\n\tFailing WCS processing.')

        t_pars, transforms = update_wcs(
            model,
            default_pa_v3=default_pa_v3,
            siaf_path=siaf_path,
            prd=prd,
            engdb_url=engdb_url,
            fgsid=fgsid,
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
    finally:
        model.close()

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


def update_wcs(model, default_pa_v3=0., default_roll_ref=0., siaf_path=None, prd=None, engdb_url=None,
               fgsid=None, tolerance=60, allow_default=False,
               reduce_func=None, **transform_kwargs):
    """Update WCS pointing information

    Given a `jwst.datamodels.JwstDataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The model to update.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the roll ref angle.

    siaf_path : str or Path-like object
        The path to the SIAF database. See `SiafDb` for more information.

    prd : str
        The PRD version from the `pysiaf` to use.
        `siaf_path` overrides this value.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    fgsid : int or None
        When in COARSE mode, the FGS to use as the guider reference.
        If None, use what is provided in telemetry.

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
        performed. In particular, FGS GUIDER data does
        not need `transforms`.
    """
    t_pars = transforms = None  # Assume telemetry is not used.

    if not prd:
        prd = model.meta.prd_software_version
    siaf_db = SiafDb(source=siaf_path, prd=prd)

    # Get model attributes
    useafter = model.meta.observation.date

    # Configure transformation parameters.
    t_pars = t_pars_from_model(
        model,
        default_pa_v3=default_pa_v3, engdb_url=engdb_url,
        tolerance=tolerance, allow_default=allow_default,
        reduce_func=reduce_func, siaf_db=siaf_db, useafter=useafter,
        **transform_kwargs
    )
    if fgsid:
        t_pars.fgsid = fgsid

    # Populate header with SIAF information.
    if t_pars.siaf is None:
        if t_pars.exp_type not in FGS_GUIDE_EXP_TYPES:
            raise ValueError('Insufficient SIAF information found in header.')
    else:
        populate_model_from_siaf(model, t_pars.siaf)

    # Calculate WCS.
    if t_pars.exp_type in FGS_GUIDE_EXP_TYPES:
        update_wcs_from_fgs_guiding(
            model, t_pars, default_roll_ref=default_roll_ref
        )
        transforms = None
    else:
        transforms = update_wcs_from_telem(model, t_pars)

    return t_pars, transforms


def update_wcs_from_fgs_guiding(model, t_pars, default_roll_ref=0.0, default_vparity=1, default_v3yangle=0.0):
    """ Update WCS pointing from header information

    For Fine Guidance guiding observations, nearly everything
    in the `wcsinfo` meta information is already populated,
    except for the PC matrix and CRVAL*. This function updates the PC
    matrix based on the rest of the `wcsinfo`.

    CRVAL* values are taken from GS_RA/GS_DEC.

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The model to update.

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.

    default_vparity : int
        The default `VIdlParity` to use and should
        be either "1" or "-1". "1" is the
        default since FGS guiding will be using the
        OSS aperture.

    default_v3yangle : float
        Default SIAF Y-angle.
    """

    logger.info('Updating WCS for Fine Guidance.')

    crpix1, crpix2, crval1, crval2, pc_matrix = calc_wcs_guiding(
        model, t_pars, default_roll_ref, default_vparity, default_v3yangle
    )

    logger.info('WCS info:'
                f'\n\tcrpix1: {crpix1} crpix2: {crpix2}'
                f'\n\tcrval1: {crval1} crval2: {crval2}'
                f'\n\tpc_matrix: {pc_matrix}')

    model.meta.wcsinfo.crpix1 = crpix1
    model.meta.wcsinfo.crpix2 = crpix2
    model.meta.wcsinfo.crval1 = crval1
    model.meta.wcsinfo.crval2 = crval2
    (
        model.meta.wcsinfo.pc1_1,
        model.meta.wcsinfo.pc1_2,
        model.meta.wcsinfo.pc2_1,
        model.meta.wcsinfo.pc2_2
    ) = pc_matrix


def update_wcs_from_telem(model, t_pars: TransformParameters):
    """Update WCS pointing information

    Given a `jwst.datamodels.JwstDataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
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

    # Setup default WCS info if actual pointing and calculations fail.
    wcsinfo = WCSRef(
        model.meta.target.ra,
        model.meta.target.dec,
        t_pars.default_pa_v3
    )
    vinfo = wcsinfo

    # Get the pointing information
    try:
        t_pars.update_pointing()
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
            t_pars.pointing = None
    else:
        logger.info('Successful read of engineering quaternions:')
        logger.info('\tPointing: %s', t_pars.pointing)

    # If pointing is available, attempt to calculate WCS information
    if t_pars.pointing is not None:
        try:
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
    model.meta.wcsinfo.roll_ref = pa_to_roll_ref(wcsinfo.pa, t_pars.siaf)
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

    # If TARG_RA/TARG_DEC still 0/0 (e.g. pure parallels with no defined target),
    # populate with RA_REF/DEC_REF values
    if (model.meta.target.ra == 0.0 and model.meta.target.dec == 0.0) and (
            'PARALLEL' in model.meta.visit.type):

        logger.warning('No target location specified for parallel observation:'
                       'copying reference point RA/Dec to TARG_RA/TARG_DEC.')
        model.meta.target.ra = model.meta.wcsinfo.ra_ref
        model.meta.target.dec = model.meta.wcsinfo.dec_ref

    return transforms


def update_s_region(model, siaf):
    """Update ``S_REGION`` sky footprint information.

    The ``S_REGION`` keyword is intended to store the spatial footprint of
    an observation using the VO standard STCS representation.

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The model to update in-place.
    siaf : namedtuple
        The ``SIAF`` tuple with values populated from the PRD database.
    """
    vertices = siaf.vertices_idl
    xvert = vertices[:4]
    yvert = vertices[4:]
    logger.info("Vertices for aperture %s: %s", model.meta.aperture.name, vertices)

    # Execute IdealToV2V3, followed by V23ToSky
    from stdatamodels.jwst.transforms.models import IdealToV2V3
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
    obstimes, wcsinfos, vinfos : [astropy.time.Time[,...]], [WCSRef[,...]], [WCSRef[,...]]
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
    """Given observatory orientation and target aperture, calculate V1 and Reference Pixel sky coordinates

    Parameters
    ----------
    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    wcsinfo, vinfo, transforms : WCSRef, WCSRef, Transforms
        A 3-tuple is returned with the WCS pointing for
        the aperture and the V1 axis, and the transformation matrices.
    """
    if t_pars.siaf is None:
        t_pars.siaf = SIAF()

    # Calculate transforms
    transforms = calc_transforms(t_pars)

    # Calculate the wcs information
    wcsinfo, vinfo = t_pars.method.calc_func(transforms)

    # That's all folks
    return wcsinfo, vinfo, transforms


def calc_wcs_tr_202111(transforms: Transforms):
    """Given observatory orientation and target aperture, calculate V1 and Reference Pixel sky coordinates

    A refactor of `calc_wcs_orig` to use the standard `calc_wcs_from_matrix` instead of the specific `calc_aperture_wcs`.

    Parameters
    ----------
    transforms : Transforms
        The transformation matrices.

    Returns
    -------
    wcsinfo, vinfo: WCSRef, WCSRef
        A 2-tuple is returned with the WCS pointing for
        the aperture and the V1 axis.
    """
    # Calculate the V1 WCS information
    vinfo = calc_wcs_from_matrix(transforms.m_eci2v)

    # Calculate the Aperture WCS
    wcsinfo = calc_wcs_from_matrix(transforms.m_eci2siaf)

    # That's all folks
    return wcsinfo, vinfo


def calc_transforms(t_pars: TransformParameters):
    """Calculate transforms  which determine reference point celestial WCS

    This implements Eq. 3 from Technical Report JWST-STScI-003222, SM-12. Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as [follows.]

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


def calc_transforms_coarse_tr_202111(t_pars: TransformParameters):
    """Modified COARSE calculation

    This implements Eq. 45 from Technical Report JWST-STScI-003222, SM-12. Rev. C, 2021-11
    From Section 4:

    In COARSE mode the measured attitude of the J-frame of the spacecraft is
    determined by the star tracker and inertial gyroscopes attitude
    measurements and is converted to an estimated guide star inertial attitude
    using the equations in section 3.2. The V-frame attitude then is determined
    using the equation below.

    One modification from the TR is the calculation of M_eci2siaf. The transformation includes
    the rotation from ICS to Ideal.

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
    quaternion pointing to aperture ra/dec/roll information
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
    logger.info('Calculating transforms using TR 202111 COARSE Tracking with SIAF modification method...')
    t_pars.method = Methods.COARSE_TR_202111

    # Choose the FGS to use.
    # Default to using FGS1 if not specified and FGS1 is not the science instrument.
    fgsid = t_pars.fgsid
    if t_pars.detector is not None:
        detector = t_pars.detector.lower()
        if detector in ['guider1', 'guider2']:
            fgsid = 1
            if detector == 'guider1':
                fgsid = 2
            logger.info(f'COARSE mode using detector {detector} implies use of FGS{fgsid}')
    if fgsid not in FGSIDS:
        fgsid = 1
    t_pars.fgsid = fgsid
    logger.info('Using FGS%s.', t_pars.fgsid)

    # Determine the M_eci_to_gs matrix. Since this is a full train, the matrix
    # is returned as part of the full Transforms object. Many of the required
    # matrices are already determined as part of this calculation.
    t = calc_m_eci2gs(t_pars)

    # Determine the M_fgsx_to_v matrix
    siaf = t_pars.siaf_db.get_wcs(FGSId2Aper[t_pars.fgsid])
    t.m_v2fgsx = calc_v2siaf_matrix(siaf)

    # Determine M_eci_to_v frame.
    t.m_eci2v = np.linalg.multi_dot([np.transpose(t.m_v2fgsx), np.transpose(t.m_fgsx2gs), M_idl2ics, t.m_eci2gs])
    logger.debug('M_eci2v: %s', t.m_eci2v)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate full transformation
    t.m_eci2siaf = np.linalg.multi_dot([M_ics2idl, t.m_v2siaf, t.m_eci2v])
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_track_tr_202111(t_pars: TransformParameters):
    """Calculate transforms for TRACK/FINEGUIDE guiding

    This implements Eq. 46 from Technical Report JWST-STScI-003222, SM-12, Rev. C,  2021-11
    From Section 5:

    Under guide star control the guide star position is measured relative to
    the V-frame. The V3 position angle at the guide star is derived from the
    measured J-frame attitude. Then the corrected guide star catalog position
    is used to determine the inertial V-frame attitude on the sky.

    One modification from the TR is the calculation of M_eci2siaf. The transformation includes
    the rotation from ICS to Ideal.

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
    quaternion pointing to aperture ra/dec/roll information
    is given by the following formula. Each term is a 3x3 matrix:

        M_eci_to_siaf =       # Complete transformation
            M_v_to_siaf    *  # V to SIAF
            M_eci_to_v        # ECI to V

        where

            M_eci_to_v = Conversion of the attitude to a DCM

    """
    logger.info('Calculating transforms using TR 202111 TRACK/FINEGUIDE Tracking method...')
    t_pars.method = Methods.TRACK_TR_202111
    t = Transforms(override=t_pars.override_transforms)  # Shorthand the resultant transforms

    # Check on telemetry for FGS ID. If invalid, use either user-specified or default to 1.
    fgsid = t_pars.pointing.fgsid
    if fgsid not in FGSIDS:
        logger.warning(f'Method {t_pars.method} requires a valid FGS ID in telementry.'
                       '\nHowever telemetry reports an invalid id of {fgsid}')
        if t_pars.fgsid in FGSIDS:
            fgsid = t_pars.fgsid
            logger.warning(f'Using user-specified ID of {fgsid}')
        else:
            fgsid = 1
            logger.warning(f'Using FGS{fgsid} as the default for the guiding FGS')
    t_pars.fgsid = fgsid

    # Determine V3PA@GS
    v3pags = calc_v3pags(t_pars)
    t_pars.guide_star_wcs = WCSRef(t_pars.guide_star_wcs.ra, t_pars.guide_star_wcs.dec, v3pags)

    # Transform the guide star location in ideal detector coordinates to the telescope/V23 frame.
    gs_pos_v23 = trans_fgs2v(t_pars.fgsid, t_pars.pointing.gs_position, t_pars.siaf_db)

    # Calculate the M_eci2v matrix. This is the attitude matrix of the observatory
    # relative to the guide star.
    t.m_eci2v = calc_attitude_matrix(t_pars.guide_star_wcs, v3pags, gs_pos_v23)

    # Calculate the SIAF transform matrix
    t.m_v2siaf = calc_v2siaf_matrix(t_pars.siaf)

    # Calculate the full ECI to SIAF transform matrix
    t.m_eci2siaf = np.linalg.multi_dot([M_ics2idl, t.m_v2siaf, t.m_eci2v])
    logger.debug('m_eci2siaf: %s', t.m_eci2siaf)

    return t


def calc_transforms_ops_tr_202111(t_pars: TransformParameters):
    """Calculate transforms in OPS using TR 2021-11

    This implements the ECI-to-SIAF transformation from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    The actual implementation depends on the guide star mode, represented by the header keyword PCS_MODE.
    For COARSE or NONE, the method COARSE is used.
    For TRACK or FINEGUIDE, the method TRACK is used.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations
    """
    method = method_from_pcs_mode(t_pars.pcs_mode)
    return method.func(t_pars)


def calc_gs2gsapp(m_eci2gsics, jwst_velocity):
    """Calculate the Velocity Aberration correction

    This implements Eq. 40 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2.5:

    The velocity aberration correction is applied in the direction of the guide
    star. The matrix that translates from ECI to the apparent guide star ICS
    frame is M_(ECI→GSAppICS), where the GS Apparent position vector is along
    the z-axis in the guide star ICS frame.

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

    # Eq. 35: Guide star position vector
    uz = np.array([0., 0., 1.])
    u_gseci = np.dot(np.transpose(m_eci2gsics), uz)

    # Eq. 36: Compute the apparent shift due to velocity aberration.
    try:
        scale_factor, u_gseci_app = compute_va_effects_vector(*velocity, u_gseci)
    except TypeError:
        logger.warning('Failure in computing velocity aberration. Returning identity matrix.')
        logger.warning('Exception: %s', sys.exc_info())
        return np.identity(3)

    # Eq. 39: Rotate from ICS into the guide star frame.
    u_gs_app = np.dot(m_eci2gsics, u_gseci_app)

    # Eq. 40: Compute the M_gs2gsapp matrix
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
    r1 = dcm(ra, dec, yangle_ra)

    r2 = np.array([
        [cos(v2) * cos(v3), -sin(v2), -cos(v2) * sin(v3)],
        [sin(v2) * cos(v3), cos(v2), -sin(v2) * sin(v3)],
        [sin(v3), 0., cos(v3)]
    ])

    # Final transformation
    m = np.dot(r2, r1)

    logger.debug('attitude DCM: %s', m)
    return m


def calc_wcs_from_matrix(m):
    """Calculate the WCS information from a DCM.

    Parameters
    ----------
    m : np.array((3, 3))
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


def calc_eci2j_matrix(q):
    """Calculate ECI to J-frame matrix from quaternions

    This implements Eq. 24 from Technical Report JWST-STScI-003222, SM-12. Rev. C, 2021-11
    From Section 3.2.1:

    The M_(ECI→J) DCM is derived from the spacecraft Attitude Control System
    (ACS) attitude quaternion telemetry using the transformation in SE-20,
    Appendix B to transform the attitude quaternion into a DCM.

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

    This implements Eq. 25 from Technical Report JWST-STScI-003222, SM-12. Rev. C, 2021-11
    From Section 3.2.2:

    The M_(J→FGS1ICS) DCM is derived from the transpose of the SC ACS telemetry

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


def calc_sifov2v_matrix():
    """Calculate the SI-FOV to V-Frame matrix

    This is currently defined as the inverse Euler rotation
    about an angle of 7.8 arcmin. Here returns the pre-calculate
    matrix.
    """
    return SIFOV2V_DEFAULT


def calc_v2siaf_matrix(siaf):
    """Calculate the SIAF transformation matrix

    This implements Eq. 12 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
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
    mat = dcm(v2 * A2R, v3 * A2R, v3idlyang * D2R)
    pmat = np.array([[0., vparity, 0.],
                     [0., 0., 1.],
                     [1., 0., 0.]])

    transform = np.dot(pmat, mat)

    logger.debug('transform: %s', transform)
    return transform


def calc_position_angle(point, ref):
    """Calculate position angle from reference to point

    Algorithm implemented is from JWST Technical Report JWST-STScI-001550, SM-12,
    2017-11-08, Rev A., Section 5.2, page 29, final equation:

    tan(pa) = cos(dec_r) * sin(ra_r - ra_p) / (sin(dec_r)cos(dec_p) - cos(dec_r)sin(dec_p)cos(ra_r-ra_p))

    where
        pa : position angle
        *_r : reference
        *_p : point

    Typically the reference is the V3 RA/DEC and point is the object RA/DEC

    Parameters
    ----------
    point : WCSRef
        The POINT wcs parameters, in radians

    ref : WCSRef
        The TARGET wcs parameters, in radians

    Returns
    -------
    point_pa : float
      The POINT position angle, in radians
    """
    y = cos(ref.dec) * sin(ref.ra - point.ra)
    x = sin(ref.dec) * cos(point.dec) - \
        cos(ref.dec) * sin(point.dec) * cos((ref.ra - point.ra))
    point_pa = np.arctan2(y, x)
    if point_pa < 0:
        point_pa += PI2
    if point_pa >= PI2:
        point_pa -= PI2

    logger.debug('Given reference: %s, point: %s, then PA: %s', ref, point, point_pa)
    return point_pa


def get_pointing(obsstart, obsend, mnemonics_to_read=TRACK_TR_202111_MNEMONICS,
                 engdb_url=None, tolerance=60, reduce_func=None):
    """
    Get telescope pointing engineering data.

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

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

    mnemonics = get_mnemonics(obsstart, obsend, mnemonics_to_read=mnemonics_to_read,
                              tolerance=tolerance, engdb_url=engdb_url)
    reduced = reduce_func(mnemonics_to_read, mnemonics)

    logger.log(DEBUG_FULL, 'Mnemonics found:')
    logger.log(DEBUG_FULL, '%s', mnemonics)
    logger.info('Reduced set of pointings:')
    logger.info('%s', reduced)

    return reduced


def vector_to_angle(v):
    """Returns tuple of spherical angles from unit direction Vector

    This implements Eq. 10 & 11 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as:

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
    return alpha, delta


def angle_to_vector(alpha, delta):
    """Convert spherical angles to unit vector

    This implements Eq. 9 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3:

    Parameters
    ----------
    alpha, delta : float
        Spherical angles in radians

    Returns
    -------
    v : [float, float, float]
        Unit vector
    """
    v0 = cos(delta) * cos(alpha)
    v1 = cos(delta) * sin(alpha)
    v2 = sin(delta)

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


def get_mnemonics(obsstart, obsend, tolerance, mnemonics_to_read=TRACK_TR_202111_MNEMONICS, engdb_url=None):
    """Retrieve pointing mnemonics from the engineering database

    Parameters
    ----------
    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to fetch. key is the mnemonic and
        value is whether it is required to be found.

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

    # Construct the mnemonic values structure.
    mnemonics = {
        mnemonic: None
        for mnemonic in mnemonics_to_read
    }

    # Retrieve the mnemonics from the engineering database.
    # Check for whether the bracket values are used and
    # within tolerance.
    for mnemonic in mnemonics:
        try:
            mnemonics[mnemonic] = engdb.get_values(
                mnemonic, obsstart, obsend,
                time_format='mjd', include_obstime=True,
                include_bracket_values=False
            )
        except Exception as exception:
            raise ValueError(f'Cannot retrieve {mnemonic} from engineering.') from exception

        # If more than two points exist, throw off the bracket values.
        # Else, ensure the bracket values are within the allowed time.
        if len(mnemonics[mnemonic]) < 2:
            logger.warning('Mnemonic %s has no telemetry within the observation time.', mnemonic)
            logger.warning('Attempting to use bracket values within %s seconds', tolerance)

            mnemonics[mnemonic] = engdb.get_values(
                mnemonic, obsstart, obsend,
                time_format='mjd', include_obstime=True,
                include_bracket_values=True
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


def all_pointings(mnemonics_to_read, mnemonics):
    """V1 of making pointings

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

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

        # Fill out the matrices
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

        gs_position = None
        if all(k in mnemonics for k in ('SA_ZFGGSPOSX', 'SA_ZFGGSPOSY')):
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
        raise ValueError('No non-zero quaternion found.')

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


def first_pointing(mnemonics_to_read, mnemonics):
    """Return first pointing

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing : Pointing
        First pointing.

    """
    pointings = all_pointings(mnemonics_to_read, mnemonics)
    return pointings[0]


def pointing_from_average(mnemonics_to_read, mnemonics):
    """Determine single pointing from average of available pointings

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    pointing : Pointing
        Pointing from average.

    """
    # Get average observation time.
    times = [
        eng_param.obstime.unix
        for key in mnemonics
        for eng_param in mnemonics[key]
        if eng_param.obstime.unix != 0.0
    ]
    if len(times) > 0:
        obstime = Time(np.average(times), format='unix')
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
        if mnemonics_to_read[mnemonic]:
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

    gs_position = None
    if all(k in mnemonic_averages for k in ('SA_ZFGGSPOSX', 'SA_ZFGGSPOSY')):
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


def fill_mnemonics_chronologically_table(mnemonics, filled_only=True):
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
    filled_by_time : astropy.table.Table
    """
    filled = fill_mnemonics_chronologically(mnemonics, filled_only=filled_only)

    names = [mnemonic for mnemonic in mnemonics]
    names = ['time'] + names
    time_idx = 0

    values = [[] for _ in names]

    for time in filled:
        values[time_idx].append(time)
        for mnemonic in filled[time]:
            idx = names.index(mnemonic)
            values[idx].append(filled[time][mnemonic].value)

    t = Table(values, names=names)

    return t


def calc_estimated_gs_wcs(t_pars: TransformParameters):
    """Calculate the estimated guide star RA/DEC/Y-angle

    This implements Eq. 18, 19, 20 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2:

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    gs_wcs : WCSRef
        Estimated RA, Dec, and Y-angle. All in degrees.
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

    This implements Eq. 21 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    v3pags : float
        The V3 Position Angle at the Guide Star, in degrees

    Notes
    -----
    Modification for `jwst` release post-1.11.3: The commanded position of the guide
    star is always relative to FGS1. Hence, the aperture to use for the SIAF
    transformation is always FGS1.
    """

    # Determine Guides Star estimated WCS information.
    gs_wcs = calc_estimated_gs_wcs(t_pars)

    # Retrieve the Ideal Y-angle for FGS1
    fgs_siaf = t_pars.siaf_db.get_wcs(FGSId2Aper[1], useafter=t_pars.useafter)

    # Calculate V3PAGS
    v3pags = gs_wcs.pa - fgs_siaf.v3yangle

    logger.debug('v3pags: %s', v3pags)
    return v3pags


def calc_m_eci2gs(t_pars: TransformParameters):
    """Calculate the M_eci2gs matrix as per TR presented in 2021-07

    This implements Eq. 16 & 17 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
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
            M_fgs1ics_to_gsics     *
            M_j_to_fgs1ics         *
            M_eci_to_j

        where

            M_eci_to_gs = ECI to Guide Star Ideal Frame
            M_gsics_to_gsappics = Velocity Aberration correction
            M_fgs1ics_to_gsics = Convert from the FGS1 ICS frame to Guide Star ICS frame
            M_j_to_fgs1ics = Convert from J frame to FGS1 ICS frame
            M_eci_to_j = ECI (quaternion) to J-frame

    Modification for `jwst` release post-1.11.3: The commanded position of the guide
    star is always relative to FGS1. Hence, the aperture to use is always FGS1.
    The formulae above have been modified appropriately.
    However, in the code, note that the transformations go to FGS1, but then
    is suddenly referred to thereafter as FGSX. The assumption to make is that X is always 1,
    for FGS1.
    """

    # Initial state of the transforms
    t = Transforms(override=t_pars.override_transforms)

    t.m_eci2j = calc_eci2j_matrix(t_pars.pointing.q)
    t.m_j2fgs1 = calc_j2fgs1_matrix(t_pars.pointing.j2fgs_matrix, t_pars.j2fgs_transpose)
    t.m_fgsx2gs = calc_m_fgsx2gs(t_pars.pointing.gs_commanded)

    # Apply the Velocity Aberration. To do so, the M_eci2gsics matrix must be created. This
    # is used to calculate the aberration matrix.
    # Also, since the aberration is to be removed, the velocity is negated.
    m_eci2gsics = np.linalg.multi_dot([t.m_fgsx2gs, t.m_j2fgs1, t.m_eci2j])
    logger.debug('m_eci2gsics: %s', m_eci2gsics)
    t.m_gs2gsapp = calc_gs2gsapp(m_eci2gsics, t_pars.jwst_velocity)

    # Put it all together
    t.m_eci2gs = np.linalg.multi_dot([M_ics2idl, t.m_gs2gsapp, m_eci2gsics])
    logger.debug('m_eci2gs: %s', t.m_eci2gs)

    # That's all folks
    return t


def calc_m_fgs12fgsx(fgsid, siaf_db):
    """Calculate the FGS1 to FGSx matrix

    This implements Eq. 27 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
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
    fgsid : [1, 2]
        The id of the FGS in use.

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

    This implements Eq. 29 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2.4.

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

    This implements Eq. 30 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2.4.

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
    fgsid : [1, 2]
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

    This implements Eq. 6 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as...

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


def pa_to_roll_ref(pa: float, siaf: SIAF):
    """Calculate Roll from the position angle of the given aperture.

    Parameters
    ----------
    pa : float
        Position angle of the aperture, in degrees.

    siaf : SIAF
        The SIAF of the aperturn

    Returns
    -------
    roll_ref : float
        The roll reference, in degrees
    """
    return pa - siaf.v3yangle


def t_pars_from_model(model, **t_pars_kwargs):
    """Initialize TransformParameters from a DataModel

    Parameters
    ----------
    model : DataModel
        Data model to initialize from.

    t_pars_kwargs : dict
        Keyword arguments used to initialize the TransformParameters object
        before reading from the model meta information.

    Returns
    -------
    t_par : TransformParameters
        The initialized parameters.
    """
    t_pars = TransformParameters(**t_pars_kwargs)

    # Retrieve SIAF information
    if t_pars.siaf is None:
        siaf = None
        useafter = None
        if t_pars.siaf_db is not None:
            aperture_name = model.meta.aperture.name.upper()
            useafter = model.meta.observation.date
            if aperture_name != "UNKNOWN":
                logger.info("Updating WCS for aperture %s", aperture_name)

                # Special case. With aperture MIRIM_TAMRS, the siaf definition is
                # for the subarray of interest. However, the whole detector is
                # read out. Hence, need to convert pixel coordinates to be detector-based.
                to_detector = False
                if aperture_name == 'MIRIM_TAMRS':
                    to_detector = True
                siaf = t_pars.siaf_db.get_wcs(aperture_name, to_detector=to_detector, useafter=useafter)
        t_pars.siaf = siaf
        t_pars.useafter = useafter
    logger.debug('SIAF: %s', t_pars.siaf)

    # Instrument details
    t_pars.detector = model.meta.instrument.detector
    try:
        exp_type = model.meta.exposure.type.lower()
    except AttributeError:
        exp_type = None
    t_pars.exp_type = exp_type

    # observation parameters
    if t_pars.exp_type in FGS_GUIDE_EXP_TYPES:
        t_pars.obsstart = Time(model.meta.observation.date_beg, format='isot').mjd
        t_pars.obsend = Time(model.meta.observation.date_end, format='isot').mjd
    else:
        t_pars.obsstart = model.meta.exposure.start_time
        t_pars.obsend = model.meta.exposure.end_time
    logger.debug('Observation time: %s - %s', t_pars.obsstart, t_pars.obsend)

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

    # Set the transform and WCS calculation method.
    t_pars.method = method_from_pcs_mode(t_pars.pcs_mode)

    # Set pointing reduction function if not already set.
    if not t_pars.reduce_func:
        t_pars.reduce_func = get_reduce_func_from_exptype(t_pars.exp_type)

    return t_pars


def dcm(alpha, delta, angle):
    """Construct the Direction Cosine Matrix (DCM)

    Typical usage is passing of (RA, DEC, PositionAngle).
    All values must be in radians.

    Parameters
    ----------
    alpha : float
        First coordinate in radians.

    delta : float
        Second coordinate in radians.

    angle : float
        Position angle in radians.

    Returns
    -------
    dcm : nd.array((3, 3))
        The 3x3 direction cosine matrix
    """
    dcm = np.array(
        [[cos(delta) * cos(alpha),
          cos(delta) * sin(alpha),
          sin(delta)],
         [-cos(angle) * sin(alpha) + sin(angle) * sin(delta) * cos(alpha),
          cos(angle) * cos(alpha) + sin(angle) * sin(delta) * sin(alpha),
          -sin(angle) * cos(delta)],
         [-sin(angle) * sin(alpha) - cos(angle) * sin(delta) * cos(alpha),
          sin(angle) * cos(alpha) - cos(angle) * sin(delta) * sin(alpha),
          cos(angle) * cos(delta)]])

    return dcm


# Determine calculation method from tracking mode.
def method_from_pcs_mode(pcs_mode):
    """Determine transform/wcs calculation method from PCS_MODE

    Pointing Control System Mode (PCS_MODE) contains the string representing
    which mode the JWST tracking system is in. The orientation calculation
    changes depending on the mode in use.

    Parameters
    ----------
    pcs_mode : str
        The PCS mode in use.

    Returns
    -------
    method : Methods
        The orientation calculation method to use.

    Raises
    ------
    ValueError
        If `pcs_mode` does not uniquely define the method to use.
    """
    if pcs_mode is None or pcs_mode in ['NONE', 'COARSE']:
        return Methods.COARSE_TR_202111
    elif pcs_mode in ['FINEGUIDE', 'MOVING', 'TRACK']:
        return Methods.TRACK_TR_202111
    else:
        raise ValueError(
            f'Invalid PCS_MODE: {pcs_mode}. Should be one of ["NONE", "COARSE", "FINEGUIDE", "MOVING", "TRACK"]'
        )


def get_reduce_func_from_exptype(exp_type):
    """Determine preferred pointing reduction based on exposure type

    Parameters
    ----------
    exp_type : str
        The exposure type.

    Returns
    -------
    reduce_func : func
        The preferred reduction function.
    """
    if exp_type in FGS_ACQ_EXP_TYPES:
        reduce_func = functools.partial(gs_position_acq, exp_type=exp_type)
    elif exp_type in FGS_GUIDED_EXP_TYPES:
        reduce_func = gs_position_fgtrack
    elif exp_type in FGS_ID_EXP_TYPES:
        reduce_func = gs_position_id
    else:
        reduce_func = pointing_from_average

    return reduce_func


def gs_position_acq(mnemonics_to_read, mnemonics, exp_type='fgs_acq1'):
    """Get the guide star position from guide star telemetry for FGS

    The ACQ1 and ACQ2 exposures share nearly the same time range of mnemonics
    with ACQ2 being slightly longer. As such, the information needed by
    both types are interleaved in the engineering.
    ACQ1 needs the first three values of the IFGS_ACQ_*POSG mnemonics.
    ACQ2 needs the next five values of the IFGS_ACQ_*POSG mnemonics.

    The corresponding values of the acquisition box location and size
    are also different between ACQ1 and ACQ2. The values to retrieve
    need be coordinated with the corresponding *POSG values.

    The *POSG values are in arcseconds and are relative to the full
    FGS array. They need to be converted to pixels and shifted
    to be with respect to the actual acquisition box in the exposure.

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    exp_type: str
        The exposure type being dealt with. Either
        'fgs_acq1' or 'fgs_acq2'

    Returns
    =======
    gs_position : GuideStarPosition
        The guide star position
    """
    exp_type = exp_type.lower()
    if exp_type not in FGS_ACQ_EXP_TYPES:
        raise ValueError(f'exp_type {exp_type} not one of {FGS_ACQ_EXP_TYPES}')

    ordered = fill_mnemonics_chronologically_table(mnemonics)
    logger.debug('%s available mnemonics:', exp_type)
    logger.debug('%s', ordered)
    valid = ordered[ordered['IFGS_ACQ_XPOSG'] != 0.0]
    if len(valid) < FGS_ACQ_MINVALUES[exp_type]:
        raise ValueError(
            f'exp_type {exp_type} IFGS_ACQ_XPOSG only has {len(valid)} locations. Requires {FGS_ACQ_MINVALUES[exp_type]}'
        )

    # Setup parameters depending on ACQ1 vs ACQ2
    if exp_type == 'fgs_acq1':
        subarray = valid[FGS_ACQ_WINDOW_INDEX[exp_type]]
        posg_slice = FGS_ACQ_SLICES[exp_type]
    else:
        subarray = valid[FGS_ACQ_WINDOW_INDEX[exp_type]]
        posg_slice = FGS_ACQ_SLICES[exp_type]

    position = (np.average(valid['IFGS_ACQ_XPOSG'][posg_slice]),
                np.average(valid['IFGS_ACQ_YPOSG'][posg_slice]))
    corner = (subarray['IFGS_ACQ_DETXCOR'], subarray['IFGS_ACQ_DETYCOR'])
    size = (subarray['IFGS_ACQ_DETXSIZ'], subarray['IFGS_ACQ_DETYSIZ'])
    gs_position = GuideStarPosition(position=position, corner=corner, size=size)

    return gs_position


def gs_position_fgtrack(mnemonics_to_read, mnemonics):
    """Get the guide star position from guide star telemetry for FGS FINEGUIDE/TRACK

    For FGS FINEGUIDE and TRACK modes, the position of the guide star is
    given by mnemonics IFGS_TFGGS_[X|Y] in arcseconds (Ideal system).
    The values are generally valid throughout the whole length of the observation.
    The average is used as the result. Ignore values where both coordinates are
    exactly zero.

    The position is relative to the full detector. However, the science data is cropped
    to an 8x8 box. The location of the box is defined by the IFGS_TFGDET_* mnemonics.
    These do not change during the duration of the exposure, so the first values are used.

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    gs_position : GuideStarPosition
        The guide star position
    """
    # Remove the zero positions.
    ordered = fill_mnemonics_chronologically_table(mnemonics)
    logger.debug('Guide Star track available mnemonics:')
    logger.debug('%s', ordered)
    valid_flags = (ordered['IFGS_TFGGS_X'] != 0.0) | (ordered['IFGS_TFGGS_Y'] != 0.0)
    valid = ordered[valid_flags]
    if len(valid) < 1:
        raise ValueError('Fine guide track mode, no valid engineering is found.')

    # Get the positions
    position = (np.average(valid['IFGS_TFGGS_X']),
                np.average(valid['IFGS_TFGGS_Y']))
    corner = (valid['IFGS_TFGDET_XCOR'][0], valid['IFGS_TFGDET_YCOR'][0])
    size = (valid['IFGS_TFGDET_XSIZ'][0], valid['IFGS_TFGDET_YSIZ'][0])
    gs_position = GuideStarPosition(position=position, corner=corner, size=size)

    return gs_position


def gs_position_id(mnemonics_to_read, mnemonics):
    """Get the guide star position from guide star telemetry for FGS ID

    For FGS ID, the position of the guide star is given by mnemonics
    IFGS_ID_[X|Y] in arcseconds (Ideal system). This mode is when the desired guide star is identified.
    As such, there is no position to report. However, when completed, the position is then reported for approximately
    five seconds after the end of the exposure. The input mnemonic array needs to have this. The first non-zero
    report is used.

    There is a box defined by the IFGS_ID_DET* mnemonics.

    Parameters
    ==========
    mnemonics_to_read: {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic

    Returns
    =======
    gs_position : GuideStarPosition
        The guide star position

    """
    # Remove the zero positions.
    ordered = fill_mnemonics_chronologically_table(mnemonics)
    logger.debug('Guide STar ID mode available mnemonics:')
    logger.debug('%s', ordered)
    valid_flags = (ordered['IFGS_ID_XPOSG'] != 0.0) | (ordered['IFGS_ID_YPOSG'] != 0.0)
    valid = ordered[valid_flags]
    if len(valid) < 1:
        raise ValueError('Guide Star ID mode, no valid engineering is found.')

    # Get the positions
    position = (valid['IFGS_ID_XPOSG'][0], valid['IFGS_ID_YPOSG'][0])
    corner = (valid['IFGS_ID_DETXCOR'][0], valid['IFGS_ID_DETYCOR'][0])
    size = (valid['IFGS_ID_DETXSIZ'][0], valid['IFGS_ID_DETYSIZ'][0])
    gs_position = GuideStarPosition(position=position, corner=corner, size=size)

    return gs_position


def gs_ideal_to_subarray(gs_position, aperture, flip=False):
    """Calculate pixel position for the guide star in the acquisition subarray

    Parameters
    ----------
    gs_position : GuideStarPosition
        The guide star position telemetry.

    aperture : pysiaf.Aperture
        The aperture in use.

    flip : bool
        Institute the magic flip. Necessary for ACQ1/ACQ2 exposures.

    Returns
    -------
    x, y : float
        The pixel position relative to the subarray
    """
    position_pixel = aperture.idl_to_det(*gs_position.position)
    position_subarray = (position_pixel[0] - gs_position.corner[0],
                         position_pixel[1] - gs_position.corner[1])

    # The magic flip. For FGS1, both axes must be reversed. For FGS2, only the X-axis changed.
    x, y = position_subarray
    if flip:
        if aperture.AperName.startswith('FGS1'):
            x = gs_position.size[0] - position_subarray[0]
            y = gs_position.size[1] - position_subarray[1]
        else:
            x = gs_position.size[0] - position_subarray[0]
            y = position_subarray[1]

    return x, y


def calc_wcs_guiding(model, t_pars, default_roll_ref=0.0, default_vparity=1, default_v3yangle=0.0):
    """Calculate WCS info for FGS guiding

    For Fine Guidance guiding observations, nearly everything
    in the `wcsinfo` meta information is already populated,
    except for the PC matrix and CRVAL*. This function updates the PC
    matrix based on the rest of the `wcsinfo`.

    CRVAL* values are taken from GS_RA/GS_DEC.

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The model to update.

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the V3 position angle.

    default_vparity : int
        The default `VIdlParity` to use and should
        be either "1" or "-1". "1" is the
        default since FGS guiding will be using the
        OSS aperture.

    default_v3yangle : float
        Default SIAF Y-angle.

    Returns
    -------
    crpix1, crpix2, crval1, crval2, pc1_1, pc1_2, pc2_1, pc2_2 : float
        The WCS info.
    """
    # Determine reference pixel

    # Retrieve the appropriate mnemonics that represent the X/Y position of guide star
    # in the image.
    if t_pars.exp_type in ['fgs_acq1', 'fgs_acq2']:
        mnemonics_to_read = FGS_ACQ_MNEMONICS
    elif t_pars.exp_type in ['fgs_fineguide', 'fgs_track']:
        mnemonics_to_read = FGS_GUIDED_MNEMONICS
    elif t_pars.exp_type in FGS_ID_EXP_TYPES:
        mnemonics_to_read = FGS_ID_MNEMONICS
    else:
        raise ValueError(f'Exposure type {t_pars.exp_type} cannot be processed as an FGS product.')

    # Getting the timing to extract from the engineering database is complicated by
    # the fact that the observation times and the processing of the guide star
    # acquisition do not always exactly match. Extend the end time by two seconds.
    #
    # For ID modes, the mnemonics are valid only after the exposure is completed.
    # The time to examine is within 10 seconds after the end of exposure
    obsstart = t_pars.obsstart
    obsend = t_pars.obsend
    if t_pars.exp_type in FGS_ID_EXP_TYPES:
        obsstart = obsend
        obsend = (Time(obsend, format='mjd') + (10 * U.second)).mjd
    elif t_pars.exp_type in FGS_ACQ_EXP_TYPES:
        obsend = (Time(obsend, format='mjd') + (2 * U.second)).mjd

    gs_position = get_pointing(obsstart, obsend,
                               mnemonics_to_read=mnemonics_to_read,
                               engdb_url=t_pars.engdb_url,
                               tolerance=t_pars.tolerance, reduce_func=t_pars.reduce_func)

    crpix1 = crpix2 = None
    apername = f'FGS{t_pars.detector[-1]}_FULL_OSS'
    aperture = t_pars.siaf_db.get_aperture(apername, t_pars.useafter)
    crpix1, crpix2 = gs_ideal_to_subarray(gs_position, aperture, flip=True)

    # Determine PC matrix

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

    pc_matrix = calc_rotation_matrix(roll_ref, np.deg2rad(v3i_yang), vparity=vparity)

    # Determine reference sky values
    crval1 = model.meta.guidestar.gs_ra
    crval2 = model.meta.guidestar.gs_dec

    return crpix1, crpix2, crval1, crval2, pc_matrix
