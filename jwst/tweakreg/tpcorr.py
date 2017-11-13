"""
A module that provides SimpleWCS class - a class that simplifies working with
GWCS objects, in particular as related to WCS transformations sky<->pix
that either include or do not include SIP. Another major purpose
for this class is to allow easy manipulation of WCS parameters related to
standard FITS WCS (CRPIX, CDELT, PC, CRVAL, LONPOLE).

.. warning::
    This class is intended mostly for the internal use by `tweakreg`. This
    class is intended to provide some sort of control in GWCS which, at this
    moment almost completely lacks any kind of standartization. In addition,
    it provides workarounds to some limitations/bugs present in GWCS such as
    the one described here: https://github.com/spacetelescope/gwcs/issues/46

    The API of this class may change in the future as GWCS evolves and bugs
    get fixed and and therefore this class should not be used in external code.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

# STDLIB
import logging
from copy import deepcopy
import math

# THIRD PARTY
import numpy as np
from astropy import units as u
from astropy.modeling import Model, Parameter
from astropy.modeling.models import (
    Shift, Scale, RotateNative2Celestial, AffineTransformation2D
)
from astropy.modeling.rotations import EulerAngleRotation
import gwcs

# LOCAL
from . import __version__
from . import __vdate__


__all__ = ['V23ToTanPlane', 'TanPlaneToV23',
           'IncompatibleCorrections', 'ImageWCS', 'TPCorr', 'rot_mat3D']

__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def rot_mat3D(angle, axis):
    cs = math.cos(angle)
    sn = math.sin(angle)
    axisv = np.array(axis * [0.0] + [1.0] + (2 - axis) * [0.0],
                     dtype=np.float)
    mat2D = np.array([[cs, sn], [-sn, cs]], dtype=np.float)
    return np.insert(np.insert(mat2D, axis, [0.0, 0.0], 1), axis, axisv, 0)


class IncompatibleCorrections(Exception):
    """
    An exception class used to report cases when two or more tangent plane
    corrections cannot be combined into a single one.
    """
    pass


class TPCorr(Model):
    """
    Apply ``V2ref``, ``V3ref``, and ``roll`` to input angles and project
    the point from the tangent plane onto a celestial sphere.

    Parameters
    ----------
    v2ref : float
        V2 position of the reference point in degrees. Default is 0 degrees.

    v3ref : float
        V3 position of the reference point in degrees. Default is 0 degrees.

    roll : float
        Roll angle in degrees. Default is 0 degrees.

    """

    v2ref = Parameter(default=0.0)
    v3ref = Parameter(default=0.0)
    roll = Parameter(default=0.0)
    matrix = Parameter(default=[[1.0, 0.0], [0.0, 1.0]])
    shift = Parameter(default=[0.0, 0.0])

    inputs = ('v2', 'v3')
    outputs = ('v2c', 'v3c')

    #input_units_strict = False
    #input_units_allow_dimensionless = True
    separable = False
    standard_broadcasting = False

    # Radius of the generating sphere. This sets the circumference to 360 deg
    # so that arc length is measured in deg.
    r0 = 3600.0 * np.rad2deg(1.0)

    def __init__(self, v2ref=v2ref.default, v3ref=v3ref.default,
                 roll=roll.default, matrix=matrix.default,
                 shift=shift.default, **kwargs):
        super(TPCorr, self).__init__(
            v2ref=v2ref, v3ref=v3ref, roll=roll, matrix=matrix,
            shift=shift, **kwargs
        )

    @property
    def input_units(self):
        return {'v2': None, 'v3': None}

    @property
    def return_units(self):
        return {'v2c': None, 'v3c': None}

    @matrix.validator
    def matrix(self, value):
        """Validates that the input matrix is a 2x2 2D array."""

        if np.shape(value) != (2, 2):
            raise InputParameterError(
                "Expected transformation matrix to be a 2x2 array")

    @shift.validator
    def shift(self, value):
        """
        Validates that the shift vector is a 2D vector.  This allows
        either a "row" vector.
        """

        if not (np.ndim(value) == 1 and np.shape(value) == (2,)):
            raise InputParameterError(
                "Expected 'shift' to be a 2 element row vector."
            )

    @staticmethod
    def cartesian2spherical(x, y, z):
        """
        Convert cartesian coordinates to spherical coordinates (in deg).
        """
        h = np.hypot(x, y)
        alpha  = np.rad2deg(np.arctan2(y, x))
        delta = np.rad2deg(np.arctan2(z, h))
        return alpha, delta

    @staticmethod
    def spherical2cartesian(alpha, delta):
        """
        Convert spherical coordinates (in deg) to cartesian.
        """
        alpha = np.deg2rad(alpha)
        delta = np.deg2rad(delta)
        x = np.cos(alpha) * np.cos(delta)
        y = np.cos(delta) * np.sin(alpha)
        z = np.sin(delta)
        return x, y, z

    def v2v3_to_tanp(self, v2, v3):
        (v2, v3), format_info = self.prepare_inputs(v2, v3)

        # convert spherical coordinates to cartesian assuming unit sphere:
        xyz = self.spherical2cartesian(v2, v3)

        # build Euler rotation matrices:
        rotm = [
            rot_mat3D(np.deg2rad(alpha), axis)
            for axis, alpha in enumerate(
                [self.roll.value, self.v3ref.value, self.v2ref.value]
            )
        ]
        euler_rot = np.linalg.multi_dot(rotm)

        # rotate cartezian coordinates:
        zr, xr, yr = np.dot(euler_rot, xyz)

        # project points onto the tanject plane
        # (tangent to a sphere of radius r0):
        xt = self.__class__.r0 * xr / zr
        yt = self.__class__.r0 * yr / zr
        zt = np.full_like(xt, self.__class__.r0)

        # apply corrections:
        # NOTE: order of transforms may need to be swapped depending on the
        #       how shifts are defined.
        xt -= self.shift[0][0]
        yt -= self.shift[0][1]
        xt, yt = np.dot(self.matrix[0], (xt, yt))

        return self.prepare_outputs(format_info, xt, yt)

    def tanp_to_v2v3(self, xt, yt):
        (xt, yt), format_info = self.prepare_inputs(xt, yt)
        zt = np.full_like(xt, self.__class__.r0)

        # build Euler rotation matrices:
        rotm = [
            rot_mat3D(np.deg2rad(alpha), axis)
            for axis, alpha in enumerate(
                [self.roll.value, self.v3ref.value, self.v2ref.value]
            )
        ]
        inv_euler_rot = np.linalg.inv(np.linalg.multi_dot(rotm))

        # "unrotate" cartezian coordinates back to their original
        # v2ref, v3ref, and roll "positions":
        zcr, xcr, ycr = np.dot(inv_euler_rot, (zt, xt, yt))

        # convert cartesian to spherical coordinates:
        v2c, v3c = self.cartesian2spherical(zcr, xcr, ycr)

        return self.prepare_outputs(format_info, v2c, v3c)

    def evaluate(self, v2, v3, v2ref, v3ref, roll, matrix, shift):
        (v2, v3), format_info = self.prepare_inputs(v2, v3)

        # convert spherical coordinates to cartesian assuming unit sphere:
        xyz = self.spherical2cartesian(v2, v3)

        # build Euler rotation matrices:
        rotm = [rot_mat3D(np.deg2rad(alpha), axis)
                for axis, alpha in enumerate([roll[0], v3ref[0], v2ref[0]])]
        euler_rot = np.linalg.multi_dot(rotm)
        inv_euler_rot = np.linalg.inv(euler_rot)

        # rotate cartezian coordinates:
        zr, xr, yr = np.dot(euler_rot, xyz)

        # project points onto the tanject plane
        # (tangent to a sphere of radius r0):
        xt = self.__class__.r0 * xr / zr
        yt = self.__class__.r0 * yr / zr
        zt = np.full_like(xt, self.__class__.r0)

        # apply corrections:
        # NOTE: order of transforms may need to be swapped depending on the
        #       how shifts are defined.
        xt -= shift[0][0]
        yt -= shift[0][1]
        xt, yt = np.dot(matrix[0], (xt, yt))

        # "unrotate" cartezian coordinates back to their original
        # v2ref, v3ref, and roll "positions":
        zcr, xcr, ycr = np.dot(inv_euler_rot, (zt, xt, yt))

        # convert cartesian to spherical coordinates:
        v2c, v3c = self.cartesian2spherical(zcr, xcr, ycr)

        return self.prepare_outputs(format_info, v2c, v3c)

    @property
    def inverse(self):
        ishift = -np.dot(self.matrix.value, self.shift.value)
        imatrix = np.linalg.inv(self.matrix.value)
        return TPCorr(v2ref=self.v2ref.value, v3ref=self.v3ref.value,
                      roll=self.roll.value, matrix=imatrix, shift=ishift)

    @classmethod
    def combine(cls, t2, t1):
        """
        Combine transformation ``t2`` with another transformation (``t1``)
        *previously applied* to the coordinates. That is,
        transformation ``t2`` is assumed to *follow* (=applied after) the
        transformation provided by the argument ``t1``.

        """
        if not isinstance(t1, TPCorr):
            raise IncompatibleCorrections(
                "Tangent plane correction 't1' is not a TPCorr instance."
            )

        if not isinstance(t2, TPCorr):
            raise IncompatibleCorrections(
                "Tangent plane correction 't2' is not a TPCorr instance."
            )

        eps_v2 = 10.0 * np.finfo(t2.v2ref.value).eps
        eps_v3 = 10.0 * np.finfo(t2.v3ref.value).eps
        eps_roll = 10.0 * np.finfo(t2.roll.value).eps
        if not (np.isclose(t2.v2ref.value, t1.v2ref.value, rtol=eps_v2) and
                np.isclose(t2.v3ref.value, t1.v3ref.value, rtol=eps_v3) and
                np.isclose(t2.roll.value, t1.roll.value, rtol=eps_roll)):
            raise IncompatibleCorrections(
                "Only combining of correction transformations within the same "
                "tangent plane is supported."
            )

        t1m = t1.matrix.value
        it1m = np.linalg.inv(t1m)
        shift = t1.shift.value + np.dot(it1m, t2.shift.value)
        matrix = np.dot(t2.matrix.value, t1m)

        name = t1.name if t2.name is None else t2.name

        return cls(v2ref=t2.v2ref.value, v3ref=t2.v3ref.value,
                   roll=t2.roll.value, matrix=matrix, shift=shift, name=name)


class ImageWCS(object):

    def __init__(self, wcs, v2_ref, v3_ref, roll_ref, ra_ref, dec_ref):
        if not self.check_wcs_structure(wcs):
            raise ValueError("Unsupported WCS structure.")

        self._ra_ref = ra_ref
        self._dec_ref = dec_ref
        self._v2_ref = v2_ref
        self._v3_ref = v3_ref
        self._roll_ref = roll_ref

        # perform additional check that if tangent plane correction is already
        # present in the WCS pipeline, it is of TPCorr class and that
        # its parameters are consistent with reference angles:
        frms = [f[0] for f in wcs.pipeline]
        if 'v2v3corr' in frms:
            self._v23name = 'v2v3corr'
            self._tpcorr = deepcopy(wcs.pipeline[frms.index('v2v3corr')-1][1])
            self._default_tpcorr = None
            if not isinstance(self._tpcorr, TPCorr):
                raise ValueError("Unsupported tangent-plance correction type.")

            # check that transformation parameters are consistent with
            # reference angles:
            v2ref = self._tpcorr.v2ref.value
            v3ref = self._tpcorr.v3ref.value
            roll = self._tpcorr.roll.value
            eps_v2 = 10.0 * np.finfo(v2_ref).eps
            eps_v3 = 10.0 * np.finfo(v3_ref).eps
            eps_roll = 10.0 * np.finfo(roll_ref).eps
            if not (np.isclose(v2_ref, v2ref, rtol=eps_v2) and
                    np.isclose(v3_ref, v3ref, rtol=eps_v3) and
                    np.isclose(roll_ref, roll, rtol=eps_roll)):
                raise ValueError(
                    "WCS/TPCorr parameters 'v2ref', 'v3ref', and/or 'roll' "
                    "differ from the corresponding reference values."
                )

        else:
            self._v23name = 'v2v3'
            self._tpcorr = None
            self._default_tpcorr = TPCorr(
                v2ref=v2_ref, v3ref=v3_ref,
                roll=roll_ref,
                name='tangent-plane linear correction'
            )


        self._owcs = wcs
        self._wcs = deepcopy(wcs)
        self._update_transformations()

    def _update_transformations(self):
        # define transformations from detector/world coordinates to
        # the tangent plane:
        detname = self._wcs.pipeline[0][0]
        worldname = self._wcs.pipeline[-1][0]

        self._world_to_v23 = self._wcs.get_transform(worldname, self._v23name)
        self._v23_to_world = self._wcs.get_transform(self._v23name, worldname)
        self._det_to_v23 = self._wcs.get_transform(detname, self._v23name)
        self._v23_to_det = self._wcs.get_transform(self._v23name, detname)

        self._det_to_world = self._wcs.__call__
        self._world_to_det = self._wcs.invert

    @property
    def ref_angles(self):
        wcsinfo = {
            'ra_ref': self._ra_ref,
            'dec_ref': self._dec_ref,
            'v2_ref': self._v2_ref,
            'v3_ref': self._v3_ref,
            'roll_ref': self._roll_ref
        }
        return wcsinfo

    @property
    def wcs(self):
        return self._wcs

    @property
    def original_wcs(self):
        return self._owcs

    def copy(self):
        return deepcopy(self)

    def set_correction(self, matrix=[[1, 0], [0, 1]], shift=[0, 0]):
        frms = [f[0] for f in self._wcs.pipeline]

        # if original WCS did not have tangent-plane corrections, create
        # new correction and add it to the WCs pipeline:
        if self._tpcorr is None:
            self._tpcorr = TPCorr(
                v2ref=self._v2_ref, v3ref=self._v3_ref,
                roll=self._roll_ref, matrix=matrix, shift=shift,
                name='tangent-plane linear correction'
            )
            idx_v2v3 = frms.index(self._v23name)
            pipeline = deepcopy(self._wcs.pipeline)
            pf, pt = pipeline[idx_v2v3]
            pipeline[idx_v2v3] = (pf, deepcopy(self._tpcorr))
            pipeline.insert(idx_v2v3 + 1, ('v2v3corr', pt))
            self._wcs = gwcs.WCS(pipeline, name=self._owcs.name)
            self._v23name = 'v2v3corr'

        else:
            # combine old and new corrections into a single one and replace
            # old transformation with the combined correction transformation:
            tpcorr2 = self._tpcorr.__class__(
                v2ref=self._tpcorr.v2ref, v3ref=self._tpcorr.v3ref,
                roll=self._tpcorr.roll, matrix=matrix, shift=shift,
                name='tangent-plane linear correction'
            )

            self._tpcorr = tpcorr2.combine(tpcorr2, self._tpcorr)

            idx_v2v3 = frms.index(self._v23name)
            pipeline = deepcopy(self._wcs.pipeline)
            pipeline[idx_v2v3 - 1] = (pipeline[idx_v2v3 - 1][0],
                                      deepcopy(self._tpcorr))
            self._wcs = gwcs.WCS(pipeline, name=self._owcs.name)

        # reset definitions of the transformations from detector/world
        # coordinates to the tangent plane:
        self._update_transformations()

    def check_wcs_structure(self, wcs):
        if wcs is None or wcs.pipeline is None:
            return False

        frms = [f[0] for f in wcs.pipeline]
        nframes = len(frms)
        if nframes < 3:
            return False

        if frms.count(frms[0]) > 1 or frms.count(frms[-1]) > 1:
            return False

        if frms.count('v2v3') != 1:
            return False

        idx_v2v3 = frms.index('v2v3')
        if idx_v2v3 == 0 or idx_v2v3 == (nframes - 1):
            return False

        nv2v3corr = frms.count('v2v3corr')
        if nv2v3corr == 0:
            return True
        elif nv2v3corr > 1:
            return False

        idx_v2v3corr = frms.index('v2v3corr')
        if idx_v2v3corr != (idx_v2v3 + 1) or idx_v2v3corr == (nframes - 1):
            return False

        return True

    def det_to_world(self, x, y):
        ra, dec = self._det_to_world(x, y)
        return ra, dec

    def world_to_det(self, ra, dec):
        x, y = self._world_to_det(ra, dec)
        return x, y

    def det_to_tanp(self, x, y):
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = self._det_to_v23(x, y)
        x, y = tpc.v2v3_to_tanp(v2, v3)
        return x, y

    def tanp_to_det(self, x, y):
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = tpc.tanp_to_v2v3(x, y)
        x, y = self._v23_to_det(v2, v3)
        return x, y

    def world_to_tanp(self, ra, dec):
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = self._world_to_v23(ra, dec)
        x, y = tpc.v2v3_to_tanp(v2, v3)
        return x, y

    def tanp_to_world(self, x, y):
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = tpc.tanp_to_v2v3(x, y)
        ra, dec = self._v23_to_world(v2, v3)
        return ra, dec

