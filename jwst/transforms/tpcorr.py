"""
A module that provides `TPCorr` class - a `~astropy.modeling.Model` derived
class that applies linear tangent-plane corrections to V2V3 coordinates of
JWST instrument's WCS.

:Authors: Mihai Cara (contact: help@stsci.edu)


"""

# STDLIB
import logging
import math

# THIRD PARTY
import numpy as np
from astropy.modeling import Model, Parameter, InputParameterError

# LOCAL
from . import __version__


__all__ = ['IncompatibleCorrections', 'rot_mat3D', 'TPCorr']

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

    # input_units_strict = False
    # input_units_allow_dimensionless = True
    _separable = False
    standard_broadcasting = False


    r0 = 3600.0 * np.rad2deg(1.0)
    """
    Radius of the generating sphere. This sets the circumference to 360 deg
    so that arc length is measured in deg.
    """

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
        """
        Validates that the input matrix is a 2x2 2D array.

        """

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
        alpha = 3600.0 * np.rad2deg(np.arctan2(y, x))
        delta = 3600.0 * np.rad2deg(np.arctan2(z, h))
        return alpha, delta

    @staticmethod
    def spherical2cartesian(alpha, delta):
        """
        Convert spherical coordinates (in deg) to cartesian.

        """
        alpha = np.deg2rad(alpha / 3600.0)
        delta = np.deg2rad(delta / 3600.0)
        x = np.cos(alpha) * np.cos(delta)
        y = np.cos(delta) * np.sin(alpha)
        z = np.sin(delta)
        return x, y, z

    def v2v3_to_tanp(self, v2, v3):
        """
        Converts V2V3 spherical coordinates to tangent plane coordinates.

        """
        (v2, v3), format_info = self.prepare_inputs(v2, v3)

        # convert spherical coordinates to cartesian assuming unit sphere:
        xyz = self.spherical2cartesian(v2.ravel(), v3.ravel())

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

        # apply corrections:
        # NOTE: order of transforms may need to be swapped depending on the
        #       how shifts are defined.
        xt -= self.shift[0][0]
        yt -= self.shift[0][1]
        xt, yt = np.dot(self.matrix[0], (xt, yt))

        return self.prepare_outputs(format_info, xt.reshape(v2.shape),
                                    yt.reshape(v3.shape))

    def tanp_to_v2v3(self, xt, yt):
        """
        Converts tangent plane coordinates to V2V3 spherical coordinates.

        """
        (xt, yt), format_info = self.prepare_inputs(xt, yt)
        zt = np.full_like(xt, self.__class__.r0)

        # undo corrections:
        xt, yt = np.dot(np.linalg.inv(self.matrix[0]), (xt, yt))
        xt += self.shift[0][0]
        yt += self.shift[0][1]

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
        zcr, xcr, ycr = np.dot(inv_euler_rot, (zt.ravel(), xt.ravel(),
                                               yt.ravel()))

        # convert cartesian to spherical coordinates:
        v2c, v3c = self.cartesian2spherical(zcr, xcr, ycr)

        return self.prepare_outputs(format_info, v2c.reshape(xt.shape),
                                    v3c.reshape(yt.shape))

    def evaluate(self, v2, v3, v2ref, v3ref, roll, matrix, shift):
        """
        Evaluate the model on some input variables.

        """
        (v2, v3), format_info = self.prepare_inputs(v2, v3)

        # convert spherical coordinates to cartesian assuming unit sphere:
        xyz = self.spherical2cartesian(v2.ravel(), v3.ravel())

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

        return self.prepare_outputs(format_info, v2c.reshape(v2.shape),
                                    v3c.reshape(v3.shape))

    @property
    def inverse(self):
        """
        Returns a new `TPCorr` instance which performs the inverse
        transformation of the transformation defined for this `TPCorr` model.

        """
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
