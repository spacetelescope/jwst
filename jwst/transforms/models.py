"""Models used by jwst_pipeline.assign_wcs.

Some of these may go in astropy.modeling in the future.
"""
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

import math
import numpy as np
from astropy.modeling.core import Model
from astropy.modeling.parameters import Parameter, InputParameterError
from astropy.modeling.models import Polynomial2D


__all__ = ['AngleFromGratingEquation', 'WavelengthFromGratingEquation',
           'NRSZCoord', 'Unitless2DirCos', 'DirCos2Unitless',
           'Rotation3DToGWA', 'Gwa2Slit', 'Slit2Msa', 'slitid_to_slit',
           'slit_to_slitid', 'Snell', 'RefractionIndex', 'Logical',
           'NirissSOSSModel']


# Number of shutters per quadrant
N_SHUTTERS_QUADRANT = 62415


class RefractionIndex(Model):
    """Calculate a refraction index given lamda and a reference file.

    Referenced in NTN-2014-004/ESA-JWST-TN-20930
    Requires  reference file information

    Parameters
    ----------
    temp: float
        the instrument temperature in kelvins
    pressure: float
        local pressure in atmospheres
    pref: float
        the reference pressure
    kcoef: list of floats
        three items which are fit coefficients
    lcoef: list of floats
        three items which are fit coefficients
    tref: float
        the reference temperature in kelvins
    """

    fittable = False
    separable = False

    inputs = ("lam",)
    outputs = ("n",)

    temp = Parameter()
    tref = Parameter()
    pref = Parameter()
    pressure = Parameter()
    kcoef = Parameter()
    lcoef = Parameter()
    tcoef = Parameter()

    def __init__(self, temp, tref, pref, pressure, kcoef, lcoef, tcoef,
                 name=None):
        super(RefractionIndex, self).__init__(temp=temp, tref=tref,
                                              pref=pref, pressure=pressure,
                                              kcoef=kcoef, lcoef=lcoef,
                                              tcoef=tcoef, name=name)

    def evaluate(self, lam, temp, tref, pref, pressure, kcoef, lcoef, tcoef):
        """Calculate and retrun the refraction index."""

        # scale input wavelength to air at reference temp and
        # pressure (T=35 K, P=1 atm)
        # convert the wavelength to microns
        lam = lam * 1e6
        KtoC = 273.15  # kelvin to celcius conversion

        nref = 1. + (6432.8 + 2949810. * lam**2 /
                     (146.0 * lam**2 - 1.) + 25540.0 * lam**2 /
                     (41.0 * lam**2 - 1.)) * 1e-8

        nair_obs = 1.0 + (nref - 1.0) * pressure / \
            (1.0 + (temp - KtoC - 15.) * 2.4785e-3)

        nair_ref = 1.0 + (nref - 1.0) * pref / \
            (1.0 + tref - KtoC)

        lamrel = lam * nair_obs / nair_ref
        nrel = np.sqrt(1.0 + kcoef[0] * lamrel**2 /
                       (lamrel**2 - lcoef[0]) +
                       kcoef[1] * lamrel**2 / (lamrel**2 - lcoef[1]) +
                       kcoef[2] * lamrel**2 / (lamrel**2 - lcoef[2]))
        nabs_ref = nrel * nair_ref

        # compute the absolute index of the glass
        delt = temp - tref
        delnabs = 0.5 * (nrel**2 - 1.) / nrel * \
            (tcoef[0] * delt + tcoef[1] * delt**2 + tcoef[2] * delt**2 +
             (tcoef[3] * delt + tcoef[4] * delt**2) / (lamrel**2 - tcoef[5]**2))
        nabs_obs = nabs_ref + delnabs
        return nabs_obs / nair_obs


class Snell(Model):
    """Computes the Prism Snell refraction through a surface.

    Parameters
    ----------
    n: float
        refraction index as calculated for a given wavelenth
    """

    inputs = ("xin", "yin", "zin")
    outputs = ("xout", "yout", "zout")

    n = Parameter(default=1.0)

    def __init__(self, n=n, name=None):
        super(Snell, self).__init__(n=n, name=name)

    def evaluate(self, x, y, z, n):
        """Compute Snell's refraction law from the front surface."""

        xout = x / n
        yout = y / n
        zout = np.sqrt(1.0 - xout**2 - yout**2)
        return xout, yout, zout

    def inverse(self, x, y, z, n):
        """Compute Snell's refraction law from the back surface."""

        xout = x * n
        yout = y * n
        zout = np.sqrt(1.0 - xout**2 - yout**2)
        return xout, yout, zout


class NRSChromaticCorrection(Polynomial2D):

    def __init__(self, degree, **coeffs):
        super(NRSChromaticCorrection, self).__init__(degree, **coeffs)

    def evaluate(self, x, y, lam, *coeffs):
        """For each input multiply the distortion coefficients by the computed lambda.
        """
        coeffs *= lam
        return super(NRSChromaticCorrection, self).evaluate(x, y, *coeffs)


class AngleFromGratingEquation(Model):
    """
    Grating Equation Model. Computes the diffracted/refracted angle.

    Parameters
    ----------
    groove_density : int
        Grating ruling density.
    order : int
        Spectral order.
    """

    separable = False

    inputs = ("lam", "alpha_in", "beta_in", "z")
    outputs = ("alpha_out", "beta_out", "zout")

    groove_density = Parameter()
    order = Parameter(default=-1)

    def evaluate(self, lam, alpha_in, beta_in, z, groove_density, order):
        if alpha_in.shape != beta_in.shape != z.shape:
            raise ValueError("Expected input arrays to have the same shape")
        orig_shape = alpha_in.shape or (1,)
        xout = -alpha_in - groove_density * order * lam
        yout = - beta_in
        zout = np.sqrt(1 - xout**2 - yout**2)
        xout.shape = yout.shape = zout.shape = orig_shape
        return xout, yout, zout


class WavelengthFromGratingEquation(Model):
    """Grating Equation Model. Computes the wavelength.

    Parameters
    ----------
    groove_density : int
        Grating ruling density.
    order : int
        Spectral order.
    """

    separable = False

    inputs = ("alpha_in", "alpha_out")
    outputs = ("lam",)

    groove_density = Parameter()
    order = Parameter(default=1)

    def evaluate(self, alpha_in, alpha_out, groove_density, order):
        return -(alpha_in + alpha_out) / (groove_density * order)


class NRSZCoord(Model):
    """
    Class to compute the z coordinate through the NIRSPEC grating wheel.

    """
    separable = False

    inputs = ("x", "y")
    outputs = ("z",)

    def evaluate(self, x, y):
        return np.sqrt(1 - (x**2 + y**2))


class Unitless2DirCos(Model):
    """
    Vector to directional cosines.
    """
    separable = False

    inputs = ('x', 'y')
    outputs = ('x', 'y', 'z')

    def evaluate(self, x, y):
        vabs = np.sqrt(1. + x**2 + y**2)
        cosa = x / vabs
        cosb = y / vabs
        cosc = 1. / vabs
        return cosa, cosb, cosc

    def inverse(self):
        return DirCos2Unitless()


class DirCos2Unitless(Model):
    """
    Directional Cosines to vector.
    """
    separable = False

    inputs = ('x', 'y', 'z')
    outputs = ('x', 'y')

    def evaluate(self, x, y, z):

        return x / z, y / z

    def inverse(self):
        return Unitless2DirCos()


class Rotation3DToGWA(Model):
    separable = False

    """
    Perform a 3D rotation given an angle in degrees.

    Positive angles represent a counter-clockwise rotation and vice-versa.

    Parameters
    ----------
    angles : array-like
        Angles of rotation in deg in the order of axes_order.
    axes_order : str
        A sequence of 'x', 'y', 'z' corresponding of axis of rotation/
    """
    standard_broadcasting = False

    inputs = ('x', 'y', 'z')
    outputs = ('x', 'y', 'z')

    angles = Parameter(getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, angles, axes_order, name=None):
        if len(angles) != len(axes_order):
            raise InputParameterError(
                "Number of angles must equal number of axes in axes_order.")

        self.axes = ['x', 'y', 'z']
        unrecognized = set(axes_order).difference(self.axes)
        if unrecognized:
            raise ValueError("Unrecognized axis label {0}; "
                             "should be one of {1} ".format(unrecognized,
                                                            self.axes))
        self.axes_order = axes_order

        self._func_map = {'x': self._xrot,
                          'y': self._yrot,
                          'z': self._zrot
                          }
        super(Rotation3DToGWA, self).__init__(angles, name=name)

    @property
    def inverse(self):
        """Inverse rotation."""
        angles = self.angles.value[::-1] * -1
        return self.__class__(angles, self.axes_order[::-1])

    def _xrot(self, x, y, z, theta):
        xout = x
        yout = y * np.cos(theta) + z * np.sin(theta)
        zout = np.sqrt(1 - xout ** 2 - yout ** 2)
        return [xout, yout, zout]

    def _yrot(self, x, y, z, theta):
        xout = x * np.cos(theta) - z * np.sin(theta)
        yout = y
        zout = np.sqrt(1 - xout ** 2 - yout ** 2)
        return [xout, yout, zout]

    def _zrot(self, x, y, z, theta):
        xout = x * np.cos(theta) + y * np.sin(theta)
        yout = -x * np.sin(theta) + y * np.cos(theta)
        zout = np.sqrt(1 - xout ** 2 - yout ** 2)
        return [xout, yout, zout]

    def evaluate(self, x, y, z, angles):
        """
        Apply the rotation to a set of 3D Cartesian coordinates.

        """

        if x.shape != y.shape != z.shape:
            raise ValueError("Expected input arrays to have the same shape")

        #  Note: If the original shape was () (an array scalar) convert to a
        #  1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)
        for ang, ax in zip(angles[0], self.axes_order):
            x, y, z = self._func_map[ax](x, y, z, theta=ang)
        x.shape = y.shape = z.shape = orig_shape

        return x, y, z


class Rotation3D(Model):

    separable = False
    """
    Perform a 3D rotation given an angle in degrees.
    Positive angles represent a counter-clockwise rotation and vice-versa.
    Parameters
    ----------
    angles : array-like
        Angles of rotation in deg in the order of axes_order.
    axes_order : str
        A sequence of 'x', 'y', 'z' corresponding of axis of rotation/
    """
    standard_broadcasting = False
    inputs = ('x', 'y', 'z')
    outputs = ('x', 'y', 'z')
    angles = Parameter(getter=np.rad2deg, setter=np.deg2rad)

    def __init__(self, angles, axes_order, name=None):
        self.axes = ['x', 'y', 'z']
        unrecognized = set(axes_order).difference(self.axes)
        if unrecognized:
            raise ValueError("Unrecognized axis label {0}; "
                             "should be one of {1} ".format(unrecognized,
                                                            self.axes))
        self.axes_order = axes_order
        if len(angles) != len(axes_order):
            raise ValueError("The number of angles {0} should match the number \
                              of axes {1}.".format(len(angles),
                                                   len(axes_order)))
        super(Rotation3D, self).__init__(angles, name=name)

    @property
    def inverse(self):
        """Inverse rotation."""
        angles = self.angles.value[::-1] * -1
        return self.__class__(angles, self.axes_order[::-1])

    def _compute_matrix(self, angles, axes_order):
        if len(angles) != len(axes_order):
            raise InputParameterError(
                "Number of angles must equal number of axes in axes_order.")
        matrices = []
        for angle, axis in zip(angles, axes_order):
            matrix = np.zeros((3, 3), dtype=np.float)
            if axis == 'x':
                mat = self._rotation_matrix_from_angle(angle)
                matrix[0, 0] = 1
                matrix[1:, 1:] = mat
            elif axis == 'y':
                mat = self._rotation_matrix_from_angle(-angle)
                matrix[1, 1] = 1
                matrix[0, 0] = mat[0, 0]
                matrix[0, 2] = mat[0, 1]
                matrix[2, 0] = mat[1, 0]
                matrix[2, 2] = mat[1, 1]
            elif axis == 'z':
                mat = self._rotation_matrix_from_angle(angle)
                matrix[2, 2] = 1
                matrix[:2, :2] = mat
            else:
                raise ValueError("Expected axes_order to be a combination \
                        of characters 'x', 'y' and 'z', got {0}".format(
                                     set(axes_order).difference(self.axes)))
            matrices.append(matrix)
        if len(angles) == 1:
            return matrix
        elif len(matrices) == 2:
            return np.dot(matrices[1], matrices[0])
        else:
            prod = np.dot(matrices[1], matrices[0])
            for m in matrices[2:]:
                prod = np.dot(m, prod)
            return prod

    def _rotation_matrix_from_angle(self, angle):
        """
        Clockwise rotation matrix.
        """
        return np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]])

    def evaluate(self, x, y, z, angles):
        """
        Apply the rotation to a set of 3D Cartesian coordinates.
        """
        if x.shape != y.shape != z.shape:
            raise ValueError("Expected input arrays to have the same shape")
        # Note: If the original shape was () (an array scalar) convert to a
        # 1-element 1-D array on output for consistency with most other models
        orig_shape = x.shape or (1,)
        inarr = np.array([x.flatten(), y.flatten(), z.flatten()])
        result = np.dot(self._compute_matrix(angles[0], self.axes_order), inarr)
        x, y, z = result[0], result[1], result[2]
        x.shape = y.shape = z.shape = orig_shape
        return x, y, z


class LRSWavelength(Model):

    standard_broadcasting = False

    linear = False
    fittable = False

    inputs = ('x', 'y')
    outputs = ('lambda',)

    def __init__(self, wavetable, zero_point, name=None):
        self._wavetable = wavetable
        self._zero_point = zero_point
        super(LRSWavelength, self).__init__(name=name)

    @property
    def wavetable(self):
        return self._wavetable

    @property
    def zero_point(self):
        return self._zero_point

    def evaluate(self, x, y):
        slitsize = 1.00076751
        imx, imy = self.zero_point
        dx = x - imx
        dy = y - imy
        if x.shape != y.shape:
            raise ValueError("Inputs have different shape.")
        x0 = self._wavetable[:, 3]
        y0 = self._wavetable[:, 4]
        x1 = self._wavetable[:, 5]
        #y1 = self._wavetable[:, 6]
        wave = self._wavetable[:, 2]

        diff0 = (dy - y0[0])
        ind = np.abs(np.asarray(diff0 / slitsize, dtype=np.int))
        condition = np.logical_and(dy < y0[0], dy > y0[-1])  #, dx>x0, dx<x1)
        xyind = condition.nonzero()
        wavelength = np.zeros(condition.shape)
        wavelength += np.nan
        wavelength[xyind] = wave[ind[xyind]]
        wavelength = wavelength.flatten()

        wavelength[(dx[xyind] < x0[ind[xyind]]).nonzero()[0]] = np.nan
        wavelength[(dx[xyind] > x1[ind[xyind]]).nonzero()[0]] = np.nan
        wavelength.shape = condition.shape

        return wavelength


def slitid_to_slit(open_slits_id):
    """
    A slit_id is a tuple of (quadrant_number, slit_number)
    Internally a slit is represented as a number
    slit = quadrant_number * N_SHUTTERS_QUADRANT + slit_number
    """
    return open_slits_id[:, 0] * N_SHUTTERS_QUADRANT + open_slits_id[:, 1]


def slit_to_slitid(slits):
    """
    Return the slitid for the slits.
    """
    slits = np.asarray(slits)
    return np.array(zip(*divmod(slits, N_SHUTTERS_QUADRANT)))


class Gwa2Slit(Model):

    inputs = ('angle1', 'angle2', 'angle3', 'quadrant', 'slitid')
    outputs = ('x_slit', 'y_slit', 'lam', 'quadrant', 'slitid')

    def __init__(self, models):
        self.slits = slit_to_slitid(list(models.keys()))
        self.models = models
        super(Gwa2Slit, self).__init__()

    def evaluate(self, quadrant, slitid, x, y, z):
        slit = int(slitid_to_slit(np.array([quadrant, slitid]).T)[0])
        return (quadrant, slitid) + self.models[slit](x, y, z)


class Slit2Msa(Model):

    inputs = ('quadrant', 'slitid', 'x_slit', 'y_slit', 'lam')
    outputs = ('x_msa', 'y_msa', 'lam')

    def __init__(self, models):
        super(Slit2Msa, self).__init__()
        self.slits = slit_to_slitid(list(models.keys()))
        self.models = models

    def evaluate(self, quadrant, slitid, x, y, lam):
        slit = int(slitid_to_slit(np.array([quadrant, slitid]).T)[0])
        return self.models[slit](x, y) + (lam,)


class NirissSOSSModel(Model):

    inputs = ('x', 'y', 'spectral_order')
    outputs = ('ra', 'dec', 'lam')

    def __init__(self, spectral_orders, models):
        super(NirissSOSSModel, self).__init__()
        self.spectral_orders = np.asarray(spectral_orders)
        self.models = dict(zip(spectral_orders, models))

    def evaluate(self, x, y, spectral_order):
        return self.models[spectral_order](x, y)


class Logical(Model):
    """
    Substitute values in an array where the condition is evaluated to True.

    Similar to numpy's where function.

    Parameters
    ----------
    condition : str
        A string representing the logical, one of GT, LT, NE, EQ
    compareto : float, ndarray
        A number to compare to using the condition
        If ndarray then the input array, compareto and value should have the
        same shape.
    value : float, ndarray
        Value to substitute where condition is True.
    """
    inputs = ('x', )
    outputs = ('x', )

    conditions = {'GT': np.greater,
                  'LT': np.less,
                  'EQ': np.equal,
                  'NE': np.not_equal
                  }

    def __init__(self, condition, compareto, value, **kwargs):
        self.condition = condition.upper()
        self.compareto = compareto
        self.value = value
        super(Logical, self).__init__(**kwargs)

    def evaluate(self, x):
        x = x.copy()
        cond_result = self.conditions[self.condition](x, self.compareto)
        if isinstance(self.compareto, np.ndarray):
            x[cond_result] = self.value[cond_result]
        else:
            x[cond_result] = self.value
        return x

    def __repr__(self):
        return "{0}(condition={1}, compareto={2}, value={3})".format(self.__class__.__name__,
                                                                     self.condition, self.compareto,
                                                                     self.value)
