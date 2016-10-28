"""Models used by jwst_pipeline.assign_wcs.

Some of these may go in astropy.modeling in the future.
"""
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function
import math
from collections import namedtuple
import numpy as np
from astropy.modeling.core import Model
from astropy.modeling.parameters import Parameter, InputParameterError
from astropy.modeling.models import Polynomial2D
from astropy.utils import isiterable

__all__ = ['AngleFromGratingEquation', 'WavelengthFromGratingEquation',
           'NRSZCoord', 'Unitless2DirCos', 'DirCos2Unitless',
           'Rotation3DToGWA', 'Gwa2Slit', 'Slit2Msa',
           'Snell', 'RefractionIndex', 'Logical',
           'NirissSOSSModel', 'V23ToSky', 'Slit']


# Number of shutters per quadrant
N_SHUTTERS_QUADRANT = 62415

# Nirspec slit definition
Slit = namedtuple('Slit', ["name", "shutter_id", "xcen", "ycen",
                           "ymin", "ymax", "quadrant", "source_id", "nshutters",
                           "source_name", "source_alias", "catalog_id", "stellarity",
                           "source_xpos", "source_ypos"])
Slit.__new__.__defaults__= ("", 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, "", "", "", 0.0, 0.0, 0.0)


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
        return self.__class__(angles, axes_order=self.axes_order[::-1])

    @staticmethod
    def _compute_matrix(angles, axes_order):
        if len(angles) != len(axes_order):
            raise InputParameterError(
                "Number of angles must equal number of axes in axes_order.")
        matrices = []
        for angle, axis in zip(angles, axes_order):
            matrix = np.zeros((3, 3), dtype=np.float)
            if axis == 'x':
                mat = Rotation3D.rotation_matrix_from_angle(angle)
                matrix[0, 0] = 1
                matrix[1:, 1:] = mat
            elif axis == 'y':
                mat = Rotation3D.rotation_matrix_from_angle(-angle)
                matrix[1, 1] = 1
                matrix[0, 0] = mat[0, 0]
                matrix[0, 2] = mat[0, 1]
                matrix[2, 0] = mat[1, 0]
                matrix[2, 2] = mat[1, 1]
            elif axis == 'z':
                mat = Rotation3D.rotation_matrix_from_angle(angle)
                matrix[2, 2] = 1
                matrix[:2, :2] = mat
            else:
                raise ValueError("Expected axes_order to be a combination \
                        of characters 'x', 'y' and 'z', got {0}".format(
                                     set(axes_order).difference(['x', 'y', 'z'])))
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

    @staticmethod
    def rotation_matrix_from_angle(angle):
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


class Gwa2Slit(Model):
    """
    NIRSpec GWA to slit transform.

    Parameters
    ----------
    slits : list
        open slits
        a slit is a namedtupe
        Slit("name", "shutter_id", "xcen", "ycen", "ymin", "ymax",
             "quadrant", "source_id", "nshutters")
    models : list
        an instance of `~astropy.modeling.core.Model`
    """
    inputs = ('name', 'angle1', 'angle2', 'angle3')
    outputs = ('name', 'x_slit', 'y_slit', 'lam')

    def __init__(self, slits, models):
        if isiterable(slits[0]):
            self._slits = [tuple(s) for s in slits]
            self.slit_ids = [s[0] for s in self._slits]
        else:
            self._slits = list(slits)
            self.slit_ids = self._slits

        self.models = models
        super(Gwa2Slit, self).__init__()

    @property
    def slits(self):
        if isiterable(self._slits[0]):
            return [Slit(*row) for row in self._slits]
        else:
            return self.slit_ids

    def get_model(self, name):
        index = self.slit_ids.index(name)
        return self.models[index]

    def evaluate(self, name, x, y, z):
        index = self.slit_ids.index(name)
        return (name, ) + self.models[index](x, y, z)


class Slit2Msa(Model):
    """
    NIRSpec slit to MSA position transform.

    Parameters
    ----------
    slits : list
        open slits
        a slit is a namedtupe
        Slit("name", "shutter_id", "xcen", "ycen", "ymin", "ymax",
             "quadrant", "source_id", "nshutters")
    models : list
        an instance of `~astropy.modeling.core.Model`
    """

    inputs = ('name', 'x_slit', 'y_slit', 'lam')
    outputs = ('x_msa', 'y_msa', 'lam')

    def __init__(self, slits, models):
        super(Slit2Msa, self).__init__()
        if isiterable(slits[0]):
            self._slits = [tuple(s) for s in slits]
            self.slit_ids = [s[0] for s in self._slits]
        else:
            self._slits = list(slits)
            self.slit_ids = self._slits
        self.models = models

    @property
    def slits(self):
        if isiterable(self._slits[0]):
            return [Slit(*row) for row in self._slits]
        else:
            return self.slit_ids

    def get_model(self, name):
        index = self.slit_ids.index(name)
        return self.models[index]

    def evaluate(self, name, x, y, lam):
        index = self.slit_ids.index(name)
        return self.models[index](x, y) + (lam,)


class NirissSOSSModel(Model):
    """
    NIRISS SOSS wavelength solution.

    Parameters
    ----------
    spectral_orders : list of int
        Spectral orders for which there is a wavelength solution.
    models : list of `~astropy.modeling.core.Model`
        A list of transforms representing the wavelength solution for
        each order in spectral orders. It should match the order in ``spectral_orders``.
    """

    inputs = ('x', 'y', 'spectral_order')
    outputs = ('ra', 'dec', 'lam')

    def __init__(self, spectral_orders, models):
        super(NirissSOSSModel, self).__init__()
        self.spectral_orders = spectral_orders
        self.models = dict(zip(spectral_orders, models))

    def get_model(self, spectral_order):
        return self.models[spectral_order]

    def evaluate(self, x, y, spectral_order):

        # The spectral_order variable is coming in as an array/list of one element.
        # So, we are going to just take the 0'th element and use that as the index.
        try:
            order_number = int(spectral_order[0])
        except Exception as e:
            raise ValueError('Spectral order is not between 1 and 3, {}'.format(spectral_order))

        return self.models[order_number](x, y)


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


class V23ToSky(Rotation3D):
    """
    Transform from V2V3 to a standard coordinate system.

    Parameters
    ----------
    angles : list
        A sequence of angles (in deg).
    axes_order : str
        A sequence of characters ('x', 'y', or 'z') corresponding to the
        axis of rotation and matching the order in ``angles``.
    """

    inputs = ("v2", "v3")
    outputs = ("ra", "dec")

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
        return np.array([x, y, z])

    @staticmethod
    def cartesian2spherical(x, y, z):
        """
        Convert cartesian coordinates to spherical coordinates (in deg).
        """
        h = np.hypot(x, y)
        alpha  = np.rad2deg(np.arctan2(y, x))
        delta = np.rad2deg(np.arctan2(z, h))
        return alpha, delta

    def evaluate(self, v2, v3, angles):
        x, y, z = self.spherical2cartesian(v2, v3)
        x1, y1, z1 = super(V23ToSky, self).evaluate(x, y, z, angles)
        return self.cartesian2spherical(x1, y1, z1)

    def __call__(self, v2, v3):
        from itertools import chain
        inputs, format_info = self.prepare_inputs(v2, v3)
        parameters = self._param_sets(raw=True)

        outputs = self.evaluate(*chain([v2, v3], parameters))

        if self.n_outputs == 1:
            outputs = (outputs,)

        return self.prepare_outputs(format_info, *outputs)
