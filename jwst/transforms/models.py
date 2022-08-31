# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Models used by the JWST pipeline.

The models are written using the astropy.modeling framework.
Since they are specific to JWST, the models and their ASDF schemas
are kept here separately from astropy. An ASDF extension for this package is
registered with ASDF through entry points.
"""


import math
from collections import namedtuple
import numpy as np
from astropy.modeling.core import Model
from astropy.modeling.parameters import Parameter, InputParameterError
from astropy.modeling.models import (Rotation2D, Identity, Mapping, Tabular1D, Const1D)
from astropy.modeling.models import math as astmath
from astropy.utils import isiterable


__all__ = ['Gwa2Slit', 'Slit2Msa', 'Logical', 'NirissSOSSModel', 'Slit',
           'NIRCAMForwardRowGrismDispersion', 'NIRCAMForwardColumnGrismDispersion',
           'NIRCAMBackwardGrismDispersion', 'MIRI_AB2Slice', 'GrismObject',
           'NIRISSForwardRowGrismDispersion', 'NIRISSForwardColumnGrismDispersion',
           'NIRISSBackwardGrismDispersion', 'V2V3ToIdeal', 'IdealToV2V3',
           'RefractionIndexFromPrism', 'Snell', 'V23ToSky', 'Rotation3DToGWA',
           'AngleFromGratingEquation', 'WavelengthFromGratingEquation']


N_SHUTTERS_QUADRANT = 62415
""" Number of shutters per quadrant in the NIRSPEC MSA shutter array"""


Slit = namedtuple('Slit', ["name", "shutter_id", "dither_position", "xcen", "ycen",
                           "ymin", "ymax", "quadrant", "source_id", "shutter_state",
                           "source_name", "source_alias", "stellarity",
                           "source_xpos", "source_ypos", "source_ra", "source_dec"])
""" Nirspec Slit structure definition"""


Slit.__new__.__defaults__ = ("", 0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, "", "", "",
                             0.0, 0.0, 0.0, 0.0, 0.0)


class GrismObject(namedtuple('GrismObject', ("sid",
                                             "order_bounding",
                                             "sky_centroid",
                                             "partial_order",
                                             "waverange",
                                             "sky_bbox_ll",
                                             "sky_bbox_lr",
                                             "sky_bbox_ur",
                                             "sky_bbox_ul",
                                             "xcentroid",
                                             "ycentroid",
                                             "is_extended",
                                             "isophotal_abmag",
                                             ), rename=False)):
    """ Grism Objects identified from a direct image catalog and segment map.

    Parameters
    ----------
    sid : int
        source identified
    xcentroid : float
        x center of object in pixels
    ycentroid : float
        y center of object in pixels
    order_bounding : dict{order: tuple}
        Contains the object x,y bounding locations on the image
        keyed on spectral order
    partial_order : bool
        True if the order is only partially contained on the image
    waverange : list
        wavelength range for the order
    is_extended : bool
        True if marked as extended in the source catalog
    isophotal_abmag : float
        isophotal_abmag from the source catalog
    sky_centroid: `~astropy.coordinates.SkyCoord`
        ra and dec of the center of the object
    sky_bbox_ll : `~astropy.coordinates.SkyCoord`
        Lower left corner of the minimum bounding box
    sky_bbox_lr : `~astropy.coordinates.SkyCoord`
        Lower right corder of the minimum bounding box
    sky_bbox_ul : `~astropy.coordinates.SkyCoord`
        Upper left corner of the minimum bounding box
    sky_bbox_ur : `~astropy.coordinates.SkyCoord`
        Upper right corner of the minimum bounding box

    Notes
    -----
    The object bounding box is computed from the segmentation map,
    using the min and max wavelength for each of the orders that
    are available. The order_bounding member is a dictionary of
    bounding boxes for the object keyed by order

    ra and dec are the sky ra and dec of the center of the object as measured
    from the non-dispersed image.

    order_bounding is stored as a lookup dictionary per order and contains
    the object x,y bounding location on the grism image
    GrismObject(order_bounding={"+1":((xmin,xmax),(ymin,ymax)),"+2":((2,3),(2,3))})

    """
    __slots__ = ()  # prevent instance dictionary for lower memory

    def __new__(cls,
                sid=None,
                order_bounding={},
                sky_centroid=None,
                partial_order={},
                waverange=None,
                sky_bbox_ll=None,
                sky_bbox_lr=None,
                sky_bbox_ur=None,
                sky_bbox_ul=None,
                xcentroid=None,
                ycentroid=None,
                is_extended=None,
                isophotal_abmag=None):

        return super(GrismObject, cls).__new__(cls,
                                               sid=sid,
                                               order_bounding=order_bounding,
                                               sky_centroid=sky_centroid,
                                               partial_order=partial_order,
                                               waverange=waverange,
                                               sky_bbox_ll=sky_bbox_ll,
                                               sky_bbox_lr=sky_bbox_lr,
                                               sky_bbox_ur=sky_bbox_ur,
                                               sky_bbox_ul=sky_bbox_ul,
                                               xcentroid=xcentroid,
                                               ycentroid=ycentroid,
                                               is_extended=is_extended,
                                               isophotal_abmag=isophotal_abmag)

    def __str__(self):
        """Return a pretty print for the object information."""
        return ("id: {0}\n"
                "order_bounding {1}\n"
                "sky_centroid: {2}\n"
                "sky_bbox_ll: {3}\n"
                "sky_bbox_lr: {4}\n"
                "sky_bbox_ur: {5}\n"
                "sky_bbox_ul:{6}\n"
                "xcentroid: {7}\n"
                "ycentroid: {8}\n"
                "partial_order: {9}\n"
                "waverange: {10}\n"
                "is_extended: {11}\n"
                "isophotal_abmag: {12}\n"
                .format(self.sid,
                        str(self.order_bounding),
                        str(self.sky_centroid),
                        str(self.sky_bbox_ll),
                        str(self.sky_bbox_lr),
                        str(self.sky_bbox_ur),
                        str(self.sky_bbox_ul),
                        self.xcentroid,
                        self.ycentroid,
                        str(self.partial_order),
                        str(self.waverange),
                        str(self.is_extended),
                        str(self.isophotal_abmag)))


class MIRI_AB2Slice(Model):
    """
    MIRI MRS alpha, beta to slice transform

    Parameters
    ----------
    beta_zero : float
    beta_del : float
    """
    standard_broadcasting = False
    _separable = False
    fittable = False

    n_inputs = 1
    n_outputs = 1

    beta_zero = Parameter('beta_zero', default=0)
    """ Beta_zero parameter"""
    beta_del = Parameter('beta_del', default=1)
    """ Beta_del parameter"""
    channel = Parameter("channel", default=1)
    """ MIRI MRS channel: one of 1, 2, 3, 4"""

    def __init__(self, beta_zero=beta_zero, beta_del=beta_del, channel=channel, **kwargs):
        super().__init__(beta_zero=beta_zero, beta_del=beta_del,
                         channel=channel, **kwargs)
        self.inputs = ("beta",)
        """ "beta": the beta angle """
        self.outputs = ("slice",)
        """ "slice": Slice number"""

    @staticmethod
    def evaluate(beta, beta_zero, beta_del, channel):
        s = channel * 100 + (beta - beta_zero) / beta_del + 1
        return _toindex(s)


class RefractionIndexFromPrism(Model):
    """
    Compute the refraction index of a prism (NIRSpec).

    Parameters
    ----------
    prism_angle : float
        Prism angle in deg.

    """
    standard_broadcasting = False
    _separable = False

    n_inputs = 3
    n_outputs = 1

    prism_angle = Parameter(setter=np.deg2rad, getter=np.rad2deg)

    def __init__(self, prism_angle, name=None):
        super(RefractionIndexFromPrism, self).__init__(prism_angle=prism_angle, name=name)
        self.inputs = ("alpha_in", "beta_in", "alpha_out",)
        self.outputs = ("n",)

    def evaluate(self, alpha_in, beta_in, alpha_out, prism_angle):
        sangle = (math.sin(prism_angle))
        cangle = (math.cos(prism_angle))
        nsq = ((alpha_out + alpha_in * (1 - 2 * sangle ** 2)) / (2 * sangle * cangle)) ** 2 + \
            alpha_in ** 2 + beta_in ** 2
        return np.sqrt(nsq)


class Gwa2Slit(Model):
    """
    NIRSpec GWA to slit transform.

    Parameters
    ----------
    slits : list
        A list of open slits.
        A slit is a namedtuple of type `~jwst.transforms.models.Slit`
        Slit("name", "shutter_id", "dither_position", "xcen", "ycen", "ymin", "ymax",
        "quadrant", "source_id", "shutter_state", "source_name",
        "source_alias", "stellarity", "source_xpos", "source_ypos",
        "source_ra", "source_dec"])
    models : list
        List of models (`~astropy.modeling.core.Model`) corresponding to the
        list of slits.
    """
    _separable = False

    n_inputs = 4
    n_outputs = 4

    def __init__(self, slits, models):
        if isiterable(slits[0]):
            self._slits = [tuple(s) for s in slits]
            self.slit_ids = [s[0] for s in self._slits]
        else:
            self._slits = list(slits)
            self.slit_ids = self._slits

        self.models = models
        super(Gwa2Slit, self).__init__()
        self.inputs = ('name', 'angle1', 'angle2', 'angle3')
        """ Name of the slit and the three angle coordinates at the GWA going from detector to sky."""
        self.outputs = ('name', 'x_slit', 'y_slit', 'lam')
        """ Name of the slit, x and y coordinates within the virtual slit and wavelength."""

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
    Transform from Nirspec ``slit_frame`` to ``msa_frame``.

    Parameters
    ----------
    slits : list
        A list of open slits.
        A slit is a namedtuple, `~jwst.transforms.models.Slit`
        Slit("name", "shutter_id", "dither_position", "xcen", "ycen", "ymin", "ymax",
        "quadrant", "source_id", "shutter_state", "source_name",
        "source_alias", "stellarity", "source_xpos", "source_ypos",
        "source_ra", "source_dec")
    models : list
        List of models (`~astropy.modeling.core.Model`) corresponding to the
        list of slits.
    """
    _separable = False

    n_inputs = 3
    n_outputs = 2

    def __init__(self, slits, models):
        super(Slit2Msa, self).__init__()
        self.inputs = ('name', 'x_slit', 'y_slit')
        """ Name of the slit, x and y coordinates within the virtual slit."""
        self.outputs = ('x_msa', 'y_msa')
        """ x and y coordinates in the MSA frame."""
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

    def evaluate(self, name, x, y):
        index = self.slit_ids.index(name)
        return self.models[index](x, y)


class NirissSOSSModel(Model):
    """
    NIRISS SOSS wavelength solution model.

    Parameters
    ----------
    spectral_orders : list of int
        Spectral orders for which there is a wavelength solution.
    models : list of `~astropy.modeling.core.Model`
        A list of transforms representing the wavelength solution for
        each order in spectral orders. It should match the order in
        ``spectral_orders``.
    """

    _separable = False

    n_inputs = 3
    n_outputs = 3

    def __init__(self, spectral_orders, models):
        super(NirissSOSSModel, self).__init__()
        self.inputs = ('x', 'y', 'spectral_order')
        """ x and y pixel coordinates and spectral order"""
        self.outputs = ('ra', 'dec', 'lam')
        """ RA and DEC coordinates and wavelength"""

        self.spectral_orders = spectral_orders
        self.models = dict(zip(spectral_orders, models))

    def get_model(self, spectral_order):
        return self.models[spectral_order]

    def evaluate(self, x, y, spectral_order):

        # The spectral_order variable is coming in as an array/list of one element.
        # So, we are going to just take the 0'th element and use that as the index.
        try:
            order_number = int(spectral_order[0])
        except Exception:
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
    n_inputs = 1
    n_outputs = 1
    _separable = False

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
        self.inputs = ('x', )
        self.outputs = ('x', )

    def evaluate(self, x):
        x = x.copy()
        m = ~np.isnan(x)
        m_ind = np.flatnonzero(m)
        if isinstance(self.compareto, np.ndarray):
            cond = self.conditions[self.condition](x[m], self.compareto[m])
            x.flat[m_ind[cond]] = self.value.flat[m_ind[cond]]
        else:
            cond = self.conditions[self.condition](x[m], self.compareto)
            x.flat[m_ind[cond]] = self.value
        return x

    def __repr__(self):
        txt = "{0}(condition={1}, compareto={2}, value={3})"
        return txt.format(self.__class__.__name__, self.condition,
                          self.compareto, self.value)


class IdealToV2V3(Model):
    """
    Performs the transform from Ideal to telescope V2,V3 coordinate system.
    The two systems have the same origin: V2_REF, V3_REF.

    Note: This model has no schema implemented - add schema if needed.
    """
    _separable = False
    n_inputs = 2
    n_outputs = 2

    v3idlyangle = Parameter()  # in deg
    v2ref = Parameter()  # in arcsec
    v3ref = Parameter()  # in arcsec
    vparity = Parameter()

    def __init__(self, v3idlyangle, v2ref, v3ref, vparity, name='idl2V', **kwargs):
        super(IdealToV2V3, self).__init__(v3idlyangle=v3idlyangle, v2ref=v2ref,
                                          v3ref=v3ref, vparity=vparity, name=name,
                                          **kwargs)
        self.inputs = ('xidl', 'yidl')
        """ x and y coordinates in the telescope Ideal frame."""
        self.outputs = ('v2', 'v3')
        """ coordinates in the telescope (V2,V3) frame."""

    @staticmethod
    def evaluate(xidl, yidl, v3idlyangle, v2ref, v3ref, vparity):
        """
        Parameters
        ----------
        xidl, yidl : ndarray-like
            Coordinates in Ideal System [in arcsec]
        v3idlyangle : float
            Angle between Ideal Y-axis and V3 [ in deg]
        v2ref, v3ref : ndarray-like
            Coordinates in V2, V3 [in arcsec]
        vparity : int
            Parity.

        Returns
        -------
        v2, v3 : ndarray-like
            Coordinates in the (V2, V3) telescope system [in arcsec].

        """
        v3idlyangle = np.deg2rad(v3idlyangle)

        v2 = v2ref + vparity * xidl * np.cos(v3idlyangle) + yidl * np.sin(v3idlyangle)
        v3 = v3ref - vparity * xidl * np.sin(v3idlyangle) + yidl * np.cos(v3idlyangle)
        return v2, v3

    def inverse(self):
        return V2V3ToIdeal(self.v3idlyangle, self.v2ref, self.v3ref, self.vparity)


class V2V3ToIdeal(Model):
    """
    Performs the transform from telescope V2,V3 to Ideal coordinate system.
    The two systems have the same origin - V2_REF, V3_REF.

    Note: This model has no schema implemented - add if needed.
    """
    _separable = False

    n_inputs = 2
    n_outputs = 2

    v3idlyangle = Parameter()  # in deg
    v2ref = Parameter()  # in arcsec
    v3ref = Parameter()  # in arcsec
    vparity = Parameter()

    def __init__(self, v3idlyangle, v2ref, v3ref, vparity, name='V2idl', **kwargs):
        super(V2V3ToIdeal, self).__init__(v3idlyangle=v3idlyangle, v2ref=v2ref,
                                          v3ref=v3ref, vparity=vparity, name=name,
                                          **kwargs)
        self.inputs = ('v2', 'v3')
        """ ('v2', 'v3'): coordinates in the telescope (V2,V3) frame."""
        self.outputs = ('xidl', 'yidl')
        """ ('xidl', 'yidl'): x and y coordinates in the telescope Ideal frame."""

    @staticmethod
    def evaluate(v2, v3, v3idlyangle, v2ref, v3ref, vparity):
        """
        Parameters
        ----------
        xidl, yidl : ndarray-like
            Coordinates in Ideal System [in arcsec]
        v3idlyangle : float
            Angle between Ideal Y-axis and V3 [ in deg]
        v2ref, v3ref : ndarray-like
            Coordinates in V2, V3 [in arcsec]
        vparity : int
            Parity.

        Returns
        -------
        xidl, yidl : ndarray-like
            Coordinates in the Ideal telescope system [in arcsec].

        """
        v3idlyangle = np.deg2rad(v3idlyangle)

        xidl = vparity * ((v2 - v2ref) * np.cos(v3idlyangle) -
                          (v3 - v3ref) * np.sin(v3idlyangle))
        yidl = ((v2 - v2ref) * np.sin(v3idlyangle) +
                (v3 - v3ref) * np.cos(v3idlyangle))

        return xidl, yidl

    def inverse(self):
        return IdealToV2V3(self.v3idlyangle, self.v2ref, self.v3ref, self.vparity)


def _toindex(value):
    """
    Convert value to an int or an int array.

    Input coordinates converted to integers
    corresponding to the center of the pixel.
    The convention is that the center of the pixel is
    (0, 0), while the lower left corner is (-0.5, -0.5).

    Examples
    --------
    >>> from jwst.transforms.models import _toindex
    >>> _toindex(np.array([-0.5, 0.49999]))
    array([0, 0])
    >>> _toindex(np.array([0.5, 1.49999]))
    array([1, 1])
    >>> _toindex(np.array([1.5, 2.49999]))
    array([2, 2])
    """
    indx = np.asarray(np.floor(np.asarray(value) + 0.5), dtype=int)
    return indx


class NIRCAMForwardRowGrismDispersion(Model):
    """Return the transform from grism to image for the given spectral order.

    Parameters
    ----------
    orders : list [int]
        List of orders which are available

    lmodels : list [astropy.modeling.Model]
        List of models which govern the wavelength solutions for each order

    xmodels : list [astropy.modeling.Model]
        List of models which govern the x solutions for each order

    ymodels : list [astropy.modeling.Model]
        List of models which givern the y solutions for each order

    Returns
    -------
    x, y, wavelength, order in the grism image for the pixel at x0,y0 that was
    specified as input using the input delta pix for the specified order

    Notes
    -----
    The evaluation here is linear currently because higher orders have not yet been
    defined for NIRCAM (NIRCAM polynomials currently do not have any field
    dependence)
    """
    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 5
    n_outputs = 4

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, name=None, meta=None):
        self.orders = orders
        self.lmodels = lmodels
        self.xmodels = xmodels
        self.ymodels = ymodels
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        meta = {"orders": orders}  # informational for users
        if name is None:
            name = 'nircam_forward_row_grism_dispersion'
        super(NIRCAMForwardRowGrismDispersion, self).__init__(name=name,
                                                              meta=meta)
        self.inputs = ("x", "y", "x0", "y0", "order")
        self.outputs = ("x", "y", "wavelength", "order")

    def evaluate(self, x, y, x0, y0, order):
        """Return the transform from grism to image for the given spectral order.

        Parameters
        ----------
        x : float
            input x pixel
        y : float
            input y pixel
        x0 : float
            input x-center of object
        y0 : float
            input y-center of object
        order : int
            the spectral order to use
        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]

        # inputs are x, y, x0, y0, order

        tmodel = astmath.SubtractUfunc() | xmodel
        model = Mapping((0, 2, 0, 2, 2, 3, 4)) | (tmodel | ymodel) & (tmodel | lmodel) & Identity(3) | \
            Mapping((2, 3, 0, 1, 4)) | Identity(1) & astmath.AddUfunc() & Identity(2) | Mapping((0, 1, 2, 3), n_inputs=4)

        return model(x, y, x0, y0, order)


class NIRCAMForwardColumnGrismDispersion(Model):
    """Return the transform from grism to image for the given spectral order.

    Parameters
    ----------
    orders : list [int]
        List of orders which are available

    lmodels : list [astropy.modeling.Model]
        List of models which govern the wavelength solutions

    xmodels : list [astropy.modeling.Model]
        List of models which govern the x solutions

    ymodels : list [astropy.modeling.Model]
        List of models which givern the y solutions

    Returns
    -------
    x, y, lam, order in the grism image for the pixel at x0,y0 that was
    specified as input using the input delta pix for the specified order

    Notes
    -----
    The evaluation here is linear because higher orders have not yet been
    defined for NIRCAM (NIRCAM polynomials currently do not have any field
    dependence)
    """
    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 5
    n_outputs = 4

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, name=None, meta=None):
        self.orders = orders
        self.lmodels = lmodels
        self.xmodels = xmodels
        self.ymodels = ymodels
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        meta = {"orders": orders}  # informational for users
        if name is None:
            name = 'nircam_forward_column_grism_dispersion'
        super(NIRCAMForwardColumnGrismDispersion, self).__init__(name=name,
                                                                 meta=meta)
        self.inputs = ("x", "y", "x0", "y0", "order")
        self.outputs = ("x", "y", "wavelength", "order")

    def evaluate(self, x, y, x0, y0, order):
        """Return the transform from grism to image for the given spectral order.

        Parameters
        ----------
        x : float
            input x pixel
        y : float
            input y pixel
        x0 : float
            input x-center of object
        y0 : float
            input y-center of object
        order : int
            the spectral order to use
        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]

        # inputs are x, y, x0, y0, order
        tmodel = astmath.SubtractUfunc() | ymodel
        dx = tmodel | xmodel
        wavelength = tmodel | lmodel
        model = Mapping((1, 3, 1, 3, 2, 3, 4)) | \
            dx & wavelength & Identity(3) |\
            Mapping((0, 2, 3, 1, 4)) | astmath.AddUfunc() & Identity(3)

        return model(x, y, x0, y0, order)


class NIRCAMBackwardGrismDispersion(Model):
    """Return the valid pixel(s) and wavelengths given center x,y and lam

    Parameters
    ----------
    orders : list [int]
        List of orders which are available

    lmodels : list [astropy.modeling.Model]
        List of models which govern the wavelength solutions

    xmodels : list [astropy.modeling.Model]
        List of models which govern the x solutions

    ymodels : list [astropy.modeling.Model]
        List of models which givern the y solutions

    Returns
    -------
    x, y, lam, order in the grism image for the pixel at x0,y0 that was
    specified as input using the wavelength l for the specified order

    Notes
    -----
    The evaluation here is linear because higher orders have not yet been defined for NIRCAM
    (NIRCAM polynomials currently do not have any field dependence)
    """
    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 4
    n_outputs = 5

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, name=None, meta=None):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.lmodels = lmodels
        self.xmodels = xmodels
        self.ymodels = ymodels
        self.orders = orders
        meta = {"orders": orders}
        if name is None:
            name = "nircam_backward_grism_dispersion"
        super(NIRCAMBackwardGrismDispersion, self).__init__(name=name,
                                                            meta=meta)
        self.inputs = ("x", "y", "wavelength", "order")
        self.outputs = ("x", "y", "x0", "y0", "order")

    def evaluate(self, x, y, wavelength, order):
        """Return the transform from image to grism for the given spectral order.

        Parameters
        ----------
        x : float
            input x pixel
        y : float
            input y pixel
        wavelength : float
            input wavelength in angstroms
        order : int
            specifies the spectral order
        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        if (wavelength < 0).any():
            raise ValueError("wavelength should be greater than zero")

        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]

        dx = lmodel | xmodel
        dy = lmodel | ymodel
        model = Mapping((0, 2, 1, 2, 0, 1, 3)) | \
            ((Identity(1) & dx) | astmath.AddUfunc()) & \
            ((Identity(1) & dy) | astmath.AddUfunc()) & Identity(3)

        return model(x, y, wavelength, order)


class NIRISSBackwardGrismDispersion(Model):
    """This model calculates the dispersion extent of NIRISS pixels.

    The dispersion is relative to the input x,y for a given wavelength.

    Parameters
    ----------
    xmodels : list[tuple]
        The list of tuple(models) for the polynomial model in x
    ymodels : list[tuple]
        The list of tuple(models) for the polynomial model in y
    lmodels : list
        The list of models for the polynomial model in l
    orders : list
        The list of orders which are available to the model
    theta : float
        Angle [deg] - defines the NIRISS filter wheel position

    Notes
    -----
    Given the x,y, wave, order as known on the direct image,
    it returns the tuple of x, y, wave, order for that wave in the dispersed image.

    This model needs to be generalized, at the moment it satisfies the
    2t x 6(xy)th order polynomial currently used by NIRISS.

    """

    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 4
    n_outputs = 5

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, theta=0.0, name=None, meta=None):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.xmodels = xmodels
        self.ymodels = ymodels
        self.lmodels = lmodels
        self.orders = orders
        self.theta = theta
        meta = {"orders": orders}
        if name is None:
            name = 'niriss_backward_grism_dispersion'
        super(NIRISSBackwardGrismDispersion, self).__init__(name=name,
                                                            meta=meta)
        self.inputs = ("x", "y", "wavelength", "order")
        self.outputs = ("x", "y", "x0", "y0", "order")

    def evaluate(self, x, y, wavelength, order):
        """Return the valid pixel(s) and wavelengths given center x,y and lam

        Parameters
        ----------
        wavelength : float
            Input wavelength you want to know about, will be converted to float
        x : float
            Input x location
        y : float
            Input y location
        wavelength : float
            Wavelength to disperse
        order : list
            The order to use


        Returns
        -------
        x, y, wavelength, order in the grism image for the pixel at x,y that was
        specified as input using the wavelength and order specified

        Notes
        -----
        There's spatial dependence for NIRISS so the forward transform
        depends on x,y as well as the filter wheel rotation. Theta is
        usu. taken to be the different between fwcpos_ref in the specwcs
        reference file and fwcpos from the input image.

        """
        if (wavelength < 0).any():
            raise ValueError("Wavelength should be greater than zero")
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        t = self.lmodels[iorder](wavelength)
        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]

        dx = xmodel[0](x, y) + t * xmodel[1](x, y) + t**2 * xmodel[2](x, y)
        dy = ymodel[0](x, y) + t * ymodel[1](x, y) + t**2 * ymodel[2](x, y)

        # rotate by theta
        if self.theta != 0.0:
            rotate = Rotation2D(self.theta)
            dx, dy = rotate(dx, dy)

        return x + dx, y + dy, x, y, order


class NIRISSForwardRowGrismDispersion(Model):
    """This model calculates the wavelengths of vertically dispersed NIRISS grism data.

    The dispersion polynomial is relative to the input x,y pixels
    in the direct image for a given wavelength.

    Parameters
    ----------
    orders : list
        The list of orders which are available to the model
    xmodels : list[tuples]
        The list of tuple(models) for the polynomial model in x
    ymodels : list[tuples]
        The list of tuple(models) for the polynomial model in y
    lmodels : list
        The list of models for the polynomial model in l
    theta : float
        Angle [deg] - defines the NIRISS filter wheel position

    Notes
    -----
    Given the x,y, source location as known on the dispersed image, as well as order,
    it returns the tuple of x,y,wavelength,order.

    This model needs to be generalized, at the moment it satisfies the
    2t x 6(xy)th order polynomial currently used by NIRISS.

    """

    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 5
    n_outputs = 4

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, theta=0., name=None, meta=None):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.xmodels = xmodels
        self.ymodels = ymodels
        self.lmodels = lmodels
        self.theta = theta
        self.orders = orders
        meta = {"orders": orders}
        if name is None:
            name = 'niriss_forward_row_grism_dispersion'
        super(NIRISSForwardRowGrismDispersion, self).__init__(name=name,
                                                              meta=meta)
        # starts with the backwards pixel and calculates the forward pixel
        self.inputs = ("x", "y", "x0", "y0", "order")
        self.outputs = ("x", "y", "wavelength", "order")

    def evaluate(self, x, y, x0, y0, order):
        """Return the valid pixel(s) and wavelengths given center x,y and lam

        Parameters
        ----------
        x0: int,float,list
            Source object x-center

        y0: int,float,list
            Source object y-center

        x :  int,float,list
            Input x location

        y :  int,float,list
            Input y location

        order : int
            Spectral order to use


        Returns
        -------
        x, y, lambda, order, theta,  in the direct image for the pixel that was
        specified as input using the wavelength l and spectral order

        Notes
        -----
        There's spatial dependence for NIRISS as well as dependence on the
        filter wheel rotation during the exposure.

        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        # The next two lines are to get around the fact that
        # modeling.standard_broadcasting=False does not work.
        x00 = x0.flatten()[0]
        y00 = y0.flatten()[0]

        t = np.linspace(0, 1, 10)  # sample t
        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]

        dx = xmodel[0](x00, y00) + t * xmodel[1](x00, y00) + t**2 * xmodel[2](x00, y00)
        dy = ymodel[0](x00, y00) + t * ymodel[1](x00, y00) + t**2 * ymodel[2](x00, y00)

        if self.theta != 0.0:
            rotate = Rotation2D(self.theta)
            dx, dy = rotate(dx, dy)

        so = np.argsort(dx)
        tab = Tabular1D(dx[so], t[so], bounds_error=False, fill_value=None)

        dxr = astmath.SubtractUfunc()
        wavelength = dxr | tab | lmodel
        model = Mapping((2, 3, 0, 2, 4)) | Const1D(x00) & Const1D(y00) & wavelength & Const1D(order)
        # returns x0, y0, lam, order
        return model(x, y, x0, y0, order)


class NIRISSForwardColumnGrismDispersion(Model):
    """This model calculates the wavelengths for horizontally dispersed NIRISS grism data.

    The dispersion polynomial is relative to the input x,y pixels
    in the direct image for a given wavelength.

    Parameters
    ----------
    orders : list
        The list of orders which are available to the model.
    xmodels : list[tuple]
        The list of tuple(models) for the polynomial model in x
    ymodels : list[tuple]
        The list of tuple(models) for the polynomial model in y
    lmodels : list
        The list of models for the polynomial model in wavelength.
    theta : float
        Angle [deg] - defines the NIRISS filter wheel position

    Notes
    -----
    Given the x,y, source location, order, it returns the tuple of
    x,y,wavelength,order on the dispersed image. It also requires
    FWCPOS from the image header, this is the filter wheel position
    in degrees.

    """

    standard_broadcasting = False
    _separable = False
    fittable = False
    linear = False

    n_inputs = 5
    n_outputs = 4

    def __init__(self, orders, lmodels=None, xmodels=None,
                 ymodels=None, theta=None, name=None, meta=None):
        self._order_mapping = {int(k): v for v, k in enumerate(orders)}
        self.xmodels = xmodels
        self.ymodels = ymodels
        self.lmodels = lmodels
        self.orders = orders
        self.theta = theta
        meta = {"orders": orders}
        if name is None:
            name = 'niriss_forward_column_grism_dispersion'
        super(NIRISSForwardColumnGrismDispersion, self).__init__(name=name,
                                                                 meta=meta)
        # starts with the backwards pixel and calculates the forward pixel
        self.inputs = ("x", "y", "x0", "y0", "order")
        self.outputs = ("x", "y", "wavelength", "order")

    def evaluate(self, x, y, x0, y0, order):
        """Return the valid pixel(s) and wavelengths given center x,y and lam

        Parameters
        ----------
        x0: int,float
            Source object x-center
        y0: int,float
            Source object y-center
        x :  int,float
            Input x location
        y :  int,float
            Input y location
        order : int
            Spectral order to use
        theta : float
            input filter wheel rotation angle in degrees

        Returns
        -------
        x, y, lambda, order,  in the direct image for the pixel that was
        specified as input using the wavelength l and spectral order

        Notes
        -----
        There's spatial dependence for NIRISS as well as rotation for the filter wheel

        """
        try:
            iorder = self._order_mapping[int(order.flatten()[0])]
        except KeyError:
            raise ValueError("Specified order is not available")

        # The next two lines are to get around the fact that
        # modeling.standard_broadcasting=False does not work.
        x00 = x0.flatten()[0]
        y00 = y0.flatten()[0]

        t = np.linspace(0, 1, 10)
        xmodel = self.xmodels[iorder]
        ymodel = self.ymodels[iorder]
        lmodel = self.lmodels[iorder]
        dx = xmodel[0](x00, y00) + t * xmodel[1](x00, y00) + t**2 * xmodel[2](x00, y00)
        dy = ymodel[0](x00, y00) + t * ymodel[1](x00, y00) + t**2 * ymodel[2](x00, y00)

        if self.theta != 0.0:
            rotate = Rotation2D(self.theta)
            dx, dy = rotate(dx, dy)
        so = np.argsort(dy)
        tab = Tabular1D(dy[so], t[so], bounds_error=False, fill_value=None)
        dyr = astmath.SubtractUfunc()
        wavelength = dyr | tab | lmodel
        model = Mapping((2, 3, 1, 3, 4)) | Const1D(x00) & Const1D(y00) & wavelength & Const1D(order)
        return model(x, y, x0, y0, order)


class Rotation3DToGWA(Model):
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
    _separable = False

    separable = False

    n_inputs = 3
    n_outputs = 3

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
        self.inputs = ('x', 'y', 'z')
        self.outputs = ('x', 'y', 'z')

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


class Snell(Model):
    """
    Apply transforms, including Snell law, through the NIRSpec prism.
    Parameters
    ----------
    angle : float
        Prism angle in deg.
    kcoef : list
        K coefficients in Sellmeier equation.
    lcoef : list
        L coefficients in Sellmeier equation.
    tcoef : list
        Thermal coefficients of glass.
    tref : float
        Reference temperature in K.
    pref : float
        Reference pressure in ATM.
    temperature : float
        System temperature during observation in K
    pressure : float
        System pressure during observation in ATM.
    """

    standard_broadcasting = False
    _separable = False

    n_inputs = 4
    n_outputs = 3

    def __init__(self, angle, kcoef, lcoef, tcoef, tref, pref,
                 temperature, pressure, name=None):
        self.prism_angle = angle
        self.kcoef = np.array(kcoef, dtype=float)
        self.lcoef = np.array(lcoef, dtype=float)
        self.tcoef = np.array(tcoef, dtype=float)
        self.tref = tref
        self.pref = pref
        self.temp = temperature
        self.pressure = pref
        super(Snell, self).__init__(name=name)
        self.inputs = ("lam", "alpha_in", "beta_in", "zin")
        self.outputs = ("alpha_out", "beta_out", "zout")

    @staticmethod
    def compute_refraction_index(lam, temp, tref, pref, pressure, kcoef, lcoef, tcoef):
        """Calculate and return the refraction index."""

        # Convert to microns
        lam = np.asarray(lam * 1e6)
        KtoC = 273.15  # kelvin to celsius conversion
        temp -= KtoC
        tref -= KtoC
        delt = temp - tref

        K1, K2, K3 = kcoef
        L1, L2, L3 = lcoef
        D0, D1, D2, E0, E1, lam_tk = tcoef

        if delt < 20:
            n = np.sqrt(1. +
                        K1 * lam**2 / (lam**2 - L1) +
                        K2 * lam**2 / (lam**2 - L2) +
                        K3 * lam**2 / (lam**2 - L3)
                        )
        else:
            # Derive the refractive index of air at the reference temperature and pressure
            # and at the operational system's temperature and pressure.
            nref = 1. + (6432.8 + 2949810. * lam**2 /
                         (146.0 * lam**2 - 1.) + (5540.0 * lam**2) /
                         (41.0 * lam**2 - 1.)) * 1e-8

            # T should be in C, P should be in ATM
            nair_obs = 1.0 + ((nref - 1.0) * pressure) / (1.0 + (temp - 15.) * 3.4785e-3)
            nair_ref = 1.0 + ((nref - 1.0) * pref) / (1.0 + (tref - 15) * 3.4785e-3)

            # Compute the relative index of the glass at Tref and Pref using Sellmeier equation I.
            lamrel = lam * nair_obs / nair_ref

            nrel = np.sqrt(1. +
                           K1 * lamrel ** 2 / (lamrel ** 2 - L1) +
                           K2 * lamrel ** 2 / (lamrel ** 2 - L2) +
                           K3 * lamrel ** 2 / (lamrel ** 2 - L3)
                           )
            # Convert the relative index of refraction at the reference temperature and pressure
            # to absolute.
            nabs_ref = nrel * nair_ref

            # Compute the absolute index of the glass
            delnabs = (0.5 * (nrel ** 2 - 1.) / nrel) * \
                (D0 * delt + D1 * delt ** 2 + D2 * delt ** 3 +
                 (E0 * delt + E1 * delt ** 2) / (lamrel ** 2 - lam_tk ** 2))
            nabs_obs = nabs_ref + delnabs

            # Define the relative index at the system's operating T and P.
            n = nabs_obs / nair_obs
        return n

    def evaluate(self, lam, alpha_in, beta_in, zin):
        """Go through the prism"""
        n = self.compute_refraction_index(lam, self.temp, self.tref, self.pref, self.pressure,
                                          self.kcoef, self.lcoef, self.tcoef)
        # Apply Snell's law through front surface, eq 5.3.3 II
        xout = alpha_in / n
        yout = beta_in / n
        zout = np.sqrt(1.0 - xout**2 - yout**2)

        # Go to back surface frame # eq 5.3.3 III
        y_rotation = Rotation3DToGWA([self.prism_angle], "y")
        xout, yout, zout = y_rotation(xout, yout, zout)

        # Reflection on back surface
        xout = -1 * xout
        yout = -1 * yout

        # Back to front surface
        y_rotation = Rotation3DToGWA([-self.prism_angle], "y")
        xout, yout, zout = y_rotation(xout, yout, zout)

        # Snell's refraction law through front surface
        xout = xout * n
        yout = yout * n
        zout = np.sqrt(1.0 - xout**2 - yout**2)
        return xout, yout, zout


class AngleFromGratingEquation(Model):
    """
    Solve the 3D Grating Dispersion Law for the refracted angle.

    Parameters
    ----------
    groove_density : int
        Grating ruling density.
    order : int
        Spectral order.
    """

    _separable = False
    n_inputs = 4
    n_outputs = 3

    groove_density = Parameter()
    """ Grating ruling density."""

    order = Parameter(default=-1)
    """ Spectral order."""

    def __init__(self, groove_density, order, **kwargs):
        super().__init__(groove_density=groove_density, order=order, **kwargs)
        self.inputs = ("lam", "alpha_in", "beta_in", "z")
        """ Wavelength and 3 angle coordinates going into the grating."""

        self.outputs = ("alpha_out", "beta_out", "zout")
        """ Three angles coming out of the grating. """

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
    """
    Solve the 3D Grating Dispersion Law for the wavelength.

    Parameters
    ----------
    groove_density : int
        Grating ruling density.
    order : int
        Spectral order.
    """

    _separable = False
    n_inputs = 3
    n_outputs = 1

    groove_density = Parameter()
    """ Grating ruling density."""
    order = Parameter(default=1)
    """ Spectral order."""

    def __init__(self, groove_density, order, **kwargs):
        super().__init__(groove_density=groove_density, order=order, **kwargs)
        self.inputs = ("alpha_in", "beta_in", "alpha_out")
        """ three angle - alpha_in and beta_in going into the grating and alpha_out coming out of the grating."""
        self.outputs = ("lam",)
        """ Wavelength."""

    def evaluate(self, alpha_in, beta_in, alpha_out, groove_density, order):
        # beta_in is not used in this equation but is here because it's
        # needed for the prism computation. Currently these two computations
        # need to have the same interface.
        return -(alpha_in + alpha_out) / (groove_density * order)


class Rotation3D(Model):
    """
    Perform a series of rotations about different axis in 3D space.

    Positive angles represent a counter-clockwise rotation.

    Parameters
    ----------
    angles : array-like
        Angles of rotation in deg in the order of axes_order.
    axes_order : str
        A sequence of 'x', 'y', 'z' corresponding of axis of rotation.
    """
    standard_broadcasting = False
    _separable = False

    n_inputs = 3
    n_outputs = 3

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
        self.inputs = ('x', 'y', 'z')
        self.outputs = ('x', 'y', 'z')

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
            matrix = np.zeros((3, 3), dtype=float)
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
                raise ValueError(f"""Expected axes_order to be a combination
                                 of characters 'x', 'y' and 'z',
                                 got {set(axes_order).difference(['x', 'y', 'z'])}""")
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


class V23ToSky(Rotation3D):
    """
    Transform from V2V3 to a standard coordinate system (ICRS).

    Parameters
    ----------
    angles : list
        A sequence of angles (in deg).
        The angles are [-V2_REF, V3_REF, -ROLL_REF, -DEC_REF, RA_REF].
    axes_order : str
        A sequence of characters ('x', 'y', or 'z') corresponding to the
        axis of rotation and matching the order in ``angles``.
        The axes are "zyxyz".
    """

    _separable = False

    n_inputs = 2
    n_outputs = 2

    def __init__(self, angles, axes_order, name=None):
        super(V23ToSky, self).__init__(angles, axes_order=axes_order, name=name)
        self._inputs = ("v2", "v3")
        """ ("v2", "v3"): Coordinates in the (V2, V3) telescope frame."""
        self._outputs = ("ra", "dec")
        """ ("ra", "dec"): RA, DEC coordinates in ICRS."""

    @property
    def inputs(self):
        return self._inputs

    @inputs.setter
    def inputs(self, val):
        self._inputs = val

    @property
    def outputs(self):
        return self._outputs

    @outputs.setter
    def outputs(self, val):
        self._outputs = val

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
        alpha = np.rad2deg(np.arctan2(y, x))
        delta = np.rad2deg(np.arctan2(z, h))
        return alpha, delta

    def evaluate(self, v2, v3, angles):
        x, y, z = self.spherical2cartesian(v2, v3)
        x1, y1, z1 = super(V23ToSky, self).evaluate(x, y, z, angles)
        ra, dec = self.cartesian2spherical(x1, y1, z1)

        return ra, dec

    def __call__(self, v2, v3, **kwargs):
        from itertools import chain
        inputs, format_info = self.prepare_inputs(v2, v3)
        parameters = self._param_sets(raw=True)

        outputs = self.evaluate(*chain(inputs, parameters))
        if self.n_outputs == 1:
            outputs = (outputs,)

        return self.prepare_outputs(format_info, *outputs)


class Unitless2DirCos(Model):
    """
    Transform a vector to directional cosines.
    """
    _separable = False
    n_inputs = 2
    n_outputs = 3

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inputs = ('x', 'y')
        self.outputs = ('x', 'y', 'z')

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
    Transform directional cosines to vector.
    """
    _separable = False
    n_inputs = 3
    n_outputs = 2

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inputs = ('x', 'y', 'z')
        self.outputs = ('x', 'y')

    def evaluate(self, x, y, z):

        return x / z, y / z

    def inverse(self):
        return Unitless2DirCos()
