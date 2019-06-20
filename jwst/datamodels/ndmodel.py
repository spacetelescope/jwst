"""
Subclass of NDDataBase to support DataModel compatibility with NDData
"""

import os.path
import numpy as np
import collections

from astropy.units import Quantity
from astropy.nddata import nddata_base

from . import util
from . import filetype
from . import properties

#---------------------------------------
# astropy.io.registry compatibility
#---------------------------------------

def identify(origin, path, fileobj, *args, **kwargs):
    """
    Identify if file is a DataModel for astropy.io.registry
    """
    if fileobj:
        file_type = filetype.check(fileobj)
    elif path:
        if os.path.isfile(path):
            file_type = filetype.check(path)
        else:
            file_type = path.lower().split(".")[1]
    else:
        file_type = None

    flag = file_type and (file_type == "asdf" or file_type == "fits")
    return flag

def read(data, *args, **kwargs):
    """
    Astropy.io registry compatibility function to wrap util.open
    """

    # Translate keyword arguments to those expected by ImageModel
    xargs = {}
    if kwargs.get("mask"):
        xargs["dq"] = kwargs["mask"]

    uncertainty = kwargs.get("uncertainty")
    if uncertainty:
        if isinstance(uncertainty, Quantity):
            uncertainty_type = uncertainty.unit
            uncertainty = uncertainty.data
        else:
            uncertainty_type = None
        xargs["err"] = uncertainty
    else:
        uncertainty_type = None

    if hasattr(data, 'mask') and hasattr(data, 'data'):
        xargs["dq"] = data.mask
        data = data.data

    if isinstance(data, Quantity):
        unit = data.unit
        data = data.value
    else:
        unit = kwargs.get("unit")

    # Create the model using the transformed arguments
    model = util.open(data, **xargs)

    # Add attributes passed as keyword arguments to model
    if unit:
        model.meta.bunit_data = unit

    wcs = kwargs.get("wcs")
    if wcs:
        model.set_fits_wcs(wcs)

    if uncertainty_type:
        model.meta.bunit_err = uncertainty_type

    return model

def write(data, path, *args, **kwargs):
    """
    Astropy.io registry compatabilty function to wrap datamodel.savw
    """
    from .model_base import DataModel

    if not isinstance(data, DataModel):
        model = DataModel(data)
    else:
        model = data

    if isinstance(path, str):
        model.save(path, *args, **kwargs)
    else:
        raise ValueError("Path to write DataModel was not found")

#---------------------------------------
# Astropy NDData compatibility
#---------------------------------------

class NDModel(nddata_base.NDDataBase):
    def my_attribute(self, attr):
        """
        Test if attribute is part of the NDData interface
        """
        properties = frozenset(("data", "mask", "unit", "wcs", "unceratainty"))
        return attr in properties

    @property
    def data(self):
        """
        Read the stored dataset.
        """
        primary_array_name = self.get_primary_array_name()
        if primary_array_name:
            primary_array = self.__getattr__(primary_array_name)
        else:
            raise AttributeError("No attribute 'data'")
        return primary_array

    @data.setter
    def data(self, value):
        """
        Write the stored dataset.
        """
        primary_array_name = self.get_primary_array_name()
        if not primary_array_name:
            primary_array_name = 'data'
        properties.ObjectNode.__setattr__(self, primary_array_name, value)

    @property
    def mask(self):
        """
        Read the mask for the dataset.
        """
        return self.__getattr__('dq')

    @mask.setter
    def mask(self, value):
        """
        Write the mask for the dataset.
        """
        properties.ObjectNode.__setattr__(self, 'dq', value)

    @property
    def unit(self):
        """
        Read the units for the dataset.
        """
        try:
            val = self.meta.bunit_data
        except AttributeError:
            val = None
        return val

    @unit.setter
    def unit(self, value):
        """
        Write the units for the dataset.
        """
        self.meta.bunit_data = value

    @property
    def wcs(self):
        """
        Read the world coordinate system (WCS) for the dataset.
        """
        return self.get_fits_wcs()

    @wcs.setter
    def wcs(self, value):
        """
        Write the world coordinate system (WCS) to the dataset.
        """
        return self.set_fits_wcs(value)

    @property
    def meta(self):
        """
        Read additional meta information about the dataset.
        """
        return self.__getattr__('meta')


    @property
    def uncertainty(self):
        """
        Read the uncertainty in the dataset.
        """
        err = self.err
        try:
            val = self.meta.bunit_err
        except AttributeError:
            val = None
        return Uncertainty(err, uncertainty_type=val)

    @uncertainty.setter
    def uncertainty(self, value):
        """
        Write the uncertainty in the dataset.
        """
        properties.ObjectNode.__setattr__(self, 'err', value)
        if hasattr(value, 'uncertainty_type'):
            self.meta.bunit_err = value.uncertainty_type

#---------------------------------------------
# The following classes provide support
# for the NDData interface to Datamodels
#---------------------------------------------

class MetaNode(properties.ObjectNode, collections.abc.MutableMapping):
    """
    NDData compatibility class for meta node
    """
    def __init__(self, name, instance, schema, ctx):
        properties.ObjectNode.__init__(self, name, instance, schema, ctx)

    def _find(self, path):
        if not path:
            return self

        cursor = self._instance
        schema = self._schema
        for attr in path:
            try:
                cursor = cursor[attr]
            except KeyError:
                raise KeyError("'%s'" % '.'.join(path))
            schema = properties._get_schema_for_property(schema, attr)

        key = '.'.join(path)
        return properties._make_node(key, cursor, schema, self._ctx)

    def __delitem__(self, key):
        path = key.split('.')
        parent = self._find(path[:-1])
        try:
            parent.__delattr__(path[-1])
        except KeyError:
            raise KeyError("'%s'" % key)

    def __getitem__(self, key):
        path = key.split('.')
        return self._find(path)

    def __len__(self):
        def recurse(val):
            n = 0
            for subval in val.values():
                if isinstance(subval, dict):
                    n += recurse(subval)
                else:
                    n += 1
            return n

        return recurse(self._instance)

    def __setitem__(self, key, value):
        path = key.split('.')
        parent = self._find(path[:-1])
        try:
            parent.__setattr__(path[-1], value)
        except KeyError:
            raise KeyError("'%s'" % key)

class Uncertainty(np.ndarray):
    """
    Subclass ndarray to include an additional property, uncertainty_type
    """
    def __new__(cls, err, uncertainty_type=None):
        # info on how to subclass np.ndarray is at
        # https://docs.scipy.org/doc/numpy/user/basics.subclassing.html
        # this code is taken from there
        obj = np.asarray(err).view(cls)
        obj.uncertainty_type = uncertainty_type
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.uncertainty_type = getattr(obj, 'uncertainty_type', None)
