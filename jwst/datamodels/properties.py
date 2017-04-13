# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

import copy
import numpy as np
import jsonschema
import collections
import warnings

from astropy.extern import six
from astropy.utils.compat.misc import override__dir__

from asdf import schema
from asdf import yamlutil
from asdf.tags.core import ndarray
from asdf import tagged

from . import util

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())

__all__ = ['ObjectNode', 'ListNode']

def _cast(val, schema):
    val = _unmake_node(val)
    if val is not None:
        if 'datatype' in schema:
            # Handle lazy array
            if isinstance(val, ndarray.NDArrayType):
                val = val._make_array()
            val = util.gentle_asarray(
                val, ndarray.asdf_datatype_to_numpy_dtype(schema['datatype']))
        if 'ndim' in schema and len(val.shape) != schema['ndim']:
            raise ValueError(
                "Array has wrong number of dimensions.  Expected {0}, got {1}".format(
                    schema['ndim'], len(val.shape)))
        if 'max_ndim' in schema and len(val.shape) > schema['max_ndim']:
            raise ValueError(
                "Array has wrong number of dimensions.  Expected <= {0}, got {1}".format(
                    schema['max_ndim'], len(val.shape)))
        tag = schema.get('tag')
        if tag is not None:
            val = tagged.tag_object(tag, val)

    return val


def _make_default_array(attr, schema, ctx):
    dtype = schema.get('datatype')
    if dtype is not None:
        dtype = ndarray.asdf_datatype_to_numpy_dtype(dtype)
    ndim = schema.get('ndim', schema.get('max_ndim'))
    default = schema.get('default', None)
    primary_array_name = ctx.get_primary_array_name()

    if attr == primary_array_name:
        if ctx.shape is not None:
            shape = ctx.shape
        elif ndim is not None:
            shape = tuple([0] * ndim)
        else:
            shape = (0,)
    else:
        if dtype.names is not None:
            if ndim is None:
                shape = (0,)
            else:
                shape = tuple([0] * ndim)
            default = None
        else:
            has_primary_array_shape = False
            if primary_array_name is not None:
                primary_array = getattr(ctx, primary_array_name, None)
                has_primary_array_shape = primary_array is not None

            if has_primary_array_shape:
                if ndim is None:
                    shape = primary_array.shape
                else:
                    shape = primary_array.shape[-ndim:]
            elif ndim is None:
                shape = (0,)
            else:
                shape = tuple([0] * ndim)

    array = np.empty(shape, dtype=dtype)
    if default is not None:
        array[...] = default
    return array


def _make_default(attr, schema, ctx):
    if 'max_ndim' in schema or 'ndim' in schema or 'datatype' in schema:
        return _make_default_array(attr, schema, ctx)
    elif 'default' in schema:
        return schema['default']
    elif (schema.get('type') == 'object' or
          schema.get('allOf') or
          schema.get('anyOf')):
        return {}
    elif schema.get('type') == 'array':
        return []
    return copy.deepcopy(schema.get('default'))


def _make_node(instance, schema, ctx):
    if isinstance(instance, dict):
        return ObjectNode(instance, schema, ctx)
    elif isinstance(instance, list):
        return ListNode(instance, schema, ctx)
    else:
        return instance


def _unmake_node(obj):
    if isinstance(obj, Node):
        return obj.instance
    return obj


def _get_schema_for_property(schema, attr):
    subschema = schema.get('properties', {}).get(attr, None)
    if subschema is not None:
        return subschema
    for combiner in ['allOf', 'anyOf']:
        for subschema in schema.get(combiner, []):
            subsubschema = _get_schema_for_property(subschema, attr)
            if subsubschema != {}:
                return subsubschema
    return {}


def _get_schema_for_index(schema, i):
    items = schema.get('items', {})
    if isinstance(items, list):
        if i >= len(items):
            return {}
        else:
            return items[i]
    else:
        return items

class ValidationWarning(Warning):
    pass

class Node(object):
    def __init__(self, instance, schema, ctx):
        self._instance = instance
        self._schema = schema
        self._schema['$schema'] = 'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema'
        self._ctx = ctx

    def _report(self, value, attr=None):
        if value is None:
            if attr is None:
                msg = "is not valid list operation"
            else:
                msgfmt = "{0} is not valid to delete"
                msg = msgfmt.format(attr)
        else:
            if isinstance(value, six.string_types):
                value = "'{0}'".format(value)
            else:
                value = str(value)
                if len(value) > 55:
                    value = value[:56] + " ..."

            if attr is None:
                msgfmt = "{0} is not valid"
                msg = msgfmt.format(value)
            else:
                msgfmt = "{0} is not valid in {1}"
                msg = msgfmt.format(value, attr)
                
        try:
            filename = self._ctx.meta.filename
        except AttributeError:
            filename = None
        
        if filename is not None:
            msg = "In {0} {1}".format(filename, msg)

        if self._ctx._pass_invalid_values:
            warnings.warn(msg, ValidationWarning)
        else:
            raise jsonschema.ValidationError(msg)


    def _validate(self):
        instance = yamlutil.custom_tree_to_tagged_tree(
            self._instance, self._ctx._asdf)
        try:
            schema.validate(instance, schema=self._schema)
            valid = True
        except jsonschema.ValidationError:
            valid = False     
        return valid

    @property
    def instance(self):
        return self._instance

class ObjectNode(Node):
    @override__dir__
    def __dir__(self):
        return list(six.iterkeys(self._schema.get('properties', {})))

    def __eq__(self, other):
        if isinstance(other, ObjectNode):
            return self._instance == other._instance
        else:
            return self._instance == other

    def __getattr__(self, attr):
        if attr.startswith('_'):
            raise AttributeError('No attribute {0}'.format(attr))

        schema = _get_schema_for_property(self._schema, attr)

        try:
            val = self._instance[attr]
        except KeyError:
            if schema == {}:
                raise AttributeError("No attribute '{0}'".format(attr))
            val = _make_default(attr, schema, self._ctx)
            self._instance[attr] = val

        if isinstance(val, dict):
            # Meta is special cased to support NDData interface
            if attr == 'meta':
                node = MetaNode(val, schema, self._ctx)
            else:
                node = ObjectNode(val, schema, self._ctx)
        elif isinstance(val, list):
            node = ListNode(val, schema, self._ctx)
        else:
            node = val

        return node

    def __setattr__(self, attr, val):
        if attr.startswith('_'):
            self.__dict__[attr] = val
        else:
            schema = _get_schema_for_property(self._schema, attr)
            if val is None:
                val = _make_default(attr, schema, self._ctx)
            val = _cast(val, schema)
            old_val = self._instance.get(attr, None)

            self._instance[attr] = val
            if not self._validate():
                # Revert the change
                if old_val is None:
                    del self._instance[attr]
                else:
                    self._instance[attr] = old_val

                self._report(val, attr)
                
    def __delattr__(self, attr):
        if attr.startswith('_'):
            del self.__dict__[attr]
        else:
            old_val = self._instance.get(attr, None)

            try:
                del self._instance[attr]
            except KeyError:
                raise AttributeError(
                    "Attribute '{0}' missing".format(attr))
            if not self._validate():
                # Revert the change
                if old_val is not None:
                    self._instance[attr] = old_val
                        
                self._report(None, attr)

    def __hasattr__(self, attr):
        return (attr in self._instance or
                _find_property(self._schema, attr))

class ListNode(Node):
    def __cast(self, other):
        if isinstance(other, ListNode):
            return other._instance
        return other

    def __repr__(self):
        return repr(self._instance)

    def __eq__(self, other):
        return self._instance == self.__cast(other)

    def __ne__(self, other):
        return self._instance != self.__cast(other)

    def __contains__(self, item):
        return item in self._instance

    def __len__(self):
        return len(self._instance)

    def __getitem__(self, i):
        schema = _get_schema_for_index(self._schema, i)
        return _make_node(self._instance[i], schema, self._ctx)

    def __setitem__(self, i, val):
        schema = _get_schema_for_index(self._schema, i)
        self._instance[i] = _cast(val, schema)
        if not self._validate():
            self._report(val)

    def __delitem__(self, i):
        del self._instance[i]
        if not self._validate():
            self._report(None)

    def __getslice__(self, i, j):
        if isinstance(self._schema['items'], list):
            r = range(*(slice(i, j).indices(len(self._instance))))
            schema_parts = [
                _get_schema_for_index(self._schema, x) for x in r
            ]
        else:
            schema_parts = self._schema['items']
        schema = {'type': 'array', 'items': schema_parts}
        return _make_node(self._instance[i:j], schema, self._ctx)

    def __setslice__(self, i, j, other):
        parts = _unmake_node(other)
        parts = [_cast(x, _get_schema_for_index(self._schema, k))
                 for (k, x) in enumerate(parts)]
        self._instance[i:j] = _unmake_node(other)
        if not self._validate():
            self._report(None)

    def __delslice__(self, i, j):
        del self._instance[i:j]
        if not self._validate():
            self._report(None)
            
    def append(self, item):
        schema = _get_schema_for_index(self._schema, len(self._instance))
        self._instance.append(_cast(item, schema))
        if not self._validate():
            self._report(item)

    def insert(self, i, item):
        schema = _get_schema_for_index(self._schema, i)
        self._instance.insert(i, _cast(item, schema))
        if not self._validate():
            self._report(item)

    def pop(self, i=-1):
        schema = _get_schema_for_index(self._schema, 0)
        x = self._instance.pop(i)
        if not self._validate():
            self._report(None)
        return _make_node(x, schema, self._ctx)

    def remove(self, item):
        self._instance.remove(item)
        if not self._validate():
            self._report(None)

    def count(self, item):
        return self._instance.count(item)

    def index(self, item):
        return self._instance.index(item)

    def reverse(self):
        self._instance.reverse()
        if not self._validate():
            self._report(None)

    def sort(self, *args, **kwargs):
        self._instance.sort(*args, **kwargs)
        if not self._validate():
            self._report(None)
            
    def extend(self, other):
        for part in _unmake_node(other):
            self.append(part)

    def item(self, **kwargs):
        assert isinstance(self._schema['items'], dict)
        obj = ObjectNode(kwargs, self._schema['items'], self._ctx)
        if not obj._validate():
            self._report(None)
        return obj

def put_value(path, value, tree):
    """
    Put a value at the given path into tree, replacing it if it is
    already present.

    Parameters
    ----------
    path : list of str or int
        The path to the element.

    value : any
        The value to place

    tree : JSON object tree
    """
    cursor = tree
    for i in range(len(path) - 1):
        part = path[i]
        if isinstance(part, int):
            while len(cursor) <= part:
                cursor.append({})
            cursor = cursor[part]
        else:
            if isinstance(path[i + 1], int) or path[i + 1] == 'items':
                cursor = cursor.setdefault(part, [])
            else:
                cursor = cursor.setdefault(part, {})

    if isinstance(path[-1], int):
        while len(cursor) <= path[-1]:
            cursor.append({})
    cursor[path[-1]] = value


def merge_tree(a, b):
    """
    Merge elements from tree `b` into tree `a`.
    """
    def recurse(a, b):
        if isinstance(b, dict):
            if not isinstance(a, dict):
                return copy.deepcopy(b)
            for key, val in six.iteritems(b):
                a[key] = recurse(a.get(key), val)
            return a
        return copy.deepcopy(b)

    recurse(a, b)
    return a

#---------------------------------------------
# The following classes provide support
# for the NDData interface to Datamodels
#---------------------------------------------

class MetaNode(ObjectNode, collections.MutableMapping):
    """
    NDData compatibility class for meta node
    """
    def __init__(self, instance, schema, ctx):
        ObjectNode.__init__(self, instance, schema, ctx)

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
            schema = _get_schema_for_property(schema, attr)
            
        return _make_node(cursor, schema, self._ctx)

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
    
    def __iter__(self):
        return MetaNodeIterator(self)

    def __len__(self):
        def recurse(val):
            n = 0
            for subval in six.itervalues(val):
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

class MetaNodeIterator(six.Iterator):
    """
    An iterator for the meta node which flattens the hierachical structure
    """
    def __init__(self, node):
        self.key_stack = []
        self.iter_stack = [six.iteritems(node._instance)]
    
    def __iter__(self):
        return self
    
    def __next__(self):
        while self.iter_stack:
            try:
                key, val = six.next(self.iter_stack[-1])
            except StopIteration:
                self.iter_stack.pop()
                if self.iter_stack:
                    self.key_stack.pop()
                continue
                
            if isinstance(val, dict):
                self.key_stack.append(key)
                self.iter_stack.append(six.iteritems(val))
            else:
                return '.'.join(self.key_stack + [key])
                
        raise StopIteration
