# -*- coding: utf-8 -*-
"""
Data model class heirarchy
"""
from __future__ import absolute_import, unicode_literals, division, print_function

import copy
import datetime
import inspect
import os
import sys
import warnings

import numpy as np

from astropy.extern import six
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

from asdf import AsdfFile
from asdf import yamlutil
from asdf import schema as asdf_schema
from asdf import extension as asdf_extension

from . import fits_support
from . import properties
from . import schema as mschema

from .extension import BaseExtension
from jwst.transforms.jwextension import JWSTExtension
from gwcs.extension import GWCSExtension


jwst_extensions = [GWCSExtension(), JWSTExtension(), BaseExtension()]

class DataModel(properties.ObjectNode):
    """
    Base class of all of the data models.
    """
    schema_url = "core.schema.yaml"

    def __init__(self, init=None, schema=None, extensions=None,
                 pass_invalid_values=False):
        """
        Parameters
        ----------
        init : shape tuple, file path, file object, astropy.io.fits.HDUList, numpy array, None

            - None: A default data model with no shape

            - shape tuple: Initialize with empty data of the given
              shape

            - file path: Initialize from the given file (FITS or ASDF)

            - readable file object: Initialize from the given file
              object

            - ``astropy.io.fits.HDUList``: Initialize from the given
              `~astropy.io.fits.HDUList`.

            - A numpy array: Used to initialize the data array

            - dict: The object model tree for the data model

        schema : tree of objects representing a JSON schema, or string naming a schema, optional
            The schema to use to understand the elements on the model.
            If not provided, the schema associated with this class
            will be used.

        extensions: classes extending the standard set of extensions, optional.
            If an extension is defined, the prefix used should be 'url'.

        pass_invalid_values: If True, values that do not validate the schema can
            be read and written, but with a warning message
        """
        filename = os.path.abspath(inspect.getfile(self.__class__))
        base_url = os.path.join(
            os.path.dirname(filename), 'schemas', '')

        if extensions is None:
            extensions = jwst_extensions[:]
        else:
            extensions.extend(jwst_extensions)
        self._extensions = extensions

        if schema is None:
            schema_path = os.path.join(base_url, self.schema_url)
            extension_list = asdf_extension.AsdfExtensionList(self._extensions)
            schema = asdf_schema.load_schema(schema_path, 
                resolver=extension_list.url_mapping, resolve_references=True)

        self._schema = mschema.flatten_combiners(schema)

        if "PASS_INVALID_VALUES" in os.environ:
            pass_invalid_values = os.environ["PASS_INVALID_VALUES"]
            try:
                self._pass_invalid_values = bool(int(pass_invalid_values))
            except ValueError:
                self._pass_invalid_values = False
        else:
            self._pass_invalid_values = pass_invalid_values

        self._files_to_close = []
        is_array = False
        is_shape = False
        shape = None
        if init is None:
            asdf = AsdfFile(extensions=extensions)
        elif isinstance(init, dict):
            asdf = AsdfFile(init, extensions=extensions)
        elif isinstance(init, np.ndarray):
            asdf = AsdfFile(extensions=extensions)
            shape = init.shape
            is_array = True
        elif isinstance(init, self.__class__):
            instance = copy.deepcopy(init._instance)
            self._schema = init._schema
            self._shape = init._shape
            self._asdf = AsdfFile(instance, extensions=self._extensions)
            self._instance = instance
            self._ctx = self
            self.__class__ = init.__class__
            return
        elif isinstance(init, DataModel):
            raise TypeError(
                "Passed in {0!r} is not of the expected subclass {1!r}".format(
                    init.__class__.__name__, self.__class__.__name__))
        elif isinstance(init, AsdfFile):
            asdf = init
        elif isinstance(init, tuple):
            for item in init:
                if not isinstance(item, int):
                    raise ValueError("shape must be a tuple of ints")
            shape = init
            asdf = AsdfFile()
            is_shape = True
        elif isinstance(init, fits.HDUList):
            asdf = fits_support.from_fits(init, self._schema,
                                          extensions=self._extensions,
                                          validate=False,
                                          pass_invalid_values=self._pass_invalid_values)
        elif isinstance(init, six.string_types):
            if isinstance(init, bytes):
                init = init.decode(sys.getfilesystemencoding())
            try:
                hdulist = fits.open(init)
            except IOError:
                try:
                    asdf = AsdfFile.open(init, extensions=self._extensions)
                    # TODO: Add json support
                except ValueError:
                    raise IOError(
                        "File does not appear to be a FITS or ASDF file.")
            else:
                asdf = fits_support.from_fits(hdulist, self._schema,
                                              extensions=self._extensions,
                                              validate=False,
                                              pass_invalid_values=self._pass_invalid_values)
                self._files_to_close.append(hdulist)
        else:
            raise ValueError(
                "Can't initialize datamodel using {0}".format(str(type(init))))

        self._shape = shape
        self._instance = asdf.tree
        self._asdf = asdf
        self._ctx = self

        # if the input model doesn't have a date set, use the current date/time
        if self.meta.date is None:
            self.meta.date = Time(datetime.datetime.now())
            if hasattr(self.meta.date, 'value'):
                self.meta.date.format = 'isot'
                self.meta.date = str(self.meta.date.value)

        # if the input is from a file, set the filename attribute
        if isinstance(init, six.string_types):
            self.meta.filename = os.path.basename(init)

        if is_array:
            primary_array_name = self.get_primary_array_name()
            if primary_array_name is None:
                raise TypeError(
                    "Array passed to DataModel.__init__, but model has "
                    "no primary array in its schema")
            setattr(self, primary_array_name, init)

        if is_shape:
            getattr(self, self.get_primary_array_name())

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        for fd in self._files_to_close:
            if fd is not None:
                fd.close()

    def copy(self, memo=None):
        """
        Returns a deep copy of this model.
        """
        result = self.__class__(
            init=copy.deepcopy(self._instance, memo=memo),
            schema=self._schema,
            extensions=self._extensions)
        result._shape = self._shape
        return result

    __copy__ = __deepcopy__ = copy

    def get_primary_array_name(self):
        """
        Returns the name "primary" array for this model, which
        controls the size of other arrays that are implicitly created.
        This is intended to be overridden in the subclasses if the
        primary array's name is not "data".
        """
        return 'data'

    def on_save(self, path=None):
        """
        This is a hook that is called just before saving the file.
        It can be used, for example, to update values in the metadata
        that are based on the content of the data.

        Override it in the subclass to make it do something, but don't
        forget to "chain up" to the base class, since it does things
        there, too.

        Parameters
        ----------
        path : str
            The path to the file that we're about to save to.
        """
        if isinstance(path, six.string_types):
            self.meta.filename = os.path.basename(path)

        self.meta.date = Time(datetime.datetime.now())
        self.meta.date.format = 'isot'
        self.meta.date = self.meta.date.value

    def save(self, path, *args, **kwargs):
        """
        Save to either a FITS or ASDF file, depending on the path.

        Parameters
        ----------
        path : string
        """
        base, ext = os.path.splitext(path)
        if isinstance(ext, bytes):
            ext = ext.decode(sys.getfilesystemencoding())

        # TODO: Support gzip-compressed fits
        if ext == '.fits':
            kwargs.setdefault('clobber', True)
            self.to_fits(path, *args, **kwargs)
        elif ext == '.asdf':
            self.to_asdf(path, *args, **kwargs)
        else:
            raise ValueError("unknown filetype {0}".format(ext))

    @classmethod
    def from_asdf(cls, init, schema=None):
        """
        Load a data model from a ASDF file.

        Parameters
        ----------
        init : file path, file object, asdf.AsdfFile object
            - file path: Initialize from the given file
            - readable file object: Initialize from the given file object
            - asdf.AsdfFile: Initialize from the given
              `~asdf.AsdfFile`.

        schema :
            Same as for `__init__`

        Returns
        -------
        model : DataModel instance
        """
        return cls(init, schema=schema, extensions=self._extensions)

    def to_asdf(self, init, *args, **kwargs):
        """
        Write a DataModel to a ASDF file.

        Parameters
        ----------
        init : file path or file object

        *args, **kwargs
            Any additional arguments are passed along to
            `asdf.AsdfFile.write_to`.
        """
        self.on_save(init)

        AsdfFile(self._instance, extensions=self._extensions).write_to(init, *args, **kwargs)

    @classmethod
    def from_fits(cls, init, schema=None):
        """
        Load a model from a FITS file.

        Parameters
        ----------
        init : file path, file object, astropy.io.fits.HDUList
            - file path: Initialize from the given file
            - readable file object: Initialize from the given file object
            - astropy.io.fits.HDUList: Initialize from the given
              `~astropy.io.fits.HDUList`.

        schema :
            Same as for `__init__`

        Returns
        -------
        model : DataModel instance
        """
        return cls(init, schema=schema)

    def to_fits(self, init, *args, **kwargs):
        """
        Write a DataModel to a FITS file.

        Parameters
        ----------
        init : file path or file object

        *args, **kwargs
            Any additional arguments are passed along to
            `astropy.io.fits.writeto`.
        """
        self.on_save(init)

        with fits_support.to_fits(self._instance, self._schema,
                                  extensions=self._extensions) as ff:
            ff.write_to(init, *args, **kwargs)

    @property
    def shape(self):
        if self._shape is None:
            if self.get_primary_array_name() in self._instance:
                return getattr(self, self.get_primary_array_name()).shape
            else:
                return None
        return self._shape

    def __setattr__(self, attr, value):
        if attr == 'shape':
            object.__setattr__(self, attr, value)
        else:
            super(DataModel, self).__setattr__(attr, value)

    def extend_schema(self, new_schema):
        """
        Extend the model's schema using the given schema, by combining
        it in an "allOf" array.

        Parameters
        ----------
        new_schema : schema tree
        """
        schema = {'allOf': [self._schema, new_schema]}
        self._schema = mschema.flatten_combiners(schema)
        self._validate()
        return self

    def add_schema_entry(self, position, new_schema):
        """
        Extend the model's schema by placing the given new_schema at
        the given dot-separated position in the tree.

        Parameters
        ----------
        position : str

        new_schema : schema tree
        """
        parts = position.split('.')
        schema = new_schema
        for part in parts[::-1]:
            schema = {'type': 'object', 'properties': {part: schema}}
        return self.extend_schema(schema)

    # return_result retained for backward compatibility
    def find_fits_keyword(self, keyword, return_result=True):
        """
        Utility function to find a reference to a FITS keyword in this
        model's schema.  This is intended for interactive use, and not
        for use within library code.

        Parameters
        ----------
        keyword : str
            A FITS keyword name

        Returns
        -------
        locations : list of str

            If `return_result` is `True, a list of the locations in
            the schema where this FITS keyword is used.  Each element
            is a dot-separated path.

        Example
        -------
        >>> model.find_fits_keyword('DATE-OBS')
        ['observation.date']
        """
        from . import schema
        return schema.find_fits_keyword(self.schema, keyword)

    def search_schema(self, substring):
        """
        Utility function to search the metadata schema for a
        particular phrase.

        This is intended for interactive use, and not for use within
        library code.

        The searching is case insensitive.

        Parameters
        ----------
        substring : str
            The substring to search for.

        Returns
        -------
        locations : list of tuples
        """
        from . import schema
        return schema.search_schema(self.schema, substring)

    def __getitem__(self, key):
        """
        Get a metadata value using a dotted name.
        """
        assert isinstance(key, six.string_types)
        meta = self
        for part in key.split('.'):
            try:
                meta = getattr(meta, part)
            except AttributeError:
                raise KeyError(repr(key))
        return meta

    def get_item_as_json_value(self, key):
        """
        Equivalent to __getitem__, except returns the value as a JSON
        basic type, rather than an arbitrary Python type.
        """
        assert isinstance(key, six.string_types)
        meta = self
        parts = key.split('.')
        for part in parts:
            try:
                meta = getattr(meta, part)
            except AttributeError:
                raise KeyError(repr(key))
        return yamlutil.custom_tree_to_tagged_tree(meta, self._instance)

    def __setitem__(self, key, value):
        """
        Set a metadata value using a dotted name.
        """
        assert isinstance(key, six.string_types)
        meta = self
        parts = key.split('.')
        for part in parts[:-1]:
            try:
                part = int(part)
            except ValueError:
                try:
                    meta = getattr(meta, part)
                except AttributeError:
                    raise KeyError(repr(key))
            else:
                meta = meta[part]

        part = parts[-1]
        try:
            part = int(part)
        except ValueError:
            setattr(meta, part, value)
        else:
            meta[part] = value

    def iteritems(self):
        """
        Iterates over all of the schema items in a flat way.

        Each element is a pair (`key`, `value`).  Each `key` is a
        dot-separated name.  For example, the schema element
        `meta.observation.date` will end up in the result as::

            ("meta.observation.date": "2012-04-22T03:22:05.432")
        """
        def recurse(tree, path=[]):
            if isinstance(tree, dict):
                for key, val in six.iteritems(tree):
                    for x in recurse(val, path + [key]):
                        yield x
            elif isinstance(tree, (list, tuple)):
                for i, val in enumerate(tree):
                    for x in recurse(val, path + [i]):
                        yield x
            elif tree is not None:
                yield (str('.'.join(six.text_type(x) for x in path)), tree)

        for x in recurse(self._instance):
            yield x

    # We are just going to define the items to return the iteritems
    items = iteritems


    def iterkeys(self):
        """
        Iterates over all of the schema keys in a flat way.

        Each result of the iterator is a `key`.  Each `key` is a
        dot-separated name.  For example, the schema element
        `meta.observation.date` will end up in the result as the
        string `"meta.observation.date"`.
        """
        for key, val in self.iteritems():
            yield key

    if six.PY3:
        keys = iterkeys
    else:
        def keys(self):
            """
            Gets all of the schema keys in a flat way.

            Each result of the iterator is a `key`.  Each `key` is a
            dot-separated name.  For example, the schema element
            `meta.observation.date` will end up in the result as the
            string `"meta.observation.date"`.
            """
            return list(self.iterkeys())

    def itervalues(self):
        """
        Iterates over all of the schema values in a flat way.
        """
        for key, val in self.iteritems():
            yield val

    if six.PY3:
        values = itervalues
    else:
        def values(self):
            """
            Gets all of the schema values in a flat way.
            """
            return list(self.itervalues())

    def update(self, d, only=''):
        """
        Updates this model with the metadata elements from another model.

        Parameters
        ----------
        d : model or dictionary-like object
            The model to copy the metadata elements from. Can also be a
            dictionary or dictionary of dictionaries or lists.
        only: only update the named hdu from extra_fits, e.g.
            only='PRIMARY'. Can either be a list of hdu names
            or a single string. If left blank, update all the hdus.
        """
        def hdu_keywords_from_data(d, path, hdu_keywords):
            # Walk tree and add paths to keywords to hdu keywords
            if isinstance(d, dict):
                for key, val in d.items():
                    if len(path) > 0 or key != 'extra_fits':
                        hdu_keywords_from_data(val, path + [key], hdu_keywords)
            elif isinstance(d, list):
                for key, val in enumerate(d):
                    hdu_keywords_from_data(val, path + [key], hdu_keywords)
            elif isinstance(d, np.ndarray):
                # skip data arrays
                pass
            else:
                hdu_keywords.append(path)

        def hdu_keywords_from_schema(subschema, path, combiner, ctx, recurse):
            # Add path to keyword to hdu_keywords if in list of hdu names
            if 'fits_keyword' in subschema:
                fits_hdu = subschema.get('fits_hdu', 'PRIMARY')
                if fits_hdu in hdu_names:
                    ctx.append(path)

        def hdu_names_from_schema(subschema, path, combiner, ctx, recurse):
            # Build a set of hdu names from the schema
            hdu_name = subschema.get('fits_hdu')
            if hdu_name:
                hdu_names.add(hdu_name)

        def included(cursor, part):
            # Test if part is in the cursor
            if isinstance(part, int):
                return part >= 0 and part < len(cursor)
            else:
                return part in cursor

        def set_hdu_keyword(this_cursor, that_cursor, path):
            # Copy an element pointed to by path from that to this
            part = path.pop(0)
            if not included(that_cursor, part):
                return
            if len(path) == 0:
                this_cursor[part] = copy.deepcopy(that_cursor[part])
            else:
                that_cursor = that_cursor[part]
                if not included(this_cursor, part):
                    if isinstance(path[0], int):
                        if isinstance(part, int):
                            this_cursor.append([])
                        else:
                            this_cursor[part] = []
                    else:
                        if isinstance(part, int):
                            this_cursor.append({})
                        else:
                            this_cursor[part] = {}
                this_cursor = this_cursor[part]
                set_hdu_keyword(this_cursor, that_cursor, path)

        # Get the list of hdu names from the model so that updates
        # are limited to those hdus

        if only:
            if isinstance(only, six.string_types):
                hdu_names = set([only])
            else:
                hdu_names = set(list(only))
        else:
            hdu_names = set(['PRIMARY'])
            mschema.walk_schema(self._schema, hdu_names_from_schema, hdu_names)

        # Get the paths to all the keywords that will be updated from

        hdu_keywords = []
        if isinstance(d, DataModel):
            schema = d._schema
            d = d._instance
            mschema.walk_schema(schema, hdu_keywords_from_schema, hdu_keywords)
        else:
            path = []
            hdu_keywords_from_data(d, path, hdu_keywords)

        # Perform the updates to the keywords mentioned in the schema

        for path in hdu_keywords:
            set_hdu_keyword(self._instance, d, path)

        # Perform updates to extra_fits area of a model

        for hdu_name in hdu_names:
            path = ['extra_fits', hdu_name, 'header']
            set_hdu_keyword(self._instance, d, path)

    def to_flat_dict(self, include_arrays=True):
        """
        Returns a dictionary of all of the schema items as a flat dictionary.

        Each dictionary key is a dot-separated name.  For example, the
        schema element `meta.observation.date` will end up in the
        dictionary as::

            { "meta.observation.date": "2012-04-22T03:22:05.432" }

        """
        def convert_val(val):
            if isinstance(val, datetime.datetime):
                return val.isoformat()
            elif isinstance(val, Time):
                return str(val)
            return val

        if include_arrays:
            return dict((key, convert_val(val)) for (key, val) in self.iteritems())
        else:
            return dict((key, convert_val(val)) for (key, val) in self.iteritems()
                        if not isinstance(val, np.ndarray))

    @property
    def schema(self):
        return self._schema

    def get_fileext(self):
        return 'fits'

    # TODO: This is just here for backward compatibility
    @property
    def _extra_fits(self):
        return self.extra_fits

    # TODO: For backward compatibility
    def get_section(self, name):
        return getattr(self, name)

    @property
    def history(self):
        return self._instance.setdefault('history', [])

    @history.setter
    def history(self, v):
        self._instance['history'] = v

    def get_fits_wcs(self, hdu_name='PRIMARY', key=' '):
        """
        Get a `astropy.wcs.WCS` object created from the FITS WCS
        information in the model.

        Note that modifying the returned WCS object will not modify
        the data in this model.  To update the model, use `set_fits_wcs`.

        Parameters
        ----------
        hdu_name : str, optional
            The name of the HDU to get the WCS from.  This must use
            named HDU's, not numerical order HDUs.  To get the primary HDU,
            pass ``'PRIMARY'`` (default).

        key : str, optional
            The name of a particular WCS transform to use.  This may
            be either ``' '`` or ``'A'``-``'Z'`` and corresponds to
            the ``"a"`` part of the ``CTYPEia`` cards.  *key* may only
            be provided if *header* is also provided.

        Returns
        -------
        wcs : `astropy.wcs.WCS` or `pywcs.WCS` object
            The type will depend on what libraries are installed on
            this system.
        """
        extensions = self._asdf._extensions
        ff = fits_support.to_fits(self._instance, self._schema,
                                  extensions=extensions)
        hdu = fits_support.get_hdu(ff._hdulist, hdu_name)
        header = hdu.header

        return WCS(header, key=key, relax=True, fix=True)

    def set_fits_wcs(self, wcs, hdu_name='PRIMARY'):
        """
        Sets the FITS WCS information on the model using the given
        `astropy.wcs.WCS` object.

        Note that the "key" of the WCS is stored in the WCS object
        itself, so it can not be set as a parameter to this method.

        Parameters
        ----------
        wcs : `astropy.wcs.WCS` or `pywcs.WCS` object
            The object containing FITS WCS information

        hdu_name : str, optional
            The name of the HDU to set the WCS from.  This must use
            named HDU's, not numerical order HDUs.  To set the primary
            HDU, pass ``'PRIMARY'`` (default).
        """
        header = wcs.to_header()
        if hdu_name == 'PRIMARY':
            hdu = fits.PrimaryHDU(header=header)
        else:
            hdu = fits.ImageHDU(name=hdu_name, header=header)
        hdulist = fits.HDUList([hdu])

        ff = fits_support.from_fits(hdulist, self._schema, validate=False)

        self._instance = properties.merge_tree(self._instance, ff.tree)
