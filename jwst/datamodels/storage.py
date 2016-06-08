"""
This model defines the concept of a pluggable storage backend for a data model.

While the same data model class offers the same front end to accessing
and changing its values (both data arrays and metadata values),
different storage backends can do the actual work of storing those
values.

There are currently 2 storage backends implemented:

  - TreeStorage: Stores values in something akin to a JSON or YAML
    tree.

  - FitsStorage: Stores values in a `astropy.io.fits.HDUList`
    hierarchy.  This includes all of the nice things about
    ``astropy.io.fits``, such as mmap'ing of array data.
"""
from __future__ import absolute_import, unicode_literals, division, print_function


class Storage(object):
    def close(self):
        raise NotImplementedError()

    def __get_array_section__(self, prop, obj, key):
        raise NotImplementedError()

    def __get_shape__(self):
        raise NotImplementedError()

    def __get__(self, prop, obj):
        raise NotImplementedError()

    def __set__(self, prop, obj, val):
        raise NotImplementedError()

    def __delete__(self, prop, obj):
        raise NotImplementedError()

    def validate(self, model):
        raise NotImplementedError()

    def extract_extra_elements(self, model):
        pass


class TreeStorage(Storage):
    def __init__(self, tree=None):
        if tree is None:
            tree = dict()
        self._tree = tree
        self._shape = None
        self._history = []

    def close(self):
        if hasattr(self, '_tree'):
            del self._tree

    def __get_array_section__(self, prop, obj, key):
        return self.__get__(prop, obj)[key]

    def __get_shape__(self):
        return self._shape

    def exists(self, prop, obj):
        cursor = self.tree
        if not prop.is_ad_hoc:
            for part in prop.path[:-1]:
                cursor = cursor.setdefault(part, {})
        else:
            for part in prop.path[:-1]:
                try:
                    cursor = cursor[part]
                except KeyError:
                    return False

        if len(prop.path):
            last_part = prop.path[-1]
            return last_part in cursor
        else:
            return True

    def __get__(self, prop, obj):
        val = self.__get_internal__(prop, obj)
        if prop.type == 'array':
            from . import schema
            items = prop.schema.get('items')
            return schema.ValidatingList(items, prop.name, val)
        return val

    def __get_internal__(self, prop, obj):
        cursor = self.tree
        if not prop.is_ad_hoc:
            for part in prop.path[:-1]:
                cursor = cursor.setdefault(part, {})
        else:
            for part in prop.path[:-1]:
                try:
                    cursor = cursor[part]
                except KeyError:
                    raise AttributeError('.'.join(prop.path))

        if len(prop.path):
            last_part = prop.path[-1]
            if last_part in cursor:
                return cursor[last_part]
            elif prop.is_data():
                val = prop._make_default(obj)
                cursor[last_part] = val
                return val
        else:
            return cursor
        raise AttributeError('.'.join(prop.path))

    def __set__(self, prop, obj, val):
        cursor = self.tree
        for part in prop.path[:-1]:
            cursor = cursor.setdefault(part, {})
        if len(prop.path):
            last_part = prop.path[-1]
            cursor[last_part] = val
        else:
            cursor = val

    def __delete__(self, prop, obj):
        cursor = self.tree
        for part in prop.path[:-1]:
            try:
                cursor = cursor[part]
            except KeyError:
                return
        if len(prop.path):
            last_part = prop.path[-1]
            if last_part in cursor:
                del cursor[last_part]
        else:
            self._tree = {}

    def validate(self, model):
        from . import schema
        schema.validate(self.tree, model.schema)

    def get_fits_header(self, model, hdu_name='PRIMARY'):
        """
        Generates a FITS header for a given FITS hdu.

        Parameters
        ----------
        model : ModelBase object

        hdu_name : str, optional
            The name of the HDU to get the WCS from.  This must use
            named HDU's, not numerical order HDUs.  To get the primary HDU,
            pass ``'PRIMARY'`` (default).

        Returns
        -------
        header : `~astropy.io.fits.Header` object
        """
        from astropy.io import fits
        from . import schema

        elements = schema.get_elements_for_fits_hdu(
            model.schema, hdu_name=hdu_name)

        header = fits.Header()

        for keyword, path in elements.items():
            val = model.get_item_as_json_value(path)
            if isinstance(val, list):
                for subval in val:
                    header[keyword] = subval
            else:
                header[keyword] = val

        return header

    @property
    def history(self):
        return self._history

    @history.setter
    def history(self, v):
        assert isinstance(v, list)
        self._history = v

    @property
    def tree(self):
        if not hasattr(self, '_tree'):
            raise IOError("Storage has been closed.")
        return self._tree

    def to_tree(self):
        return self.tree


class HasStorage(object):
    def __init__(self, storage=None):
        if storage is None:
            storage = TreeStorage()
        assert isinstance(storage, Storage)
        self._storage = storage
        self._owns_storage = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if hasattr(self, '_parent'):
            del self._parent
        if hasattr(self, '_storage') and self._owns_storage:
            self._storage.close()

    @property
    def storage(self):
        return self._storage
