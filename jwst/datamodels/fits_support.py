# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

import datetime
import os
import re
import warnings

import numpy as np

import jsonschema

from astropy.extern import six
from astropy.io import fits
from astropy import time

from asdf import fits_embed
from asdf import resolver
from asdf import schema as asdf_schema
from asdf.tags.core import ndarray, HistoryEntry
from asdf import treeutil
from asdf.util import HashableDict
from asdf import tagged

from jsonschema import validators

from . import properties
from . import schema as mschema
from . import util


__all__ = ['to_fits', 'from_fits', 'fits_hdu_name', 'get_hdu']


_builtin_regexes = [
    '', 'NAXIS[0-9]{0,3}', 'BITPIX', 'XTENSION', 'PCOUNT', 'GCOUNT',
    'EXTEND', 'BSCALE', 'BZERO', 'BLANK', 'DATAMAX', 'DATAMIN',
    'EXTNAME', 'EXTVER', 'EXTLEVEL', 'GROUPS', 'PYTPE[0-9]',
    'PSCAL[0-9]', 'PZERO[0-9]', 'SIMPLE', 'TFIELDS',
    'TBCOL[0-9]{1,3}', 'TFORM[0-9]{1,3}', 'TTYPE[0-9]{1,3}',
    'TUNIT[0-9]{1,3}', 'TSCAL[0-9]{1,3}', 'TZERO[0-9]{1,3}',
    'TNULL[0-9]{1,3}', 'TDISP[0-9]{1,3}', 'HISTORY'
    ]


_builtin_regex = re.compile(
    '|'.join('(^{0}$)'.format(x) for x in _builtin_regexes))


def _is_builtin_fits_keyword(key):
    """
    Returns `True` if the given `key` is a built-in FITS keyword, i.e.
    a keyword that is managed by ``astropy.io.fits`` and we wouldn't
    want to propagate through the `_extra_fits` mechanism.
    """
    return _builtin_regex.match(key) is not None


_keyword_indices = [
    ('nnn', 1000, None),
    ('nn', 100, None),
    ('n', 10, None),
    ('s', 27, ' ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    ]


def _get_indexed_keyword(keyword, i):
    for (sub, max, r) in _keyword_indices:
        if sub in keyword:
            if i >= max:
                raise ValueError(
                    "Too many entries for given keyword '{0}'".format(keyword))
            if r is None:
                val = six.text_type(i)
            else:
                val = r[i]
            keyword = keyword.replace(sub, val)

    return keyword


if six.PY3:
    def fits_hdu_name(name):
        """
        Returns a FITS hdu name in the correct form for the current
        version of Python.
        """
        if isinstance(name, bytes):
            return name.decode('ascii')
        return name
else:
    def fits_hdu_name(name):
        """
        Returns a FITS hdu name in the correct form for the current
        version of Python.
        """
        if isinstance(name, six.text_type):
            return name.encode('ascii')
        return name


def _get_hdu_name(schema):
    hdu_name = schema.get('fits_hdu')
    if hdu_name in (None, 'PRIMARY'):
        hdu_name = 0
    else:
        hdu_name = fits_hdu_name(hdu_name)
    return hdu_name


def _make_new_hdu(hdulist, value, hdu_name, index=None):
    try:
        defs = fits.ColDefs(value) # just to test if value is a table
        hdu = fits.BinTableHDU(value, name=hdu_name)
    except TypeError:
        hdu = fits.ImageHDU(value, name=hdu_name)
    if index is not None:
        hdu.ver = index + 1
    hdulist.append(hdu)
    return hdu


def _get_hdu_pair(hdu_name, index=None):
    if index is None:
        pair = hdu_name
    else:
        pair = (hdu_name, index + 1)
    return pair


def get_hdu(hdulist, hdu_name, index=None):
    pair = _get_hdu_pair(hdu_name, index=index)
    try:
        hdu = hdulist[pair]
    except (KeyError, IndexError, AttributeError):
        try:
            if index is None:
                hdu = hdulist[(pair, 1)]
            elif index == 0:
                hdu = hdulist[pair[0]]
            else:
                raise
        except (KeyError, IndexError, AttributeError):
            raise AttributeError(
                "Property missing because FITS file has no "
                "'{0!r}' HDU".format(
                    pair))

    if index is not None:
        if hdu.header.get('EXTVER', 1) != index + 1:
            raise AttributeError(
                "Property missing because FITS file has no "
                "{0!r} HDU".format(
                    pair))

    return hdu


def _make_hdu(hdulist, hdu_name, index=None, hdu_type=None, value=None):
    if hdu_type is None:
        if hdu_name in (0, 'PRIMARY'):
            hdu_type = fits.PrimaryHDU
        else:
            hdu_type = fits.ImageHDU
    if hdu_type == fits.PrimaryHDU:
        hdu = hdu_type(value)
    else:
        hdu = hdu_type(value, name=hdu_name)
    if index is not None:
        hdu.ver = index + 1
    hdulist.append(hdu)
    return hdu


def _get_or_make_hdu(hdulist, hdu_name, index=None, hdu_type=None, value=None):
    try:
        hdu = get_hdu(hdulist, hdu_name, index=index)
    except AttributeError:
        hdu = _make_hdu(hdulist, hdu_name, index=index, hdu_type=hdu_type,
                        value=value)
    else:
        if hdu_type is not None and not isinstance(hdu, hdu_type):
            new_hdu = _make_hdu(hdulist, hdu_name, index=index,
                                hdu_type=hdu_type, value=value)
            for key, val in six.iteritems(hdu.header):
                if not _is_builtin_fits_keyword(key):
                    new_hdu.header[key] = val
            hdulist.remove(hdu)
            hdu = new_hdu
        elif value is not None:
            hdu.data = value
    return hdu


def _assert_non_primary_hdu(hdu_name):
    if hdu_name in (None, 0, 'PRIMARY'):
        raise ValueError(
            "Schema for data property does not specify a non-primary hdu name")


##############################################################################
# WRITER


def _fits_comment_section_handler(validator, properties, instance, schema):
    if not validator.is_type(instance, "object"):
        return

    title = schema.get('title')
    if title is not None:
        current_comment_stack = validator.comment_stack
        current_comment_stack.append(util.ensure_ascii(title))

    for property, subschema in six.iteritems(properties):
        if property in instance:
            for error in validator.descend(
                instance[property],
                subschema,
                path=property,
                schema_path=property,
            ):
                yield error

    if title is not None:
        current_comment_stack.pop(-1)


def _fits_element_writer(validator, fits_keyword, instance, schema):
    if schema.get('type', 'object') == 'array':
        raise ValueError("'fits_keyword' not valid with type of 'array'")

    hdu_name = _get_hdu_name(schema)
    index = getattr(validator, 'sequence_index', None)
    hdu = _get_or_make_hdu(validator.hdulist, hdu_name, index=index)

    for comment in validator.comment_stack:
        hdu.header.append((' ', ''), end=True)
        hdu.header.append((' ', comment), end=True)
        hdu.header.append((' ', ''), end=True)
    validator.comment_stack = []

    comment = util.ensure_ascii(util.get_short_doc(schema))
    instance = util.ensure_ascii(instance)

    if fits_keyword in ('COMMENT', 'HISTORY'):
        for item in instance:
            hdu.header[fits_keyword] = util.ensure_ascii(item)
    elif fits_keyword in hdu.header:
        hdu.header[fits_keyword] = (instance, comment)
    else:
        hdu.header.append((fits_keyword, instance, comment), end=True)


def _fits_array_writer(validator, _, instance, schema):
    if instance is None:
        return

    instance = np.asanyarray(instance)

    if not len(instance.shape):
        return

    if 'ndim' in schema:
        ndarray.validate_ndim(validator, schema['ndim'], instance, schema)
    if 'max_ndim' in schema:
        ndarray.validate_max_ndim(validator, schema['max_ndim'], instance, schema)
    if 'dtype' in schema:
        ndarray.validate_dtype(validator, schema['dtype'], instance, schema)

    hdu_name = _get_hdu_name(schema)
    _assert_non_primary_hdu(hdu_name)
    if instance.dtype.names is not None:
        hdu_type = fits.BinTableHDU
    else:
        hdu_type = fits.ImageHDU
    index = getattr(validator, 'sequence_index', None)

    hdu = _get_or_make_hdu(validator.hdulist, hdu_name, index=index, hdu_type=hdu_type)

    hdu.data = instance


# This is copied from jsonschema._validators and modified to keep track
# of the index of the item we've recursed into.
def _fits_item_recurse(validator, items, instance, schema):
    if not validator.is_type(instance, "array"):
        return

    if validator.is_type(items, "object"):
        for index, item in enumerate(instance):
            validator.sequence_index = index
            for error in validator.descend(item, items, path=index):
                yield error
    else:
        # We don't do the index trick on "tuple validated" sequences
        for (index, item), subschema in zip(enumerate(instance), items):
            for error in validator.descend(
                item, subschema, path=index, schema_path=index,
            ):
                yield error


def _fits_type(validator, items, instance, schema):
    if instance in ('N/A', '#TODO', '', None):
        return
    return validators._validators.type_draft4(validator, items, instance, schema)


FITS_VALIDATORS = HashableDict(asdf_schema.YAML_VALIDATORS)


FITS_VALIDATORS.update({
    'fits_keyword': _fits_element_writer,
    'ndim': _fits_array_writer,
    'max_ndim': _fits_array_writer,
    'datatype': _fits_array_writer,
    'items': _fits_item_recurse,
    'properties': _fits_comment_section_handler,
    'type': _fits_type
})


META_SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'metaschema'))


FITS_SCHEMA_URL_MAPPING = resolver.Resolver(
    [
        ('http://stsci.edu/schemas/fits-schema/',
         'file://' + META_SCHEMA_PATH + '/{url_suffix}.yaml')
    ] + resolver.DEFAULT_URL_MAPPING, 'url')


def _save_from_schema(hdulist, tree, schema):
    def convert_datetimes(node, json_id):
        if isinstance(node, datetime.datetime):
            node = time.Time(node)
        if isinstance(node, time.Time):
            node = six.text_type(time.Time(node, format='iso'))
        return node
    tree = treeutil.walk_and_modify(tree, convert_datetimes)

    validator = asdf_schema.get_validator(
        schema, None, FITS_VALIDATORS, FITS_SCHEMA_URL_MAPPING)

    validator.hdulist = hdulist
    # TODO: Handle comment stack on per-hdu-basis
    validator.comment_stack = []
    # Tag the tree values first, the validator requires it
    _tag_values(tree, schema)
    # This actually kicks off the saving
    validator.validate(tree, _schema=schema)


def _save_extra_fits(hdulist, tree):
    # Handle _extra_fits
    for hdu_name, parts in six.iteritems(tree.get('extra_fits', {})):
        hdu_name = fits_hdu_name(hdu_name)
        if 'data' in parts:
            hdu = _make_new_hdu(hdulist, parts['data'], hdu_name)
        if 'header' in parts:
            hdu = _get_or_make_hdu(hdulist, hdu_name)
            for key, val, comment in parts['header']:
                if _is_builtin_fits_keyword(key):
                    continue
                hdu.header.append((key, val, comment), end=True)


def _save_history(hdulist, tree):
    history = tree.get('history', [])
    for i in range(len(history)):
        # There is no guarantee the user has added proper HistoryEntry records
        if not isinstance(history[i], HistoryEntry):
            if isinstance(history[i], dict):
                history[i] = HistoryEntry(history[i])
            else:
                history[i] = HistoryEntry({'description': str(history[i])})
        hdulist[0].header['HISTORY'] = history[i]['description']


def _tag_values(tree, schema):
    # Replace tag value in tree with tagged versions

    def included(cursor, part):
        if isinstance(part, int):
            return part > 0 and part < len(cursor)
        else:
            return part in cursor

    def callback(subschema, path, combiner, ctx, recurse):
        tag = subschema.get('tag')
        if tag is not None:
            cursor = tree
            for part in path[:-1]:
                if included(cursor, part):
                    cursor = cursor[part]
                else:
                    return
            part = path[-1]
            if included(cursor, part):
                cursor[part] = tagged.tag_object(tag, cursor[part])

    mschema.walk_schema(schema, callback)


def to_fits(tree, schema, extensions=None):
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())

    _save_from_schema(hdulist, tree, schema)
    _save_extra_fits(hdulist, tree)
    _save_history(hdulist, tree)

    asdf = fits_embed.AsdfInFits(hdulist, tree, extensions=extensions)
    return asdf


##############################################################################
# READER


def _fits_keyword_loader(hdulist, fits_keyword, schema, hdu_index, known_keywords):
    hdu_name = _get_hdu_name(schema)
    try:
        hdu = get_hdu(hdulist, hdu_name, hdu_index)
    except AttributeError:
        return None

    try:
        val = hdu.header[fits_keyword]
    except KeyError:
        return None

    tag = schema.get('tag')
    if tag is not None:
        val = tagged.tag_object(tag, val)

    known_keywords.setdefault(hdu, set()).add(fits_keyword)

    return val


def _fits_array_loader(hdulist, schema, hdu_index, known_datas):
    hdu_name = _get_hdu_name(schema)
    _assert_non_primary_hdu(hdu_name)
    try:
        hdu = get_hdu(hdulist, hdu_name, hdu_index)
    except AttributeError:
        return None

    known_datas.add(hdu)

    data = hdu.data
    data = properties._cast(data, schema)
    return data


def _schema_has_fits_hdu(schema):
    has_fits_hdu = [False]

    for node in treeutil.iter_tree(schema):
        if isinstance(node, dict) and 'fits_hdu' in node:
            has_fits_hdu[0] = True

    return has_fits_hdu[0]


def _load_from_schema(hdulist, schema, tree, validate=True,
                      pass_invalid_values=False):
    known_keywords = {}
    known_datas = set()

    def callback(schema, path, combiner, ctx, recurse):
        result = None

        if 'fits_keyword' in schema:
            fits_keyword = schema['fits_keyword']
            result = _fits_keyword_loader(
                hdulist, fits_keyword, schema,
                ctx.get('hdu_index'), known_keywords)
            if result is not None:
                temp_schema = {
                    '$schema':
                    'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema'}
                temp_schema.update(schema)
                try:
                    asdf_schema.validate(result, schema=temp_schema)
                except jsonschema.ValidationError:
                    if validate:
                        raise
                    else:
                        msgfmt = "'{0}' is not valid in keyword '{1}'"
                        warnings.warn(msgfmt.format(result, fits_keyword))
                        if pass_invalid_values:
                            properties.put_value(path, result, tree)
                else:
                    properties.put_value(path, result, tree)

        elif 'fits_hdu' in schema and (
                'max_ndim' in schema or 'ndim' in schema or 'datatype' in schema):
            result = _fits_array_loader(
                hdulist, schema, ctx.get('hdu_index'), known_datas)
            if result is not None:
                temp_schema = {
                    '$schema':
                    'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema'}
                temp_schema.update(schema)
                asdf_schema.validate(result, schema=temp_schema)
                properties.put_value(path, result, tree)

        if schema.get('type') == 'array':
            has_fits_hdu = _schema_has_fits_hdu(schema)
            if has_fits_hdu:
                for i in range(len(hdulist)):
                    recurse(schema['items'],
                            path + [i],
                            combiner,
                            {'hdu_index': i})
                return True

    mschema.walk_schema(schema, callback)
    return known_keywords, known_datas


def _load_extra_fits(hdulist, known_keywords, known_datas, tree):
    # Handle _extra_fits
    for hdu in hdulist:
        known = known_keywords.get(hdu, set())

        cards = []
        for key, val, comment in hdu.header.cards:
            if not (_is_builtin_fits_keyword(key) or
                    key in known):
                cards.append([key, val, comment])

        if len(cards):
            properties.put_value(
                ['extra_fits', hdu.name, 'header'], cards, tree)

        if hdu not in known_datas:
            if hdu.data is not None:
                properties.put_value(
                    ['extra_fits', hdu.name, 'data'], hdu.data, tree)


def _load_history(hdulist, tree):
    try:
        hdu = get_hdu(hdulist, 0)
    except AttributeError:
        return

    header = hdu.header
    if 'HISTORY' not in header:
        return

    history = tree['history'] = []

    for entry in header['HISTORY']:
        history.append(HistoryEntry({'description': entry}))


def from_fits(hdulist, schema, extensions=None, validate=True,
              pass_invalid_values=False):
    ff = fits_embed.AsdfInFits.open(hdulist, extensions=extensions)

    known_keywords, known_datas = _load_from_schema(
        hdulist, schema, ff.tree, validate,
        pass_invalid_values=pass_invalid_values)
    _load_extra_fits(hdulist, known_keywords, known_datas, ff.tree)
    _load_history(hdulist, ff.tree)

    return ff
