import datetime
import hashlib
import os
from pkg_resources import parse_version
import re
import warnings

import numpy as np
from astropy.io import fits
from astropy import time
from astropy.utils.exceptions import AstropyWarning
import asdf
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
from . import validate

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())


__all__ = ['to_fits', 'from_fits', 'fits_hdu_name', 'get_hdu']


_ASDF_GE_2_6 = parse_version(asdf.__version__) >= parse_version('2.6')


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

# Key where the FITS hash is stored in the ASDF tree
FITS_HASH_KEY = '_fits_hash'


def _get_indexed_keyword(keyword, i):
    for (sub, max, r) in _keyword_indices:
        if sub in keyword:
            if i >= max:
                raise ValueError(
                    "Too many entries for given keyword '{0}'".format(keyword))
            if r is None:
                val = str(i)
            else:
                val = r[i]
            keyword = keyword.replace(sub, val)

    return keyword


def fits_hdu_name(name):
    """
    Returns a FITS hdu name in the correct form for the current
    version of Python.
    """
    if isinstance(name, bytes):
        return name.decode('ascii')
    return name


def _get_hdu_name(schema):
    hdu_name = schema.get('fits_hdu')
    if hdu_name in (None, 'PRIMARY'):
        hdu_name = 0
    else:
        hdu_name = fits_hdu_name(hdu_name)
    return hdu_name


def _get_hdu_type(hdu_name, schema=None, value=None):
    hdu_type = None
    if hdu_name in (0, 'PRIMARY'):
        hdu_type = fits.PrimaryHDU
    elif schema is not None:
        dtype = ndarray.asdf_datatype_to_numpy_dtype(schema['datatype'])
        if dtype.fields is not None:
            hdu_type = fits.BinTableHDU
    elif value is not None:
        if hasattr(value, 'dtype') and value.dtype.names is not None:
            hdu_type = fits.BinTableHDU
    return hdu_type


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
            if isinstance(pair, str):
                hdu = hdulist[(pair, 1)]
            elif isinstance(pair, tuple) and index == 0:
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
        hdu_type = _get_hdu_type(hdu_name, value=value)
        if hdu_type is None:
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
            for key, val in hdu.header.items():
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

    for property, subschema in properties.items():
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
        raise ValueError("'fits_keyword' is not valid with type of 'array'")

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
    index = getattr(validator, 'sequence_index', 0)

    hdu_type = _get_hdu_type(hdu_name, schema=schema, value=instance)
    hdu = _get_or_make_hdu(validator.hdulist, hdu_name,
                           index=index, hdu_type=hdu_type)

    hdu.data = instance
    hdu.ver = index + 1


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
    return validators.Draft4Validator.VALIDATORS["type"](validator, items, instance, schema)


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
            node = str(time.Time(node, format='iso'))
        return node
    tree = treeutil.walk_and_modify(tree, convert_datetimes)

    if _ASDF_GE_2_6:
        kwargs = {"_visit_repeat_nodes": True}
    else:
        kwargs = {}

    validator = asdf_schema.get_validator(
        schema, None, FITS_VALIDATORS, FITS_SCHEMA_URL_MAPPING, **kwargs)

    validator.hdulist = hdulist
    # TODO: Handle comment stack on per-hdu-basis
    validator.comment_stack = []
    # This actually kicks off the saving
    validator.validate(tree, _schema=schema)


def _save_extra_fits(hdulist, tree):
    # Handle _extra_fits
    for hdu_name, parts in tree.get('extra_fits', {}).items():
        hdu_name = fits_hdu_name(hdu_name)
        if 'data' in parts:
            hdu_type = _get_hdu_type(hdu_name, value=parts['data'])
            hdu = _get_or_make_hdu(hdulist, hdu_name, hdu_type=hdu_type,
                                   value=parts['data'])
        if 'header' in parts:
            hdu = _get_or_make_hdu(hdulist, hdu_name)
            for key, val, comment in parts['header']:
                if _is_builtin_fits_keyword(key):
                    continue
                hdu.header.append((key, val, comment), end=True)


def _save_history(hdulist, tree):
    if 'history' not in tree:
        return

    # Support the older way of representing ASDF history entries
    if isinstance(tree['history'], list):
        history = tree['history']
    else:
        history = tree['history'].get('entries', [])

    for i in range(len(history)):
        # There is no guarantee the user has added proper HistoryEntry records
        if not isinstance(history[i], HistoryEntry):
            if isinstance(history[i], dict):
                history[i] = HistoryEntry(history[i])
            else:
                history[i] = HistoryEntry({'description': str(history[i])})
        hdulist[0].header['HISTORY'] = history[i]['description']


def to_fits(tree, schema):
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())

    _save_from_schema(hdulist, tree, schema)
    _save_extra_fits(hdulist, tree)
    _save_history(hdulist, tree)

    # Store the FITS hash in the tree
    tree[FITS_HASH_KEY] = fits_hash(hdulist)

    asdf = fits_embed.AsdfInFits(hdulist, tree)
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
    return from_fits_hdu(hdu, schema)


def _schema_has_fits_hdu(schema):
    has_fits_hdu = [False]

    for node in treeutil.iter_tree(schema):
        if isinstance(node, dict) and 'fits_hdu' in node:
            has_fits_hdu[0] = True

    return has_fits_hdu[0]


def _load_from_schema(hdulist, schema, tree, context, skip_fits_update=False):
    known_keywords = {}
    known_datas = set()

    # Check if there are any table HDU's. If not, this whole process
    # can be skipped.
    if skip_fits_update:
        if not any(isinstance(hdu, fits.BinTableHDU) for hdu in hdulist if hdu.name != 'ASDF'):
            log.debug('Skipping FITS updating completely.')
            return known_keywords, known_datas
        log.debug('Skipping FITS keyword updating except for BinTableHDU and its associated header keywords.')

    # Determine maximum EXTVER that could be used in finding named HDU's.
    # This is needed to constrain the loop over HDU's when resolving arrays.
    max_extver = max(hdu.ver for hdu in hdulist) if len(hdulist) else 0

    def callback(schema, path, combiner, ctx, recurse):
        result = None
        if not skip_fits_update and 'fits_keyword' in schema:
            fits_keyword = schema['fits_keyword']
            result = _fits_keyword_loader(
                hdulist, fits_keyword, schema,
                ctx.get('hdu_index'), known_keywords)

            if result is None:
                validate.value_change(path, result, schema,
                                      context._pass_invalid_values,
                                      context._strict_validation)
            else:
                if validate.value_change(path, result, schema,
                                         context._pass_invalid_values,
                                         context._strict_validation):
                    properties.put_value(path, result, tree)

        elif 'fits_hdu' in schema and (
                'max_ndim' in schema or 'ndim' in schema or 'datatype' in schema):
            result = _fits_array_loader(
                hdulist, schema, ctx.get('hdu_index'), known_datas)

            if result is None:
                validate.value_change(path, result, schema,
                                      context._pass_invalid_values,
                                      context._strict_validation)
            else:
                if validate.value_change(path, result, schema,
                                         context._pass_invalid_values,
                                         context._strict_validation):
                    properties.put_value(path, result, tree)

        if schema.get('type') == 'array':
            has_fits_hdu = _schema_has_fits_hdu(schema)
            if has_fits_hdu:
                for i in range(max_extver):
                    recurse(schema['items'],
                            path + [i],
                            combiner,
                            {'hdu_index': i})
                return True

    mschema.walk_schema(schema, callback)
    return known_keywords, known_datas


def _load_extra_fits(hdulist, known_keywords, known_datas, tree):
    # Remove any extra_fits from tree
    if 'extra_fits' in tree:
        del tree['extra_fits']

    # Add header keywords and data not in schema to extra_fits
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
            if hdu.name.lower() != 'asdf':
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

    history = tree['history'] = {'entries': []}

    for entry in header['HISTORY']:
        history['entries'].append(HistoryEntry({'description': entry}))


def from_fits(hdulist, schema, context, skip_fits_update=None, **kwargs):
    """Read model information from a FITS HDU list

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The FITS HDUList

    schema : dict
        The schema defining the ASDF > FITS_KEYWORD, FITS_HDU mapping.

    context: DataModel
        The `DataModel` to update

    skip_fits_update : bool or None
        When `False`, models opened from FITS files will proceed
        and load the FITS header values into the model.
        When `True` and the FITS file has an ASDF extension, the
        loading/validation of the FITS header will be skipped, loading
        the model only from the ASDF extension.
        When `None`, the value is taken from the environmental SKIP_FITS_UPDATE.
        Otherwise, the default is `False`
    """
    try:
        ff = from_fits_asdf(hdulist, **kwargs)
    except Exception as exc:
        raise exc.__class__("ERROR loading embedded ASDF: " + str(exc)) from exc

    # Determine whether skipping the FITS loading can be done.
    skip_fits_update = _verify_skip_fits_update(
        skip_fits_update, hdulist, ff, context
    )

    known_keywords, known_datas = _load_from_schema(
        hdulist, schema, ff.tree, context, skip_fits_update=skip_fits_update
    )
    if not skip_fits_update:
        _load_extra_fits(hdulist, known_keywords, known_datas, ff.tree)

    _load_history(hdulist, ff.tree)

    return ff


def from_fits_asdf(hdulist,
                   ignore_version_mismatch=True,
                   ignore_unrecognized_tag=False,
                   **kwargs):
    """
    Wrap asdf call to extract optional arguments
    """
    ignore_missing_extensions = kwargs.pop('ignore_missing_extensions')
    return fits_embed.AsdfInFits.open(hdulist,
                                      ignore_version_mismatch=ignore_version_mismatch,
                                      ignore_unrecognized_tag=ignore_unrecognized_tag,
                                      ignore_missing_extensions=ignore_missing_extensions)


def from_fits_hdu(hdu, schema):
    """
    Read the data from a fits hdu into a numpy ndarray
    """
    data = hdu.data

    # Save the column listeners for possible restoration
    if hasattr(data, '_coldefs'):
        listeners = data._coldefs._listeners
    else:
        listeners = None

    # Cast array to type mentioned in schema
    data = properties._cast(data, schema)

    # Casting a table loses the column listeners, so restore them
    if listeners is not None:
        data._coldefs._listeners = listeners

    return data


def _verify_skip_fits_update(skip_fits_update, hdulist, asdf_struct, context):
    """Ensure all conditions for skipping FITS updating are true

    Returns True if either 1) the FITS hash in the asdf structure matches the input
    FITS structure. Or 2) skipping has been explicitly asked for in `skip_fits_update`.

    Parameters
    ----------
    skip_fits_update : bool
        Regardless of FIT hash check, attempt to skip if requested.

    hdulist : astropy.io.fits.HDUList
        The input FITS information

    asdf_struct : asdf.ASDFFile
        The associated ASDF structure

    context : DataModel
        The DataModel being built.

    Returns
    -------
    skip_fits_update : bool
        All conditions are satisfied for skipping FITS updating.
    """
    if skip_fits_update is None:
        skip_fits_update = util.get_envar_as_boolean('SKIP_FITS_UPDATE', None)

    # If skipping has been explicitly disallowed, indicate as such.
    if skip_fits_update is False:
        return False

    # Skipping has either been requested or has been left to be determined automatically.
    # Continue checking conditions necessary for skipping.

    # Need an already existing ASDF. If not, cannot skip.
    if not len(asdf_struct.tree):
        log.debug('No ASDF information found. Cannot skip updating from FITS headers.')
        return False

    # Ensure model types match
    hdulist_class = util._class_from_model_type(hdulist)
    if hdulist_class is None:
        log.debug('Cannot determine model of the FITS file.'
                  ' Cannot skip updating from FITS headers.')
        return False
    if not isinstance(context, hdulist_class):
        log.debug(f'Input model {hdulist_class} does not match the'
                  f' requested model {type(context)}.'
                  ' Cannot skip updating from FITS headers.')
        return False

    # Check for FITS hash and compare to current. If equal, automatically skip.
    if asdf_struct.tree.get(FITS_HASH_KEY, None) is not None:
        if asdf_struct.tree[FITS_HASH_KEY] == fits_hash(hdulist):
            log.debug('FITS hash matches. Skipping FITS updating.')
            return True

    # If skip only if explicitly requested.
    return False if skip_fits_update is None else True


def fits_hash(hdulist):
    """Calculate a hash based on all HDU headers

    Uses basic SHA-256 hash to calculate.

    Parameters
    ----------
    hdulist : astropy.fits.HDUList
        The FITS structure.

    Returns
    -------
    fits_hash : str
        The hash of all HDU headers.
    """
    fits_hash = hashlib.sha256()

    # Ignore FITS header warnings, such as "Card is too long".
    # Such issues are inconsequential to hash calculation.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        fits_hash.update(''.join(
            str(hdu.header)
            for hdu in hdulist
            if hdu.name != 'ASDF').encode()
        )
    return fits_hash.hexdigest()
